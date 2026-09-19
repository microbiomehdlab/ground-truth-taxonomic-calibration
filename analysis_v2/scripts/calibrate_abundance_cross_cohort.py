#!/usr/bin/env python3
"""Learn spike-induced abundance inflation in training cohorts and correct held-out cohorts.

The correction is deliberately conservative: only reproducible positive
departures from the ideal dilution reference are subtracted, and direct target
features are protected. Outputs retain the disease-model input schema.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import math
from collections import defaultdict
from pathlib import Path

MANIFEST_KEY = ("cohort", "study", "analysis_population", "sample_id", "condition",
                "target_label", "assembly_arm", "profiler", "profile_id",
                "baseline_profile_id")
MANIFEST_REQUIRED = set(MANIFEST_KEY) | {
    "spike_fraction_target", "dose_level", "source_profile", "target_feature",
    "age", "sex", "bmi", "include", "exclusion_reason",
}
ENDPOINT_REQUIRED = set(MANIFEST_KEY) | {
    "spike_fraction_total", "spike_fraction_target", "source_baseline_profile",
    "source_profile",
}
ABUNDANCE_REQUIRED = {"profiler", "source_profile", "feature", "abundance_fraction"}


def num(value: str) -> float:
    x = float(value)
    if not math.isfinite(x):
        raise ValueError(f"non-finite number: {value}")
    return x


def render(value: float) -> str:
    return format(value, ".17g")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def table(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"empty table: {path}")
        return list(reader), list(reader.fieldnames)


def key(row: dict[str, str], fields=MANIFEST_KEY) -> tuple[str, ...]:
    return tuple(row[x] for x in fields)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--profile-manifest", required=True, type=Path)
    p.add_argument("--abundance-long", required=True, type=Path)
    p.add_argument("--paired-endpoints", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--validation-cohort", action="append",
                   help="Held-out cohort; repeat for several. Default: every cohort.")
    p.add_argument("--min-prevalence", type=float, default=0.10)
    p.add_argument("--min-training-samples", type=int, default=10)
    p.add_argument("--min-positive-residual-fraction", type=float, default=0.60)
    args = p.parse_args()
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise SystemExit("[ERROR] OUTDIR must be new or empty")
    if not 0 <= args.min_prevalence <= 1 or not 0 <= args.min_positive_residual_fraction <= 1:
        raise SystemExit("[ERROR] invalid fraction threshold")
    if args.min_training_samples < 2:
        raise SystemExit("[ERROR] min-training-samples must be at least two")

    try:
        manifest, manifest_fields = table(args.profile_manifest)
        endpoints, endpoint_fields = table(args.paired_endpoints)
        abundance_rows, abundance_fields = table(args.abundance_long)
        if not MANIFEST_REQUIRED <= set(manifest_fields):
            raise ValueError("profile manifest columns are incomplete")
        if not ENDPOINT_REQUIRED <= set(endpoint_fields):
            raise ValueError("paired endpoint columns are incomplete")
        if not ABUNDANCE_REQUIRED <= set(abundance_fields):
            raise ValueError("abundance columns are incomplete")
        manifest = [row for row in manifest if row["include"] == "1"]
        if not manifest:
            raise ValueError("no included profile contexts")
        if len({key(row) + (row["dose_level"],) for row in manifest}) != len(manifest):
            raise ValueError("duplicate included manifest context")
        endpoint_map = {key(row): row for row in endpoints}
        if len(endpoint_map) != len(endpoints):
            raise ValueError("duplicate paired endpoint context")

        profiles: dict[tuple[str, str], dict[str, float]] = defaultdict(dict)
        for row in abundance_rows:
            value = num(row["abundance_fraction"])
            if value < 0 or value > 1.00001:
                raise ValueError("abundance outside [0,1]")
            token = (row["profiler"], row["source_profile"])
            if row["feature"] in profiles[token]:
                raise ValueError(f"duplicate abundance feature: {token} {row['feature']}")
            profiles[token][row["feature"]] = value
        del abundance_rows

        # Attach exact endpoint fractions and baseline paths to positive contexts.
        positive = []
        for row in manifest:
            f = num(row["spike_fraction_target"])
            if row["dose_level"] == "baseline":
                if f != 0:
                    raise ValueError("baseline has positive target fraction")
                continue
            endpoint = endpoint_map.get(key(row))
            if endpoint is None:
                raise ValueError(f"positive context lacks endpoint: {key(row)}")
            F, endpoint_f = num(endpoint["spike_fraction_total"]), num(endpoint["spike_fraction_target"])
            if not 0 < endpoint_f <= F < 1 or abs(endpoint_f - f) > 1e-12:
                raise ValueError("invalid or inconsistent endpoint fractions")
            if endpoint["source_profile"] != row["source_profile"]:
                raise ValueError("manifest/endpoint perturbed path mismatch")
            attached = dict(row, _F=F, _f=f, _baseline=endpoint["source_baseline_profile"])
            positive.append(attached)
        if not positive:
            raise ValueError("no positive-dose contexts")

        cohorts = sorted({row["cohort"] for row in manifest})
        validation = args.validation_cohort or cohorts
        if len(cohorts) < 2 or any(x not in cohorts for x in validation):
            raise ValueError("cross-cohort calibration requires valid held-out cohorts")

        # A community source represents one joint mixture but occurs once per
        # component target in the manifest. Collapse it before learning/applying.
        obs_groups = defaultdict(list)
        for row in positive:
            obs_key = (row["cohort"], row["study"], row["analysis_population"],
                       row["sample_id"], row["condition"], row["assembly_arm"],
                       row["profiler"], row["source_profile"], row["_baseline"])
            obs_groups[obs_key].append(row)
        observations = []
        for obs_key, rows in obs_groups.items():
            Fs = {render(row["_F"]) for row in rows}
            if len(Fs) != 1:
                raise ValueError("one profile has conflicting total implanted fractions")
            population = rows[0]["analysis_population"]
            community = population == "community"
            if not community and len({row["target_label"] for row in rows}) != 1:
                raise ValueError("independent profile maps to multiple targets")
            perturbation = "COMMUNITY_MIXTURE" if community else rows[0]["target_label"]
            predictor = rows[0]["_F"] if community else rows[0]["_f"]
            observations.append({
                "key": obs_key, "rows": rows, "cohort": rows[0]["cohort"],
                "population": population, "profiler": rows[0]["profiler"],
                "sample_id": rows[0]["sample_id"], "source": rows[0]["source_profile"],
                "baseline": rows[0]["_baseline"], "F": rows[0]["_F"],
                "predictor": predictor, "perturbation": perturbation,
                "direct": {row["target_feature"] for row in rows},
            })

        # Baseline prevalence universes prevent observed-only residual selection.
        baseline_records = {}
        for row in manifest:
            if row["dose_level"] != "baseline":
                continue
            token = (row["cohort"], row["analysis_population"], row["profiler"],
                     row["sample_id"], row["source_profile"])
            baseline_records[token] = profiles.get((row["profiler"], row["source_profile"]), {})

        args.outdir.mkdir(parents=True)
        top_ledger = []
        for test in validation:
            train = tuple(x for x in cohorts if x != test)
            label = "_".join(train) + "__to__" + test
            root = args.outdir / "transfers" / label
            root.mkdir(parents=True)

            universes = defaultdict(set)
            baseline_by_group = defaultdict(list)
            for (cohort, population, profiler, sample, source), values in baseline_records.items():
                if cohort in train:
                    baseline_by_group[(population, profiler)].append(values)
            for group, sample_profiles in baseline_by_group.items():
                counts = defaultdict(int)
                for values in sample_profiles:
                    for feature, value in values.items():
                        if value > 0:
                            counts[feature] += 1
                n = len(sample_profiles)
                universes[group] = {feature for feature, count in counts.items()
                                    if count / n >= args.min_prevalence}
            if not universes:
                raise ValueError(f"no training baseline universes for {label}")

            accum = defaultdict(lambda: {"num": 0.0, "den": 0.0, "n": 0,
                                         "positive": 0, "samples": set()})
            for obs in observations:
                if obs["cohort"] not in train:
                    continue
                baseline_values = profiles.get((obs["profiler"], obs["baseline"]), {})
                observed_values = profiles.get((obs["profiler"], obs["source"]), {})
                x = obs["predictor"]
                for feature in universes[(obs["population"], obs["profiler"])] - obs["direct"]:
                    expected = (1 - obs["F"]) * baseline_values.get(feature, 0.0)
                    delta = observed_values.get(feature, 0.0) - expected
                    item = accum[(obs["population"], obs["profiler"], obs["perturbation"], feature)]
                    item["num"] += x * delta; item["den"] += x * x; item["n"] += 1
                    item["positive"] += delta > 0
                    item["samples"].add((obs["cohort"], obs["sample_id"]))
            coefficients = {}
            coefficient_rows = []
            for ckey, item in sorted(accum.items()):
                slope = item["num"] / item["den"] if item["den"] else 0.0
                positive_fraction = item["positive"] / item["n"]
                active = (slope > 0 and len(item["samples"]) >= args.min_training_samples and
                          positive_fraction >= args.min_positive_residual_fraction)
                coefficients[ckey] = slope if active else 0.0
                coefficient_rows.append({
                    "training_cohorts": ",".join(train), "validation_cohort": test,
                    "analysis_population": ckey[0], "profiler": ckey[1],
                    "perturbation": ckey[2], "feature": ckey[3],
                    "slope_fraction_per_implanted_fraction": render(slope),
                    "observations": item["n"], "training_samples": len(item["samples"]),
                    "positive_residual_fraction": render(positive_fraction),
                    "correction_active": int(active),
                })

            # Generate held-out uncorrected and corrected model inputs.
            test_manifest = [dict(row) for row in manifest if row["cohort"] == test]
            corrected_source = {}
            correction_metrics = []
            corrected_profiles = {}
            for obs in observations:
                if obs["cohort"] != test:
                    continue
                source_key = (obs["population"], obs["profiler"], obs["source"])
                virtual = str((root / "virtual_profiles" /
                               (hashlib.sha256("\r".join(source_key).encode()).hexdigest() + ".profile")).resolve())
                corrected_source[source_key] = virtual
                baseline_values = profiles.get((obs["profiler"], obs["baseline"]), {})
                observed_values = profiles.get((obs["profiler"], obs["source"]), {})
                corrected = {}
                before_abs = []; after_abs = []; altered = 0; improved = 0
                evaluation_universe = universes.get((obs["population"], obs["profiler"]), set())
                for feature, abundance in observed_values.items():
                    if feature in obs["direct"]:
                        value = abundance
                    else:
                        ckey = (obs["population"], obs["profiler"], obs["perturbation"], feature)
                        predicted = max(0.0, coefficients.get(ckey, 0.0) * obs["predictor"])
                        value = max(0.0, abundance - predicted)
                        altered += abs(value - abundance) > 1e-15
                    if value > 0:
                        corrected[feature] = value
                corrected_profiles[(obs["profiler"], virtual)] = corrected
                for feature in evaluation_universe - obs["direct"]:
                    expected = (1 - obs["F"]) * baseline_values.get(feature, 0.0)
                    before = abs(observed_values.get(feature, 0.0) - expected)
                    after = abs(corrected.get(feature, 0.0) - expected)
                    before_abs.append(before); after_abs.append(after); improved += after < before
                correction_metrics.append({
                    "training_cohorts": ",".join(train), "validation_cohort": test,
                    "analysis_population": obs["population"], "profiler": obs["profiler"],
                    "sample_id": obs["sample_id"], "source_profile": obs["source"],
                    "perturbation": obs["perturbation"], "spike_fraction_total": render(obs["F"]),
                    "features_evaluated": len(before_abs), "features_altered": altered,
                    "features_improved": improved,
                    "mae_before": render(sum(before_abs) / len(before_abs)) if before_abs else "NA",
                    "mae_after": render(sum(after_abs) / len(after_abs)) if after_abs else "NA",
                    "relative_mae_restoration": render(1 - sum(after_abs) / sum(before_abs)) if sum(before_abs) else "NA",
                    "direct_targets_protected": len(obs["direct"]),
                })

            corrected_manifest = []
            for row in test_manifest:
                record = dict(row)
                if row["dose_level"] != "baseline":
                    source_key = (row["analysis_population"], row["profiler"], row["source_profile"])
                    if source_key not in corrected_source:
                        raise ValueError("held-out positive profile lacks correction mapping")
                    record["source_profile"] = corrected_source[source_key]
                corrected_manifest.append(record)

            baseline_tokens = {(row["profiler"], row["source_profile"]) for row in test_manifest
                               if row["dose_level"] == "baseline"}
            uncorrected_tokens = {(row["profiler"], row["source_profile"]) for row in test_manifest}
            def materialize(tokens):
                result = []
                for profiler, source in sorted(tokens):
                    for feature, value in sorted(profiles.get((profiler, source), {}).items()):
                        result.append({"profiler": profiler, "source_profile": source,
                                       "feature": feature, "abundance_fraction": render(value)})
                return result
            uncorrected_abundance = materialize(uncorrected_tokens)
            corrected_abundance = materialize(baseline_tokens)
            for (profiler, source), values in sorted(corrected_profiles.items()):
                corrected_abundance.extend({"profiler": profiler, "source_profile": source,
                                            "feature": feature, "abundance_fraction": render(value)}
                                           for feature, value in sorted(values.items()))

            def write(name: str, fields: list[str] | tuple[str, ...], rows) -> Path:
                path = root / name
                with path.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                    writer.writeheader(); writer.writerows(rows)
                return path

            files = []
            files.append(write("calibration_coefficients.tsv", list(coefficient_rows[0]), coefficient_rows))
            files.append(write("abundance_restoration_metrics.tsv", list(correction_metrics[0]), correction_metrics))
            files.append(write("uncorrected_profile_manifest.tsv", manifest_fields, test_manifest))
            files.append(write("corrected_profile_manifest.tsv", manifest_fields, corrected_manifest))
            files.append(write("uncorrected_abundance_long.tsv", abundance_fields, uncorrected_abundance))
            files.append(write("corrected_abundance_long.tsv", abundance_fields, corrected_abundance))
            settings = root / "calibration_settings.tsv"
            settings.write_text(
                "field\tvalue\n" +
                f"training_cohorts\t{','.join(train)}\nvalidation_cohort\t{test}\n" +
                f"min_prevalence\t{render(args.min_prevalence)}\n" +
                f"min_training_samples\t{args.min_training_samples}\n" +
                f"min_positive_residual_fraction\t{render(args.min_positive_residual_fraction)}\n" +
                "negative_or_unstable_correction\tzero\ndirect_target_policy\tprotected\n" +
                "community_policy\tjoint_mixture_total_fraction\nstatus\tDEVELOPMENT_ONLY\n",
                encoding="utf-8")
            files.append(settings)
            note = root / "DEVELOPMENT_ONLY.txt"
            note.write_text("status=DEVELOPMENT_ONLY\nheldout_cohort_used_for_training=NO\n", encoding="utf-8")
            files.append(note)
            checksum = root / "calibration.sha256"
            checksum.write_text("".join(f"{sha256(path)}  {path.resolve()}\n"
                                        for path in [args.profile_manifest, args.abundance_long,
                                                     args.paired_endpoints] + files), encoding="utf-8")
            (root / "SUCCESS").write_text("status=PASS\n", encoding="utf-8")
            top_ledger.append({"transfer": label, "training_cohorts": ",".join(train),
                               "validation_cohort": test, "directory": str(root.resolve()),
                               "active_coefficients": sum(int(row["correction_active"])
                                                          for row in coefficient_rows),
                               "validation_profiles": len(corrected_profiles)})

        ledger = args.outdir / "transfer_ledger.tsv"
        with ledger.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(top_ledger[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(top_ledger)
        (args.outdir / "DEVELOPMENT_ONLY.txt").write_text(
            "status=DEVELOPMENT_ONLY\nautomatic_feature_removal=NO\n", encoding="utf-8")
        (args.outdir / "SUCCESS").write_text(f"transfers={len(top_ledger)}\nstatus=PASS\n", encoding="utf-8")
        print(f"[PASS] Cross-cohort abundance calibration inputs: {args.outdir}")
    except (ValueError, OSError, KeyError) as error:
        raise SystemExit(f"[ERROR] {error}") from error


if __name__ == "__main__":
    main()
