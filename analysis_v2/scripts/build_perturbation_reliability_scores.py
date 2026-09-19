#!/usr/bin/env python3
"""Build descriptive, target-excluded reliability scores for native biomarkers."""
from __future__ import annotations

import argparse
import csv
import hashlib
import math
import statistics
from collections import defaultdict
from pathlib import Path

GROUP = ("cohort", "analysis_population", "profiler", "contrast", "feature")
REQUIRED = set(GROUP) | {
    "target_label", "assembly_arm", "spike_fraction_target", "q_threshold",
    "feature_role", "baseline_called", "dose_called", "baseline_effect",
    "dose_effect", "effect_sign_changed",
}


def number(value: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite value: {value}")
    return result


def render(value: float) -> str:
    return format(value, ".17g")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def median(values: list[float]) -> float:
    return statistics.median(values)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transition-ledger", required=True, type=Path)
    parser.add_argument("--artifact-scores", type=Path)
    parser.add_argument("--disease-results", type=Path,
                        help="Full disease-model table for external replication labels")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--contrast", default="CRC_vs_Control")
    parser.add_argument("--assembly-arm", default="original")
    parser.add_argument("--q-threshold", type=float, default=0.05)
    parser.add_argument("--stability-threshold", type=float, default=0.80)
    args = parser.parse_args()
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise SystemExit("[ERROR] OUTDIR must be new or empty")
    if not 0 < args.q_threshold < 1 or not 0 <= args.stability_threshold <= 1:
        raise SystemExit("[ERROR] invalid threshold")

    grouped: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    with args.transition_ledger.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not REQUIRED <= set(reader.fieldnames or []):
            raise SystemExit("[ERROR] transition ledger columns are incomplete")
        seen = set()
        for row in reader:
            if row["contrast"] != args.contrast or row["assembly_arm"] != args.assembly_arm:
                continue
            if abs(number(row["q_threshold"]) - args.q_threshold) > 1e-12:
                continue
            if row["feature_role"] != "bystander" or number(row["spike_fraction_target"]) <= 0:
                continue
            unique = tuple(row[x] for x in GROUP) + (
                row["target_label"], render(number(row["spike_fraction_target"]))
            )
            if unique in seen:
                raise SystemExit(f"[ERROR] duplicate stress-test row: {unique}")
            seen.add(unique)
            grouped[tuple(row[x] for x in GROUP)].append(row)
    if not grouped:
        raise SystemExit("[ERROR] no eligible bystander rows")

    # Artificial off-target recurrence is learned separately. It is one score
    # component, never a deletion rule and never uses held-out target identity.
    artifact: dict[tuple[str, str, str, str], float] = {}
    if args.artifact_scores:
        with args.artifact_scores.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            needed = {"cohort", "analysis_population", "profiler", "feature", "artifact_score"}
            if not needed <= set(reader.fieldnames or []):
                raise SystemExit("[ERROR] artifact-score columns are incomplete")
            for row in reader:
                key = tuple(row[x] for x in ("cohort", "analysis_population", "profiler", "feature"))
                value = number(row["artifact_score"])
                if not 0 <= value <= 1:
                    raise SystemExit("[ERROR] artifact score outside [0,1]")
                if key in artifact and abs(artifact[key] - value) > 1e-12:
                    raise SystemExit(f"[ERROR] conflicting artifact score: {key}")
                artifact[key] = value

    # Baseline summaries from the sparse transition ledger are sufficient for
    # scoring, but not for declaring that a feature was tested in another cohort.
    baseline = {}
    for key, rows in grouped.items():
        baseline[key] = {
            "effect": median([number(row["baseline_effect"]) for row in rows]),
            "called_rate": sum(int(row["baseline_called"]) for row in rows) / len(rows),
        }
    cohorts = sorted({key[0] for key in grouped})
    external_baseline = baseline
    if args.disease_results:
        full = defaultdict(lambda: {"effects": [], "called": []})
        needed = set(GROUP) | {"dose_level", "assembly_arm", "effect", "q_value", "include"}
        with args.disease_results.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not needed <= set(reader.fieldnames or []):
                raise SystemExit("[ERROR] full disease-result columns are incomplete")
            for row in reader:
                if (row["include"] != "1" or row["dose_level"] != "baseline" or
                        row["contrast"] != args.contrast or
                        row["assembly_arm"] != args.assembly_arm):
                    continue
                key = tuple(row[x] for x in GROUP)
                full[key]["effects"].append(number(row["effect"]))
                full[key]["called"].append(number(row["q_value"]) <= args.q_threshold)
        external_baseline = {
            key: {"effect": median(values["effects"]),
                  "called_rate": sum(values["called"]) / len(values["called"])}
            for key, values in full.items() if values["effects"]
        }
        if not external_baseline:
            raise SystemExit("[ERROR] no eligible full baseline disease results")

    output = []
    for key, rows in sorted(grouped.items()):
        cohort, population, profiler, contrast, feature = key
        # Reliability is defined for baseline biomarkers, not for all tested taxa.
        called_rows = [row for row in rows if int(row["baseline_called"]) == 1]
        if not called_rows:
            continue
        persistence = sum(int(row["dose_called"]) for row in called_rows) / len(called_rows)
        direction = 1 - sum(int(row["effect_sign_changed"]) for row in called_rows) / len(called_rows)
        fidelities = []
        for row in called_rows:
            before = number(row["baseline_effect"])
            after = number(row["dose_effect"])
            # Bounded symmetric closeness: 1 is unchanged, 0 is maximal change.
            denominator = abs(before) + abs(after)
            fidelities.append(1 - abs(after - before) / denominator if denominator else 1.0)
        effect_fidelity = sum(fidelities) / len(fidelities)
        artifact_key = (cohort, population, profiler, feature)
        artifact_score = artifact.get(artifact_key)
        components = [persistence, direction, effect_fidelity]
        if artifact_score is not None:
            components.append(1 - artifact_score)
        reliability = sum(components) / len(components)

        source_effect = baseline[key]["effect"]
        external = []
        for other in cohorts:
            if other == cohort:
                continue
            other_key = (other, population, profiler, contrast, feature)
            if other_key not in external_baseline:
                continue
            other_summary = external_baseline[other_key]
            external.append((other, other_summary["called_rate"], other_summary["effect"]))
        significant_same = [x for x in external if x[1] >= 0.5 and x[2] * source_effect > 0]
        significant_opposite = [x for x in external if x[1] >= 0.5 and x[2] * source_effect < 0]
        if significant_same:
            replication = "DIRECTIONALLY_REPLICATED"
        elif significant_opposite:
            replication = "DIRECTIONALLY_DISCORDANT"
        elif external:
            replication = "NOT_SIGNIFICANT_ELSEWHERE"
        else:
            replication = ("NOT_AVAILABLE_IN_SPARSE_LEDGER" if not args.disease_results
                           else "NOT_IN_OTHER_COHORT_UNIVERSE")
        stable = reliability >= args.stability_threshold
        if replication == "DIRECTIONALLY_REPLICATED":
            tier = "replicated_stable" if stable else "replicated_sensitive"
        elif replication == "DIRECTIONALLY_DISCORDANT":
            tier = "directionally_discordant"
        else:
            tier = "cohort_specific_stable" if stable else "cohort_specific_sensitive"
        output.append({
            "cohort": cohort, "analysis_population": population, "profiler": profiler,
            "contrast": contrast, "feature": feature,
            "baseline_effect": render(source_effect),
            "eligible_stress_tests": len(called_rows),
            "eligible_implanted_targets": len({row["target_label"] for row in called_rows}),
            "perturbation_persistence": render(persistence),
            "direction_stability": render(direction),
            "effect_fidelity": render(effect_fidelity),
            "artificial_artifact_score": "NA" if artifact_score is None else render(artifact_score),
            "artifact_resistance": "NA" if artifact_score is None else render(1 - artifact_score),
            "score_components": len(components),
            "perturbation_reliability_score": render(reliability),
            "external_cohorts_evaluable": len(external),
            "external_replication_status": replication,
            "evidence_tier": tier,
        })
    if not output:
        raise SystemExit("[ERROR] no baseline biomarkers were eligible")

    args.outdir.mkdir(parents=True, exist_ok=True)
    score_path = args.outdir / "taxon_reliability_scores.tsv"
    fields = tuple(output[0])
    with score_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)

    summary_path = args.outdir / "evidence_tier_summary.tsv"
    counts = defaultdict(int)
    for row in output:
        counts[(row["cohort"], row["analysis_population"], row["profiler"], row["evidence_tier"])] += 1
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        fields2 = ("cohort", "analysis_population", "profiler", "evidence_tier", "biomarkers")
        writer = csv.DictWriter(handle, fieldnames=fields2, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for key2, count in sorted(counts.items()):
            writer.writerow(dict(zip(fields2[:-1], key2), biomarkers=count))

    settings = args.outdir / "reliability_score_settings.tsv"
    settings.write_text(
        "field\tvalue\n"
        f"contrast\t{args.contrast}\nassembly_arm\t{args.assembly_arm}\n"
        f"q_threshold\t{render(args.q_threshold)}\nstability_threshold\t{render(args.stability_threshold)}\n"
        "component_weighting\tequal\n"
        "interpretation\trobustness_annotation_not_false_positive_probability\n"
        "replication_in_score\tNO\nstatus\tDEVELOPMENT_ONLY\n",
        encoding="utf-8",
    )
    note = args.outdir / "DEVELOPMENT_ONLY.txt"
    note.write_text(
        "status=DEVELOPMENT_ONLY\nautomatic_feature_removal=NO\n"
        "score_is_false_positive_probability=NO\nheldout_validation_required=YES\n",
        encoding="utf-8",
    )
    manifest = args.outdir / "reliability_scores.sha256"
    inputs = [args.transition_ledger]
    if args.artifact_scores:
        inputs.append(args.artifact_scores)
    if args.disease_results:
        inputs.append(args.disease_results)
    outputs = [score_path, summary_path, settings, note]
    manifest.write_text(
        "".join(f"{digest(path)}  {path.resolve()}\n" for path in inputs + outputs),
        encoding="utf-8",
    )
    (args.outdir / "SUCCESS").write_text(
        f"biomarkers={len(output)}\ncohorts={len(cohorts)}\nstatus=PASS\n", encoding="utf-8"
    )
    print(f"[PASS] Perturbation reliability scores: {args.outdir}")


if __name__ == "__main__":
    main()
