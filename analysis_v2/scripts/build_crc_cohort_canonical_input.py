#!/usr/bin/env python3
"""Build a complete, fail-closed CRC-cohort canonical v2 evidence table."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from decimal import Decimal
from pathlib import Path

from build_assembly_sensitivity_input import (
    HEADER, SCHEMA, bracken_abundance, metaphlan_abundance, read_rows,
    require_columns, unique_file,
)


PROFILERS = ("kraken2_bracken", "metaphlan4")


def fraction_tag(raw: str) -> str:
    text = ("{:.6f}".format(float(raw))).rstrip("0").rstrip(".")
    return "f" + text.replace(".", "p")


def allocate(total: int, panel: list[dict[str, str]]) -> dict[str, int]:
    """Reproduce the frozen community allocator exactly."""
    weights = [Decimal(row.get("weight") or "1") for row in panel]
    total_weight = sum(weights, Decimal(0))
    shares = [Decimal(total) * weight / total_weight for weight in weights]
    values = [int(share.to_integral_value()) for share in shares[:-1]]
    values.append(total - sum(values))
    if values[-1] < 0:
        values = [int(share) for share in shares]
        order = sorted(range(len(shares)), key=lambda i: (-(shares[i] - int(shares[i])), i))
        for index in order[: total - sum(values)]:
            values[index] += 1
    if sum(values) != total or min(values) < 0:
        raise ValueError("community allocation failed")
    return {row["label"]: value for row, value in zip(panel, values)}


class Profiles:
    def __init__(self) -> None:
        self.cache: dict[tuple[Path, str, str], float] = {}

    def abundance(self, path: Path, profiler: str, alias: str) -> float:
        key = (path, profiler, alias)
        if key not in self.cache:
            reader = bracken_abundance if profiler == "kraken2_bracken" else metaphlan_abundance
            self.cache[key] = reader(path, alias)
        return self.cache[key]


def make_row(meta: dict[str, str], population: str, target: dict[str, str],
             profiler: str, profile_id: str, baseline_id: str, total_fraction: float,
             target_fraction: float, pairs: int, abundance: float, profile: Path,
             design: str, cohort: str) -> dict[str, str]:
    unit = "fraction_total_reads" if profiler == "kraken2_bracken" else "relative_abundance_pct"
    abundance_fraction = abundance if profiler == "kraken2_bracken" else abundance / 100.0
    return {
        "schema_version": SCHEMA, "cohort": cohort, "study": meta["study"],
        "sample_id": meta["sample_id"], "condition": meta["condition"],
        "analysis_population": population, "target_label": target["label"],
        "target_taxon": target["taxon_name"], "assembly_arm": "original",
        "profiler": profiler, "profile_id": profile_id,
        "baseline_profile_id": baseline_id,
        "spike_fraction_total": format(total_fraction, ".17g"),
        "spike_fraction_target": format(target_fraction, ".17g"),
        "implanted_read_pairs_target": str(pairs),
        "native_abundance": format(abundance, ".17g"), "native_unit": unit,
        "abundance_fraction": format(abundance_fraction, ".17g"),
        "detected_native_nonzero": "1" if abundance_fraction > 0 else "0",
        "source_profile": str(profile.resolve()), "source_design": design,
        "include": "1", "exclusion_reason": "",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--independent-manifest", required=True, type=Path)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--spike-panel", required=True, type=Path)
    parser.add_argument("--aliases", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--cohort", choices=("yachida", "feng", "zeller"), default="yachida")
    parser.add_argument("--expected-samples", type=int, default=201)
    parser.add_argument("--expected-independent", type=int, default=30)
    parser.add_argument("--expected-independent-doses", type=int, default=6)
    parser.add_argument("--expected-community-doses", type=int, default=7)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        manifest = read_rows(args.manifest)
        independent_manifest = read_rows(args.independent_manifest)
        panel = read_rows(args.spike_panel)
        aliases = read_rows(args.aliases, delimiter=",")
        require_columns(manifest, {"sample_id"}, "manifest")
        for row in manifest:
            row["study"] = row.get("study") or row.get("Study") or ""
            row["condition"] = row.get("condition") or row.get("Target_Condition") or ""
            if not row["study"] or not row["condition"]:
                raise ValueError("manifest lacks study/condition values")
        require_columns(independent_manifest, {"sample_id"}, "independent manifest")
        require_columns(panel, {"label", "taxon_name", "weight"}, "spike panel")
        require_columns(aliases, {"canonical", "alias", "tool"}, "aliases")
        if len(manifest) != args.expected_samples:
            raise ValueError(f"expected {args.expected_samples} samples; found {len(manifest)}")
        independent = {row["sample_id"] for row in independent_manifest}
        if len(independent) != args.expected_independent:
            raise ValueError(f"expected {args.expected_independent} independent samples; found {len(independent)}")
        all_samples = {row["sample_id"] for row in manifest}
        if len(all_samples) != len(manifest) or not independent.issubset(all_samples):
            raise ValueError("sample manifests are duplicated or not nested")
        alias_map = {(row["canonical"], row["tool"]): row["alias"] for row in aliases}
        profiles = Profiles()
        output_rows: list[dict[str, str]] = []

        for meta in manifest:
            sample = meta["sample_id"]
            sample_root = args.results_root / meta["study"] / sample
            baseline_root = sample_root / "profiles" / "baseline"
            community_design = sample_root / "spike_design" / "community" / "CRCpanel.tsv"
            community_rows = read_rows(community_design)
            require_columns(community_rows, {"sample_id", "fraction", "R", "N_total", "f_hat"}, str(community_design))
            community_rows = [row for row in community_rows if row["sample_id"] == sample]
            if len(community_rows) != args.expected_community_doses:
                raise ValueError(f"expected {args.expected_community_doses} community doses in {community_design}; found {len(community_rows)}")

            for target in panel:
                for profiler in PROFILERS:
                    alias = alias_map.get((target["taxon_name"], profiler))
                    if not alias:
                        raise ValueError(f"missing alias for {target['taxon_name']} {profiler}")
                    baseline_profile = unique_file(baseline_root, sample + (".bracken.S.tsv" if profiler == "kraken2_bracken" else ".metaphlan.tsv"))
                    baseline_abundance = profiles.abundance(baseline_profile, profiler, alias)
                    output_rows.append(make_row(meta, "community", target, profiler, sample, sample,
                                                0.0, 0.0, 0, baseline_abundance,
                                                baseline_profile, "BASELINE", args.cohort))
                    for design in sorted(community_rows, key=lambda row: float(row["fraction"])):
                        profile_id = f"{sample}_CRCpanel_{fraction_tag(design['fraction'])}"
                        suffix = ".bracken.S.tsv" if profiler == "kraken2_bracken" else ".metaphlan.tsv"
                        profile = unique_file(sample_root / "profiles" / "community", profile_id + suffix)
                        abundance = profiles.abundance(profile, profiler, alias)
                        total_pairs = int(design["N_total"])
                        target_pairs = allocate(total_pairs, panel)[target["label"]]
                        denominator = int(design["R"]) + total_pairs
                        output_rows.append(make_row(
                            meta, "community", target, profiler, profile_id, sample,
                            float(design["f_hat"]), target_pairs / denominator, target_pairs,
                            abundance, profile, str(community_design.resolve()), args.cohort))

                if sample not in independent:
                    continue
                independent_design = sample_root / "spike_design" / "independent" / f"{target['label']}.tsv"
                independent_rows = read_rows(independent_design)
                require_columns(independent_rows, {"sample_id", "label", "fraction", "N_inserted", "f_hat"}, str(independent_design))
                independent_rows = [row for row in independent_rows if row["sample_id"] == sample and row["label"] == target["label"]]
                if len(independent_rows) != args.expected_independent_doses:
                    raise ValueError(f"expected {args.expected_independent_doses} independent doses in {independent_design}; found {len(independent_rows)}")
                for profiler in PROFILERS:
                    alias = alias_map[(target["taxon_name"], profiler)]
                    baseline_profile = unique_file(baseline_root, sample + (".bracken.S.tsv" if profiler == "kraken2_bracken" else ".metaphlan.tsv"))
                    baseline_abundance = profiles.abundance(baseline_profile, profiler, alias)
                    output_rows.append(make_row(meta, "independent", target, profiler, sample, sample,
                                                0.0, 0.0, 0, baseline_abundance,
                                                baseline_profile, "BASELINE", args.cohort))
                    for design in sorted(independent_rows, key=lambda row: float(row["fraction"])):
                        profile_id = f"{sample}_{target['label']}_{fraction_tag(design['fraction'])}"
                        suffix = ".bracken.S.tsv" if profiler == "kraken2_bracken" else ".metaphlan.tsv"
                        profile = unique_file(sample_root / "profiles" / "independent" / target["label"], profile_id + suffix)
                        abundance = profiles.abundance(profile, profiler, alias)
                        achieved = float(design["f_hat"])
                        output_rows.append(make_row(
                            meta, "independent", target, profiler, profile_id, sample,
                            achieved, achieved, int(design["N_inserted"]), abundance,
                            profile, str(independent_design.resolve()), args.cohort))

        expected_rows = args.expected_samples * len(panel) * len(PROFILERS) * (args.expected_community_doses + 1)
        expected_rows += args.expected_independent * len(panel) * len(PROFILERS) * (args.expected_independent_doses + 1)
        if len(output_rows) != expected_rows:
            raise ValueError(f"expected {expected_rows} canonical rows; built {len(output_rows)}")
        args.outdir.mkdir(parents=True, exist_ok=True)
        output = args.outdir / "canonical_input.tsv"
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=HEADER, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(output_rows)
        validator = Path(__file__).with_name("validate_canonical_input.py")
        subprocess.run([sys.executable, str(validator), "--input", str(output),
                        "--outdir", str(args.outdir / "validation")], check=True)
        print(f"[PASS] Complete {args.cohort} canonical input built: {len(output_rows)} rows")
        print(f"[INFO] {output}")
    except (ValueError, OSError, KeyError, subprocess.CalledProcessError) as error:
        raise SystemExit(f"[ERROR] {error}") from error


if __name__ == "__main__":
    main()
