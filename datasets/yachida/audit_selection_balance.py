#!/usr/bin/env python3
"""Describe and freeze Yachida pilot and nested-subset selection provenance."""
from __future__ import annotations

import argparse
import csv
import hashlib
import pathlib
import statistics
from collections import Counter


CONDITIONS = ("Control", "Adenoma", "CRC")


def read_table(path: pathlib.Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames or []
    if not rows:
        raise SystemExit(f"[ERROR] Empty manifest: {path}")
    return rows, fields


def digest(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * probability
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-manifest", required=True, type=pathlib.Path)
    parser.add_argument("--independent-manifest", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    args = parser.parse_args()
    pilot, pilot_fields = read_table(args.pilot_manifest)
    independent, independent_fields = read_table(args.independent_manifest)
    required = {"sample_id", "Target_Condition", "age", "sex", "matching_age_bin", "matching_stratum"}
    missing = required.difference(pilot_fields)
    if missing:
        raise SystemExit(f"[ERROR] Pilot manifest lacks fields: {sorted(missing)}")
    independent_required = {"sample_id", "Target_Condition", "selection_hash", "selection_rank", "selection_seed"}
    missing = independent_required.difference(independent_fields)
    if missing:
        raise SystemExit(f"[ERROR] Independent manifest lacks fields: {sorted(missing)}")
    by_id = {row["sample_id"]: row for row in pilot}
    if len(by_id) != len(pilot):
        raise SystemExit("[ERROR] Pilot sample IDs are not unique")
    independent_ids = [row["sample_id"] for row in independent]
    if len(set(independent_ids)) != len(independent_ids) or not set(independent_ids) <= set(by_id):
        raise SystemExit("[ERROR] Independent sample IDs are duplicate or absent from the pilot")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    frozen = args.output_dir / "independent_selection_provenance.tsv"
    fields = [
        "sample_id", "Target_Condition", "age", "sex", "matching_age_bin", "matching_stratum",
        "pilot_selection_rank", "pilot_selection_hash", "pilot_selection_seed", "pilot_matching_reference",
        "independent_selection_rank", "independent_selection_hash", "independent_selection_seed",
    ]
    with frozen.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for nested in independent:
            parent = by_id[nested["sample_id"]]
            writer.writerow({
                "sample_id": nested["sample_id"],
                "Target_Condition": parent["Target_Condition"],
                "age": parent["age"],
                "sex": parent["sex"],
                "matching_age_bin": parent["matching_age_bin"],
                "matching_stratum": parent["matching_stratum"],
                "pilot_selection_rank": parent.get("selection_rank", ""),
                "pilot_selection_hash": parent.get("selection_hash", ""),
                "pilot_selection_seed": parent.get("selection_seed", ""),
                "pilot_matching_reference": parent.get("matching_reference", ""),
                "independent_selection_rank": nested["selection_rank"],
                "independent_selection_hash": nested["selection_hash"],
                "independent_selection_seed": nested["selection_seed"],
            })

    summary = args.output_dir / "selection_balance.tsv"
    with summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("set", "condition", "variable", "level", "n", "proportion", "value"))
        for set_name, rows in (("pilot_67", pilot), ("independent_10", independent)):
            for condition in CONDITIONS:
                group = [by_id[row["sample_id"]] if set_name == "independent_10" else row
                         for row in rows if row["Target_Condition"] == condition]
                if not group:
                    raise SystemExit(f"[ERROR] No {set_name}/{condition} samples")
                for variable in ("sex", "matching_age_bin", "matching_stratum"):
                    counts = Counter(row[variable] for row in group)
                    for level in sorted(counts):
                        writer.writerow((set_name, condition, variable, level, counts[level],
                                         f"{counts[level] / len(group):.6f}", ""))
                ages = [float(row["age"]) for row in group]
                age_values = {
                    "mean": statistics.fmean(ages), "median": statistics.median(ages),
                    "q1": percentile(ages, 0.25), "q3": percentile(ages, 0.75),
                    "minimum": min(ages), "maximum": max(ages),
                }
                for level, value in age_values.items():
                    writer.writerow((set_name, condition, "age", level, len(ages), "", f"{value:.6f}"))

    methods = args.output_dir / "selection_methods.txt"
    pilot_seed = sorted({row.get("selection_seed", "") for row in pilot})
    independent_seed = sorted({row.get("selection_seed", "") for row in independent})
    methods.write_text(
        "pilot_rule\tAll 67 adenomas retained; 67 controls and 67 CRC samples selected to exactly match adenoma sex-by-age-bin counts\n"
        "pilot_age_bins\t<50, 50-59, 60-69, >=70\n"
        f"pilot_selection_seed\t{','.join(pilot_seed)}\n"
        "pilot_tie_breaking\tSHA-256 stable sample hash, then sample identifier\n"
        "independent_rule\tTen samples selected without replacement within each condition from the frozen pilot\n"
        f"independent_selection_seed\t{','.join(independent_seed)}\n"
        "independent_tie_breaking\tSHA-256 stable sample hash, then sample identifier\n"
        "outcome_use\tNo profiling, recovery, or biomarker result entered either selection\n",
        encoding="utf-8",
    )
    checksum = args.output_dir / "SHA256SUMS"
    checksum.write_text("".join(f"{digest(path)}  {path.name}\n" for path in (frozen, summary, methods)), encoding="utf-8")
    print(f"[PASS] Selection provenance and descriptive balance written under: {args.output_dir}")


if __name__ == "__main__":
    main()
