#!/usr/bin/env python3
"""Validate and describe a deterministic condition-balanced sample subset."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import pathlib
import statistics


def read_table(path: pathlib.Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        return rows, list(reader.fieldnames or [])


def stable_hash(seed: str, sample: str) -> str:
    return hashlib.sha256(
        f"stable-selection-v1\0{seed}\0{sample}".encode("utf-8")
    ).hexdigest()


def numeric(values: list[str]) -> list[float]:
    result = []
    for value in values:
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(number):
            result.append(number)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eligible-manifest", required=True, type=pathlib.Path)
    parser.add_argument("--selection", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    parser.add_argument("--id-column", default="Name")
    parser.add_argument("--condition-column", default="Study condition")
    parser.add_argument("--age-column", default="Age")
    parser.add_argument("--sex-column", default="Sex")
    parser.add_argument("--bmi-column", default="BMI")
    parser.add_argument("--expected-condition", action="append", dest="conditions")
    parser.add_argument("--per-condition", type=int, default=10)
    args = parser.parse_args()

    expected_conditions = args.conditions or ["Control", "Adenoma", "CRC"]
    eligible, eligible_fields = read_table(args.eligible_manifest)
    selected, selected_fields = read_table(args.selection)
    required_eligible = {
        args.id_column, args.condition_column, args.age_column, args.sex_column,
    }
    required_selected = required_eligible | {
        "selection_rank", "selection_hash", "selection_seed",
    }
    if missing := required_eligible.difference(eligible_fields):
        raise SystemExit(f"[ERROR] eligible manifest lacks fields: {sorted(missing)}")
    if missing := required_selected.difference(selected_fields):
        raise SystemExit(f"[ERROR] selection lacks fields: {sorted(missing)}")

    eligible_by_id: dict[str, dict[str, str]] = {}
    for row in eligible:
        sample = row[args.id_column]
        if not sample or sample in eligible_by_id:
            raise SystemExit(f"[ERROR] blank or duplicate eligible sample ID: {sample!r}")
        eligible_by_id[sample] = row

    selected_ids = [row[args.id_column] for row in selected]
    if len(selected_ids) != len(set(selected_ids)):
        raise SystemExit("[ERROR] selected sample IDs are not unique")
    if not set(selected_ids) <= set(eligible_by_id):
        raise SystemExit("[ERROR] selection is not nested within the eligible manifest")

    observed_conditions = {row[args.condition_column] for row in selected}
    if observed_conditions != set(expected_conditions):
        raise SystemExit(
            f"[ERROR] selected conditions {sorted(observed_conditions)} != "
            f"expected {sorted(expected_conditions)}"
        )

    counts = {condition: 0 for condition in expected_conditions}
    seeds = set()
    for row in selected:
        sample = row[args.id_column]
        condition = row[args.condition_column]
        counts[condition] += 1
        seeds.add(row["selection_seed"])
        if row["selection_hash"] != stable_hash(row["selection_seed"], sample):
            raise SystemExit(f"[ERROR] invalid selection hash: {sample}")
        source = eligible_by_id[sample]
        if source[args.condition_column] != condition:
            raise SystemExit(f"[ERROR] condition changed between manifests: {sample}")
    if any(count != args.per_condition for count in counts.values()):
        raise SystemExit(f"[ERROR] unexpected selected counts: {counts}")
    if len(seeds) != 1 or "" in seeds:
        raise SystemExit(f"[ERROR] expected one non-empty selection seed, found: {seeds}")

    for condition in expected_conditions:
        rows = [row for row in selected if row[args.condition_column] == condition]
        observed_ranks = sorted(int(row["selection_rank"]) for row in rows)
        if observed_ranks != list(range(1, args.per_condition + 1)):
            raise SystemExit(f"[ERROR] invalid ranks for {condition}: {observed_ranks}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / "selection_balance.tsv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        fields = ["condition", "selected_n", "variable", "level", "count", "mean", "minimum", "maximum", "missing"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for condition in expected_conditions:
            rows = [row for row in selected if row[args.condition_column] == condition]
            for variable, column in (("sex", args.sex_column),):
                levels: dict[str, int] = {}
                for row in rows:
                    level = row.get(column, "") or "MISSING"
                    levels[level] = levels.get(level, 0) + 1
                for level in sorted(levels):
                    writer.writerow({"condition": condition, "selected_n": len(rows), "variable": variable,
                                     "level": level, "count": levels[level], "mean": "", "minimum": "",
                                     "maximum": "", "missing": levels.get("MISSING", 0)})
            for variable, column in (("age", args.age_column), ("bmi", args.bmi_column)):
                raw = [row.get(column, "") for row in rows]
                values = numeric(raw)
                writer.writerow({"condition": condition, "selected_n": len(rows), "variable": variable,
                                 "level": "continuous", "count": len(values),
                                 "mean": f"{statistics.fmean(values):.6g}" if values else "",
                                 "minimum": f"{min(values):.6g}" if values else "",
                                 "maximum": f"{max(values):.6g}" if values else "",
                                 "missing": len(raw) - len(values)})

    methods_path = args.output_dir / "selection_methods.txt"
    methods_path.write_text(
        "eligible_rule\tall rows in the frozen eligible cohort manifest\n"
        f"independent_rule\t{args.per_condition} samples selected without replacement within each condition\n"
        f"conditions\t{','.join(expected_conditions)}\n"
        f"selection_seed\t{next(iter(seeds))}\n"
        "tie_breaking\tSHA-256 stable sample hash, then sample identifier\n"
        "selection_inputs\tsample identifier and diagnostic condition only\n"
        "outcome_use\tno profiler, recovery, spike, or biomarker outcome entered selection\n"
        "demographic_use\tage, sex, and BMI summarized after selection; not used to replace samples\n",
        encoding="utf-8",
    )
    digest = hashlib.sha256(summary_path.read_bytes()).hexdigest()
    (args.output_dir / "selection_balance.tsv.sha256").write_text(
        f"{digest}  {summary_path.name}\n", encoding="utf-8"
    )
    print(f"[PASS] deterministic independent selection validated: {counts}")
    print(f"[INFO] Audit: {args.output_dir}")


if __name__ == "__main__":
    main()
