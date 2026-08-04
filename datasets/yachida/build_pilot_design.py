#!/usr/bin/env python3
"""Build the deterministic 67/67/67 Yachida pilot matched to adenomas."""

from __future__ import annotations

import argparse
import csv
import hashlib
import pathlib
from collections import Counter, defaultdict


CONDITIONS = ("Control", "Adenoma", "CRC")


def age_bin(raw: str) -> str:
    try:
        age = float(raw)
    except ValueError:
        return "Unknown"
    if age < 50:
        return "<50"
    if age < 60:
        return "50-59"
    if age < 70:
        return "60-69"
    return ">=70"


def stable_hash(seed: str, sample: str) -> str:
    return hashlib.sha256(f"yachida-pilot-v1\0{seed}\0{sample}".encode("utf-8")).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--selection-seed", default="ground-truth-taxonomic-calibration-yachida-pilot-v1")
    parser.add_argument("--reference-condition", default="Adenoma", choices=CONDITIONS)
    parser.add_argument("--per-condition", type=int, default=67)
    args = parser.parse_args()

    with args.manifest.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames or []
    required = {"sample_id", "Target_Condition", "age", "sex"}
    missing = required.difference(fields)
    if missing:
        raise SystemExit(f"[ERROR] manifest lacks columns: {sorted(missing)}")
    if len({row["sample_id"] for row in rows}) != len(rows):
        raise SystemExit("[ERROR] sample_id values must be unique")

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        condition = row["Target_Condition"]
        if condition not in CONDITIONS:
            raise SystemExit(f"[ERROR] unrecognized condition: {condition!r}")
        row = dict(row)
        row["matching_age_bin"] = age_bin(row["age"])
        row["matching_sex"] = row["sex"].strip().title() or "Unknown"
        row["matching_stratum"] = f"{row['matching_sex']}|{row['matching_age_bin']}"
        row["selection_hash"] = stable_hash(args.selection_seed, row["sample_id"])
        grouped[condition].append(row)

    reference = grouped[args.reference_condition]
    if len(reference) != args.per_condition:
        raise SystemExit(
            f"[ERROR] {args.reference_condition} must contain exactly {args.per_condition} samples; "
            f"found {len(reference)}"
        )
    quotas = Counter(row["matching_stratum"] for row in reference)

    selected = []
    for condition in CONDITIONS:
        candidates = grouped[condition]
        by_stratum: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in candidates:
            by_stratum[row["matching_stratum"]].append(row)
        chosen = []
        for stratum in sorted(quotas):
            available = sorted(by_stratum[stratum], key=lambda row: (row["selection_hash"], row["sample_id"]))
            needed = quotas[stratum]
            if len(available) < needed:
                raise SystemExit(
                    f"[ERROR] {condition}/{stratum}: need {needed} samples to match "
                    f"{args.reference_condition}, found {len(available)}"
                )
            chosen.extend(available[:needed])
        chosen.sort(key=lambda row: (row["selection_hash"], row["sample_id"]))
        if len(chosen) != args.per_condition:
            raise SystemExit(f"[ERROR] internal selection error for {condition}: {len(chosen)}")
        for rank, row in enumerate(chosen, 1):
            row["selection_rank"] = str(rank)
            row["selection_seed"] = args.selection_seed
            row["matching_reference"] = args.reference_condition
            selected.append(row)

    for condition in CONDITIONS:
        observed = Counter(row["matching_stratum"] for row in selected if row["Target_Condition"] == condition)
        if observed != quotas:
            raise SystemExit(f"[ERROR] stratum mismatch for {condition}: {observed} != {quotas}")

    output_fields = fields + [
        "matching_age_bin", "matching_sex", "matching_stratum", "selection_hash",
        "selection_rank", "selection_seed", "matching_reference",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected)
    digest = hashlib.sha256(args.output.read_bytes()).hexdigest()
    args.output.with_suffix(args.output.suffix + ".sha256").write_text(
        f"{digest}  {args.output.name}\n", encoding="utf-8"
    )
    print(f"[OK] selected 201 samples: 67 Control, 67 Adenoma, 67 CRC -> {args.output}")
    print(f"[OK] All conditions exactly match {args.reference_condition} sex x age-bin counts: {dict(sorted(quotas.items()))}")


if __name__ == "__main__":
    main()
