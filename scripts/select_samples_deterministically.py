#!/usr/bin/env python3
"""Select a balanced subset using stable accession hashes, not row order."""

from __future__ import annotations

import argparse
import csv
import hashlib
import pathlib


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--per-condition", type=int, required=True)
    parser.add_argument("--selection-seed", default="ground-truth-taxonomic-calibration-yachida-v1")
    parser.add_argument("--id-column", default="sample_id")
    parser.add_argument("--condition-column", default="Target_Condition")
    args = parser.parse_args()

    with args.manifest.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames or []
    for field in (args.id_column, args.condition_column):
        if field not in fields:
            raise SystemExit(f"[ERROR] missing column {field!r}")

    seen = set()
    grouped: dict[str, list[tuple[str, dict[str, str]]]] = {}
    for row in rows:
        sample = row[args.id_column]
        if not sample or sample in seen:
            raise SystemExit(f"[ERROR] blank or duplicate sample ID: {sample!r}")
        seen.add(sample)
        key = hashlib.sha256(
            f"stable-selection-v1\0{args.selection_seed}\0{sample}".encode("utf-8")
        ).hexdigest()
        grouped.setdefault(row[args.condition_column], []).append((key, row))

    selected = []
    for condition in sorted(grouped):
        candidates = sorted(grouped[condition], key=lambda item: (item[0], item[1][args.id_column]))
        if len(candidates) < args.per_condition:
            raise SystemExit(f"[ERROR] {condition}: need {args.per_condition}, found {len(candidates)}")
        for rank, (key, row) in enumerate(candidates[: args.per_condition], 1):
            selected.append({**row, "selection_rank": str(rank), "selection_hash": key,
                             "selection_seed": args.selection_seed})

    output_fields = fields + ["selection_rank", "selection_hash", "selection_seed"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected)
    print(f"[OK] selected {len(selected)} samples across {len(grouped)} conditions -> {args.output}")


if __name__ == "__main__":
    main()
