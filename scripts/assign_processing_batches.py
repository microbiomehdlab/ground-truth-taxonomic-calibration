#!/usr/bin/env python3
"""Assign a frozen selected cohort to deterministic, balanced small batches."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import pathlib
from collections import defaultdict, deque


def key(seed: str, sample: str) -> str:
    return hashlib.sha256(f"processing-batch-v1\0{seed}\0{sample}".encode("utf-8")).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--batch-dir", type=pathlib.Path, required=True)
    parser.add_argument("--max-batch-size", type=int, default=10)
    parser.add_argument("--batch-seed", default="ground-truth-taxonomic-calibration-yachida-batches-v1")
    args = parser.parse_args()
    if args.max_batch_size < 1:
        raise SystemExit("[ERROR] max batch size must be positive")

    with args.manifest.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames or []
    required = {"sample_id", "Target_Condition"}
    missing = required.difference(fields)
    if missing:
        raise SystemExit(f"[ERROR] manifest lacks columns: {sorted(missing)}")

    grouped: dict[str, deque[dict[str, str]]] = {}
    temporary: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        row = dict(row)
        row["batch_hash"] = key(args.batch_seed, row["sample_id"])
        temporary[row["Target_Condition"]].append(row)
    for condition, members in temporary.items():
        grouped[condition] = deque(sorted(members, key=lambda row: (row["batch_hash"], row["sample_id"])))

    ordered = []
    condition_order = [condition for condition in ("Control", "Adenoma", "CRC") if condition in grouped]
    condition_order += sorted(set(grouped).difference(condition_order))
    while any(grouped[condition] for condition in condition_order):
        for condition in condition_order:
            if grouped[condition]:
                ordered.append(grouped[condition].popleft())

    n_batches = math.ceil(len(ordered) / args.max_batch_size)
    base, extra = divmod(len(ordered), n_batches)
    sizes = [base + (index < extra) for index in range(n_batches)]
    if max(sizes, default=0) > args.max_batch_size:
        raise SystemExit("[ERROR] internal batch-size error")

    output_fields = fields + ["batch_hash", "processing_order", "batch_id", "batch_position", "batch_size", "batch_seed"]
    assigned = []
    offset = 0
    for batch_index, size in enumerate(sizes, 1):
        batch_id = f"batch_{batch_index:03d}"
        batch_rows = ordered[offset: offset + size]
        offset += size
        for position, row in enumerate(batch_rows, 1):
            assigned.append({**row, "processing_order": str(len(assigned) + 1), "batch_id": batch_id,
                             "batch_position": str(position), "batch_size": str(size),
                             "batch_seed": args.batch_seed})

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.batch_dir.mkdir(parents=True, exist_ok=True)
    stale = list(args.batch_dir.glob("batch_*.tsv"))
    if stale:
        raise SystemExit(f"[ERROR] batch directory already contains manifests; use an empty directory: {args.batch_dir}")
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(assigned)
    for batch_index in range(1, n_batches + 1):
        batch_id = f"batch_{batch_index:03d}"
        batch_rows = [row for row in assigned if row["batch_id"] == batch_id]
        path = args.batch_dir / f"{batch_id}.tsv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(batch_rows)
    checksum_paths = [args.output, *sorted(args.batch_dir.glob("batch_*.tsv"))]
    with (args.batch_dir / "SHA256SUMS").open("w", encoding="utf-8") as handle:
        for path in checksum_paths:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            relative = os.path.relpath(path.resolve(), args.batch_dir.resolve())
            handle.write(f"{digest}  {relative}\n")
    print(f"[OK] assigned {len(assigned)} samples to {n_batches} batches (sizes {min(sizes)}-{max(sizes)})")
    print(f"[OK] master manifest: {args.output}")
    print(f"[OK] per-batch manifests: {args.batch_dir}")


if __name__ == "__main__":
    main()
