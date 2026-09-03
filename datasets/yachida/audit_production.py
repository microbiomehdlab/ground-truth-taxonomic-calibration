#!/usr/bin/env python3
"""Rehash all Yachida receipts and seal the complete frozen production cohort."""
from __future__ import annotations

import argparse
import csv
import hashlib
import pathlib
from collections import Counter, defaultdict


def digest(path: pathlib.Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def verify_receipt(receipt: pathlib.Path, scratch: pathlib.Path) -> int:
    if not receipt.is_file():
        raise SystemExit(f"[ERROR] Missing receipt: {receipt}")
    with receipt.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or set(rows[0]) != {"path", "sha256", "bytes"}:
        raise SystemExit(f"[ERROR] Invalid or empty receipt: {receipt}")
    scratch = scratch.resolve()
    for row in rows:
        path = pathlib.Path(row["path"]).resolve()
        if path == scratch or scratch in path.parents:
            raise SystemExit(f"[ERROR] Persistent output is inside scratch: {path}")
        if not path.is_file() or not path.stat().st_size:
            raise SystemExit(f"[ERROR] Missing or empty retained output: {path}")
        if path.stat().st_size != int(row["bytes"]) or digest(path) != row["sha256"]:
            raise SystemExit(f"[ERROR] Retained-output mismatch: {path}")
    return len(rows)


def write_table(path: pathlib.Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    path.with_suffix(path.suffix + ".sha256").write_text(
        f"{digest(path)}  {path.name}\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=pathlib.Path)
    parser.add_argument("--scratch-root", required=True, type=pathlib.Path)
    parser.add_argument("--state-dir", required=True, type=pathlib.Path)
    parser.add_argument("--expected-samples", type=int, default=201)
    args = parser.parse_args()

    with args.manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"sample_id", "Target_Condition", "batch_id"}
    if not rows or not required.issubset(rows[0]):
        raise SystemExit("[ERROR] Manifest is empty or lacks required columns")
    if len(rows) != args.expected_samples:
        raise SystemExit(f"[ERROR] Expected {args.expected_samples} samples; observed {len(rows)}")
    ids = [row["sample_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise SystemExit("[ERROR] Manifest contains duplicate sample IDs")
    counts = Counter(row["Target_Condition"] for row in rows)
    if args.expected_samples == 201 and counts != {"Control": 67, "Adenoma": 67, "CRC": 67}:
        raise SystemExit(f"[ERROR] Expected frozen 67/67/67 design; observed {dict(counts)}")

    state = args.state_dir.resolve()
    sample_state = state / "samples"
    completion: list[dict[str, str]] = []
    by_batch: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        sample = row["sample_id"]
        marker = sample_state / f"{sample}.verified"
        if not marker.is_file() or not marker.stat().st_size:
            raise SystemExit(f"[ERROR] Missing verified marker: {sample}")
        outputs = verify_receipt(
            sample_state / f"{sample}.retained_outputs.tsv",
            args.scratch_root.resolve() / sample,
        )
        record = {
            "sample_id": sample,
            "condition": row["Target_Condition"],
            "batch_id": row["batch_id"],
            "retained_outputs": str(outputs),
            "status": "PASS",
        }
        completion.append(record)
        by_batch[row["batch_id"]].append(record)
        print(f"[OK] {sample}: {outputs} retained outputs")

    fields = ["sample_id", "condition", "batch_id", "retained_outputs", "status"]
    for batch, batch_rows in sorted(by_batch.items()):
        root = state / "batches" / batch
        write_table(root / "batch_completion.tsv", fields, batch_rows)
        (root / "SUCCESS").write_text(
            f"batch_id\t{batch}\nsamples\t{len(batch_rows)}\nstatus\tPASS\n",
            encoding="utf-8",
        )

    seal = state / "production_seal"
    write_table(seal / "dataset_completion.tsv", fields, completion)
    manifest_copy = seal / "pilot_batched.tsv"
    manifest_copy.write_bytes(args.manifest.read_bytes())
    (seal / "SUCCESS").write_text(
        "dataset\tYachidaS_2019\n"
        f"samples\t{len(completion)}\n"
        f"batches\t{len(by_batch)}\n"
        "conditions\tControl=67;Adenoma=67;CRC=67\n"
        "status\tPASS\n",
        encoding="utf-8",
    )
    sealed = [
        seal / "dataset_completion.tsv", seal / "dataset_completion.tsv.sha256",
        manifest_copy, seal / "SUCCESS",
    ]
    (seal / "production_seal.sha256").write_text(
        "".join(f"{digest(path)}  {path.name}\n" for path in sealed), encoding="utf-8"
    )
    print(f"[PASS] Yachida production sealed: {len(completion)} samples, {len(by_batch)} batches")
    print(f"[INFO] Seal: {seal}")


if __name__ == "__main__":
    main()
