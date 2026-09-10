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


def read_key_value_table(path: pathlib.Path) -> dict[str, str]:
    if not path.is_file() or not path.stat().st_size:
        raise SystemExit(f"[ERROR] Missing sample-completion table: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or rows[0] != ["field", "value"]:
        raise SystemExit(f"[ERROR] Invalid sample-completion table: {path}")
    values = {row[0]: row[1] for row in rows[1:] if len(row) == 2}
    if len(values) != len(rows) - 1:
        raise SystemExit(f"[ERROR] Duplicate or malformed sample-completion fields: {path}")
    return values


def count_success(root: pathlib.Path) -> int:
    return sum(1 for path in root.rglob("SUCCESS") if path.is_file()) if root.is_dir() else 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=pathlib.Path)
    parser.add_argument("--scratch-root", required=True, type=pathlib.Path)
    parser.add_argument("--results-root", required=True, type=pathlib.Path)
    parser.add_argument("--state-dir", required=True, type=pathlib.Path)
    parser.add_argument("--independent-manifest", required=True, type=pathlib.Path)
    parser.add_argument("--expected-samples", type=int, default=201)
    args = parser.parse_args()

    state = args.state_dir.resolve()
    seal = state / "production_seal"
    seal.mkdir(parents=True, exist_ok=True)
    # Invalidate publication authority before doing any potentially long audit.
    # A failed, cancelled, or interrupted rerun must not leave an older seal
    # looking current.
    for marker in (seal / "SUCCESS", seal / "production_seal.sha256"):
        marker.unlink(missing_ok=True)
    in_progress = seal / "AUDIT_IN_PROGRESS"
    in_progress.write_text(
        "status\tIN_PROGRESS\n"
        f"manifest\t{args.manifest.resolve()}\n",
        encoding="utf-8",
    )

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

    with args.independent_manifest.open(newline="", encoding="utf-8") as handle:
        independent_rows = list(csv.DictReader(handle, delimiter="\t"))
    if not independent_rows or "sample_id" not in independent_rows[0]:
        raise SystemExit("[ERROR] Independent manifest is empty or lacks sample_id")
    independent = [row["sample_id"] for row in independent_rows]
    if len(independent) != len(set(independent)) or not set(independent).issubset(ids):
        raise SystemExit("[ERROR] Independent manifest is duplicated or not nested in production manifest")
    if args.expected_samples == 201 and len(independent) != 30:
        raise SystemExit(f"[ERROR] Expected 30 independent-subset samples; observed {len(independent)}")
    independent_set = set(independent)

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
        sample_root = args.results_root.resolve() / "YachidaS_2019" / sample
        completion_values = read_key_value_table(sample_root / "sample_completion.tsv")
        baseline_profiles = count_success(sample_root / "profiles" / "baseline")
        independent_profiles = count_success(sample_root / "profiles" / "independent")
        community_profiles = count_success(sample_root / "profiles" / "community")
        expected_independent = 60 if sample in independent_set else 0
        expected_total = 1 + expected_independent + 7
        observed_total = baseline_profiles + independent_profiles + community_profiles
        expected_completion = {
            "sample_id": sample,
            "condition": row["Target_Condition"],
            "independent_subset": "1" if sample in independent_set else "0",
            "expected_profiles": str(expected_total),
            "observed_profiles": str(expected_total),
            "community_design_rows": "7",
            "independent_design_rows": str(expected_independent),
        }
        mismatches = {
            field: (completion_values.get(field), expected)
            for field, expected in expected_completion.items()
            if completion_values.get(field) != expected
        }
        if mismatches:
            raise SystemExit(f"[ERROR] Sample-completion mismatch for {sample}: {mismatches}")
        if (baseline_profiles, independent_profiles, community_profiles, observed_total) != (
            1, expected_independent, 7, expected_total
        ):
            raise SystemExit(
                f"[ERROR] Profile topology mismatch for {sample}: baseline={baseline_profiles}, "
                f"independent={independent_profiles}, community={community_profiles}, "
                f"total={observed_total}"
            )
        record = {
            "sample_id": sample,
            "condition": row["Target_Condition"],
            "batch_id": row["batch_id"],
            "retained_outputs": str(outputs),
            "baseline_profiles": str(baseline_profiles),
            "independent_profiles": str(independent_profiles),
            "community_profiles": str(community_profiles),
            "observed_profiles": str(observed_total),
            "status": "PASS",
        }
        completion.append(record)
        by_batch[row["batch_id"]].append(record)
        print(f"[OK] {sample}: {outputs} retained outputs")

    fields = [
        "sample_id", "condition", "batch_id", "retained_outputs",
        "baseline_profiles", "independent_profiles", "community_profiles",
        "observed_profiles", "status",
    ]
    for batch, batch_rows in sorted(by_batch.items()):
        root = state / "batches" / batch
        write_table(root / "batch_completion.tsv", fields, batch_rows)
        (root / "SUCCESS").write_text(
            f"batch_id\t{batch}\nsamples\t{len(batch_rows)}\nstatus\tPASS\n",
            encoding="utf-8",
        )

    write_table(seal / "dataset_completion.tsv", fields, completion)
    manifest_copy = seal / "pilot_batched.tsv"
    manifest_copy.write_bytes(args.manifest.read_bytes())
    independent_copy = seal / "independent_10_per_condition.tsv"
    independent_copy.write_bytes(args.independent_manifest.read_bytes())
    (seal / "SUCCESS").write_text(
        "dataset\tYachidaS_2019\n"
        f"samples\t{len(completion)}\n"
        f"batches\t{len(by_batch)}\n"
        f"independent_subset\t{len(independent)}\n"
        f"profiles\tbaseline={len(completion)};independent={len(independent) * 60};community={len(completion) * 7}\n"
        "conditions\tControl=67;Adenoma=67;CRC=67\n"
        "status\tPASS\n",
        encoding="utf-8",
    )
    sealed = [
        seal / "dataset_completion.tsv", seal / "dataset_completion.tsv.sha256",
        manifest_copy, independent_copy, seal / "SUCCESS",
    ]
    checksum = seal / "production_seal.sha256"
    checksum_temporary = seal / "production_seal.sha256.tmp"
    checksum_temporary.write_text(
        "".join(f"{digest(path)}  {path.name}\n" for path in sealed), encoding="utf-8"
    )
    checksum_temporary.replace(checksum)
    in_progress.unlink()
    print(f"[PASS] Yachida production sealed: {len(completion)} samples, {len(by_batch)} batches")
    print(f"[INFO] Seal: {seal}")


if __name__ == "__main__":
    main()
