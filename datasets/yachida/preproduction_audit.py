#!/usr/bin/env python3
"""Audit frozen manifests, spike design, pools, and integer allocations before production."""
from __future__ import annotations

import argparse
import csv
import hashlib
import pathlib
import subprocess
from collections import Counter
from decimal import Decimal


def table(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"[ERROR] Empty table: {path}")
    return rows


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_unique(rows: list[dict[str, str]], field: str, source: pathlib.Path) -> set[str]:
    values = [row[field] for row in rows]
    duplicates = sorted(value for value, count in Counter(values).items() if count > 1)
    if duplicates:
        raise SystemExit(f"[ERROR] Duplicate {field} values in {source}: {duplicates[:5]}")
    return set(values)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-manifest", required=True, type=pathlib.Path)
    parser.add_argument("--independent-manifest", required=True, type=pathlib.Path)
    parser.add_argument("--panel", required=True, type=pathlib.Path)
    parser.add_argument("--community", required=True, type=pathlib.Path)
    parser.add_argument("--pools-dir", required=True, type=pathlib.Path)
    parser.add_argument("--allocator", required=True, type=pathlib.Path)
    parser.add_argument("--paired-fastq-validator", required=True, type=pathlib.Path)
    parser.add_argument("--expected-pilot-per-condition", type=int, default=67)
    parser.add_argument("--expected-independent-per-condition", type=int, default=10)
    parser.add_argument("--verify-pool-contents", action="store_true")
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    pilot = table(args.pilot_manifest)
    independent = table(args.independent_manifest)
    panel = table(args.panel)
    community = table(args.community)
    pilot_ids = require_unique(pilot, "sample_id", args.pilot_manifest)
    independent_ids = require_unique(independent, "sample_id", args.independent_manifest)
    if not independent_ids <= pilot_ids:
        raise SystemExit("[ERROR] Independent-spike samples are not a subset of the pilot cohort")
    pilot_counts = Counter(row["Target_Condition"] for row in pilot)
    independent_counts = Counter(row["Target_Condition"] for row in independent)
    expected_conditions = {"Control", "Adenoma", "CRC"}
    if set(pilot_counts) != expected_conditions or set(pilot_counts.values()) != {args.expected_pilot_per_condition}:
        raise SystemExit(f"[ERROR] Unexpected pilot condition counts: {dict(pilot_counts)}")
    if set(independent_counts) != expected_conditions or set(independent_counts.values()) != {args.expected_independent_per_condition}:
        raise SystemExit(f"[ERROR] Unexpected independent condition counts: {dict(independent_counts)}")

    panel_labels = require_unique(panel, "label", args.panel)
    community_labels = require_unique(community, "label", args.community)
    if panel_labels != community_labels or len(panel_labels) != 10:
        raise SystemExit("[ERROR] Panel and community must contain the same ten unique labels")
    weights = {row["label"]: Decimal(row.get("weight", "1") or "1") for row in community}
    if any(weight <= 0 for weight in weights.values()):
        raise SystemExit("[ERROR] Community weights must be positive")

    pair_index = args.pools_dir / "pool_pair_counts.tsv"
    pair_seal = args.pools_dir / "pool_pair_counts.tsv.sha256"
    pool_manifest = args.pools_dir / "pool_files.sha256"
    pool_manifest_seal = args.pools_dir / "pool_files.sha256.sha256"
    for path in (pair_index, pair_seal, pool_manifest, pool_manifest_seal):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"[ERROR] Missing finalized-pool provenance: {path}")
    expected_pair_sha = pair_seal.read_text(encoding="utf-8").split()[0]
    if sha256(pair_index) != expected_pair_sha:
        raise SystemExit("[ERROR] Pool pair-count index checksum mismatch")
    expected_manifest_sha = pool_manifest_seal.read_text(encoding="utf-8").split()[0]
    if sha256(pool_manifest) != expected_manifest_sha:
        raise SystemExit("[ERROR] Pool checksum-manifest seal mismatch")
    indexed = table(pair_index)
    indexed_labels = require_unique(indexed, "label", pair_index)
    if indexed_labels != panel_labels:
        raise SystemExit("[ERROR] Pool index labels differ from the spike panel")
    capacities: dict[str, int] = {}
    for row in indexed:
        n1, n2, total = (int(row[key]) for key in ("mate1_pairs", "mate2_pairs", "pool_pairs"))
        if n1 <= 0 or n1 != n2 or n1 != total:
            raise SystemExit(f"[ERROR] Invalid pool pair counts for {row['label']}: {row}")
        capacities[row["label"]] = total
        for mate in (1, 2):
            path = args.pools_dir / f"{row['label']}.pool_{mate}.fq"
            if not path.is_file() or path.stat().st_size == 0:
                raise SystemExit(f"[ERROR] Missing/empty pool: {path}")
    checksum_entries = [line for line in pool_manifest.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(checksum_entries) != 20:
        raise SystemExit(f"[ERROR] Expected 20 pool checksum entries; found {len(checksum_entries)}")
    if args.verify_pool_contents:
        subprocess.run(["sha256sum", "-c", str(pool_manifest)], cwd=args.pools_dir, check=True)
        pool_qc = args.output_dir / "pool_pair_integrity"
        pool_qc.mkdir(parents=True, exist_ok=True)
        for label in sorted(panel_labels):
            subprocess.run(
                [
                    str(args.paired_fastq_validator),
                    "--r1", str(args.pools_dir / f"{label}.pool_1.fq"),
                    "--r2", str(args.pools_dir / f"{label}.pool_2.fq"),
                    "--output", str(pool_qc / f"{label}.paired_fastq_integrity.tsv"),
                ],
                check=True,
            )

    members = [f"{row['label']}:{row.get('weight', '1') or '1'}" for row in community]
    for total in range(1, 1001):
        command = [str(args.allocator), "--total", str(total)]
        for member in members:
            command.extend(("--member", member))
        result = subprocess.run(command, text=True, capture_output=True, check=True)
        allocated = list(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
        counts = [int(row["n"]) for row in allocated]
        if len(counts) != 10 or sum(counts) != total or min(counts) < 0:
            raise SystemExit(f"[ERROR] Invalid equal-weight community allocation for N={total}: {counts}")

    report = args.output_dir / "preproduction_audit.tsv"
    report.write_text(
        "check\tstatus\tdetail\n"
        f"pilot_design\tPASS\t{dict(sorted(pilot_counts.items()))}\n"
        f"independent_subset\tPASS\t{dict(sorted(independent_counts.items()))}\n"
        f"target_panel\tPASS\t{len(panel_labels)} unique targets\n"
        f"pool_index\tPASS\t{len(capacities)} synchronized paired pools\n"
        f"pool_content_and_pair_integrity\t{'PASS' if args.verify_pool_contents else 'NOT_RUN'}\t{'checksums, FASTQ structure, and mate IDs verified' if args.verify_pool_contents else 'enable --verify-pool-contents for the one-time full scan'}\n"
        "community_allocation\tPASS\tdeterministic exact allocation checked for totals 1-1000\n",
        encoding="utf-8",
    )
    digest = sha256(report)
    (args.output_dir / "preproduction_audit.tsv.sha256").write_text(f"{digest}  preproduction_audit.tsv\n", encoding="utf-8")
    (args.output_dir / "SUCCESS").write_text("preproduction_audit\tPASS\n", encoding="utf-8")
    print(f"[PASS] Pre-production audit completed: {report}")


if __name__ == "__main__":
    main()
