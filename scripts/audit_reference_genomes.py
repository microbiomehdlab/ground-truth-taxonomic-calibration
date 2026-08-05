#!/usr/bin/env python3
"""Checksum target FASTAs and optionally compare them with an expected audit."""

from __future__ import annotations

import argparse
import csv
import hashlib
import pathlib


FIELDS = [
    "label", "assembly", "fasta", "bytes", "sha256",
    "normalized_sequence_sha256", "contigs", "bases",
]


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def normalized_content(path: pathlib.Path) -> tuple[str, int, int]:
    """Hash the sorted multiset of contig-sequence hashes.

    This intentionally ignores FASTA headers, line wrapping, and contig order,
    while retaining contig multiplicity and exact nucleotide content.
    """
    records: list[str] = []
    sequence: list[str] = []
    bases = 0

    def finish_record() -> None:
        nonlocal bases
        if not sequence:
            return
        joined = "".join(sequence).upper()
        records.append(hashlib.sha256(joined.encode()).hexdigest())
        bases += len(joined)
        sequence.clear()

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                finish_record()
            else:
                sequence.append("".join(line.split()))
    finish_record()
    if not records:
        raise SystemExit(f"[ERROR] FASTA contains no sequences: {path}")
    digest = hashlib.sha256("\n".join(sorted(records)).encode()).hexdigest()
    return digest, len(records), bases


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument(
        "--fasta-root", type=pathlib.Path,
        help="Override panel FASTA paths with FASTA_ROOT/<label>.fa",
    )
    parser.add_argument("--expected", type=pathlib.Path, help="Expected audit TSV to compare")
    parser.add_argument("--comparison-output", type=pathlib.Path)
    args = parser.parse_args()

    panel_rows = read_tsv(args.panel)
    if not panel_rows or not {"label", "assembly", "fasta"}.issubset(panel_rows[0]):
        raise SystemExit("[ERROR] Panel must contain label, assembly, and fasta columns")

    observed: list[dict[str, object]] = []
    for row in panel_rows:
        path = args.fasta_root / f"{row['label']}.fa" if args.fasta_root else pathlib.Path(row["fasta"])
        if not path.is_file() or not path.stat().st_size:
            raise SystemExit(f"[ERROR] Missing/empty FASTA for {row['label']}: {path}")
        normalized, contigs, bases = normalized_content(path)
        observed.append({
            "label": row["label"],
            "assembly": row["assembly"],
            "fasta": str(path),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
            "normalized_sequence_sha256": normalized,
            "contigs": contigs,
            "bases": bases,
        })
    write_tsv(args.output, observed, FIELDS)
    print(f"[OK] Wrote reference audit: {args.output}")

    if not args.expected:
        return
    if not args.comparison_output:
        raise SystemExit("[ERROR] --comparison-output is required with --expected")
    expected_rows = read_tsv(args.expected)
    expected = {(row["label"], row["assembly"]): row for row in expected_rows}
    observed_keys = {(str(row["label"]), str(row["assembly"])) for row in observed}
    if observed_keys != set(expected):
        missing = sorted(set(expected) - observed_keys)
        extra = sorted(observed_keys - set(expected))
        raise SystemExit(f"[ERROR] Audit keys differ; missing={missing}, extra={extra}")

    comparison: list[dict[str, object]] = []
    mismatch = False
    for row in observed:
        key = (str(row["label"]), str(row["assembly"]))
        exp = expected[key]
        normalized_match = str(row["normalized_sequence_sha256"]) == exp["normalized_sequence_sha256"]
        contigs_match = str(row["contigs"]) == exp["contigs"]
        bases_match = str(row["bases"]) == exp["bases"]
        biological_match = normalized_match and contigs_match and bases_match
        mismatch |= not biological_match
        comparison.append({
            "label": key[0],
            "assembly": key[1],
            "byte_sha256_match": str(row["sha256"]) == exp["sha256"],
            "normalized_sequence_match": normalized_match,
            "contig_count_match": contigs_match,
            "base_count_match": bases_match,
            "biological_content_match": biological_match,
        })
    comparison_fields = [
        "label", "assembly", "byte_sha256_match", "normalized_sequence_match",
        "contig_count_match", "base_count_match", "biological_content_match",
    ]
    write_tsv(args.comparison_output, comparison, comparison_fields)
    if mismatch:
        raise SystemExit(f"[ERROR] Reference biological-content mismatch; see {args.comparison_output}")
    print(f"[PASS] All {len(comparison)} reference genomes have identical biological sequence content")
    print(f"[OK] Wrote comparison: {args.comparison_output}")


if __name__ == "__main__":
    main()
