#!/usr/bin/env python3
"""Extract the authoritative clade-to-genome-size mapping from a MetaPhlAn database.

This reads the exact MetaPhlAn database pickle used for production profiling and
emits a versioned, checksummed mapping table. It deliberately does not assume a
fixed internal schema: it inspects the loaded object, reports what it found, and
fails closed when it cannot identify a taxonomy mapping with genome lengths.

Genome length here is the database's representative genome length for a clade.
It is NOT marker length, and marker records are never used as a size source.

This script must run where the production database is installed (the cluster).
"""

from __future__ import annotations

import argparse
import bz2
import csv
import hashlib
import pickle
import sys
from pathlib import Path

# Keys under which MetaPhlAn releases have stored the clade taxonomy mapping.
TAXONOMY_KEYS = ("taxonomy", "taxonomy_ids", "clade_taxonomy")
MIN_GENOME_BP = 100_000          # below this a value is not a bacterial genome
MAX_GENOME_BP = 50_000_000       # above this it is not a single prokaryote genome


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_database(path: Path) -> object:
    opener = bz2.open if path.suffix == ".bz2" else open
    with opener(path, "rb") as handle:
        return pickle.load(handle)


def describe(obj: object, limit: int = 20) -> list[str]:
    """Human-readable schema description, printed on failure so a human can adapt."""
    lines = [f"top-level type: {type(obj).__name__}"]
    if isinstance(obj, dict):
        for key in list(obj)[:limit]:
            value = obj[key]
            detail = f"len={len(value)}" if hasattr(value, "__len__") else repr(value)[:60]
            lines.append(f"  key {key!r}: {type(value).__name__} {detail}")
            if isinstance(value, dict) and value:
                sample_key = next(iter(value))
                lines.append(f"    example {sample_key!r} -> {value[sample_key]!r}"[:200])
    return lines


def genome_length(value: object) -> int | None:
    """Accept (taxid, length), [taxid, length], or a bare length."""
    candidate = value
    if isinstance(value, (tuple, list)):
        if not value:
            return None
        candidate = value[-1]
    if isinstance(candidate, bool) or not isinstance(candidate, (int, float)):
        return None
    number = int(candidate)
    return number if MIN_GENOME_BP <= number <= MAX_GENOME_BP else None


def tax_id(value: object) -> str:
    if isinstance(value, (tuple, list)) and len(value) >= 2:
        return str(value[0])
    return ""


def find_taxonomy(db: object) -> dict:
    if not isinstance(db, dict):
        raise SystemExit(
            "[FAIL] database is not a dict; cannot locate a taxonomy mapping\n"
            + "\n".join(describe(db)))
    for key in TAXONOMY_KEYS:
        value = db.get(key)
        if isinstance(value, dict) and value:
            usable = sum(1 for v in value.values() if genome_length(v) is not None)
            if usable:
                print(f"[INFO] using database key {key!r}: "
                      f"{usable}/{len(value)} entries carry a usable genome length")
                return value
    raise SystemExit(
        "[FAIL] no taxonomy mapping with usable genome lengths found.\n"
        f"[FAIL] tried keys: {', '.join(TAXONOMY_KEYS)}\n"
        "[FAIL] observed schema follows; adapt TAXONOMY_KEYS/genome_length "
        "deliberately and record the change in METHODS_DECISION_LOG.md\n"
        + "\n".join(describe(db)))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True,
                        help="MetaPhlAn database pickle, e.g. mpa_vJan25_*.pkl")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--inspect-only", action="store_true",
                        help="Print the database schema and exit without writing.")
    parser.add_argument("--skip-checksum", action="store_true",
                        help="Skip hashing a very large database file.")
    args = parser.parse_args()

    if not args.database.is_file():
        raise SystemExit(f"[FAIL] database not found: {args.database}")

    db = load_database(args.database)
    if args.inspect_only:
        print("\n".join(describe(db)))
        return

    taxonomy = find_taxonomy(db)
    args.outdir.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, str]] = []
    skipped_no_length = 0
    conflicts: dict[str, set[int]] = {}
    for lineage, value in taxonomy.items():
        lineage = str(lineage)
        length = genome_length(value)
        if length is None:
            skipped_no_length += 1
            continue
        conflicts.setdefault(lineage, set()).add(length)
        leaf = lineage.split("|")[-1]
        records.append({
            "clade_lineage": lineage,
            "clade_leaf": leaf,
            "rank": leaf.split("__", 1)[0] if "__" in leaf else "unknown",
            "ncbi_tax_id": tax_id(value),
            "genome_size_bp": str(length),
        })

    # A single lineage resolving to two different lengths is unresolvable here;
    # never pick one silently.
    ambiguous = sorted(k for k, v in conflicts.items() if len(v) > 1)
    if ambiguous:
        raise SystemExit(
            f"[FAIL] {len(ambiguous)} lineage(s) map to conflicting genome lengths, "
            f"first: {ambiguous[0]}")
    if not records:
        raise SystemExit("[FAIL] no clade produced a usable genome length")

    records.sort(key=lambda r: r["clade_lineage"])
    mapping = args.outdir / "metaphlan_genome_sizes.tsv"
    fields = ["clade_lineage", "clade_leaf", "rank", "ncbi_tax_id", "genome_size_bp"]
    with mapping.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)

    by_rank: dict[str, int] = {}
    for record in records:
        by_rank[record["rank"]] = by_rank.get(record["rank"], 0) + 1

    provenance = args.outdir / "genome_size_provenance.tsv"
    lines = [
        "field\tvalue",
        f"database_path\t{args.database.resolve()}",
        f"database_filename\t{args.database.name}",
        f"database_sha256\t{'SKIPPED' if args.skip_checksum else sha256(args.database)}",
        f"database_bytes\t{args.database.stat().st_size}",
        f"extraction_command\t{' '.join(sys.argv)}",
        f"taxonomy_entries\t{len(taxonomy)}",
        f"mapped_clades\t{len(records)}",
        f"skipped_without_usable_length\t{skipped_no_length}",
        f"min_genome_size_bp\t{min(int(r['genome_size_bp']) for r in records)}",
        f"max_genome_size_bp\t{max(int(r['genome_size_bp']) for r in records)}",
    ]
    lines += [f"clades_rank_{rank}\t{count}" for rank, count in sorted(by_rank.items())]
    provenance.write_text("\n".join(lines) + "\n", encoding="utf-8")

    checksum = args.outdir / "genome_sizes.sha256"
    checksum.write_text(
        "".join(f"{sha256(p)}  {p.resolve()}\n" for p in (mapping, provenance)),
        encoding="utf-8")
    (args.outdir / "SUCCESS").write_text(
        f"mapped_clades\t{len(records)}\nstatus\tPASS\n", encoding="utf-8")
    print(f"[PASS] wrote {len(records)} clade genome sizes")
    print(f"[INFO] {mapping}")


if __name__ == "__main__":
    main()
