#!/usr/bin/env python3
"""Generate the ten-target genome-size table from the frozen spike panel.

Lengths are measured directly from the exact FASTA files referenced by
`spikes/spike_panel.tsv` -- the assemblies actually used to simulate the
implanted reads. Remembered or generic species lengths are never used when the
FASTAs are available, and there is no fallback that would silently substitute
them.

Outputs `target_genome_sizes.tsv` with provenance, checksums and SUCCESS.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

EXPECTED_LABELS = ["Bfrag", "Csym", "Dpne", "Fnuc", "Hhat",
                   "Pmic", "Pana", "Psto", "Porp", "Pint"]
FIELDS = ["target_label", "target_feature", "assembly_accession", "fasta_path",
          "fasta_sha256", "genome_size_bp"]
MIN_GENOME_BP = 100_000
MAX_GENOME_BP = 50_000_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta_length(path: Path) -> int:
    """Sum non-header sequence characters, ignoring whitespace."""
    total = 0
    saw_header = False
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                saw_header = True
                continue
            total += len(line.strip())
    if not saw_header:
        raise SystemExit(f"[FAIL] {path}: no FASTA header found")
    return total


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spike-panel", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--fasta-root", type=Path,
                        help="Base directory for relative FASTA paths in the panel.")
    parser.add_argument("--expected-labels", default=",".join(EXPECTED_LABELS),
                        help="Comma-separated labels that must be present exactly once.")
    args = parser.parse_args()

    expected = [x.strip() for x in args.expected_labels.split(",") if x.strip()]
    if not args.spike_panel.is_file() or args.spike_panel.stat().st_size == 0:
        raise SystemExit(f"[FAIL] missing or empty spike panel: {args.spike_panel}")

    rows = []
    with args.spike_panel.open(newline="", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for required in ("label", "taxon_name", "assembly", "fasta"):
            if required not in header:
                raise SystemExit(f"[FAIL] spike panel missing column {required}")
        index = {name: i for i, name in enumerate(header)}
        for line_no, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            rows.append({
                "line": line_no,
                "label": parts[index["label"]].strip(),
                "taxon_name": parts[index["taxon_name"]].strip(),
                "assembly": parts[index["assembly"]].strip(),
                "fasta": parts[index["fasta"]].strip(),
            })

    seen: dict[str, int] = {}
    for entry in rows:
        if entry["label"] in seen:
            raise SystemExit(
                f"[FAIL] duplicate target label {entry['label']!r} at lines "
                f"{seen[entry['label']]} and {entry['line']}")
        seen[entry["label"]] = entry["line"]

    unexpected = sorted(set(seen) - set(expected))
    missing = sorted(set(expected) - set(seen))
    if unexpected or missing:
        raise SystemExit(
            "[FAIL] spike-panel labels do not match the expected set; "
            f"unexpected={unexpected or 'none'} missing={missing or 'none'}")

    records = []
    for entry in rows:
        fasta = Path(entry["fasta"])
        if args.fasta_root and not fasta.is_absolute():
            fasta = args.fasta_root / fasta
        if not fasta.is_file():
            raise SystemExit(
                f"[FAIL] FASTA for {entry['label']} not found: {fasta}. Genome "
                "sizes are never substituted from memory.")
        if fasta.stat().st_size == 0:
            raise SystemExit(f"[FAIL] FASTA for {entry['label']} is empty: {fasta}")
        length = fasta_length(fasta)
        if not MIN_GENOME_BP <= length <= MAX_GENOME_BP:
            raise SystemExit(
                f"[FAIL] {entry['label']}: measured length {length} outside "
                f"[{MIN_GENOME_BP}, {MAX_GENOME_BP}]")
        records.append({
            "target_label": entry["label"],
            "target_feature": entry["taxon_name"],
            "assembly_accession": entry["assembly"],
            "fasta_path": str(fasta.resolve()),
            "fasta_sha256": sha256(fasta),
            "genome_size_bp": str(length),
        })

    records.sort(key=lambda r: r["target_label"])
    args.outdir.mkdir(parents=True, exist_ok=True)
    table = args.outdir / "target_genome_sizes.tsv"
    with table.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(FIELDS) + "\n")
        for record in records:
            handle.write("\t".join(record[f] for f in FIELDS) + "\n")

    provenance = args.outdir / "target_genome_sizes_provenance.tsv"
    provenance.write_text("\n".join([
        "field\tvalue",
        f"spike_panel\t{args.spike_panel.resolve()}",
        f"spike_panel_sha256\t{sha256(args.spike_panel)}",
        f"targets\t{len(records)}",
        f"expected_labels\t{','.join(expected)}",
        f"command\t{' '.join(sys.argv)}",
        f"min_genome_size_bp\t{min(int(r['genome_size_bp']) for r in records)}",
        f"max_genome_size_bp\t{max(int(r['genome_size_bp']) for r in records)}",
        "status\tPASS",
    ]) + "\n", encoding="utf-8")

    (args.outdir / "target_genome_sizes.sha256").write_text(
        "".join(f"{sha256(p)}  {p.resolve()}\n" for p in (table, provenance)),
        encoding="utf-8")
    (args.outdir / "SUCCESS").write_text(
        f"targets\t{len(records)}\nstatus\tPASS\n", encoding="utf-8")
    print(f"[PASS] measured {len(records)} target genome size(s)")
    print(f"[INFO] {table}")


if __name__ == "__main__":
    main()
