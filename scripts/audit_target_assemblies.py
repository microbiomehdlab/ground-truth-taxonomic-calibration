#!/usr/bin/env python3
"""Produce deterministic sequence-integrity metrics for spike target FASTAs.

This is an integrity and identity-support audit, not a taxonomic contamination
classifier. Taxonomic representation is handled by the separate UHGG/MetaPhlAn
audit and empirical assignment by the pure-pool audit.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


VALID = set("ACGTRYSWKMBDHVN")


def fail(message: str) -> "NoReturn":
    raise SystemExit(f"[ERROR] {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve(panel: Path, value: str) -> Path:
    path = Path(value)
    return (path if path.is_absolute() else panel.parent / path).resolve()


def fasta_records(path: Path):
    name = None
    pieces: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(pieces).upper()
                name = line[1:].strip()
                if not name:
                    fail(f"empty FASTA header in {path}:{number}")
                pieces = []
            else:
                if name is None:
                    fail(f"sequence before first FASTA header in {path}:{number}")
                pieces.append(re.sub(r"\s+", "", line))
    if name is not None:
        yield name, "".join(pieces).upper()


def n50(lengths: list[int]) -> int:
    threshold = sum(lengths) / 2
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= threshold:
            return length
    return 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spike-panel", required=True, type=Path)
    parser.add_argument("--ncbi-report-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    with args.spike_panel.open(newline="", encoding="utf-8-sig") as handle:
        panel = list(csv.DictReader(handle, delimiter="\t"))
    required = {"label", "taxon_name", "assembly", "fasta"}
    if not panel or not required.issubset(panel[0]):
        fail("spike panel is empty or lacks label, taxon_name, assembly, or fasta")
    labels = [row["label"].strip() for row in panel]
    if len(labels) != len(set(labels)):
        fail("spike panel contains duplicate labels")

    rows = []
    ncbi_rows = []
    ncbi_files: list[Path] = []
    all_contig_hashes: dict[str, list[tuple[str, str]]] = {}
    for target in panel:
        path = resolve(args.spike_panel, target["fasta"])
        if not path.is_file() or path.stat().st_size == 0:
            fail(f"missing or empty FASTA for {target['label']}: {path}")
        records = list(fasta_records(path))
        if not records or any(not sequence for _, sequence in records):
            fail(f"empty FASTA record for {target['label']}: {path}")
        headers = [name for name, _ in records]
        sequences = [sequence for _, sequence in records]
        invalid = Counter(base for sequence in sequences for base in sequence if base not in VALID)
        if invalid:
            fail(f"invalid FASTA symbols for {target['label']}: {dict(invalid)}")
        lengths = [len(sequence) for sequence in sequences]
        bases = sum(lengths)
        gc = sum(sequence.count("G") + sequence.count("C") for sequence in sequences)
        ambiguous = sum(sum(base not in "ACGT" for base in sequence) for sequence in sequences)
        hashes = [hashlib.sha256(sequence.encode()).hexdigest() for sequence in sequences]
        duplicates = len(hashes) - len(set(hashes))
        for header, digest in zip(headers, hashes):
            all_contig_hashes.setdefault(digest, []).append((target["label"], header))
        normalized = hashlib.sha256("\n".join(sorted(hashes)).encode()).hexdigest()
        rows.append({
            "label": target["label"], "canonical_target": target["taxon_name"],
            "assembly": target["assembly"], "fasta": str(path),
            "fasta_bytes": path.stat().st_size, "fasta_sha256": sha256(path),
            "normalized_sequence_sha256": normalized, "contigs": len(records),
            "total_bases": bases, "n50": n50(lengths),
            "longest_contig": max(lengths), "gc_pct": f"{100 * gc / bases:.4f}",
            "ambiguous_bases": ambiguous,
            "ambiguous_pct": f"{100 * ambiguous / bases:.6f}",
            "duplicate_contigs_within_assembly": duplicates,
        })

        reports = sorted(args.ncbi_report_dir.glob(
            f"{target['label']}.{target['assembly']}.assembly_data_report.jsonl"
        ))
        if len(reports) != 1:
            fail(
                f"expected exactly one NCBI assembly report for {target['label']} "
                f"{target['assembly']}; found {len(reports)} in {args.ncbi_report_dir}"
            )
        report_path = reports[0]
        with report_path.open(encoding="utf-8") as handle:
            report_lines = [line for line in handle if line.strip()]
        if len(report_lines) != 1:
            fail(f"expected one JSON record in {report_path}; found {len(report_lines)}")
        report = json.loads(report_lines[0])
        if report.get("accession") != target["assembly"]:
            fail(
                f"NCBI accession mismatch for {target['label']}: panel="
                f"{target['assembly']} report={report.get('accession')}"
            )
        assembly_info = report.get("assemblyInfo") or {}
        organism = report.get("organism") or {}
        ani = report.get("averageNucleotideIdentity") or {}
        best_ani = ani.get("bestAniMatch") or {}
        checkm = report.get("checkmInfo") or {}
        ncbi_files.append(report_path)
        catalog_path = args.ncbi_report_dir / (
            f"{target['label']}.{target['assembly']}.dataset_catalog.json"
        )
        if not catalog_path.is_file() or catalog_path.stat().st_size == 0:
            fail(f"missing NCBI dataset catalog: {catalog_path}")
        ncbi_files.append(catalog_path)
        ncbi_rows.append({
            "label": target["label"],
            "canonical_target": target["taxon_name"],
            "panel_accession": target["assembly"],
            "report_accession": report.get("accession", ""),
            "current_accession": report.get("currentAccession", ""),
            "assembly_status": assembly_info.get("assemblyStatus", ""),
            "assembly_level": assembly_info.get("assemblyLevel", ""),
            "refseq_category": assembly_info.get("refseqCategory", ""),
            "assembly_release_date": assembly_info.get("releaseDate", ""),
            "organism_name": organism.get("organismName", ""),
            "organism_taxid": organism.get("taxId", ""),
            "submitted_species": ani.get("submittedSpecies", ""),
            "taxonomy_check_status": ani.get("taxonomyCheckStatus", ""),
            "best_ani_match_species": best_ani.get("organismName", ""),
            "best_ani_match_assembly": best_ani.get("assembly", ""),
            "best_ani_pct": best_ani.get("ani", ""),
            "best_ani_assembly_coverage_pct": best_ani.get("assemblyCoverage", ""),
            "checkm_version": checkm.get("checkmVersion", ""),
            "checkm_marker_set": checkm.get("checkmMarkerSet", ""),
            "checkm_completeness_pct": checkm.get("completeness", ""),
            "checkm_contamination_pct": checkm.get("contamination", ""),
            "report_sha256": sha256(report_path),
            "catalog_sha256": sha256(catalog_path),
        })

    cross = []
    for digest, occurrences in sorted(all_contig_hashes.items()):
        target_labels = sorted({label for label, _ in occurrences})
        if len(target_labels) > 1:
            cross.append({
                "contig_sequence_sha256": digest,
                "target_count": len(target_labels),
                "targets": ";".join(target_labels),
                "occurrences": json.dumps(occurrences, separators=(",", ":")),
            })

    args.outdir.mkdir(parents=True, exist_ok=True)
    summary = args.outdir / "target_assembly_integrity.tsv"
    with summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    shared = args.outdir / "exact_cross_target_contigs.tsv"
    fields = ["contig_sequence_sha256", "target_count", "targets", "occurrences"]
    with shared.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(cross)
    ncbi = args.outdir / "target_assembly_ncbi_quality.tsv"
    with ncbi.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(ncbi_rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader(); writer.writerows(ncbi_rows)
    provenance = args.outdir / "target_assembly_integrity_PROVENANCE.txt"
    provenance.write_text(
        "\n".join([
            "Target assembly sequence-integrity audit",
            "scope=FASTA structure, sequence checksums, assembly summaries, exact duplicate contigs",
            "limitation=not a taxonomic contamination classifier",
            f"spike_panel={args.spike_panel.resolve()}",
            f"spike_panel_sha256={sha256(args.spike_panel)}",
            f"ncbi_report_dir={args.ncbi_report_dir.resolve()}",
            f"script_sha256={sha256(Path(__file__).resolve())}",
            f"targets={len(rows)}", f"exact_cross_target_contigs={len(cross)}",
            f"ncbi_provenance_files={len(ncbi_files)}",
        ]) + "\n", encoding="utf-8",
    )
    checksums = args.outdir / "target_assembly_audit_outputs.sha256"
    checksums.write_text("".join(
        f"{sha256(path)}  {path.name}\n" for path in (summary, shared, ncbi, provenance)
    ), encoding="utf-8")
    print(f"[PASS] Audited {len(rows)} target assemblies")
    print(f"[INFO] Exact contigs shared across targets: {len(cross)}")
    print(f"[INFO] Output: {args.outdir}")


if __name__ == "__main__":
    main()
