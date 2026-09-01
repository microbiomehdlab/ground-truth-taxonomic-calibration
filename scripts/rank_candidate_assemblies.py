#!/usr/bin/env python3
"""Rank same-species NCBI assemblies using prespecified sequence-quality rules.

This selection uses only NCBI identity/assembly metadata. It must be run before
profiling replacement pools so profiler outcomes cannot influence selection.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


def fail(message: str) -> "NoReturn":
    raise SystemExit(f"[ERROR] {message}")


def norm_species(value: object) -> str:
    text = str(value or "").lower().replace("[", "").replace("]", "")
    text = re.sub(r"[^a-z0-9]+", " ", text).strip()
    return " ".join(text.split()[:2])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def numeric(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def pick(mapping: dict, camel: str, snake: str, default=""):
    return mapping.get(camel, mapping.get(snake, default))


def enum_is(value: object, expected: str) -> bool:
    normalized = re.sub(r"[^a-z0-9]+", "_", str(value or "").lower()).strip("_")
    wanted = re.sub(r"[^a-z0-9]+", "_", expected.lower()).strip("_")
    return normalized == wanted or normalized.endswith("_" + wanted)


def load_jsonl(path: Path) -> list[dict]:
    records = []
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            value = json.loads(line)
            if "reports" in value and isinstance(value["reports"], list):
                records.extend(value["reports"])
            elif isinstance(value, dict):
                records.append(value)
            else:
                fail(f"unexpected JSON value in {path}:{number}")
    if not records:
        fail(f"candidate snapshot contains no records: {path}")
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-manifest", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--minimum-completeness", type=float, default=95.0)
    parser.add_argument("--maximum-contamination", type=float, default=5.0)
    parser.add_argument("--minimum-ani", type=float, default=95.0)
    parser.add_argument("--minimum-ani-coverage", type=float, default=90.0)
    args = parser.parse_args()

    with args.candidate_manifest.open(newline="", encoding="utf-8-sig") as handle:
        specs = list(csv.DictReader(handle, delimiter="\t"))
    required = {"label", "canonical_species", "current_accession", "jsonl"}
    if not specs or not required.issubset(specs[0]):
        fail("candidate manifest lacks label, canonical_species, current_accession, or jsonl")

    level_order = {"complete genome": 0, "chromosome": 1, "scaffold": 2, "contig": 3}
    category_order = {"reference genome": 0, "representative genome": 1, "": 2}
    all_rows = []
    selected = []
    no_candidate_labels = []
    snapshot_hashes = []
    for spec in specs:
        snapshot = Path(spec["jsonl"])
        if not snapshot.is_absolute():
            snapshot = (args.candidate_manifest.parent / snapshot).resolve()
        if not snapshot.is_file() or snapshot.stat().st_size == 0:
            fail(f"missing candidate snapshot: {snapshot}")
        snapshot_hashes.append((snapshot, sha256(snapshot)))
        expected = norm_species(spec["canonical_species"])
        cohort = []
        accessions = set()
        for record in load_jsonl(snapshot):
            accession = str(record.get("accession", ""))
            if not accession or accession in accessions:
                continue
            accessions.add(accession)
            info = record.get("assemblyInfo") or record.get("assembly_info") or {}
            organism = record.get("organism") or {}
            ani = record.get("averageNucleotideIdentity") or record.get("average_nucleotide_identity") or {}
            best = ani.get("bestAniMatch") or ani.get("best_ani_match") or {}
            checkm = record.get("checkmInfo") or record.get("checkm_info") or {}
            completeness = numeric(checkm.get("completeness"))
            contamination = numeric(checkm.get("contamination"))
            best_ani = numeric(best.get("ani"))
            ani_coverage = numeric(pick(best, "assemblyCoverage", "assembly_coverage"))
            submitted = norm_species(
                pick(ani, "submittedSpecies", "submitted_species")
                or pick(organism, "organismName", "organism_name")
            )
            best_species = norm_species(pick(best, "organismName", "organism_name"))
            reasons = []
            assembly_status = pick(info, "assemblyStatus", "assembly_status")
            taxonomy_status = pick(ani, "taxonomyCheckStatus", "taxonomy_check_status")
            if not enum_is(assembly_status, "current"): reasons.append("not_current")
            if submitted != expected: reasons.append("submitted_species_mismatch")
            if not enum_is(taxonomy_status, "ok"): reasons.append("taxonomy_not_ok")
            if best_species != expected: reasons.append("ani_species_mismatch")
            if completeness is None or completeness < args.minimum_completeness:
                reasons.append("completeness_below_threshold_or_missing")
            if contamination is None or contamination > args.maximum_contamination:
                reasons.append("contamination_above_threshold_or_missing")
            if best_ani is None or best_ani < args.minimum_ani:
                reasons.append("ani_below_threshold_or_missing")
            if ani_coverage is None or ani_coverage < args.minimum_ani_coverage:
                reasons.append("ani_coverage_below_threshold_or_missing")
            stats = record.get("assemblyStats") or record.get("assembly_stats") or {}
            type_material = bool(record.get("typeMaterial") or record.get("type_material"))
            assembly_level = pick(info, "assemblyLevel", "assembly_level")
            refseq_category = pick(info, "refseqCategory", "refseq_category")
            row = {
                "label": spec["label"], "canonical_species": spec["canonical_species"],
                "accession": accession, "is_original": "yes" if accession == spec["current_accession"] else "no",
                "eligible": "yes" if not reasons else "no", "exclusion_reasons": ";".join(reasons),
                "organism_name": pick(organism, "organismName", "organism_name"),
                "taxonomy_check_status": taxonomy_status,
                "best_ani_species": pick(best, "organismName", "organism_name"), "best_ani_pct": best_ani if best_ani is not None else "",
                "best_ani_coverage_pct": ani_coverage if ani_coverage is not None else "",
                "assembly_level": assembly_level, "refseq_category": refseq_category,
                "from_type_material": "yes" if type_material else "no",
                "checkm_marker_set": pick(checkm, "checkmMarkerSet", "checkm_marker_set"),
                "checkm_completeness_pct": completeness if completeness is not None else "",
                "checkm_contamination_pct": contamination if contamination is not None else "",
                "contig_n50": pick(stats, "contigN50", "contig_n50"),
                "contigs": pick(stats, "numberOfContigs", "number_of_contigs"),
                "total_sequence_length": pick(stats, "totalSequenceLength", "total_sequence_length"),
                "selection_rank": "",
            }
            row["_sort"] = (
                category_order.get(str(refseq_category).lower(), 2),
                0 if type_material else 1,
                level_order.get(str(assembly_level).lower(), 4),
                contamination if contamination is not None else float("inf"),
                -(completeness if completeness is not None else -1),
                -(ani_coverage if ani_coverage is not None else -1),
                -numeric(pick(stats, "contigN50", "contig_n50"), 0), accession,
            )
            cohort.append(row)
        eligible = sorted((row for row in cohort if row["eligible"] == "yes"), key=lambda row: row["_sort"])
        for rank, row in enumerate(eligible, 1):
            row["selection_rank"] = rank
        if eligible:
            selected.append(dict(eligible[0]))
        else:
            no_candidate_labels.append(spec["label"])
        all_rows.extend(cohort)

    for row in all_rows + selected:
        row.pop("_sort", None)
    args.outdir.mkdir(parents=True, exist_ok=True)
    fields = list(all_rows[0])
    for name, rows in (("candidate_assembly_audit.tsv", all_rows), ("selected_candidate_assemblies.tsv", selected)):
        with (args.outdir / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
    provenance = args.outdir / "candidate_assembly_selection_PROVENANCE.txt"
    provenance.write_text("\n".join([
        "Outcome-independent NCBI candidate-assembly selection",
        f"minimum_completeness={args.minimum_completeness}",
        f"maximum_contamination={args.maximum_contamination}",
        f"minimum_ani={args.minimum_ani}",
        f"minimum_ani_coverage={args.minimum_ani_coverage}",
        "ranking=reference category; type material; assembly level; contamination; completeness; ANI coverage; contig N50; accession",
        f"targets_without_eligible_candidate={';'.join(no_candidate_labels)}",
        f"manifest_sha256={sha256(args.candidate_manifest)}",
        f"script_sha256={sha256(Path(__file__).resolve())}",
        *[f"snapshot_sha256={digest}  {path}" for path, digest in snapshot_hashes],
    ]) + "\n", encoding="utf-8")
    outputs = [args.outdir / "candidate_assembly_audit.tsv", args.outdir / "selected_candidate_assemblies.tsv", provenance]
    (args.outdir / "candidate_assembly_selection.sha256").write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in outputs), encoding="utf-8"
    )
    print(f"[PASS] Ranked {len(all_rows)} assemblies across {len(specs)} targets")
    for row in selected:
        print(f"[SELECTED] {row['label']} {row['accession']}")
    for label in no_candidate_labels:
        print(f"[REVIEW] {label}: no assembly passed every frozen eligibility rule")


if __name__ == "__main__":
    main()
