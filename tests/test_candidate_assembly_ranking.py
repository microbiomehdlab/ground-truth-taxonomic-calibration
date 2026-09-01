#!/usr/bin/env python3
import csv, json, subprocess, tempfile
from pathlib import Path


def record(accession, contamination, category="", level="Contig"):
    return {
        "accession": accession,
        "assemblyInfo": {"assemblyStatus": "current", "assemblyLevel": level, "refseqCategory": category},
        "organism": {"organismName": "Species alpha"},
        "averageNucleotideIdentity": {"submittedSpecies": "Species alpha", "taxonomyCheckStatus": "OK", "bestAniMatch": {"organismName": "Species alpha", "ani": 99, "assemblyCoverage": 99}},
        "checkmInfo": {"completeness": 99, "contamination": contamination},
        "assemblyStats": {"contigN50": 1000, "numberOfContigs": 10, "totalSequenceLength": 10000},
    }


def snake_record(accession):
    return {
        "accession": accession,
        "assembly_info": {"assembly_status": "ASSEMBLY_STATUS_CURRENT", "assembly_level": "Complete Genome", "refseq_category": "reference genome"},
        "organism": {"organism_name": "Species alpha"},
        "average_nucleotide_identity": {"submitted_species": "Species alpha", "taxonomy_check_status": "TAXONOMY_CHECK_STATUS_OK", "best_ani_match": {"organism_name": "Species alpha", "ani": 99, "assembly_coverage": 99}},
        "checkm_info": {"completeness": 99, "contamination": 0.2, "checkm_marker_set": "Species alpha"},
        "assembly_stats": {"contig_n50": 2000, "number_of_contigs": 2, "total_sequence_length": 10000},
    }


def main():
    repo = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory() as name:
        root = Path(name); snapshot = root / "candidates.jsonl"
        snapshot.write_text("\n".join(json.dumps(x) for x in [record("GCF_BAD", 12), record("GCF_GOOD", 1, "reference genome")]) + "\n")
        manifest = root / "manifest.tsv"
        with manifest.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t"); writer.writerow(["label", "canonical_species", "current_accession", "jsonl"]); writer.writerow(["X", "Species alpha", "GCF_BAD", snapshot])
        out = root / "out"
        subprocess.run(["python3", str(repo / "scripts/rank_candidate_assemblies.py"), "--candidate-manifest", str(manifest), "--outdir", str(out)], check=True)
        chosen = list(csv.DictReader((out / "selected_candidate_assemblies.tsv").open(), delimiter="\t"))
        assert chosen[0]["accession"] == "GCF_GOOD"
        audit = list(csv.DictReader((out / "candidate_assembly_audit.tsv").open(), delimiter="\t"))
        assert next(row for row in audit if row["accession"] == "GCF_BAD")["eligible"] == "no"
        # A target with no eligible candidate is a review result, not a broken
        # workflow; the complete exclusion ledger must still be retained.
        snapshot.write_text(json.dumps(record("GCF_ONLY_BAD", 12)) + "\n")
        subprocess.run(["python3", str(repo / "scripts/rank_candidate_assemblies.py"), "--candidate-manifest", str(manifest), "--outdir", str(root / "none")], check=True)
        empty = list(csv.DictReader((root / "none/selected_candidate_assemblies.tsv").open(), delimiter="\t"))
        assert empty == []
        assert (root / "none/candidate_assembly_audit.tsv").stat().st_size > 0
        snapshot.write_text(json.dumps(snake_record("GCF_SNAKE")) + "\n")
        subprocess.run(["python3", str(repo / "scripts/rank_candidate_assemblies.py"), "--candidate-manifest", str(manifest), "--outdir", str(root / "snake")], check=True)
        snake = list(csv.DictReader((root / "snake/selected_candidate_assemblies.tsv").open(), delimiter="\t"))
        assert snake[0]["accession"] == "GCF_SNAKE"
        print("[PASS] Candidate assembly ranking test passed.")


if __name__ == "__main__": main()
