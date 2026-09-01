#!/usr/bin/env python3
import csv
import json
import subprocess
import tempfile
from pathlib import Path


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="target_assembly_audit.") as name:
        root = Path(name)
        panel = root / "panel.tsv"
        reports = root / "provenance"
        reports.mkdir()
        with panel.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["label", "taxon_name", "assembly", "fasta"])
            for label, sequence in (("A", "ACGTNN"), ("B", "GGGG")):
                fasta = root / f"{label}.fa"
                fasta.write_text(f">{label}_1\n{sequence}\n>{label}_2\nACGT\n")
                accession = f"GCF_{label}"
                writer.writerow([label, f"Species {label}", accession, fasta])
                report = {
                    "accession": accession, "currentAccession": accession,
                    "assemblyInfo": {"assemblyStatus": "current", "assemblyLevel": "Complete Genome"},
                    "organism": {"organismName": f"Species {label}", "taxId": 1},
                    "averageNucleotideIdentity": {
                        "submittedSpecies": f"Species {label}", "taxonomyCheckStatus": "OK",
                        "bestAniMatch": {"organismName": f"Species {label}", "ani": 99.9},
                    },
                    "checkmInfo": {"checkmVersion": "test", "completeness": 99, "contamination": 0.1},
                }
                (reports / f"{label}.{accession}.assembly_data_report.jsonl").write_text(
                    json.dumps(report) + "\n"
                )
                (reports / f"{label}.{accession}.dataset_catalog.json").write_text("{}\n")
        out = root / "out"
        subprocess.run([
            "python3", str(repo / "scripts/audit_target_assemblies.py"),
            "--spike-panel", str(panel), "--ncbi-report-dir", str(reports),
            "--outdir", str(out),
        ], check=True)
        rows = list(csv.DictReader((out / "target_assembly_integrity.tsv").open(), delimiter="\t"))
        assert len(rows) == 2
        assert rows[0]["ambiguous_bases"] == "2"
        shared = list(csv.DictReader((out / "exact_cross_target_contigs.tsv").open(), delimiter="\t"))
        assert len(shared) == 1 and shared[0]["targets"] == "A;B"
        ncbi = list(csv.DictReader((out / "target_assembly_ncbi_quality.tsv").open(), delimiter="\t"))
        assert len(ncbi) == 2 and all(row["taxonomy_check_status"] == "OK" for row in ncbi)
        assert (out / "target_assembly_audit_outputs.sha256").stat().st_size > 0
        print("[PASS] Target assembly audit smoke test passed.")


if __name__ == "__main__":
    main()
