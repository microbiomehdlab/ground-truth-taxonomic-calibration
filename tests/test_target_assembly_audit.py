#!/usr/bin/env python3
import csv
import subprocess
import tempfile
from pathlib import Path


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="target_assembly_audit.") as name:
        root = Path(name)
        panel = root / "panel.tsv"
        with panel.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["label", "taxon_name", "assembly", "fasta"])
            for label, sequence in (("A", "ACGTNN"), ("B", "GGGG")):
                fasta = root / f"{label}.fa"
                fasta.write_text(f">{label}_1\n{sequence}\n>{label}_2\nACGT\n")
                writer.writerow([label, f"Species {label}", f"GCF_{label}", fasta])
        out = root / "out"
        subprocess.run([
            "python3", str(repo / "scripts/audit_target_assemblies.py"),
            "--spike-panel", str(panel), "--outdir", str(out),
        ], check=True)
        rows = list(csv.DictReader((out / "target_assembly_integrity.tsv").open(), delimiter="\t"))
        assert len(rows) == 2
        assert rows[0]["ambiguous_bases"] == "2"
        shared = list(csv.DictReader((out / "exact_cross_target_contigs.tsv").open(), delimiter="\t"))
        assert len(shared) == 1 and shared[0]["targets"] == "A;B"
        assert (out / "target_assembly_audit_outputs.sha256").stat().st_size > 0
        print("[PASS] Target assembly audit smoke test passed.")


if __name__ == "__main__":
    main()
