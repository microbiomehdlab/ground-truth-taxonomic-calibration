#!/usr/bin/env python3
"""Fixture tests for native-profile semantic auditing."""

import csv
import subprocess
import tempfile
from pathlib import Path


def main() -> None:
    repo = Path(__file__).resolve().parents[2]
    script = repo / "analysis_v2/scripts/audit_profiler_semantics.py"
    with tempfile.TemporaryDirectory(prefix="profiler_semantics.") as name:
        root = Path(name)
        bracken = root / "sample.bracken.S.tsv"
        bracken.write_text(
            "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\t"
            "added_reads\tnew_est_reads\tfraction_total_reads\n"
            "Target species\t1\tS\t60\t20\t80\t0.80\n"
            "Other species\t2\tS\t5\t5\t10\t0.10\n",
            encoding="utf-8",
        )
        metaphlan = root / "sample.metaphlan.tsv"
        metaphlan.write_text(
            "#mpa_vJan25_CHOCOPhlAnSGB_202503\n"
            "#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n"
            "k__Bacteria|p__P|g__G|s__Target_species\t1\t72.5\t\n"
            "k__Bacteria|p__P|g__G\t2\t17.5\t\n"
            "UNCLASSIFIED\t-1\t10.0\t\n",
            encoding="utf-8",
        )
        outdir = root / "audit"
        subprocess.run([
            "python3", str(script), "--bracken", str(bracken),
            "--metaphlan", str(metaphlan), "--outdir", str(outdir),
        ], check=True)
        with (outdir / "profile_semantics_audit.tsv").open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == 2
        assert rows[0]["species_total"] == "0.9"
        assert rows[0]["unclassified_total"] == ""
        assert rows[1]["reported_total"] == "100.0"
        assert rows[1]["species_total"] == "72.5"
        assert rows[1]["non_species_total"] == "17.5"
        assert rows[1]["unclassified_total"] == "10.0"
        assert (outdir / "SUCCESS").is_file()

    print("[PASS] profiler-semantics fixture test")


if __name__ == "__main__":
    main()
