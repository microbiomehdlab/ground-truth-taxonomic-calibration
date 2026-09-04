#!/usr/bin/env python3
"""Fixture test for expansion of native profiles into biomarker input."""

import csv
import subprocess
import tempfile
from pathlib import Path


def write(path, header, data):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header); writer.writerows(data)


def main():
    repo = Path(__file__).resolve().parents[2]
    with tempfile.TemporaryDirectory(prefix="biomarker_input.") as name:
        root = Path(name); canonical = root / "canonical.tsv"; outdir = root / "out"
        bracken = root / "S1.bracken.S.tsv"
        write(bracken, ["name", "taxonomy_lvl", "fraction_total_reads"],
              [["Peptostreptococcus anaerobius", "S", "0.02"], ["Other", "S", "0.98"]])
        header = ["cohort", "study", "analysis_population", "sample_id", "condition",
                  "target_label", "assembly_arm", "profiler", "profile_id",
                  "baseline_profile_id", "spike_fraction_target", "source_profile",
                  "target_taxon", "include", "exclusion_reason"]
        base = ["yachida", "Study1", "independent", "S1", "CRC", "Pana", "original",
                "kraken2_bracken"]
        write(canonical, header, [base + ["S1", "S1", "0", str(bracken),
                                          "Peptostreptococcus anaerobius", "1", ""],
                                  base + ["S1_Pana_f0p01", "S1", "0.01", str(bracken),
                                          "Peptostreptococcus anaerobius", "1", ""]])
        result = subprocess.run(["python3", str(repo / "analysis_v2/scripts/build_biomarker_abundance_input.py"),
                                 "--canonical", str(canonical), "--aliases",
                                 str(repo / "examples/spike_taxon_aliases.csv"),
                                 "--outdir", str(outdir)], text=True, capture_output=True)
        assert result.returncode == 0, result.stderr
        with (outdir / "biomarker_profile_manifest.tsv").open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert [row["dose_level"] for row in rows] == ["baseline", "dose_01"]
        assert (outdir / "SUCCESS").is_file()
    print("[PASS] biomarker abundance-input fixture")


if __name__ == "__main__":
    main()
