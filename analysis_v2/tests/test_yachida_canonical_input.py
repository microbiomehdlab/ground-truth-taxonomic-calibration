#!/usr/bin/env python3
"""Integration fixture for the complete Yachida canonical builder."""
import csv
import subprocess
import sys
import tempfile
from pathlib import Path


repo = Path(__file__).resolve().parents[2]


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    results = root / "results"; sample_root = results / "Study" / "S1"
    manifest = root / "manifest.tsv"
    independent = root / "independent.tsv"
    panel = root / "panel.tsv"; aliases = root / "aliases.csv"
    write(manifest, "sample_id\tStudy\tTarget_Condition\nS1\tStudy\tControl\n")
    write(independent, "sample_id\nS1\n")
    write(panel, "label\ttaxon_name\tweight\nTarget\tTarget species\t1\n")
    write(aliases,
          "canonical,alias,tool\n"
          "Target species,Target species,kraken2_bracken\n"
          "Target species,Target species,metaphlan4\n")
    bracken_header = "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n"
    write(sample_root / "profiles/baseline/S1/S1.bracken.S.tsv",
          bracken_header + "Target species\t1\tS\t1\t0\t1\t0.01\n")
    write(sample_root / "profiles/baseline/S1/S1.metaphlan.tsv",
          "#fixture\nk__Bacteria|s__Target_species\t1\t2.0\n")
    write(sample_root / "spike_design/community/CRCpanel.tsv",
          "sample_id\tcommunity\tfraction\tR1\tR2\tR\tN_total\tf_hat\tmember_seeds\tseed_scheme\n"
          "S1\tCRCpanel\t0.1\t100\t100\t100\t11\t0.0990990990990991\tx\tv1\n")
    write(sample_root / "profiles/community/S1_CRCpanel_f0p1/S1_CRCpanel_f0p1.bracken.S.tsv",
          bracken_header + "Target species\t1\tS\t10\t0\t10\t0.1\n")
    write(sample_root / "profiles/community/S1_CRCpanel_f0p1/S1_CRCpanel_f0p1.metaphlan.tsv",
          "#fixture\nk__Bacteria|s__Target_species\t1\t20.0\n")
    write(sample_root / "spike_design/independent/Target.tsv",
          "sample_id\tlabel\tfraction\tN_inserted\tf_hat\n"
          "S1\tTarget\t0.1\t11\t0.0990990990990991\n")
    write(sample_root / "profiles/independent/Target/S1_Target_f0p1/S1_Target_f0p1.bracken.S.tsv",
          bracken_header + "Target species\t1\tS\t10\t0\t10\t0.1\n")
    write(sample_root / "profiles/independent/Target/S1_Target_f0p1/S1_Target_f0p1.metaphlan.tsv",
          "#fixture\nk__Bacteria|s__Target_species\t1\t20.0\n")
    out = root / "out"
    command = [sys.executable, str(repo / "analysis_v2/scripts/build_yachida_canonical_input.py"),
               "--manifest", str(manifest), "--independent-manifest", str(independent),
               "--results-root", str(results), "--spike-panel", str(panel),
               "--aliases", str(aliases), "--outdir", str(out),
               "--expected-samples", "1", "--expected-independent", "1",
               "--expected-independent-doses", "1", "--expected-community-doses", "1"]
    subprocess.run(command, check=True)
    with (out / "canonical_input.tsv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 8
    positives = [row for row in rows if float(row["spike_fraction_target"]) > 0]
    assert len(positives) == 4
    assert {row["analysis_population"] for row in positives} == {"community", "independent"}
    assert all(row["implanted_read_pairs_target"] == "11" for row in positives)
    assert (out / "validation/SUCCESS").is_file()
print("[PASS] complete Yachida canonical-input fixture")
