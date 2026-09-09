#!/usr/bin/env python3
import csv
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True); path.write_text(text, encoding="utf-8")

with tempfile.TemporaryDirectory() as tmp:
    t = Path(tmp); results = t / "results"; baselines = t / "baselines"; out = t / "out"
    manifest = "sample_id\tcondition\tstudy\trun_accessions\nS1\tControl\tStudy1\tERR1;ERR2\n"
    write(t / "feng.tsv", manifest); write(t / "zeller.tsv", "sample_id\tcondition\tstudy\trun_accessions\n")
    write(t / "panel.tsv", "label\ttaxon_name\tassembly\tfasta\tweight\turl\nBfrag\tBacteroides fragilis\tA\tx\t1\t\n")
    write(t / "aliases.csv", "canonical,alias,tool,spike_label\nBacteroides fragilis,Bacteroides fragilis,kraken2_bracken,\nBacteroides fragilis,Bacteroides fragilis,metaphlan4,\n")
    bracken = "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\nBacteroides fragilis\t1\tS\t1\t0\t1\t0.1\n"
    metaphlan = "#mpa_vJan25_CHOCOPhlAnSGB_202503\n#clade_name\tNCBI_tax_id\trelative_abundance\n" + "k__Bacteria|s__Bacteroides_fragilis\t1\t10\n"
    for name in ("ERR1_Bfrag_f0p01", "ERR1_CRCpanel_f0p1"):
        write(results / name / "kraken_bracken" / f"{name}.bracken.S.tsv", bracken)
        write(results / name / "metaphlan4" / f"{name}.metaphlan.tsv", metaphlan)
    write(baselines / "ERR1" / "kraken_bracken" / "ERR1.bracken.S.tsv", bracken)
    write(baselines / "ERR1" / "metaphlan4" / "ERR1.metaphlan.tsv", metaphlan)
    # A historical spike available for only one profiler must not enter the
    # default common-profiler comparison set.
    write(results / "ERR1_Bfrag_f0p005" / "metaphlan4" /
          "ERR1_Bfrag_f0p005.metaphlan.tsv", metaphlan)
    for profiler, suffix, content in (
        ("kraken2_bracken", ".bracken.S.tsv", bracken),
        ("metaphlan4", ".metaphlan.tsv", metaphlan),
    ):
        write(baselines / "ERR2" / profiler / ("ERR2" + suffix), content)
        write(results / "ERR2_Bfrag_f0p01" / profiler /
              ("ERR2_Bfrag_f0p01" + suffix), content)
    subprocess.run([
        "python3", str(ROOT / "analysis_v2/scripts/build_legacy_native_development_input.py"),
        "--results-root", str(results), "--feng-manifest", str(t / "feng.tsv"),
        "--baseline-root", str(baselines),
        "--zeller-manifest", str(t / "zeller.tsv"), "--spike-panel", str(t / "panel.tsv"),
        "--aliases", str(t / "aliases.csv"), "--trajectory-coverage", "available",
        "--outdir", str(out)], check=True)
    with (out / "canonical_input.tsv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 8
    assert {r["analysis_population"] for r in rows} == {"independent", "community"}
    assert {r["cohort"] for r in rows} == {"feng"}
    assert "profiler_coverage=paired" in (out / "DEVELOPMENT_ONLY.txt").read_text()
    assert "incomplete_common_profiler_pair" in (out / "exclusion_ledger.tsv").read_text()
    assert "non_primary_run_for_multirun_sample" in (out / "exclusion_ledger.tsv").read_text()
    assert (out / "DEVELOPMENT_ONLY.txt").is_file()
    assert (out / "validation/SUCCESS").is_file()
print("[PASS] legacy-native development adapter fixture")
