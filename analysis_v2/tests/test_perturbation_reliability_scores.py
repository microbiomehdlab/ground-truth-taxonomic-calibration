#!/usr/bin/env python3
import csv
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
header = ["cohort", "analysis_population", "profiler", "contrast", "feature",
          "target_label", "assembly_arm", "spike_fraction_target", "q_threshold",
          "feature_role", "baseline_called", "dose_called", "baseline_effect",
          "dose_effect", "effect_sign_changed"]
rows = []
for cohort, dose_called, dose_effect in (("feng", 1, 1.8), ("zeller", 1, 1.9)):
    for target in ("A", "B"):
        rows.append([cohort, "community", "metaphlan4", "CRC_vs_Control", "marker",
                     target, "original", "0.01", "0.05", "bystander", 1,
                     dose_called, 2.0, dose_effect, 0])
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp); ledger = tmp / "ledger.tsv"; out = tmp / "out"
    with ledger.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header); writer.writerows(rows)
    subprocess.run(["python3", str(ROOT / "analysis_v2/scripts/build_perturbation_reliability_scores.py"),
                    "--transition-ledger", str(ledger), "--outdir", str(out)], check=True)
    with (out / "taxon_reliability_scores.tsv").open(newline="", encoding="utf-8") as handle:
        scores = list(csv.DictReader(handle, delimiter="\t"))
    assert len(scores) == 2
    assert all(row["external_replication_status"] == "DIRECTIONALLY_REPLICATED" for row in scores)
    assert all(float(row["perturbation_reliability_score"]) > 0.95 for row in scores)
    assert (out / "SUCCESS").is_file()
print("[PASS] perturbation-reliability score fixture")
