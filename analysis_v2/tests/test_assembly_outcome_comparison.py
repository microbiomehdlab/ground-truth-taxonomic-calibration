#!/usr/bin/env python3
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

repo = Path(__file__).resolve().parents[2]


def table(path: Path, fields: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields); writer.writerows(rows)


with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp); canonical = root / "canonical.tsv"; metrics = root / "metrics.tsv"
    cfields = ["sample_id", "target_label", "assembly_arm", "profiler", "analysis_population",
               "spike_fraction_target", "detected_native_nonzero", "include"]
    crows = []
    for sample, original, clean in (("S1", 0, 1), ("S2", 1, 1)):
        for arm, detected in (("original", original), ("clean", clean)):
            crows.append([sample, "Pana", arm, "metaphlan4", "independent", .001, detected, 1])
    table(canonical, cfields, crows)
    mfields = ["target_label", "assembly_arm", "profiler", "contrast", "spike_fraction_target",
               "q_threshold", "target_called", "target_effect", "off_target_enriched_calls", "precision"]
    table(metrics, mfields, [
        ["Pana", "original", "metaphlan4", "pooled", .001, .05, 0, 1.0, 2, 0],
        ["Pana", "clean", "metaphlan4", "pooled", .001, .05, 1, 1.5, 1, .5],
    ])
    out = root / "out"
    subprocess.run([sys.executable, str(repo / "analysis_v2/scripts/compare_assembly_sensitivity_outcomes.py"),
                    "--canonical", str(canonical), "--biomarker-metrics", str(metrics),
                    "--outdir", str(out)], check=True)
    detection = list(csv.DictReader((out / "detection_arm_comparison.tsv").open(), delimiter="\t"))
    biomarker = list(csv.DictReader((out / "biomarker_arm_comparison.tsv").open(), delimiter="\t"))
    assert detection[0]["clean_only_detections"] == "1"
    assert detection[0]["original_only_detections"] == "0"
    assert float(detection[0]["clean_minus_original_detection_rate"]) == .5
    assert float(biomarker[0]["clean_minus_original_target_effect"]) == .5
    assert float(biomarker[0]["clean_minus_original_off_target_enriched_calls"]) == -1
    assert (out / "SUCCESS").is_file()
print("[PASS] assembly outcome-comparison fixture")
