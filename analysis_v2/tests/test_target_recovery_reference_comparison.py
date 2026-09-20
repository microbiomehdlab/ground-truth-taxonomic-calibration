#!/usr/bin/env python3
import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis_v2/scripts/compare_target_recovery_references.py"
FIELDS = [
    "cohort", "study", "sample_id", "condition", "analysis_population",
    "profiler", "baseline_id", "target_feature", "dose_index",
    "nominal_target_fraction", "target_fraction_for_feature",
    "response_signal_selected", "implanted_signal_selected",
    "observed_over_expected", "absolute_relative_error", "recovery_class",
    "reference_type",
]


def write(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def row(profiler, sample, ratio, error, recovery, reference):
    return {
        "cohort": "yachida", "study": "study", "sample_id": sample,
        "condition": "CRC", "analysis_population": "community",
        "profiler": profiler, "baseline_id": sample + "_base",
        "target_feature": "Fusobacterium nucleatum", "dose_index": "1",
        "nominal_target_fraction": "0.001",
        "target_fraction_for_feature": "0.001",
        "response_signal_selected": str(ratio * 0.001),
        "implanted_signal_selected": "0.001", "observed_over_expected": str(ratio),
        "absolute_relative_error": str(error), "recovery_class": recovery,
        "reference_type": reference,
    }


class ComparisonTest(unittest.TestCase):
    def test_paired_comparison_and_bracken_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            primary = root / "primary.tsv"
            sensitivity = root / "sensitivity.tsv"
            write(primary, [
                row("kraken2_bracken", "b1", 1.0, 0.0, "Good", "read_proportional"),
                row("metaphlan4", "m1", 1.05, 0.05, "Good", "genome_equivalent"),
                row("metaphlan4", "m2", 0.8, 0.2, "Average", "genome_equivalent"),
            ])
            write(sensitivity, [
                row("kraken2_bracken", "b1", 1.0, 0.0, "Good", "read_proportional"),
                row("metaphlan4", "m1", 1.8, 0.8, "Poor / missed", "read_proportional"),
                row("metaphlan4", "m2", 0.8, 0.2, "Average", "read_proportional"),
            ])
            out = root / "out"
            done = subprocess.run(["python3", str(SCRIPT), "--primary", str(primary),
                                   "--sensitivity", str(sensitivity), "--outdir", str(out)],
                                  text=True, capture_output=True)
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            with (out / "target_recovery_reference_comparison.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            m1 = next(item for item in rows if item["sample_id"] == "m1")
            self.assertEqual(m1["class_transition"], "Poor / missed -> Good")
            self.assertEqual(m1["primary_class_improved"], "1")
            self.assertEqual(m1["primary_class_worsened"], "0")
            with (out / "target_recovery_reference_overall.tsv").open(newline="") as handle:
                overall = list(csv.DictReader(handle, delimiter="\t"))
            metaphlan = next(item for item in overall if item["profiler"] == "metaphlan4")
            self.assertEqual(metaphlan["observations"], "2")
            self.assertEqual(metaphlan["primary_improved_fraction"], "0.5")
            self.assertEqual(metaphlan["primary_worsened_fraction"], "0")
            self.assertTrue((out / "SUCCESS").is_file())

            broken = root / "broken.tsv"
            bad_rows = [
                row("kraken2_bracken", "b1", 1.2, 0.2, "Average", "read_proportional"),
                row("metaphlan4", "m1", 1.8, 0.8, "Poor / missed", "read_proportional"),
                row("metaphlan4", "m2", 0.8, 0.2, "Average", "read_proportional"),
            ]
            write(broken, bad_rows)
            failed = subprocess.run(["python3", str(SCRIPT), "--primary", str(primary),
                                     "--sensitivity", str(broken),
                                     "--outdir", str(root / "bad")],
                                    text=True, capture_output=True)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("Bracken differs", failed.stdout + failed.stderr)


if __name__ == "__main__":
    unittest.main()
