#!/usr/bin/env python3
"""Fixture test for visible cohort IQR boxes on corrected community endpoints."""
import csv
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PLOT = ROOT / "analysis_v2/scripts/plot_three_cohort_quantitative_recovery.R"
LABELS = ("Bfrag", "Csym", "Dpne", "Fnuc", "Hhat", "Pmic", "Pana", "Psto", "Porp", "Pint")


class RecoveryBoxFigure(unittest.TestCase):
    def fixture(self, path):
        rows = []
        for cohort in ("feng", "yachida", "zeller"):
            for condition in ("Control", "Adenoma", "CRC"):
                for profiler in ("kraken2_bracken", "metaphlan4"):
                    for target in LABELS:
                        for dose in (.0001, .0005, .001):
                            for i, ratio in enumerate((.5, 1, 2, 0)):
                                rows.append(dict(cohort=cohort, condition=condition,
                                    analysis_population="community", assembly_arm="original",
                                    profiler=profiler, target_label=target,
                                    spike_fraction_target=dose, sample_id="s" + str(i),
                                    reference_type=("genome_equivalent" if profiler == "metaphlan4"
                                                    else "read_proportional"),
                                    observed_abundance_fraction=dose * ratio,
                                    expected_abundance_profiler_scale=dose))
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def test_boxes_are_nonzero_and_zeros_reported(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.fixture(root / "endpoints.tsv")
            result = subprocess.run(["Rscript", str(PLOT), "--endpoints", str(root / "endpoints.tsv"),
                                     "--outdir", str(root / "figure")], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with (root / "figure/quantitative_recovery_box_summary.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 540)
            self.assertTrue(all(int(row["n_total"]) == 4 and
                                int(row["n_positive"]) == 3 and
                                int(row["n_zero"]) == 1 for row in rows))
            self.assertTrue(all(float(row["q1"]) < float(row["median"]) <
                                float(row["q3"]) for row in rows))
            for name in ("three_cohort_quantitative_recovery_boxes.pdf",
                         "three_cohort_quantitative_recovery_boxes.png", "SUCCESS"):
                self.assertGreater((root / "figure" / name).stat().st_size, 0)

    def test_wrong_reference_fails(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.fixture(root / "endpoints.tsv")
            text = (root / "endpoints.tsv").read_text()
            (root / "endpoints.tsv").write_text(text.replace("genome_equivalent", "read_proportional", 1))
            result = subprocess.run(["Rscript", str(PLOT), "--endpoints", str(root / "endpoints.tsv"),
                                     "--outdir", str(root / "figure")], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Wrong endpoint reference type", result.stderr)


if __name__ == "__main__":
    unittest.main()
