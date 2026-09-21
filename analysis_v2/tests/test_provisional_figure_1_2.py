#!/usr/bin/env python3
"""Development-only Figure 1 schematics and Figure 2 detection gates."""
import csv
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DESIGN = ROOT / "analysis_v2/scripts/build_study_design_schematics.py"
DETECTION = ROOT / "analysis_v2/scripts/plot_three_cohort_detection.R"
TARGETS = ("Bfrag", "Csym", "Dpne", "Fnuc", "Hhat", "Pmic", "Pana", "Psto", "Porp", "Pint")


class ProvisionalFigure12(unittest.TestCase):
    def test_schematics_parse_and_contain_reference_formulas(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            done = subprocess.run(["python3", str(DESIGN), "--outdir", str(root / "design")],
                                  capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            for name in ("figure_1A_study_design.svg", "figure_1B_profiler_estimands.svg"):
                ET.parse(root / "design" / name)
            text = (root / "design/figure_1B_profiler_estimands.svg").read_text()
            self.assertIn("Gₑff", text)
            self.assertIn("[(1 - F) o + qₜ] / D", text)

    def fixture(self, path):
        rows = []
        for cohort in ("feng", "yachida", "zeller"):
            for condition in ("Control", "Adenoma", "CRC"):
                for profiler in ("kraken2_bracken", "metaphlan4"):
                    for target in TARGETS:
                        for dose in (.0001, .0005, .001):
                            for sample in range(3):
                                observed = 0 if sample == 0 else dose
                                rows.append(dict(cohort=cohort, condition=condition,
                                    analysis_population="community", assembly_arm="original",
                                    profiler=profiler, target_label=target,
                                    spike_fraction_target=dose, sample_id="s" + str(sample),
                                    reference_type=("genome_equivalent" if profiler == "metaphlan4"
                                                    else "read_proportional"),
                                    observed_abundance_fraction=observed,
                                    observed_detected_native_nonzero=int(observed > 0)))
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def test_detection_is_complete_and_wilson_counts_are_correct(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.fixture(root / "endpoints.tsv")
            done = subprocess.run(["Rscript", str(DETECTION), "--endpoints",
                                   str(root / "endpoints.tsv"), "--outdir", str(root / "plot")],
                                  capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            with (root / "plot/detection_prevalence_wilson.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 540)
            self.assertTrue(all(int(r["detected"]) == 2 and int(r["samples"]) == 3 for r in rows))
            self.assertTrue(all(0 < float(r["wilson_low"]) < float(r["prevalence"]) <
                                float(r["wilson_high"]) < 1 for r in rows))
            self.assertGreater((root / "plot/three_cohort_detection_heatmap.png").stat().st_size, 0)

    def test_read_reference_metaphlan_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.fixture(root / "endpoints.tsv")
            text = (root / "endpoints.tsv").read_text()
            (root / "endpoints.tsv").write_text(text.replace("genome_equivalent", "read_proportional", 1))
            done = subprocess.run(["Rscript", str(DETECTION), "--endpoints",
                                   str(root / "endpoints.tsv"), "--outdir", str(root / "plot")],
                                  capture_output=True, text=True)
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("Wrong endpoint reference type", done.stderr)


if __name__ == "__main__":
    unittest.main()
