#!/usr/bin/env python3
"""Small real-plot fixture for development presentation panels."""
import csv
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PLOT = ROOT / "analysis_v2/scripts/plot_institute_presentation_panels.R"
RUNNER = ROOT / "analysis_v2/run_institute_presentation_panels.sh"


def write(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class PresentationPanels(unittest.TestCase):
    def test_plot_from_corrected_reference(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            baseline, endpoints = [], []
            for cohort in ("feng", "yachida", "zeller"):
                for condition in ("Control", "Adenoma", "CRC"):
                    for profiler in ("kraken2_bracken", "metaphlan4"):
                        sample = "{}_{}_{}".format(cohort, condition, profiler)
                        for target in ("Bfrag", "Fnuc"):
                            baseline.append(dict(cohort=cohort, sample_id=sample,
                                                 condition=condition, profiler=profiler,
                                                 target_label=target, abundance_fraction=.001))
                        for target in ("Bfrag", "Dpne", "Fnuc"):
                            for dose in (.0001, .0005, .001):
                                endpoints.append(dict(cohort=cohort, sample_id=sample,
                                    condition=condition, analysis_population="community",
                                    assembly_arm="original", profiler=profiler,
                                    target_label=target, spike_fraction_target=dose,
                                    reference_type="genome_equivalent" if profiler == "metaphlan4"
                                                   else "read_proportional",
                                    observed_abundance_fraction=dose,
                                    expected_abundance_profiler_scale=dose * 2 if profiler == "metaphlan4" else dose,
                                    read_proportional_reference=dose,
                                    baseline_abundance_fraction=0,
                                    spike_fraction_total=dose))
            write(root / "baseline.tsv", baseline)
            write(root / "endpoints.tsv", endpoints)
            done = subprocess.run(["Rscript", str(PLOT), "--baseline", str(root / "baseline.tsv"),
                                   "--endpoints", str(root / "endpoints.tsv"),
                                   "--outdir", str(root / "slides")], capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            for name in ("slide3_baseline_visibility", "slide4_focused_recovery",
                         "slide4_metaphlan_reference_comparison", "slide5_fnuc_detection"):
                self.assertGreater((root / "slides" / (name + ".png")).stat().st_size, 0)
            with (root / "slides/sources/slide4_metaphlan_reference_comparison.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 54)
            original = [r for r in rows if r["reference"] == "Original read-proportional"]
            corrected = [r for r in rows if r["reference"] == "Corrected genome-equivalent"]
            self.assertTrue(all(abs(float(r["median"])) < 1e-9 for r in original))
            self.assertTrue(all(abs(float(r["median"]) + 1) < 1e-9 for r in corrected))
            self.assertTrue((root / "slides/SUCCESS").is_file())

    def test_runner_syntax(self):
        done = subprocess.run(["bash", "-n", str(RUNNER)], capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stderr)


if __name__ == "__main__":
    unittest.main()
