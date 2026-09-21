#!/usr/bin/env python3
"""End-to-end gates for the corrected three-cohort recoverability figure."""
import csv
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILDER = ROOT / "analysis_v2/scripts/build_three_cohort_recoverability_figure.py"
PLOT = ROOT / "analysis_v2/scripts/plot_three_cohort_recoverability.R"
PANEL = ROOT / "spikes/spike_panel.tsv"
DOSES = (.0001, .0005, .001, .005, .01, .05)


def write(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class ThreeCohortRecoverability(unittest.TestCase):
    def fixture(self):
        with PANEL.open(newline="", encoding="utf-8") as handle:
            labels = [row["label"] for row in csv.DictReader(handle, delimiter="\t")]
        biomarker, endpoints = [], []
        for cohort in ("feng", "yachida", "zeller"):
            for condition in ("Control", "Adenoma", "CRC"):
                for profiler in ("kraken2_bracken", "metaphlan4"):
                    for target_index, label in enumerate(labels):
                        for dose_index, dose in enumerate(DOSES):
                            q = .001 if dose_index >= target_index % 6 else .5
                            biomarker.append(dict(cohort=cohort, analysis_population="independent",
                                assembly_arm="original", profiler=profiler, target_label=label,
                                contrast="spiked_vs_matched_baseline__background_" + condition,
                                spike_fraction_target=dose, q_threshold=.05,
                                target_called=int(q < .05), target_q_value=q))
                            for sample_index in range(3):
                                implanted = dose * (1.4 if profiler == "metaphlan4" else 1)
                                endpoints.append(dict(cohort=cohort, condition=condition,
                                    analysis_population="independent", assembly_arm="original",
                                    profiler=profiler, target_label=label,
                                    spike_fraction_target=dose, sample_id="S" + str(sample_index),
                                    reference_type=("genome_equivalent" if profiler == "metaphlan4"
                                                    else "read_proportional"),
                                    baseline_abundance_fraction=.001 * (target_index + 1),
                                    implanted_signal_profiler_scale=implanted,
                                    recovered_spike_signal_profiler_scale=implanted * (.8 + .1 * sample_index),
                                    expected_abundance_profiler_scale=.001 + implanted,
                                    observed_abundance_fraction=.001 + implanted * (.8 + .1 * sample_index),
                                    observed_detected_native_nonzero=1))
        return biomarker, endpoints

    def invoke(self, biomarker, endpoints, root):
        calls, paired = root / "biomarkers.tsv", root / "endpoints.tsv"
        write(calls, biomarker)
        write(paired, endpoints)
        return subprocess.run(["python3", str(BUILDER), "--biomarkers", str(calls),
                               "--endpoints", str(paired), "--panel", str(PANEL),
                               "--outdir", str(root / "source")],
                              capture_output=True, text=True)

    def test_all_panels_and_corrected_scale(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            biomarkers, endpoints = self.fixture()
            # A zero reported q-value is a numerical lower-bound case. It must
            # remain in the source but cannot set the display axis to 300.
            biomarkers[0]["target_q_value"] = 0
            # An extreme recovered/implanted ratio must remain auditable while
            # the variability scatterplot stays readable.
            endpoints[2]["recovered_spike_signal_profiler_scale"] = .002
            result = self.invoke(biomarkers, endpoints, root)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with (root / "source/minimum_biomarker_fraction.tsv").open() as handle:
                threshold = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(threshold), 180)
            self.assertEqual(threshold[0]["tested_doses"], "6")
            plot = subprocess.run(["Rscript", str(PLOT), "--source-dir", str(root / "source"),
                                   "--outdir", str(root / "figures")], capture_output=True, text=True)
            self.assertEqual(plot.returncode, 0, plot.stdout + plot.stderr)
            for stem in ("panel_A_minimum_spike_fraction", "panel_B_recovery_drivers",
                         "panel_B_driver_tertiles", "panel_C_spearman_associations",
                         "three_cohort_recoverability_combined"):
                self.assertGreater((root / "figures" / (stem + ".pdf")).stat().st_size, 0)
            with (root / "figures/panel_B_display_capped_points.tsv").open() as handle:
                capped = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(float(row["biomarker_q"]) == 0 for row in capped))
            with (root / "figures/panel_B_variability_capped_points.tsv").open() as handle:
                x_capped = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(float(row["driver_value"]) > 2 for row in x_capped))
            with (root / "figures/panel_B_driver_tertiles_source.tsv").open() as handle:
                tertiles = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(tertiles), 24)
            self.assertTrue(all(int(row["n"]) == 30 for row in tertiles))

    def test_read_scale_metaphlan_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            biomarkers, endpoints = self.fixture()
            for row in endpoints:
                if row["profiler"] == "metaphlan4":
                    row["reference_type"] = "read_proportional"
                    break
            result = self.invoke(biomarkers, endpoints, root)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("wrong row reference type", result.stderr)

    def test_missing_biomarker_dose_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            biomarkers, endpoints = self.fixture()
            result = self.invoke(biomarkers[1:], endpoints, root)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("complete three-cohort six-dose grid", result.stderr)


if __name__ == "__main__":
    unittest.main()
