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
    # summarize_target_recovery.py stamps all three: the ORIGINAL row-level
    # provenance plus the selection that produced the arm.
    "reference_type", "reference_scale", "selected_reference_type",
]

# (reference_scale, selected_reference_type) per arm.
PRIMARY_SELECTION = ("profiler_scale", "profiler_scale_primary")
SENSITIVITY_SELECTION = ("read_proportional", "read_proportional")
# Original row-level provenance per profiler; unchanged by the selection.
ROW_REFERENCE = {"kraken2_bracken": "read_proportional",
                 "metaphlan4": "genome_equivalent"}


def write(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def row(profiler, sample, ratio, error, recovery, reference, selection=None):
    scale, selected = selection or PRIMARY_SELECTION
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
        "reference_scale": scale, "selected_reference_type": selected,
    }


class ComparisonTest(unittest.TestCase):
    def test_paired_comparison_and_bracken_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            primary = root / "primary.tsv"
            sensitivity = root / "sensitivity.tsv"
            # Row provenance is a property of the SOURCE PROFILE, so a
            # MetaPhlAn row stays genome_equivalent in both arms; only the
            # selection columns differ.
            write(primary, [
                row("kraken2_bracken", "b1", 1.0, 0.0, "Good", "read_proportional"),
                row("metaphlan4", "m1", 1.05, 0.05, "Good", "genome_equivalent"),
                row("metaphlan4", "m2", 0.8, 0.2, "Average", "genome_equivalent"),
            ])
            write(sensitivity, [
                row("kraken2_bracken", "b1", 1.0, 0.0, "Good", "read_proportional",
                    SENSITIVITY_SELECTION),
                row("metaphlan4", "m1", 1.8, 0.8, "Poor / missed",
                    "genome_equivalent", SENSITIVITY_SELECTION),
                row("metaphlan4", "m2", 0.8, 0.2, "Average", "genome_equivalent",
                    SENSITIVITY_SELECTION),
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
                row("kraken2_bracken", "b1", 1.2, 0.2, "Average",
                    "read_proportional", SENSITIVITY_SELECTION),
                row("metaphlan4", "m1", 1.8, 0.8, "Poor / missed",
                    "genome_equivalent", SENSITIVITY_SELECTION),
                row("metaphlan4", "m2", 0.8, 0.2, "Average", "genome_equivalent",
                    SENSITIVITY_SELECTION),
            ]
            write(broken, bad_rows)
            failed = subprocess.run(["python3", str(SCRIPT), "--primary", str(primary),
                                     "--sensitivity", str(broken),
                                     "--outdir", str(root / "bad")],
                                    text=True, capture_output=True)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("Bracken differs", failed.stdout + failed.stderr)


class SelectedReferenceContractTest(unittest.TestCase):
    """The comparator must carry the selected estimand, not only the row type."""

    def arms(self, root, primary_rows=None, sensitivity_rows=None):
        primary = root / "primary.tsv"
        sensitivity = root / "sensitivity.tsv"
        write(primary, primary_rows if primary_rows is not None else [
            row("kraken2_bracken", "b1", 1.0, 0.0, "Good", "read_proportional"),
            row("metaphlan4", "m1", 1.05, 0.05, "Good", "genome_equivalent"),
        ])
        write(sensitivity, sensitivity_rows if sensitivity_rows is not None else [
            row("kraken2_bracken", "b1", 1.0, 0.0, "Good", "read_proportional",
                SENSITIVITY_SELECTION),
            row("metaphlan4", "m1", 1.8, 0.8, "Poor / missed", "genome_equivalent",
                SENSITIVITY_SELECTION),
        ])
        return primary, sensitivity

    def run_comparator(self, root, **kwargs):
        primary, sensitivity = self.arms(root, **kwargs)
        return subprocess.run(
            ["python3", str(SCRIPT), "--primary", str(primary),
             "--sensitivity", str(sensitivity), "--outdir", str(root / "out")],
            text=True, capture_output=True)

    def test_selected_fields_are_written_per_arm(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            done = self.run_comparator(root)
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            with (root / "out/target_recovery_reference_comparison.tsv").open(
                    newline="") as handle:
                rows = {r["profiler"]: r
                        for r in csv.DictReader(handle, delimiter="\t")}
            metaphlan = rows["metaphlan4"]
            # A genome_equivalent SOURCE row measured on the read-proportional
            # scale must be reported as the read_proportional estimand.
            self.assertEqual(metaphlan["sensitivity_row_reference_type"],
                             "genome_equivalent")
            self.assertEqual(metaphlan["sensitivity_selected_reference_type"],
                             "read_proportional")
            self.assertEqual(metaphlan["sensitivity_reference_scale"],
                             "read_proportional")
            self.assertEqual(metaphlan["sensitivity_selected_estimand"],
                             "read_proportional")
            self.assertEqual(metaphlan["primary_selected_estimand"],
                             "genome_equivalent")
            self.assertEqual(metaphlan["primary_selected_reference_type"],
                             "profiler_scale_primary")
            bracken = rows["kraken2_bracken"]
            self.assertEqual(bracken["primary_selected_estimand"],
                             "read_proportional")
            self.assertEqual(bracken["sensitivity_selected_estimand"],
                             "read_proportional")
            # The ambiguous legacy alias survives for old consumers.
            self.assertEqual(metaphlan["primary_reference_type"],
                             "genome_equivalent")
            self.assertEqual(metaphlan["sensitivity_reference_type"],
                             "genome_equivalent")

    def expect_failure(self, needle, **kwargs):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            done = self.run_comparator(root, **kwargs)
            self.assertNotEqual(done.returncode, 0,
                                f"comparator accepted bad input ({needle})")
            self.assertIn(needle, done.stdout + done.stderr)

    def test_sensitivity_claiming_genome_equivalent_selection_fails(self):
        self.expect_failure("must be summarised as", sensitivity_rows=[
            row("kraken2_bracken", "b1", 1.0, 0.0, "Good", "read_proportional",
                PRIMARY_SELECTION),
            row("metaphlan4", "m1", 1.8, 0.8, "Poor / missed", "genome_equivalent",
                PRIMARY_SELECTION),
        ])

    def test_primary_claiming_read_proportional_selection_fails(self):
        self.expect_failure("must be summarised as", primary_rows=[
            row("kraken2_bracken", "b1", 1.0, 0.0, "Good", "read_proportional",
                SENSITIVITY_SELECTION),
            row("metaphlan4", "m1", 1.05, 0.05, "Good", "genome_equivalent",
                SENSITIVITY_SELECTION),
        ])

    def test_bracken_with_inconsistent_selection_fails(self):
        self.expect_failure("unknown reference selection", primary_rows=[
            row("kraken2_bracken", "b1", 1.0, 0.0, "Good", "read_proportional",
                ("profiler_scale", "read_proportional")),
            row("metaphlan4", "m1", 1.05, 0.05, "Good", "genome_equivalent"),
        ])

    def test_bracken_with_genome_equivalent_row_type_fails(self):
        self.expect_failure("the source rows must be 'read_proportional'",
                            primary_rows=[
            row("kraken2_bracken", "b1", 1.0, 0.0, "Good", "genome_equivalent"),
            row("metaphlan4", "m1", 1.05, 0.05, "Good", "genome_equivalent"),
        ])

    def test_metaphlan_with_read_proportional_row_type_fails(self):
        self.expect_failure("the source rows must be 'genome_equivalent'",
                            primary_rows=[
            row("kraken2_bracken", "b1", 1.0, 0.0, "Good", "read_proportional"),
            row("metaphlan4", "m1", 1.05, 0.05, "Good", "read_proportional"),
        ])

    def test_blank_selection_fails(self):
        self.expect_failure("blank selected_reference_type", primary_rows=[
            row("kraken2_bracken", "b1", 1.0, 0.0, "Good", "read_proportional"),
            row("metaphlan4", "m1", 1.05, 0.05, "Good", "genome_equivalent",
                ("profiler_scale", "")),
        ])

    def test_missing_selection_columns_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            primary, sensitivity = self.arms(root)
            for path in (primary, sensitivity):
                with path.open(newline="") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                fields = [f for f in FIELDS if f != "selected_reference_type"]
                for item in rows:
                    item.pop("selected_reference_type")
                with path.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=fields,
                                            delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                    writer.writerows(rows)
            done = subprocess.run(
                ["python3", str(SCRIPT), "--primary", str(primary),
                 "--sensitivity", str(sensitivity), "--outdir", str(root / "out")],
                text=True, capture_output=True)
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("lacks columns: selected_reference_type",
                          done.stdout + done.stderr)


if __name__ == "__main__":
    unittest.main()
