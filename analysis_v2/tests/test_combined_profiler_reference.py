#!/usr/bin/env python3
"""Combined-profiler primary tables and profiler-aware reference validation.

A primary table normally holds two row-level reference types at once:
    kraken2_bracken -> read_proportional
    metaphlan4      -> genome_equivalent
That is correct, not a mixed-reference error. Validation is within profiler.
"""

from __future__ import annotations

import csv
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import duckdb

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analysis_v2/scripts"))
from reference_scale import validate_profiler_reference_pairs  # noqa: E402

RECOV = ROOT / "analysis_v2/scripts/summarize_target_recovery.py"

COLUMNS = [
    "cohort", "study", "sample_id", "condition", "analysis_population", "profiler",
    "baseline_id", "dose_index", "nominal_target_fraction",
    "target_fraction_for_feature", "feature", "baseline_abundance_fraction",
    "observed_abundance_fraction", "response_signal", "observed_detected",
    "baseline_detected", "is_direct_target", "reference_type",
    "response_signal_profiler_scale", "implanted_signal_profiler_scale",
    "expected_abundance_profiler_scale", "retained_baseline_profiler_scale",
    "response_delta_profiler_scale", "quantitative_log2_ratio_profiler_scale",
    "signed_bounded_error_profiler_scale", "absolute_bounded_error_profiler_scale",
]


def _bounded(observed, expected):
    total = observed + expected
    return (observed - expected) / total if total > 0 else 0.0


def row(profiler, reference_type, sample, feature, *, f_read, signal_read,
        f_ps=None, signal_ps=None, detected=1, direct=1):
    f_ps = f_read if f_ps is None else f_ps
    signal_ps = signal_read if signal_ps is None else signal_ps
    return {
        "cohort": "yachida", "study": "S", "sample_id": sample, "condition": "CRC",
        "analysis_population": "community", "profiler": profiler,
        "baseline_id": f"b_{sample}", "dose_index": 1,
        "nominal_target_fraction": f_read, "target_fraction_for_feature": f_read,
        "feature": feature, "baseline_abundance_fraction": 0.0,
        "observed_abundance_fraction": signal_read, "response_signal": signal_read,
        "observed_detected": detected, "baseline_detected": 0,
        "is_direct_target": direct, "reference_type": reference_type,
        "response_signal_profiler_scale": signal_ps,
        "implanted_signal_profiler_scale": f_ps,
        "expected_abundance_profiler_scale": f_ps,
        "retained_baseline_profiler_scale": 0.0,
        "response_delta_profiler_scale": signal_ps - f_ps,
        "quantitative_log2_ratio_profiler_scale": None,
        "signed_bounded_error_profiler_scale": _bounded(signal_ps, f_ps),
        "absolute_bounded_error_profiler_scale": abs(_bounded(signal_ps, f_ps)),
    }


def write_parquet(path, rows):
    con = duckdb.connect()
    con.execute("CREATE TABLE t (" + ", ".join(
        f'"{c}" DOUBLE' if c not in {
            "cohort", "study", "sample_id", "condition", "analysis_population",
            "profiler", "baseline_id", "feature", "reference_type"} else f'"{c}" VARCHAR'
        for c in COLUMNS) + ")")
    con.executemany(
        f"INSERT INTO t VALUES ({', '.join('?' for _ in COLUMNS)})",
        [[r[c] for c in COLUMNS] for r in rows])
    con.execute(f"COPY t TO '{str(path).replace(chr(39), chr(39)*2)}' (FORMAT PARQUET)")


def metrics_file(path):
    path.write_text(
        "cohort\tanalysis_population\ttarget_label\tprofiler\tcontrast\t"
        "spike_fraction_target\tq_threshold\ttarget_called\ttarget_effect\ttarget_q_value\n"
        "yachida\tcommunity\tFnuc\tkraken2_bracken\t"
        "spiked_vs_matched_baseline__pooled_conditions\t0.001\t0.05\t1\t2\t0.01\n",
        encoding="utf-8")
    return path


def run_recovery(parquet, outdir, scale, metrics):
    return subprocess.run(
        ["python3", str(RECOV), "--responses", str(parquet),
         "--biomarker-metrics", str(metrics), "--outdir", str(outdir),
         "--reference-scale", scale], text=True, capture_output=True)


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


VALID_COMBINED = [("kraken2_bracken", "read_proportional"),
                  ("metaphlan4", "genome_equivalent")]


class ValidatorTest(unittest.TestCase):
    """Rules 1-8 exercised directly against the validator."""

    def test_combined_primary_table_is_accepted(self):
        self.assertEqual(
            validate_profiler_reference_pairs(VALID_COMBINED, "profiler_scale"),
            "profiler_scale_primary")

    def test_read_proportional_accepts_the_same_combined_table(self):
        self.assertEqual(
            validate_profiler_reference_pairs(VALID_COMBINED, "read_proportional"),
            "read_proportional")

    def test_wrong_bracken_reference_fails(self):
        with self.assertRaises(SystemExit) as caught:
            validate_profiler_reference_pairs(
                [("kraken2_bracken", "genome_equivalent")], "profiler_scale")
        self.assertIn("requires 'read_proportional'", str(caught.exception))

    def test_wrong_metaphlan_reference_fails(self):
        with self.assertRaises(SystemExit) as caught:
            validate_profiler_reference_pairs(
                [("metaphlan4", "read_proportional")], "profiler_scale")
        self.assertIn("requires 'genome_equivalent'", str(caught.exception))

    def test_multiple_types_within_one_profiler_fail(self):
        with self.assertRaises(SystemExit) as caught:
            validate_profiler_reference_pairs(
                [("metaphlan4", "genome_equivalent"),
                 ("metaphlan4", "read_proportional")], "profiler_scale")
        self.assertIn("multiple reference types", str(caught.exception))

    def test_unknown_profiler_fails(self):
        with self.assertRaises(SystemExit) as caught:
            validate_profiler_reference_pairs(
                [("mystery_profiler", "read_proportional")], "profiler_scale")
        self.assertIn("unknown profiler", str(caught.exception))

    def test_blank_reference_type_fails(self):
        with self.assertRaises(SystemExit) as caught:
            validate_profiler_reference_pairs(
                [("metaphlan4", "")], "profiler_scale")
        self.assertIn("blank reference_type", str(caught.exception))


class CombinedRecoveryTest(unittest.TestCase):

    def _combined_rows(self):
        # Bracken: read fraction 0.001, response 0.001 -> perfect recovery.
        # MetaPhlAn: read fraction 0.001 but the profiler-scale implanted signal
        # is 0.002 and the profiler-scale response is 0.002. Judged correctly it
        # is Good; judged against the read fraction it would look 100% wrong.
        return [
            row("kraken2_bracken", "read_proportional", "s1", "Fnuc_feature",
                f_read=0.001, signal_read=0.001),
            row("metaphlan4", "genome_equivalent", "s2", "Fnuc_feature",
                f_read=0.001, signal_read=0.001, f_ps=0.002, signal_ps=0.002),
        ]

    def test_combined_profiler_scale_accepted_and_classes_correct(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            parquet = root / "responses.parquet"
            write_parquet(parquet, self._combined_rows())
            out = root / "out"
            done = run_recovery(parquet, out, "profiler_scale",
                                metrics_file(root / "m.tsv"))
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            rows = {r["profiler"]: r for r in
                    read_tsv(out / "target_recovery_observations.tsv")}
            self.assertEqual(set(rows), {"kraken2_bracken", "metaphlan4"})
            # Row-level provenance preserved.
            self.assertEqual(rows["kraken2_bracken"]["reference_type"],
                             "read_proportional")
            self.assertEqual(rows["metaphlan4"]["reference_type"],
                             "genome_equivalent")
            self.assertEqual(rows["metaphlan4"]["selected_reference_type"],
                             "profiler_scale_primary")
            # MetaPhlAn recovery uses the profiler-scale implanted signal.
            self.assertAlmostEqual(
                float(rows["metaphlan4"]["implanted_signal_selected"]), 0.002)
            self.assertAlmostEqual(
                float(rows["metaphlan4"]["response_signal_selected"]), 0.002)
            self.assertAlmostEqual(
                float(rows["metaphlan4"]["absolute_relative_error"]), 0.0)
            self.assertEqual(rows["metaphlan4"]["recovery_class"], "Good")
            # Had it been compared with the read fraction it would have been
            # |0.002-0.001|/0.001 = 1.0 -> Poor / missed.
            self.assertNotEqual(rows["metaphlan4"]["recovery_class"], "Poor / missed")

    def test_read_proportional_sensitivity_from_the_same_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            parquet = root / "responses.parquet"
            write_parquet(parquet, self._combined_rows())
            out = root / "out"
            done = run_recovery(parquet, out, "read_proportional",
                                metrics_file(root / "m.tsv"))
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            rows = {r["profiler"]: r for r in
                    read_tsv(out / "target_recovery_observations.tsv")}
            self.assertEqual(rows["metaphlan4"]["selected_reference_type"],
                             "read_proportional")
            # Row-level provenance is not rewritten.
            self.assertEqual(rows["metaphlan4"]["reference_type"],
                             "genome_equivalent")
            self.assertAlmostEqual(
                float(rows["metaphlan4"]["implanted_signal_selected"]), 0.001)

    def test_bracken_identical_between_scales(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            parquet = root / "responses.parquet"
            write_parquet(parquet, self._combined_rows())
            metrics = metrics_file(root / "m.tsv")
            a, b = root / "ps", root / "rp"
            self.assertEqual(run_recovery(parquet, a, "profiler_scale", metrics).returncode, 0)
            self.assertEqual(run_recovery(parquet, b, "read_proportional", metrics).returncode, 0)
            ps = {r["profiler"]: r for r in read_tsv(a / "target_recovery_observations.tsv")}
            rp = {r["profiler"]: r for r in read_tsv(b / "target_recovery_observations.tsv")}
            for field in ("observed_over_expected", "absolute_relative_error",
                          "recovery_class", "implanted_signal_selected",
                          "response_signal_selected"):
                self.assertEqual(ps["kraken2_bracken"][field],
                                 rp["kraken2_bracken"][field],
                                 f"Bracken {field} changed with the scale")
            # Detection is invariant for both profilers.
            for profiler in ("kraken2_bracken", "metaphlan4"):
                self.assertEqual(ps[profiler]["observed_detected"],
                                 rp[profiler]["observed_detected"])

    def test_wrong_bracken_reference_rejected_end_to_end(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bad = [row("kraken2_bracken", "genome_equivalent", "s1", "F",
                       f_read=0.001, signal_read=0.001)]
            parquet = root / "bad.parquet"
            write_parquet(parquet, bad)
            done = run_recovery(parquet, root / "out", "profiler_scale",
                                metrics_file(root / "m.tsv"))
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("requires 'read_proportional'", done.stdout + done.stderr)

    def test_no_silent_fallback_when_profiler_columns_absent(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            con = duckdb.connect()
            legacy = root / "legacy.parquet"
            con.execute("CREATE TABLE t AS SELECT 'yachida' AS cohort, "
                        "'kraken2_bracken' AS profiler, 'read_proportional' AS reference_type, "
                        "0.001 AS target_fraction_for_feature, 0.001 AS response_signal")
            con.execute(f"COPY t TO '{legacy}' (FORMAT PARQUET)")
            done = run_recovery(legacy, root / "out", "profiler_scale",
                                metrics_file(root / "m.tsv"))
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("needs column", done.stdout + done.stderr)


if __name__ == "__main__":
    unittest.main()
