#!/usr/bin/env python3
"""Input contract for the profiler-scale continuous model.

`mgcv` is absent locally, so model fitting cannot run here. These checks are
package-independent: they exercise the argument contract and the column
selection statically, and the fixture is written so the same file can be run
inside the analysis image where mgcv exists.

Container command:
    apptainer exec --cleanenv --pwd "$PWD" "$ANALYSIS_SIF" \
      python3 analysis_v2/tests/test_continuous_model_reference_contract.py
"""

from __future__ import annotations

import csv
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis_v2/scripts/fit_continuous_dose_response.R"

FIELDS = ["cohort", "sample_id", "condition", "analysis_population", "profiler",
          "target_label", "assembly_arm", "spike_fraction_target",
          "nominal_target_fraction", "recovered_spike_signal", "reference_type",
          "implanted_signal_profiler_scale", "recovered_spike_signal_profiler_scale"]


def mgcv_available() -> bool:
    if not shutil.which("Rscript"):
        return False
    done = subprocess.run(
        ["Rscript", "-e", 'cat(requireNamespace("mgcv", quietly=TRUE))'],
        text=True, capture_output=True)
    return done.stdout.strip().upper() == "TRUE"


def write_fixture(path: Path, *, metaphlan_reference="genome_equivalent",
                  bracken_reference="read_proportional", mixed=False):
    """Two profilers; MetaPhlAn implanted signals differ per sample via G_eff."""
    rows = []
    doses = [0.001, 0.005, 0.01]
    for sample_index, sample in enumerate(("S1", "S2", "S3"), start=1):
        for target in ("Fnuc", "Dpne"):
            for dose in doses:
                # Bracken: profiler scale equals the read scale exactly.
                rows.append({
                    "cohort": "yachida", "sample_id": f"K{sample}",
                    "condition": "CRC" if sample_index % 2 else "Control",
                    "analysis_population": "community",
                    "profiler": "kraken2_bracken", "target_label": target,
                    "assembly_arm": "original", "spike_fraction_target": dose,
                    "nominal_target_fraction": dose,
                    "recovered_spike_signal": dose * 0.98,
                    "reference_type": bracken_reference,
                    "implanted_signal_profiler_scale": dose,
                    "recovered_spike_signal_profiler_scale": dose * 0.98,
                })
                # MetaPhlAn: sample-specific G_eff makes the implanted signal
                # differ from the read fraction and differ between samples.
                scale = 1.0 + 0.10 * sample_index
                rows.append({
                    "cohort": "yachida", "sample_id": f"M{sample}",
                    "condition": "CRC" if sample_index % 2 else "Control",
                    "analysis_population": "community",
                    "profiler": "metaphlan4", "target_label": target,
                    "assembly_arm": "original", "spike_fraction_target": dose,
                    "nominal_target_fraction": dose,
                    "recovered_spike_signal": dose * 0.98,
                    "reference_type": (bracken_reference if mixed and target == "Dpne"
                                       else metaphlan_reference),
                    "implanted_signal_profiler_scale": dose * scale,
                    "recovered_spike_signal_profiler_scale": dose * scale * 0.97,
                })
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def run_model(fixture, outdir, extra=()):
    return subprocess.run(
        ["Rscript", str(SCRIPT), "--input", str(fixture), "--outdir", str(outdir),
         "--cohort", "yachida", "--population", "community",
         "--assembly-arm", "original", *map(str, extra)],
        text=True, capture_output=True)


class StaticContractTest(unittest.TestCase):
    """Package-independent: these always run."""

    def test_script_selects_profiler_scale_columns(self):
        text = SCRIPT.read_text()
        block = text[text.index("if (profiler_scale) {"):]
        self.assertIn("dat$dose_value <- as.numeric(dat$implanted_signal_profiler_scale)",
                      block, "predictor must be the profiler-scale implanted signal")
        self.assertIn(
            "dat$recovered <- as.numeric(dat$recovered_spike_signal_profiler_scale)",
            block, "response must be the profiler-scale recovered signal")

    def test_script_requires_the_profiler_scale_columns(self):
        text = SCRIPT.read_text()
        self.assertIn('"implanted_signal_profiler_scale"', text)
        self.assertIn('"recovered_spike_signal_profiler_scale"', text)
        self.assertIn('"reference_type"', text)

    def test_script_validates_within_profiler(self):
        text = SCRIPT.read_text()
        self.assertIn("kraken2_bracken = \"read_proportional\"", text)
        self.assertIn("metaphlan4 = \"genome_equivalent\"", text)
        self.assertIn("multiple reference types", text)

    def test_script_has_no_reference_scale_default(self):
        text = SCRIPT.read_text()
        self.assertIn("--reference-scale {read_proportional|profiler_scale} is required",
                      text)
        self.assertNotIn('value("--reference-scale", "profiler_scale")', text)
        self.assertNotIn('value("--reference-scale", "read_proportional")', text)

    def test_script_records_omitted_components(self):
        text = SCRIPT.read_text()
        self.assertIn("reference_scale_manifest.tsv", text)
        self.assertIn("omitted_components", text)
        self.assertIn("grouped_dose_available", text)

    def test_script_parses(self):
        if not shutil.which("Rscript"):
            self.skipTest("Rscript unavailable")
        done = subprocess.run(
            ["Rscript", "-e", f"invisible(parse('{SCRIPT}'))"],
            text=True, capture_output=True)
        self.assertEqual(done.returncode, 0, done.stderr)

    def test_bracken_fixture_rows_are_scale_identical(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = write_fixture(Path(tmp) / "f.tsv")
        for row in rows:
            if row["profiler"] != "kraken2_bracken":
                continue
            self.assertAlmostEqual(row["implanted_signal_profiler_scale"],
                                   row["spike_fraction_target"])
            self.assertAlmostEqual(row["recovered_spike_signal_profiler_scale"],
                                   row["recovered_spike_signal"])

    def test_metaphlan_fixture_signals_vary_by_sample(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = write_fixture(Path(tmp) / "f.tsv")
        mpa = [r for r in rows if r["profiler"] == "metaphlan4"]
        ratios = {round(r["implanted_signal_profiler_scale"]
                        / r["spike_fraction_target"], 6) for r in mpa}
        self.assertGreater(len(ratios), 1,
                           "G_eff must make the implanted signal sample-specific")


class ExecutionTest(unittest.TestCase):
    """Only runs where mgcv exists (the analysis image)."""

    def setUp(self):
        if not mgcv_available():
            self.skipTest("mgcv unavailable locally; run inside the analysis image")

    def test_missing_reference_scale_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_fixture(root / "f.tsv")
            done = run_model(root / "f.tsv", root / "out")
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("--reference-scale", done.stderr)

    def test_profiler_scale_fits_and_records_components(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_fixture(root / "f.tsv")
            out = root / "out"
            done = run_model(root / "f.tsv", out,
                             ["--reference-scale", "profiler_scale"])
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            manifest = (out / "reference_scale_manifest.tsv").read_text()
            self.assertIn("profiler_scale", manifest)
            self.assertRegex(manifest, r"omitted_components\t(none|categorical)")

    def test_mixed_reference_within_profiler_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_fixture(root / "f.tsv", mixed=True)
            done = run_model(root / "f.tsv", root / "out",
                             ["--reference-scale", "profiler_scale"])
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("multiple reference types", done.stderr)

    def test_wrong_reference_for_profiler_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_fixture(root / "f.tsv", metaphlan_reference="read_proportional")
            done = run_model(root / "f.tsv", root / "out",
                             ["--reference-scale", "profiler_scale"])
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("requires genome_equivalent", done.stderr)


if __name__ == "__main__":
    unittest.main()
