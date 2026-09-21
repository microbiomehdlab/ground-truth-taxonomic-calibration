#!/usr/bin/env python3
"""Gates for partial-snapshot response and biomarker figure runners."""
import csv
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESPONSE = ROOT / "analysis_v2/run_partial_snapshot_response_figures.sh"
BIOMARKER = ROOT / "analysis_v2/run_partial_snapshot_biomarker_figures.sh"
MATRIX = ROOT / "analysis_v2/scripts/build_partial_replication_matrix.py"


class PartialSnapshotFigures(unittest.TestCase):
    def test_shell_interfaces_fail_closed_without_inputs(self):
        for script in (RESPONSE, BIOMARKER):
            syntax = subprocess.run(["bash", "-n", str(script)], capture_output=True, text=True)
            self.assertEqual(syntax.returncode, 0, syntax.stderr)
            missing = subprocess.run(["env", "-i", "PATH=/usr/bin:/bin", "bash", str(script)],
                                     capture_output=True, text=True)
            self.assertNotEqual(missing.returncode, 0)
        text = RESPONSE.read_text()
        self.assertIn("--metaphlan-reference genome_equivalent", text)
        self.assertGreaterEqual(text.count("--reference-scale profiler_scale"), 2)
        self.assertIn("--cohort-validation require_holdout", text)
        self.assertIn("--report-status DEVELOPMENT_ONLY", text)

    def fixture(self, path, conflict=False):
        rows = []
        for profiler in ("kraken2_bracken", "metaphlan4"):
            for cohort in ("feng", "yachida", "zeller"):
                for feature, q, effect in (("A", .01, 1), ("B", .02, -1), ("C", .5, 1)):
                    if feature == "B" and cohort == "zeller":
                        q = .5
                    rows.append(dict(cohort=cohort, analysis_population="community",
                        assembly_arm="original", profiler=profiler, feature=feature,
                        contrast="CRC_vs_Control", spike_fraction_target=0,
                        effect=effect, q_value=q, include=1,
                        model_spec="primary_age_sex", n_samples=20))
        if conflict:
            row = rows[0].copy()
            row["effect"] = 2
            rows.append(row)
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def test_replication_denominator_and_svg(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.fixture(root / "calls.tsv")
            run = subprocess.run(["python3", str(MATRIX), "--calls", str(root / "calls.tsv"),
                                  "--outdir", str(root / "matrix")], capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            with (root / "matrix/directional_replication_matrix.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 18)
            feng_to_zeller = next(r for r in rows if r["profiler"] == "metaphlan4" and
                                  r["discovery"] == "feng" and r["validation"] == "zeller")
            self.assertEqual((feng_to_zeller["discovered"], feng_to_zeller["evaluable"],
                              feng_to_zeller["replicated"]), ("2", "2", "1"))
            ET.parse(root / "matrix/directional_replication_matrix.svg")

    def test_conflicting_repeated_baseline_fit_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.fixture(root / "calls.tsv", conflict=True)
            run = subprocess.run(["python3", str(MATRIX), "--calls", str(root / "calls.tsv"),
                                  "--outdir", str(root / "matrix")], capture_output=True, text=True)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("contradictory repeated baseline fit", run.stderr)


if __name__ == "__main__":
    unittest.main()
