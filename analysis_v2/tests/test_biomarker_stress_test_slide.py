#!/usr/bin/env python3
import csv
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis_v2/scripts/plot_biomarker_stress_test.py"


class StressTestSlide(unittest.TestCase):
    def fixture(self, root, omit=False):
        panel = root / "panel.tsv"
        panel.write_text("label\ttaxon_name\nA\tTaxon A\nB\tTaxon B\n", encoding="utf-8")
        rows = []
        for profiler in ("kraken2_bracken", "metaphlan4"):
            for target in ("A", "B"):
                if omit and profiler == "metaphlan4" and target == "B":
                    continue
                rows.append(dict(cohort="yachida", analysis_population="independent",
                                 assembly_arm="original", profiler=profiler,
                                 contrast="CRC_vs_Control", target_label=target,
                                 spike_fraction_target="0.001", q_threshold="0.05",
                                 feature="Candidate species", feature_role="bystander",
                                 baseline_called="1", dose_called="0" if target == "B" else "1",
                                 baseline_q_value="0.002", effect_sign_changed="0"))
        ledger = root / "ledger.tsv"
        with ledger.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        return panel, ledger

    def test_slide_and_counts(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            panel, ledger = self.fixture(root)
            done = subprocess.run(["python3", str(SCRIPT), "--ledger", str(ledger),
                                   "--panel", str(panel), "--outdir", str(root / "out")],
                                  capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            self.assertTrue((root / "out/SUCCESS").is_file())
            ET.parse(root / "out/biomarker_stress_test.svg")
            with (root / "out/biomarker_stress_test_cells.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({r["status"] for r in rows}, {"retained", "lost"})

    def test_incomplete_grid_fails_closed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            panel, ledger = self.fixture(root, omit=True)
            done = subprocess.run(["python3", str(SCRIPT), "--ledger", str(ledger),
                                   "--panel", str(panel), "--outdir", str(root / "out")],
                                  capture_output=True, text=True)
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("incomplete challenge grid", done.stderr)


if __name__ == "__main__":
    unittest.main()
