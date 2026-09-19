#!/usr/bin/env python3
"""Runner interfaces must state the reference scale explicitly.

Static checks that every maintained caller supplies `--reference-scale`, plus
fixture checks of the DiD runner's REFERENCE_SCALE / RESPONSE_TABLE contract.
"""

from __future__ import annotations

import os
import re
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
V2 = ROOT / "analysis_v2"
DID_RUNNER = V2 / "run_ideal_counterfactual_did.sh"

# Maintained runners that invoke the continuous model on a primary endpoint set.
CONTINUOUS_CALLERS = ["run_yachida_definitive_analysis.sh",
                      "run_crc_cohort_definitive_analysis.sh"]


def shell_files():
    return sorted(p for p in V2.glob("*.sh"))


class StaticInterfaceTest(unittest.TestCase):

    def test_every_continuous_caller_selects_profiler_scale(self):
        for name in CONTINUOUS_CALLERS:
            text = (V2 / name).read_text()
            self.assertIn("fit_continuous_dose_response.R", text, name)
            # The flag must appear inside the same backslash-continued command.
            block = re.search(
                r"fit_continuous_dose_response\.R((?:.*\\\n)*.*)", text)
            self.assertIsNotNone(block, name)
            self.assertIn("--reference-scale profiler_scale", block.group(1),
                          f"{name} does not select the primary scale")

    def test_no_shell_runner_invokes_continuous_model_without_scale(self):
        for path in shell_files():
            text = path.read_text()
            if "fit_continuous_dose_response.R" not in text:
                continue
            block = re.search(
                r"fit_continuous_dose_response\.R((?:.*\\\n)*.*)", text)
            self.assertIn("--reference-scale", block.group(1),
                          f"{path.name} omits --reference-scale")

    def test_did_runner_passes_reference_args(self):
        text = DID_RUNNER.read_text()
        self.assertIn('REFERENCE_SCALE:?', text, "REFERENCE_SCALE must be required")
        self.assertIn('RESPONSE_TABLE:?', text,
                      "RESPONSE_TABLE must be required for profiler_scale")
        self.assertIn('"${reference_args[@]}"', text,
                      "reference args must reach the R script")

    def test_no_runner_silently_defaults_to_read_proportional(self):
        for path in shell_files():
            text = path.read_text()
            if "--reference-scale" not in text:
                continue
            self.assertNotIn("--reference-scale read_proportional\n",
                             text.replace('"$REFERENCE_SCALE"', ""),
                             f"{path.name} hard-codes the sensitivity scale")


class DidRunnerContractTest(unittest.TestCase):
    """Exercise the guard clauses without reaching apptainer."""

    def _run(self, env_extra):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("disease", "endpoints.tsv", "manifest.tsv", "abundance.tsv"):
                (root / name).write_text("x", encoding="utf-8")
            env = {
                "PATH": os.environ["PATH"],
                "DISEASE_RUN": str(root / "disease"),
                "PAIRED_ENDPOINTS": str(root / "endpoints.tsv"),
                "ANALYSIS_SIF": str(root / "image.sif"),
                "OUTDIR": str(root / "out"),
                "ANALYSIS_STATUS": "DEVELOPMENT_ONLY",
                "PROFILE_MANIFEST": str(root / "manifest.tsv"),
                "ABUNDANCE_LONG": str(root / "abundance.tsv"),
            }
            env.update({k: v.replace("<ROOT>", str(root))
                        for k, v in env_extra.items()})
            return subprocess.run(["bash", str(DID_RUNNER)], text=True,
                                  capture_output=True, env=env)

    def test_missing_reference_scale_fails(self):
        done = self._run({})
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("REFERENCE_SCALE", done.stderr)

    def test_invalid_reference_scale_fails(self):
        done = self._run({"REFERENCE_SCALE": "genome_equivalent"})
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("must be profiler_scale or read_proportional", done.stderr)

    def test_profiler_scale_without_response_table_fails(self):
        done = self._run({"REFERENCE_SCALE": "profiler_scale"})
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("RESPONSE_TABLE", done.stderr)

    def test_profiler_scale_with_empty_response_table_fails(self):
        with tempfile.NamedTemporaryFile(suffix=".parquet") as empty:
            done = self._run({"REFERENCE_SCALE": "profiler_scale",
                              "RESPONSE_TABLE": empty.name})
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("missing or empty", done.stderr)

    def test_read_proportional_does_not_require_response_table(self):
        done = self._run({"REFERENCE_SCALE": "read_proportional"})
        # It proceeds past the guards and fails later on the absent image,
        # which proves no RESPONSE_TABLE requirement was applied.
        self.assertNotIn("RESPONSE_TABLE", done.stderr)


if __name__ == "__main__":
    unittest.main()
