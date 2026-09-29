#!/usr/bin/env python3
"""Fail-closed smoke tests for the common three-cohort entry point."""
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "analysis_v2/run_cohort_definitive_analysis.sh"


class SharedRunnerTests(unittest.TestCase):
    def test_each_cohort_requires_its_unified_seal_before_writing(self):
        for cohort in ("yachida", "feng", "zeller"):
            with self.subTest(cohort=cohort), tempfile.TemporaryDirectory() as temp:
                base = Path(temp)
                envfile = base / "cohort.env"
                envfile.write_text(
                    f"YACHIDA_STATE_DIR={base / 'state'}\n"
                    f"CRC_STATE_DIR={base / 'state'}\n"
                    f"PERSISTENT_RESULTS_ROOT={base / 'results'}\n",
                    encoding="utf-8",
                )
                output = base / "output"
                env = os.environ.copy()
                env.update(
                    COHORT=cohort,
                    COHORT_ENV=str(envfile),
                    ANALYSIS_SIF=str(base / "analysis.sif"),
                    ASSEMBLY_SENSITIVITY_ROOT=str(base / "assembly"),
                    RUN_ROOT=str(output),
                    PREFLIGHT_ONLY="1",
                )
                result = subprocess.run(
                    ["bash", str(RUNNER)], cwd=ROOT, env=env,
                    text=True, capture_output=True, check=False,
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("Missing unified upstream seal", result.stderr)
                self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
