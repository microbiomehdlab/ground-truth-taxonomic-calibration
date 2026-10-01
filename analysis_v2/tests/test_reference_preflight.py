#!/usr/bin/env python3
"""Reference gate regression tests; no cluster or model execution."""
from __future__ import annotations
import subprocess
import tempfile
import unittest
from pathlib import Path
from test_metaphlan_reference_runner import ROOT, HELPER, write, pair


class ReferencePreflightTests(unittest.TestCase):
    def test_maintained_runners_do_not_precreate_endpoints_and_gate_preflight(self):
        for name in ("run_cohort_definitive_analysis.sh",
                     "run_crc_cohort_definitive_analysis.sh",
                     "run_yachida_definitive_analysis.sh"):
            text = (ROOT / "analysis_v2" / name).read_text()
            mkdir_lines = [s for s in text.splitlines() if s.startswith("mkdir -p")]
            self.assertFalse(any("endpoints" in s for s in mkdir_lines), name)
            self.assertLess(text.index('if require_reference_inputs "$CANONICAL_INPUT"'),
                            text.index('"${PREFLIGHT_ONLY:-0}"'), name)

    def test_read_only_gate_validates_values_targets_and_population_coverage(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            canonical = root / "canonical.tsv"
            write(canonical, pair("metaphlan4", "1.0", "1.2"))
            targets, effective = root / "targets.tsv", root / "effective.tsv"
            valid_target = "target_label\tgenome_size_bp\nFnuc\t2180101\n"
            valid_effective = "cohort\tsample_id\teffective_genome_size_bp\nyachida\tS1\t3477000\n"
            cases = [
                (valid_target, valid_effective, True),
                (valid_target.replace("2180101", "inf"), valid_effective, False),
                (valid_target.replace("2180101", "0"), valid_effective, False),
                (valid_target.replace("Fnuc", "OTHER"), valid_effective, False),
                (valid_target, valid_effective + "yachida\tS1\t3477000\n", False),
                (valid_target, "cohort\tsample_id\tanalysis_population\teffective_genome_size_bp\nyachida\tS1\tWRONG\t3477000\n", False),
            ]
            for target_text, effective_text, ok in cases:
                targets.write_text(target_text)
                effective.write_text(effective_text)
                before = {p.name: p.read_bytes() for p in root.iterdir()}
                result = subprocess.run(
                    ["bash", "-c", 'set -euo pipefail; source "$1"; '
                     'TARGET_GENOME_SIZES="$3"; EFFECTIVE_GENOME_SIZES="$4"; '
                     'require_reference_inputs "$2"', "test", str(HELPER),
                     str(canonical), str(targets), str(effective)],
                    cwd=ROOT, text=True, capture_output=True)
                self.assertEqual(result.returncode == 0, ok, result.stderr)
                self.assertEqual(before, {p.name: p.read_bytes() for p in root.iterdir()})


if __name__ == "__main__":
    unittest.main()
