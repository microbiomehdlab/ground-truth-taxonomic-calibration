#!/usr/bin/env python3
from __future__ import annotations
import tempfile
import unittest
from pathlib import Path
from test_sealed_metaphlan_baselines import SelectionTests, tsv, prepare
from check_geff_canonical_baselines import check


class BaselineIdentityTests(unittest.TestCase):
    def test_matching_missing_wrong_and_changed_baselines(self):
        for mode in ("valid", "wrong_path", "wrong_cohort", "missing", "changed"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                state, results, profile, receipt = SelectionTests().fixture(root, "feng")
                selection = root / "selection"
                prepare("feng", state, results, selection, 1)
                canonical = root / "canonical.tsv"
                rows = [] if mode == "missing" else [[
                    "zeller" if mode == "wrong_cohort" else "feng", "metaphlan4", "0", "1", "S1",
                    str(root / "old.tsv") if mode == "wrong_path" else str(profile)]]
                tsv(canonical, ["cohort", "profiler", "spike_fraction_total", "include", "sample_id", "source_profile"], rows)
                if mode == "changed": profile.write_text("changed\n")
                if mode == "valid": check(canonical, selection, "feng", 1)
                else:
                    with self.assertRaises(ValueError): check(canonical, selection, "feng", 1)


if __name__ == "__main__": unittest.main()
