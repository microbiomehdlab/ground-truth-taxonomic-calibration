#!/usr/bin/env python3
from __future__ import annotations
import csv
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from prepare_sealed_metaphlan_baselines import prepare, digest


def tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


class SelectionTests(unittest.TestCase):
    def fixture(self, root, cohort):
        state, results = root / "state", root / "results"
        seal = state / "production_seal_v2"
        seal.mkdir(parents=True)
        profile = results / "Study" / "S1" / "profiles" / "baseline" / "S1" / "S1.metaphlan.tsv"
        profile.parent.mkdir(parents=True)
        profile.write_text("# raw fixture\n")
        tsv(seal / "production_manifest.tsv", ["sample_id", "study"], [["S1", "Study"]])
        tsv(seal / "sample_flow.tsv", ["sample_id", "status"], [["S1", "PASS"]])
        (seal / "SUCCESS").write_text("status\tPASS\ncohort\t" + cohort + "\nsamples\t1\n")
        (seal / "production_seal.sha256").write_text("".join(
            digest(seal / n) + "  " + n + "\n" for n in
            ("SUCCESS", "production_manifest.tsv", "sample_flow.tsv")))
        receipt = state / "samples" / "S1.retained_outputs.tsv"
        tsv(receipt, ["path", "sha256", "bytes"],
            [[str(profile), digest(profile), profile.stat().st_size]])
        return state, results, profile, receipt

    def test_common_selection_and_fail_closed(self):
        for cohort in ("yachida", "feng", "zeller"):
            for failure in (None, "count", "corrupt", "absent_receipt", "duplicate_profile", "seal"):
                with self.subTest(cohort=cohort, failure=failure), tempfile.TemporaryDirectory() as tmp:
                    root = Path(tmp)
                    state, results, profile, receipt = self.fixture(root, cohort)
                    out = root / "out"
                    expected = 1
                    if failure == "count": expected = 2
                    if failure == "corrupt": profile.write_text("changed bytes\n")
                    if failure == "absent_receipt": tsv(receipt, ["path", "sha256", "bytes"], [])
                    if failure == "duplicate_profile":
                        extra = profile.parent.parent / "S1.metaphlan.tsv"
                        extra.write_bytes(profile.read_bytes())
                    if failure == "seal": (state / "production_seal_v2" / "SUCCESS").write_text("bad\n")
                    if failure:
                        with self.assertRaises(ValueError): prepare(cohort, state, results, out, expected)
                        self.assertFalse(out.exists())
                    else:
                        prepare(cohort, state, results, out, expected)
                        with (out / "baseline_manifest.tsv").open() as handle:
                            rows = list(csv.DictReader(handle, delimiter="\t"))
                        self.assertEqual(len(rows), 1)
                        self.assertEqual(rows[0]["source_profile"], str(profile))
                        self.assertEqual(rows[0]["analysis_population"], "")
                        with self.assertRaises(ValueError): prepare(cohort, state, results, out, expected)


if __name__ == "__main__": unittest.main()
