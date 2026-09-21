#!/usr/bin/env python3
"""Regression gates for community biomarker-fate physical-unit handling."""
import csv
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis_v2/scripts/evaluate_disease_biomarker_propagation.py"
FIELDS = ["cohort", "study", "analysis_population", "target_label",
          "assembly_arm", "profiler", "contrast", "dose_level",
          "spike_fraction_target", "feature", "effect", "p_value", "q_value",
          "include", "exclusion_reason"]


def write_tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


class CommunityExclusionTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.root = Path(self.tmp.name)
        self.calls = self.root / "calls.tsv"
        self.aliases = self.root / "aliases.csv"
        self.panel = self.root / "panel.tsv"
        self.aliases.write_text(
            "canonical,alias,tool\nSpecies A,A_ALIAS,metaphlan4\n"
            "Species B,B_ALIAS,metaphlan4\n", encoding="utf-8")
        self.panel.write_text(
            "label\ttaxon_name\nA\tSpecies A\nB\tSpecies B\n", encoding="utf-8")

    def tearDown(self): self.tmp.cleanup()

    def rows(self):
        rows=[]
        # The same physical model result is repeated for A and B. Species A is
        # reported by its profiler alias; Species B by its canonical name. Both
        # must be excluded before calculating non-implanted biomarker fate.
        features={"A_ALIAS": (.4, .01), "Species B": (-.3, .01),
                  "Stable": (.2, .01), "Noise": (.1, .2)}
        for level, fractions in (("baseline", {"A": 0, "B": 0}),
                                 ("dose_01", {"A": .01, "B": .02})):
            for label in ("A", "B"):
                for feature,(effect,q) in features.items():
                    rows.append(dict(
                        cohort="yachida", study="S", analysis_population="community",
                        target_label=label, assembly_arm="original",
                        profiler="metaphlan4", contrast="CRC_vs_Control",
                        dose_level=level, spike_fraction_target=fractions[label],
                        feature=feature, effect=effect, p_value=q*q, q_value=q,
                        include=1, exclusion_reason=""))
        return rows

    def invoke(self, rows, name="out"):
        write_tsv(self.calls, FIELDS, rows)
        return subprocess.run(
            ["python3", str(SCRIPT), "--calls", str(self.calls),
             "--aliases", str(self.aliases), "--spike-panel", str(self.panel),
             "--outdir", str(self.root / name), "--q-thresholds", "0.05"],
            text=True, capture_output=True)

    def test_one_physical_mixture_and_all_targets_excluded(self):
        done=self.invoke(self.rows()); self.assertEqual(done.returncode,0,done.stdout+done.stderr)
        out=self.root/"out"
        with (out/"disease_biomarker_propagation_metrics.tsv").open() as handle:
            metrics=list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(metrics),1)
        row=metrics[0]
        self.assertEqual(row["target_label"],"CRCpanel")
        self.assertAlmostEqual(float(row["spike_fraction_target"]),.015)
        self.assertAlmostEqual(float(row["spike_fraction_total"]),.03)
        self.assertEqual(row["baseline_bystander_biomarkers"],"1")
        self.assertEqual(row["retained_bystanders"],"1")
        self.assertEqual(row["target_significant"],"NA")
        with (out/"disease_biomarker_transition_ledger.tsv").open() as handle:
            ledger=list(csv.DictReader(handle, delimiter="\t"))
        roles={r["feature"]:r["feature_role"] for r in ledger}
        self.assertEqual(roles["A_ALIAS"],"implanted_target")
        self.assertEqual(roles["Species B"],"implanted_target")
        self.assertEqual(roles["Stable"],"bystander")
        self.assertTrue(all(abs(float(r["spike_fraction_total"])-.03)<1e-12
                            for r in ledger))
        with (out/"disease_biomarker_propagation_summary.tsv").open() as handle:
            summary={r["metric"]:r["value"] for r in csv.DictReader(handle,delimiter="\t")}
        self.assertEqual(summary["community_repeated_rows_collapsed"],"8")
        self.assertEqual(summary["community_target_exclusion"],"all_panel_members")
        self.assertEqual(summary["feature_identity_policy"],"canonicalize_before_exclusion")

    def test_missing_panel_member_fails(self):
        rows=[r for r in self.rows() if r["target_label"] != "B"]
        done=self.invoke(rows,"missing")
        self.assertNotEqual(done.returncode,0)
        self.assertIn("does not contain every implanted target",done.stdout+done.stderr)

    def test_equal_member_dose_keeps_plot_axis_per_target(self):
        rows=self.rows()
        for row in rows:
            if row["dose_level"] == "dose_01":
                row["spike_fraction_target"] = .01
        done=self.invoke(rows,"equal_dose")
        self.assertEqual(done.returncode,0,done.stdout+done.stderr)
        with (self.root/"equal_dose"/"disease_biomarker_propagation_metrics.tsv").open() as handle:
            result=list(csv.DictReader(handle,delimiter="\t"))
        self.assertEqual(len(result),1)
        self.assertAlmostEqual(float(result[0]["spike_fraction_target"]),.01)
        self.assertAlmostEqual(float(result[0]["spike_fraction_total"]),.02)

    def test_disagreeing_repeated_fit_fails(self):
        rows=self.rows()
        for row in rows:
            if row["target_label"] == "B" and row["feature"] == "Stable":
                row["effect"] = .25; break
        done=self.invoke(rows,"conflict")
        self.assertNotEqual(done.returncode,0)
        self.assertIn("repeated fits disagree on effect",done.stdout+done.stderr)


if __name__ == "__main__": unittest.main()
