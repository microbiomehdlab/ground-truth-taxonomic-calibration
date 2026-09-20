#!/usr/bin/env python3
"""Resume safety: a SUCCESS marker is not proof that a stage may be reused.

The completed Yachida development run holds a `reference_comparison` produced
before `compare_target_recovery_references.py` emitted the explicit per-arm
reference-provenance columns. Reusing it would carry an old-schema table into
the genome-size audit, which would then fail. These tests build temporary run
roots and drive the real shell gate (`analysis_v2/lib/stage_compatibility.sh`)
and the real schema checker, never a reimplementation of either.
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LIBRARY = ROOT / "analysis_v2/lib/stage_compatibility.sh"
CHECKER = ROOT / "analysis_v2/scripts/check_stage_schema.py"
COMPARATOR = ROOT / "analysis_v2/scripts/compare_target_recovery_references.py"
AUDIT = ROOT / "analysis_v2/scripts/audit_metaphlan_genome_size_residual.py"
sys.path.insert(0, str(Path(__file__).resolve().parent))

import test_metaphlan_genome_size_residual_audit as AUDITFIX  # noqa: E402

COMPARISON_STAGE = "reference_comparison"
AUDIT_STAGE = "metaphlan_genome_size_residual_audit"

# The eight explicit provenance columns a reusable comparison must carry.
NEW_COLUMNS = [
    "primary_row_reference_type", "sensitivity_row_reference_type",
    "primary_selected_reference_type", "sensitivity_selected_reference_type",
    "primary_reference_scale", "sensitivity_reference_scale",
    "primary_selected_estimand", "sensitivity_selected_estimand",
]
NEW_METRICS = ["primary_selection", "sensitivity_selection",
               "selected_estimand_contract"]


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def gate(stage, directory, run_root, stamp="20260921T000000Z", label="stage"):
    """Invoke the real shell gate exactly as the runner does."""
    script = (
        f'set -euo pipefail\n'
        f'source "{LIBRARY}"\n'
        f'if reuse_or_quarantine {stage} "{directory}" "{run_root}" '
        f'"{stamp}" "{label}"; then echo "DECISION=REUSE"; '
        f'else echo "DECISION=REGENERATE"; fi\n')
    done = subprocess.run(["bash", "-c", script], text=True, capture_output=True,
                          cwd=str(ROOT))
    assert done.returncode == 0, done.stdout + done.stderr
    return done.stdout


class SchemaCheckerTest(unittest.TestCase):
    """The header test must be exact, not a substring or grep match."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def check(self, stage, directory):
        return subprocess.run(
            [sys.executable, str(CHECKER), "--stage", stage,
             "--directory", str(directory)], text=True, capture_output=True)

    def comparison_dir(self, columns=NEW_COLUMNS, metrics=NEW_METRICS,
                       name="reference_comparison"):
        directory = self.root / name
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "SUCCESS").write_text("status\tPASS\n", encoding="utf-8")
        write_tsv(directory / "target_recovery_reference_comparison.tsv",
                  ["cohort", "profiler", *columns],
                  [{"cohort": "yachida", "profiler": "metaphlan4",
                    **{c: "x" for c in columns}}])
        write_tsv(directory / "reference_comparison_validation.tsv",
                  ["metric", "value"],
                  [{"metric": m, "value": "v"} for m in metrics])
        return directory

    def test_current_schema_is_accepted(self):
        done = self.check(COMPARISON_STAGE, self.comparison_dir())
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)

    def test_legacy_alias_does_not_satisfy_the_explicit_column(self):
        """`primary_reference_type` must not pass for
        `primary_row_reference_type`: exactly the substring trap."""
        legacy = ["primary_reference_type", "sensitivity_reference_type"]
        done = self.check(COMPARISON_STAGE, self.comparison_dir(columns=legacy))
        self.assertEqual(done.returncode, 1)
        self.assertIn("primary_row_reference_type", done.stdout)
        self.assertIn("predates the explicit reference-provenance schema",
                      done.stdout)

    def test_each_required_column_is_individually_required(self):
        for column in NEW_COLUMNS:
            remaining = [c for c in NEW_COLUMNS if c != column]
            done = self.check(COMPARISON_STAGE,
                              self.comparison_dir(columns=remaining,
                                                  name=f"c_{column}"))
            self.assertEqual(done.returncode, 1, column)
            self.assertIn(column, done.stdout)

    def test_each_required_metric_is_individually_required(self):
        for metric in NEW_METRICS:
            remaining = [m for m in NEW_METRICS if m != metric]
            done = self.check(COMPARISON_STAGE,
                              self.comparison_dir(metrics=remaining,
                                                  name=f"m_{metric}"))
            self.assertEqual(done.returncode, 1, metric)
            self.assertIn(metric, done.stdout)
            self.assertIn("predates the current selection contract", done.stdout)

    def test_metric_present_but_blank_is_rejected(self):
        directory = self.comparison_dir(name="blankmetric")
        write_tsv(directory / "reference_comparison_validation.tsv",
                  ["metric", "value"],
                  [{"metric": m, "value": "" if m == "primary_selection" else "v"}
                   for m in NEW_METRICS])
        done = self.check(COMPARISON_STAGE, directory)
        self.assertEqual(done.returncode, 1)
        self.assertIn("primary_selection", done.stdout)

    def test_missing_or_empty_required_file_is_rejected(self):
        for name in ("SUCCESS", "target_recovery_reference_comparison.tsv",
                     "reference_comparison_validation.tsv"):
            directory = self.comparison_dir(name=f"missing_{name}")
            (directory / name).unlink()
            done = self.check(COMPARISON_STAGE, directory)
            self.assertEqual(done.returncode, 1, name)
            self.assertIn(f"missing {name}", done.stdout)

            directory = self.comparison_dir(name=f"empty_{name}")
            (directory / name).write_text("", encoding="utf-8")
            done = self.check(COMPARISON_STAGE, directory)
            self.assertEqual(done.returncode, 1, name)
            self.assertIn(f"empty {name}", done.stdout)


class ComparisonResumeTest(unittest.TestCase):
    """Old-schema comparison is quarantined; current-schema one is reused."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.run_root = Path(self._tmp.name) / "run"
        self.run_root.mkdir(parents=True)
        self.arms = self._write_recovery_arms()

    def tearDown(self):
        self._tmp.cleanup()

    def _write_recovery_arms(self):
        """Two genuine target-recovery arms, as stages 4 and 5 leave them."""
        import test_target_recovery_reference_comparison as CMPFIX

        arms = {}
        for name, selection in (("target_recovery", CMPFIX.PRIMARY_SELECTION),
                                ("target_recovery_read_sensitivity",
                                 CMPFIX.SENSITIVITY_SELECTION)):
            directory = self.run_root / name
            directory.mkdir(parents=True)
            CMPFIX.write(directory / "target_recovery_observations.tsv", [
                CMPFIX.row("kraken2_bracken", "b1", 1.0, 0.0, "Good",
                           "read_proportional", selection),
                CMPFIX.row("metaphlan4", "m1", 1.05, 0.05, "Good",
                           "genome_equivalent", selection),
            ])
            (directory / "SUCCESS").write_text("status\tPASS\n", encoding="utf-8")
            arms[name] = directory
        return arms

    def _arm_state(self):
        return {name: sorted(
            (p.name, p.stat().st_size, p.stat().st_mtime_ns)
            for p in directory.iterdir())
            for name, directory in self.arms.items()}

    def generate(self):
        """Run the real comparator, as the runner's stage 6 does."""
        done = subprocess.run(
            [sys.executable, str(COMPARATOR),
             "--primary",
             str(self.arms["target_recovery"] / "target_recovery_observations.tsv"),
             "--sensitivity",
             str(self.arms["target_recovery_read_sensitivity"]
                 / "target_recovery_observations.tsv"),
             "--outdir", str(self.run_root / "reference_comparison")],
            text=True, capture_output=True)
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        return self.run_root / "reference_comparison"

    def write_old_schema(self, drop_metrics=False, drop_columns=True):
        """A completed pre-20-September comparison: SUCCESS, old columns."""
        directory = self.run_root / "reference_comparison"
        directory.mkdir(parents=True)
        (directory / "SUCCESS").write_text(
            "status\tPASS\npaired_observations\t2\n", encoding="utf-8")
        columns = (["primary_reference_type", "sensitivity_reference_type"]
                   if drop_columns else NEW_COLUMNS)
        write_tsv(directory / "target_recovery_reference_comparison.tsv",
                  ["cohort", "profiler", *columns],
                  [{"cohort": "yachida", "profiler": "metaphlan4",
                    **{c: "genome_equivalent" for c in columns}}])
        metrics = ["paired_observations", "bracken_identical"]
        if not drop_metrics:
            metrics += NEW_METRICS
        write_tsv(directory / "reference_comparison_validation.tsv",
                  ["metric", "value"],
                  [{"metric": m, "value": "PASS"} for m in metrics])
        return directory

    # ------------------------------------------------------------------ tests
    def test_old_schema_comparison_is_quarantined_and_regenerated(self):
        stale = self.write_old_schema()
        marker = stale / "target_recovery_reference_comparison.tsv"
        before = marker.read_text()
        output = gate(COMPARISON_STAGE, stale, self.run_root,
                      label="reference comparison")
        self.assertIn("DECISION=REGENERATE", output)
        self.assertIn("stale or incompatible reference comparison", output)
        self.assertIn("primary_row_reference_type", output)
        self.assertFalse(stale.exists(), "the stale directory must be moved away")
        quarantined = list((self.run_root / "failed_attempts").iterdir())
        self.assertEqual(len(quarantined), 1)
        self.assertTrue(quarantined[0].name.startswith("reference_comparison_"))
        # Nothing is destroyed: the old evidence is kept intact.
        self.assertEqual(
            (quarantined[0] / "target_recovery_reference_comparison.tsv").read_text(),
            before)
        # And regeneration then produces a compatible directory.
        fresh = self.generate()
        with (fresh / "target_recovery_reference_comparison.tsv").open(
                newline="", encoding="utf-8") as handle:
            header = handle.readline().rstrip("\n").split("\t")
        for column in NEW_COLUMNS:
            self.assertIn(column, header)
        self.assertIn("DECISION=REUSE",
                      gate(COMPARISON_STAGE, fresh, self.run_root))

    def test_current_schema_comparison_is_reused(self):
        fresh = self.generate()
        listing = sorted(p.name for p in fresh.iterdir())
        output = gate(COMPARISON_STAGE, fresh, self.run_root,
                      label="reference comparison")
        self.assertIn("DECISION=REUSE", output)
        self.assertIn("[REUSE] reference comparison", output)
        self.assertTrue(fresh.exists())
        self.assertEqual(listing, sorted(p.name for p in fresh.iterdir()))
        self.assertFalse((self.run_root / "failed_attempts").exists())

    def test_missing_validation_contract_metrics_force_regeneration(self):
        stale = self.write_old_schema(drop_metrics=True, drop_columns=False)
        output = gate(COMPARISON_STAGE, stale, self.run_root)
        self.assertIn("DECISION=REGENERATE", output)
        self.assertIn("primary_selection", output)
        self.assertIn("predates the current selection contract", output)
        self.assertFalse(stale.exists())

    def test_recovery_arms_survive_the_migration_untouched(self):
        before = self._arm_state()
        stale = self.write_old_schema()
        gate(COMPARISON_STAGE, stale, self.run_root)
        self.generate()
        self.assertEqual(before, self._arm_state(),
                         "the target-recovery arms must not be regenerated, "
                         "modified or deleted by the comparison migration")
        for directory in self.arms.values():
            self.assertTrue((directory / "SUCCESS").is_file())
            self.assertTrue(
                (directory / "target_recovery_observations.tsv").is_file())

    def test_absent_directory_simply_regenerates(self):
        output = gate(COMPARISON_STAGE, self.run_root / "reference_comparison",
                      self.run_root)
        self.assertIn("DECISION=REGENERATE", output)
        self.assertNotIn("QUARANTINE", output)
        self.assertFalse((self.run_root / "failed_attempts").exists())

    def test_incomplete_directory_without_success_is_quarantined(self):
        partial = self.run_root / "reference_comparison"
        partial.mkdir(parents=True)
        (partial / "partial.tsv").write_text("x\n", encoding="utf-8")
        output = gate(COMPARISON_STAGE, partial, self.run_root)
        self.assertIn("DECISION=REGENERATE", output)
        self.assertIn("missing SUCCESS", output)
        self.assertFalse(partial.exists())


class AuditResumeTest(unittest.TestCase):
    """The same gate protects a completed genome-size audit directory."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.run_root = self.root / "run"
        self.run_root.mkdir(parents=True)
        self.fixture_dir = self.root / "inputs"
        self.fixture_dir.mkdir()
        self.fixture = AUDITFIX.Fixture(self.fixture_dir)

    def tearDown(self):
        self._tmp.cleanup()

    def generate(self):
        out = self.run_root / "metaphlan_genome_size_residual_audit"
        done = self.fixture.run(out, replicates=200)
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        return out

    def write_old_schema(self):
        """A completed audit from before the selected-provenance columns."""
        directory = self.run_root / "metaphlan_genome_size_residual_audit"
        directory.mkdir(parents=True)
        for name in ("SUCCESS", "metaphlan_genome_size_regression.tsv",
                     "metaphlan_genome_size_paired_change.tsv",
                     "metaphlan_genome_size_audit.sha256",
                     "DEVELOPMENT_ONLY.txt"):
            (directory / name).write_text("placeholder\n", encoding="utf-8")
        write_tsv(directory / "metaphlan_genome_size_target_level.tsv",
                  ["scope", "arm", "target_label", "reference_type"],
                  [{"scope": "pooled", "arm": "genome_equivalent_primary",
                    "target_label": "Fnuc", "reference_type": "genome_equivalent"}])
        write_tsv(directory / "metaphlan_genome_size_audit_validation.tsv",
                  ["metric", "value"],
                  [{"metric": "audit_status", "value": "PASS"},
                   {"metric": "bootstrap_seed", "value": "20260920"}])
        return directory

    def test_old_schema_audit_is_quarantined_and_regenerated(self):
        stale = self.write_old_schema()
        output = gate(AUDIT_STAGE, stale, self.run_root,
                      label="genome-size residual audit")
        self.assertIn("DECISION=REGENERATE", output)
        self.assertIn("stale or incompatible genome-size residual audit", output)
        self.assertIn("row_reference_type", output)
        self.assertFalse(stale.exists())
        quarantined = list((self.run_root / "failed_attempts").iterdir())
        self.assertEqual(len(quarantined), 1)
        self.assertTrue(quarantined[0].name.startswith(
            "metaphlan_genome_size_residual_audit_"))
        fresh = self.generate()
        self.assertIn("DECISION=REUSE", gate(AUDIT_STAGE, fresh, self.run_root))

    def test_current_schema_audit_is_reused(self):
        fresh = self.generate()
        listing = sorted(p.name for p in fresh.iterdir())
        output = gate(AUDIT_STAGE, fresh, self.run_root,
                      label="genome-size residual audit")
        self.assertIn("DECISION=REUSE", output)
        self.assertTrue(fresh.exists())
        self.assertEqual(listing, sorted(p.name for p in fresh.iterdir()))
        self.assertFalse((self.run_root / "failed_attempts").exists())

    def test_missing_audit_validation_metric_forces_regeneration(self):
        fresh = self.generate()
        rows = [r for r in read_tsv(
            fresh / "metaphlan_genome_size_audit_validation.tsv")
            if r["metric"] != "observed_cohorts"]
        write_tsv(fresh / "metaphlan_genome_size_audit_validation.tsv",
                  ["metric", "value"], rows)
        output = gate(AUDIT_STAGE, fresh, self.run_root)
        self.assertIn("DECISION=REGENERATE", output)
        self.assertIn("observed_cohorts", output)

    def test_each_required_audit_output_is_required(self):
        for name in ("SUCCESS", "metaphlan_genome_size_target_level.tsv",
                     "metaphlan_genome_size_regression.tsv",
                     "metaphlan_genome_size_paired_change.tsv",
                     "metaphlan_genome_size_audit_validation.tsv",
                     "metaphlan_genome_size_audit.sha256",
                     "DEVELOPMENT_ONLY.txt"):
            fresh = self.generate()
            (fresh / name).unlink()
            output = gate(AUDIT_STAGE, fresh, self.run_root)
            self.assertIn("DECISION=REGENERATE", output, name)
            self.assertIn(f"missing {name}", output)
            quarantine = self.run_root / "failed_attempts"
            for item in quarantine.iterdir():
                subprocess.run(["rm", "-rf", str(item)], check=True)


class EstimandContradictionTest(unittest.TestCase):
    """Requirement 10: contradictory provenance must fail, not be ignored."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def fixture(self, name):
        directory = self.root / name
        directory.mkdir(parents=True)
        return AUDITFIX.Fixture(directory)

    def test_contradictory_selected_estimand_fails(self):
        for column, wrong in (("sensitivity_selected_estimand", "genome_equivalent"),
                              ("primary_selected_estimand", "read_proportional")):
            fixture = self.fixture(f"f_{column}")
            for row in fixture.rows:
                if row["profiler"] == "metaphlan4":
                    row[column] = wrong
            done = fixture.run(self.root / f"out_{column}", replicates=50)
            self.assertNotEqual(done.returncode, 0,
                                f"audit ignored a contradictory {column}")
            message = done.stdout + done.stderr
            self.assertIn(column, message)
            self.assertIn("contradicts itself", message)

    def test_blank_selected_estimand_fails(self):
        fixture = self.fixture("blank")
        fixture.rows[0]["primary_selected_estimand"] = ""
        done = fixture.run(self.root / "out_blank", replicates=50)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("missing primary_selected_estimand",
                      done.stdout + done.stderr)

    def test_consistent_estimands_pass(self):
        fixture = self.fixture("clean")
        done = fixture.run(self.root / "out_clean", replicates=50)
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)


if __name__ == "__main__":
    unittest.main()
