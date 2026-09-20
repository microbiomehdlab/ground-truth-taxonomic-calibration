#!/usr/bin/env python3
"""Numerical tests for the MetaPhlAn genome-size residual audit.

The real script is always executed as a subprocess against a deterministic
synthetic fixture whose target-level answers are known in closed form:

    sensitivity target median log2(observed/expected) = C - log2(G_t)   slope -1
    primary     target median log2(observed/expected) = constant        slope  0

Nothing here imports the audit's own helpers, so the assertions are independent
of its implementation.
"""

from __future__ import annotations

import csv
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis_v2/scripts/audit_metaphlan_genome_size_residual.py"

# Measured spike-FASTA lengths (CLAUDE.md audited table). Ten distinct sizes.
GENOME_SIZES = {
    "Bfrag": 5_241_700, "Csym": 4_823_675, "Dpne": 1_247_407,
    "Fnuc": 2_180_101, "Hhat": 5_697_783, "Pmic": 1_661_863,
    "Pana": 2_256_756, "Psto": 1_988_044, "Porp": 2_199_119,
    "Pint": 2_871_118,
}
LABELS = tuple(sorted(GENOME_SIZES))
# Canonical spike-panel taxon names.
TAXON = {
    "Bfrag": "Bacteroides fragilis",
    "Csym": "Clostridium symbiosum",
    "Dpne": "Dialister pneumosintes",
    "Fnuc": "Fusobacterium nucleatum subsp. nucleatum",
    "Hhat": "Hungatella hathewayi",
    "Pmic": "Parvimonas micra",
    "Pana": "Peptostreptococcus anaerobius",
    "Psto": "Peptostreptococcus stomatis",
    "Porp": "Porphyromonas asaccharolytica",
    "Pint": "Prevotella intermedia",
}
# The MetaPhlAn feature actually written by the upstream comparison. Fnuc is the
# one that genuinely differs from the canonical taxon name, which is why the
# mapping must go through the frozen alias table and not through string surgery.
METAPHLAN_ALIAS = dict(TAXON, Fnuc="Fusobacterium nucleatum")

COMPARISON_FIELDS = [
    "cohort", "study", "sample_id", "condition", "analysis_population",
    "profiler", "baseline_id", "target_feature", "dose_index",
    "nominal_target_fraction",
    "primary_row_reference_type", "sensitivity_row_reference_type",
    "primary_selected_reference_type", "sensitivity_selected_reference_type",
    "primary_reference_scale", "sensitivity_reference_scale",
    "primary_selected_estimand", "sensitivity_selected_estimand",
    "primary_reference_type", "sensitivity_reference_type",
    "primary_implanted_signal", "sensitivity_implanted_signal",
    "primary_response_signal", "sensitivity_response_signal",
    "primary_observed_over_expected", "sensitivity_observed_over_expected",
    "primary_absolute_relative_error", "sensitivity_absolute_relative_error",
    "absolute_relative_error_change_primary_minus_sensitivity",
    "primary_recovery_class", "sensitivity_recovery_class",
    "class_transition", "primary_class_improved", "primary_class_worsened",
]

# Intercepts chosen so every synthetic ratio is of realistic magnitude.
SENSITIVITY_INTERCEPT_INDEPENDENT = 21.7
SENSITIVITY_INTERCEPT_COMMUNITY = 21.3
PRIMARY_LOG2 = 0.1


def g(value: float) -> str:
    return format(value, ".17g")


def write_tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def ols(points):
    """Plain OLS, written independently of the script under test."""
    n = len(points)
    mx = sum(x for x, _ in points) / n
    my = sum(y for _, y in points) / n
    sxx = sum((x - mx) ** 2 for x, _ in points)
    sxy = sum((x - mx) * (y - my) for x, y in points)
    return sxy / sxx


def comparison_row(label, population, sample, sensitivity_log2, primary_log2,
                   profiler="metaphlan4", feature=None,
                   sensitivity_ratio=None, primary_ratio=None):
    s_ratio = 2.0 ** sensitivity_log2 if sensitivity_ratio is None else sensitivity_ratio
    p_ratio = 2.0 ** primary_log2 if primary_ratio is None else primary_ratio
    return {
        "cohort": "yachida", "study": "YachidaS_2019", "sample_id": sample,
        "condition": "CRC", "analysis_population": population,
        "profiler": profiler, "baseline_id": sample + "_base",
        "target_feature": feature or METAPHLAN_ALIAS[label],
        "dose_index": "4", "nominal_target_fraction": "0.005",
        # Row provenance stays genome_equivalent in BOTH arms: it records the
        # source profile, not the estimand. The selected fields are what
        # distinguish the arms.
        "primary_row_reference_type": "genome_equivalent",
        "sensitivity_row_reference_type": "genome_equivalent",
        "primary_selected_reference_type": "profiler_scale_primary",
        "sensitivity_selected_reference_type": "read_proportional",
        "primary_reference_scale": "profiler_scale",
        "sensitivity_reference_scale": "read_proportional",
        "primary_selected_estimand": "genome_equivalent",
        "sensitivity_selected_estimand": "read_proportional",
        "primary_reference_type": "genome_equivalent",
        "sensitivity_reference_type": "genome_equivalent",
        "primary_implanted_signal": "0.005",
        "sensitivity_implanted_signal": "0.005",
        "primary_response_signal": "0.005",
        "sensitivity_response_signal": "0.005",
        "primary_observed_over_expected": g(p_ratio) if isinstance(p_ratio, float) else str(p_ratio),
        "sensitivity_observed_over_expected": g(s_ratio) if isinstance(s_ratio, float) else str(s_ratio),
        "primary_absolute_relative_error": g(abs(p_ratio - 1.0)) if isinstance(p_ratio, float) else "0.1",
        "sensitivity_absolute_relative_error": g(abs(s_ratio - 1.0)) if isinstance(s_ratio, float) else "0.1",
        "absolute_relative_error_change_primary_minus_sensitivity": "0",
        "primary_recovery_class": "Good", "sensitivity_recovery_class": "Good",
        "class_transition": "Good -> Good",
        "primary_class_improved": "0", "primary_class_worsened": "0",
    }


class Fixture:
    """A complete, valid set of audit inputs that individual tests mutate."""

    def __init__(self, root, replication=None, residual=None,
                 community_replication=None):
        self.root = Path(root)
        self.replication = replication or {}
        self.community_replication = community_replication or {}
        self.residual = residual or {}
        self.rows = []
        for label in LABELS:
            x = math.log2(GENOME_SIZES[label])
            bump = self.residual.get(label, 0.0)
            for index in range(self.replication.get(label, 1)):
                self.rows.append(comparison_row(
                    label, "independent", f"ind_{label}_{index}",
                    SENSITIVITY_INTERCEPT_INDEPENDENT - x + bump, PRIMARY_LOG2))
            for index in range(self.community_replication.get(label, 1)):
                self.rows.append(comparison_row(
                    label, "community", f"com_{label}_{index}",
                    SENSITIVITY_INTERCEPT_COMMUNITY - x + bump, PRIMARY_LOG2))
        # Bracken must be present as an input gate and must never be regressed.
        # Its ROW provenance is read_proportional in both arms.
        for row in (comparison_row("Bfrag", "community", "brk_1", 9.0, 9.0,
                                   profiler="kraken2_bracken"),
                    comparison_row("Hhat", "independent", "brk_2", -9.0, -9.0,
                                   profiler="kraken2_bracken")):
            row["primary_row_reference_type"] = "read_proportional"
            row["sensitivity_row_reference_type"] = "read_proportional"
            row["primary_reference_type"] = "read_proportional"
            row["sensitivity_reference_type"] = "read_proportional"
            row["primary_selected_estimand"] = "read_proportional"
            self.rows.append(row)

        self.panel_rows = [
            {"label": label, "taxon_name": TAXON[label], "assembly": "GCF_X",
             "fasta": f"references/genomes/{label}.fa", "weight": "1", "url": ""}
            for label in LABELS]
        self.size_rows = [
            {"target_label": label, "target_feature": TAXON[label],
             "assembly_accession": "GCF_X",
             "fasta_path": f"/p/{label}.fa", "fasta_sha256": "0" * 64,
             "genome_size_bp": str(GENOME_SIZES[label])}
            for label in LABELS]
        self.alias_rows = (
            [{"canonical": TAXON[label], "alias": METAPHLAN_ALIAS[label],
              "tool": "metaphlan4", "spike_label": ""} for label in LABELS]
            + [{"canonical": TAXON[label], "alias": TAXON[label],
                "tool": "kraken2_bracken", "spike_label": ""} for label in LABELS])
        self.validation_rows = [
            {"metric": "paired_observations", "value": str(len(self.rows))},
            {"metric": "bracken_observations", "value": "2"},
            {"metric": "metaphlan_observations", "value": str(len(self.rows) - 2)},
            {"metric": "bracken_identical", "value": "PASS"},
            {"metric": "comparison_status", "value": "DEVELOPMENT_ONLY"},
        ]

    # -- expected target-level answers, computed here and not by the script --
    def expected_sensitivity_median(self, scope, label):
        x = math.log2(GENOME_SIZES[label])
        bump = self.residual.get(label, 0.0)
        independent = [SENSITIVITY_INTERCEPT_INDEPENDENT - x + bump] * \
            self.replication.get(label, 1)
        community = [SENSITIVITY_INTERCEPT_COMMUNITY - x + bump] * \
            self.community_replication.get(label, 1)
        values = {"independent": independent, "community": community,
                  "pooled": independent + community}[scope]
        ordered = sorted(values)
        middle = len(ordered) // 2
        if len(ordered) % 2:
            return ordered[middle]
        return (ordered[middle - 1] + ordered[middle]) / 2.0

    def paths(self):
        comparison = self.root / "target_recovery_reference_comparison.tsv"
        validation = self.root / "reference_comparison_validation.tsv"
        sizes = self.root / "target_genome_sizes.tsv"
        panel = self.root / "spike_panel.tsv"
        aliases = self.root / "spike_taxon_aliases.csv"
        write_tsv(comparison, COMPARISON_FIELDS, self.rows)
        write_tsv(validation, ["metric", "value"], self.validation_rows)
        write_tsv(sizes, ["target_label", "target_feature", "assembly_accession",
                          "fasta_path", "fasta_sha256", "genome_size_bp"],
                  self.size_rows)
        write_tsv(panel, ["label", "taxon_name", "assembly", "fasta", "weight",
                          "url"], self.panel_rows)
        with aliases.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=["canonical", "alias", "tool", "spike_label"],
                lineterminator="\n")
            writer.writeheader()
            writer.writerows(self.alias_rows)
        return comparison, validation, sizes, panel, aliases

    def run(self, outdir, replicates=200, extra=()):
        comparison, validation, sizes, panel, aliases = self.paths()
        return subprocess.run(
            [sys.executable, str(SCRIPT),
             "--comparison", str(comparison),
             "--comparison-validation", str(validation),
             "--target-genome-sizes", str(sizes),
             "--spike-panel", str(panel),
             "--feature-aliases", str(aliases),
             "--outdir", str(outdir),
             "--bootstrap-replicates", str(replicates), *extra],
            text=True, capture_output=True)


class AuditTestCase(unittest.TestCase):

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def fixture(self, **kwargs):
        directory = self.root / f"fix{len(list(self.root.iterdir()))}"
        directory.mkdir()
        return Fixture(directory, **kwargs)

    def succeed(self, fixture, name="out", replicates=200):
        out = self.root / name
        done = fixture.run(out, replicates=replicates)
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertTrue((out / "SUCCESS").is_file())
        return out

    def fail(self, fixture, needle, name="bad"):
        out = self.root / name
        done = fixture.run(out)
        self.assertNotEqual(done.returncode, 0,
                            f"audit accepted bad input (expected {needle!r})")
        message = done.stdout + done.stderr
        self.assertIn(needle, message)
        self.assertFalse((out / "SUCCESS").exists())
        return message

    def regression(self, out):
        return {(r["scope"], r["arm"]): r
                for r in read_tsv(out / "metaphlan_genome_size_regression.tsv")}


# ------------------------------------------------------------------ case A
class KnownSlopeTest(AuditTestCase):

    def test_slopes_match_the_closed_form(self):
        out = self.succeed(self.fixture())
        rows = self.regression(out)
        for scope in ("independent", "community", "pooled"):
            sensitivity = rows[(scope, "read_proportional_sensitivity")]
            primary = rows[(scope, "genome_equivalent_primary")]
            self.assertEqual(int(sensitivity["n_targets"]), 10)
            self.assertAlmostEqual(float(sensitivity["slope"]), -1.0, places=12)
            self.assertEqual(sensitivity["expected_slope"], "-1")
            self.assertAlmostEqual(
                float(sensitivity["absolute_distance_from_expected"]), 0.0,
                places=12)
            self.assertAlmostEqual(float(primary["slope"]), 0.0, places=12)
            self.assertEqual(primary["expected_slope"], "0")
            self.assertAlmostEqual(
                float(primary["absolute_distance_from_expected"]), 0.0, places=12)

    def test_sensitivity_fit_quality_and_intercepts(self):
        out = self.succeed(self.fixture())
        rows = self.regression(out)
        independent = rows[("independent", "read_proportional_sensitivity")]
        community = rows[("community", "read_proportional_sensitivity")]
        self.assertAlmostEqual(float(independent["r_squared"]), 1.0, places=12)
        self.assertAlmostEqual(float(independent["pearson_correlation"]), -1.0,
                               places=12)
        self.assertAlmostEqual(float(independent["intercept"]),
                               SENSITIVITY_INTERCEPT_INDEPENDENT, places=9)
        self.assertAlmostEqual(float(community["intercept"]),
                               SENSITIVITY_INTERCEPT_COMMUNITY, places=9)

    def test_zero_variance_primary_is_explicitly_undefined_not_zero(self):
        """A perfectly flat arm must not crash and must not fake an R^2."""
        out = self.succeed(self.fixture())
        rows = self.regression(out)
        for scope in ("independent", "community", "pooled"):
            primary = rows[(scope, "genome_equivalent_primary")]
            self.assertEqual(primary["r_squared"], "",
                             "zero-variance R^2 must be blank, not 0 or NaN")
            self.assertEqual(primary["pearson_correlation"], "")
            self.assertNotEqual(primary["slope"], "")
            self.assertNotEqual(primary["intercept"], "")

    def test_defined_variance_populates_r_squared_and_correlation(self):
        """The other side of the policy: real variance yields real statistics."""
        out = self.succeed(self.fixture(residual={"Hhat": 0.4, "Dpne": -0.3}))
        row = self.regression(out)[("independent", "read_proportional_sensitivity")]
        r_squared = float(row["r_squared"])
        correlation = float(row["pearson_correlation"])
        self.assertTrue(math.isfinite(r_squared) and math.isfinite(correlation))
        self.assertAlmostEqual(correlation ** 2, r_squared, places=12)
        self.assertLess(r_squared, 1.0)

    def test_target_level_medians_match_the_closed_form(self):
        fixture = self.fixture()
        out = self.succeed(fixture)
        rows = {(r["scope"], r["arm"], r["target_label"]): r
                for r in read_tsv(out / "metaphlan_genome_size_target_level.tsv")}
        for scope in ("independent", "community", "pooled"):
            for label in LABELS:
                row = rows[(scope, "read_proportional_sensitivity", label)]
                self.assertAlmostEqual(
                    float(row["median_log2_observed_over_expected"]),
                    fixture.expected_sensitivity_median(scope, label), places=12,
                    msg=f"{scope}/{label}")
                self.assertEqual(float(row["genome_size_bp"]),
                                 GENOME_SIZES[label])
                self.assertAlmostEqual(float(row["log2_genome_size_bp"]),
                                       math.log2(GENOME_SIZES[label]), places=12)
                self.assertEqual(row["target_feature"], METAPHLAN_ALIAS[label])

    def test_paired_change_table_is_consistent(self):
        fixture = self.fixture()
        out = self.succeed(fixture)
        rows = {(r["scope"], r["target_label"]): r
                for r in read_tsv(out / "metaphlan_genome_size_paired_change.tsv")}
        for scope in ("independent", "community", "pooled"):
            for label in LABELS:
                row = rows[(scope, label)]
                sensitivity = float(row["sensitivity_median_log2_ratio"])
                primary = float(row["primary_median_log2_ratio"])
                self.assertAlmostEqual(
                    sensitivity, fixture.expected_sensitivity_median(scope, label),
                    places=12)
                self.assertAlmostEqual(primary, PRIMARY_LOG2, places=12)
                self.assertAlmostEqual(
                    float(row["primary_minus_sensitivity_log2_ratio"]),
                    primary - sensitivity, places=12)


# ------------------------------------------------------------------ case B
class ReplicationTest(AuditTestCase):
    """The implanted taxon is the unit; observation rows are repeated measures."""

    def test_duplicating_rows_does_not_change_the_slope(self):
        plain = self.regression(self.succeed(self.fixture(), name="plain"))
        replicated = self.regression(self.succeed(
            self.fixture(replication={label: 3 for label in LABELS},
                         community_replication={label: 3 for label in LABELS}),
            name="replicated"))
        for key, row in plain.items():
            self.assertEqual(row["slope"], replicated[key]["slope"], key)
            self.assertEqual(row["intercept"], replicated[key]["intercept"], key)

    def test_unbalanced_replication_cannot_reweight_the_regression(self):
        """Row-level OLS would be dominated by the heavily replicated taxon."""
        residual = {"Hhat": 0.6}
        balanced = self.fixture(residual=residual)
        skewed = self.fixture(residual=residual, replication={"Hhat": 40},
                              community_replication={"Hhat": 40})
        balanced_slope = float(self.regression(
            self.succeed(balanced, name="balanced")
        )[("independent", "read_proportional_sensitivity")]["slope"])
        skewed_slope = float(self.regression(
            self.succeed(skewed, name="skewed")
        )[("independent", "read_proportional_sensitivity")]["slope"])
        # Target-level answer, unchanged by replication.
        target_points = [(math.log2(GENOME_SIZES[label]),
                          balanced.expected_sensitivity_median("independent", label))
                         for label in LABELS]
        expected = ols(target_points)
        self.assertAlmostEqual(balanced_slope, expected, places=12)
        self.assertAlmostEqual(skewed_slope, expected, places=12)
        # The observation-level slope the audit must NOT produce.
        observation_points = []
        for label in LABELS:
            x = math.log2(GENOME_SIZES[label])
            y = SENSITIVITY_INTERCEPT_INDEPENDENT - x + residual.get(label, 0.0)
            observation_points.extend([(x, y)] * (40 if label == "Hhat" else 1))
        pseudoreplicated = ols(observation_points)
        self.assertGreater(abs(pseudoreplicated - expected), 0.05,
                           "fixture cannot distinguish the two regression units")
        self.assertNotAlmostEqual(skewed_slope, pseudoreplicated, places=3)


# ------------------------------------------------------------------ case C
class ScopeSeparationTest(AuditTestCase):

    def test_populations_are_analyzed_separately(self):
        rows = self.regression(self.succeed(self.fixture()))
        independent = float(
            rows[("independent", "read_proportional_sensitivity")]["intercept"])
        community = float(
            rows[("community", "read_proportional_sensitivity")]["intercept"])
        self.assertAlmostEqual(independent, SENSITIVITY_INTERCEPT_INDEPENDENT,
                               places=9)
        self.assertAlmostEqual(community, SENSITIVITY_INTERCEPT_COMMUNITY, places=9)
        self.assertNotAlmostEqual(independent, community, places=3)

    def test_pooled_uses_observations_not_population_medians(self):
        """Two community observations per target outvote one independent one."""
        fixture = self.fixture(community_replication={label: 2 for label in LABELS})
        rows = self.regression(self.succeed(fixture))
        pooled = float(rows[("pooled", "read_proportional_sensitivity")]["intercept"])
        # Median of [independent, community, community] is the community value.
        self.assertAlmostEqual(pooled, SENSITIVITY_INTERCEPT_COMMUNITY, places=9)
        mean_of_population_medians = (SENSITIVITY_INTERCEPT_INDEPENDENT
                                      + SENSITIVITY_INTERCEPT_COMMUNITY) / 2.0
        self.assertNotAlmostEqual(pooled, mean_of_population_medians, places=3)

    def test_pooled_target_counts_are_the_sum_of_the_populations(self):
        fixture = self.fixture(replication={label: 2 for label in LABELS})
        out = self.succeed(fixture)
        rows = {(r["scope"], r["arm"], r["target_label"]): r
                for r in read_tsv(out / "metaphlan_genome_size_target_level.tsv")}
        for label in LABELS:
            key = lambda scope: rows[(scope, "genome_equivalent_primary", label)]
            self.assertEqual(int(key("independent")["eligible_observations"]), 2)
            self.assertEqual(int(key("community")["eligible_observations"]), 1)
            self.assertEqual(int(key("pooled")["eligible_observations"]), 3)


# ------------------------------------------------------------------ case D
class ExactMappingTest(AuditTestCase):

    def test_missing_target_fails(self):
        fixture = self.fixture()
        fixture.rows = [r for r in fixture.rows
                        if r["target_feature"] != METAPHLAN_ALIAS["Pmic"]]
        self.fail(fixture, "missing=['Pmic']")

    def test_unexpected_target_fails(self):
        fixture = self.fixture()
        fixture.rows.append(comparison_row(
            "Bfrag", "community", "extra_1", 1.0, 0.1,
            feature="Escherichia coli"))
        self.fail(fixture, "does not match exactly one implanted target")

    def test_duplicate_inconsistent_genome_size_fails(self):
        fixture = self.fixture()
        fixture.size_rows.append({
            "target_label": "Bfrag", "target_feature": TAXON["Bfrag"],
            "assembly_accession": "GCF_X", "fasta_path": "/p/Bfrag2.fa",
            "fasta_sha256": "1" * 64, "genome_size_bp": "4000000"})
        self.fail(fixture, "inconsistent duplicate genome sizes")

    def test_duplicate_consistent_genome_size_is_accepted(self):
        fixture = self.fixture()
        fixture.size_rows.append(dict(fixture.size_rows[0]))
        self.succeed(fixture)

    def test_ambiguous_alias_fails(self):
        """Two implanted taxa may never collapse onto one profiler feature."""
        fixture = self.fixture()
        for row in fixture.alias_rows:
            if row["tool"] == "metaphlan4" and row["canonical"] == TAXON["Psto"]:
                row["alias"] = METAPHLAN_ALIAS["Pana"]
        self.fail(fixture, "maps to more than one target label")

    def test_conflicting_alias_rows_fail(self):
        fixture = self.fixture()
        fixture.alias_rows.append({"canonical": TAXON["Bfrag"],
                                   "alias": "Bacteroides fragilis_A",
                                   "tool": "metaphlan4", "spike_label": ""})
        self.fail(fixture, "conflicting metaphlan4 aliases")

    def test_prefix_match_is_not_accepted(self):
        fixture = self.fixture()
        for row in fixture.rows:
            if row["target_feature"] == METAPHLAN_ALIAS["Bfrag"]:
                row["target_feature"] = "Bacteroides"
        self.fail(fixture, "does not match exactly one implanted target")

    def test_extended_name_is_not_accepted(self):
        fixture = self.fixture()
        for row in fixture.rows:
            if row["target_feature"] == METAPHLAN_ALIAS["Pmic"]:
                row["target_feature"] = METAPHLAN_ALIAS["Pmic"] + " subsp. x"
        self.fail(fixture, "does not match exactly one implanted target")

    def test_case_difference_is_not_accepted(self):
        fixture = self.fixture()
        for row in fixture.rows:
            if row["target_feature"] == METAPHLAN_ALIAS["Csym"]:
                row["target_feature"] = METAPHLAN_ALIAS["Csym"].upper()
        self.fail(fixture, "does not match exactly one implanted target")

    def test_one_label_with_two_features_fails(self):
        """Both the canonical name and the alias resolve to Fnuc; using both
        in one comparison would double-count the taxon."""
        fixture = self.fixture()
        fixture.rows.append(comparison_row(
            "Fnuc", "community", "dup_feature_1", 1.0, PRIMARY_LOG2,
            feature=TAXON["Fnuc"]))
        self.fail(fixture, "maps to multiple MetaPhlAn features")

    def test_genome_size_table_must_agree_with_the_panel(self):
        fixture = self.fixture()
        fixture.size_rows[0]["target_feature"] = "Something else"
        self.fail(fixture, "does not exactly equal the spike-panel taxon_name")

    def test_implausible_genome_size_fails(self):
        fixture = self.fixture()
        for row in fixture.size_rows:
            if row["target_label"] == "Hhat":
                row["genome_size_bp"] = "90000000"
        self.fail(fixture, "plausibility bounds")

    def test_nonnumeric_genome_size_fails(self):
        fixture = self.fixture()
        for row in fixture.size_rows:
            if row["target_label"] == "Hhat":
                row["genome_size_bp"] = "unknown"
        self.fail(fixture, "nonnumeric genome_size_bp")

    def test_unexpected_panel_label_fails(self):
        fixture = self.fixture()
        fixture.panel_rows.append({
            "label": "Ecoli", "taxon_name": "Escherichia coli",
            "assembly": "GCF_Y", "fasta": "references/genomes/Ecoli.fa",
            "weight": "1", "url": ""})
        self.fail(fixture, "spike panel labels do not match the frozen set")


# ------------------------------------------------------------------ case E
class InvalidRatioTest(AuditTestCase):

    def _corrupt(self, values):
        """Give Pana one bad observation per supplied value, plus good ones."""
        fixture = self.fixture(replication={"Pana": 1 + len(values)},
                               community_replication={"Pana": 1})
        replaced = 0
        for row in fixture.rows:
            if (row["target_feature"] == METAPHLAN_ALIAS["Pana"]
                    and row["analysis_population"] == "independent"
                    and replaced < len(values)):
                row["sensitivity_observed_over_expected"] = values[replaced]
                replaced += 1
        self.assertEqual(replaced, len(values))
        return fixture

    def test_invalid_ratios_are_excluded_without_pseudocounts(self):
        fixture = self._corrupt(["0", "-0.5", "nan", "inf"])
        out = self.succeed(fixture)
        rows = {(r["scope"], r["arm"], r["target_label"]): r
                for r in read_tsv(out / "metaphlan_genome_size_target_level.tsv")}
        row = rows[("independent", "read_proportional_sensitivity", "Pana")]
        self.assertEqual(int(row["eligible_observations"]), 1)
        self.assertEqual(int(row["excluded_observations"]), 4)
        reasons = dict(part.split("=") for part in row["exclusion_reasons"].split(";"))
        self.assertEqual(reasons["read_proportional_sensitivity_nonpositive"], "2")
        self.assertEqual(reasons["read_proportional_sensitivity_nonfinite"], "2")
        # The surviving observation carries the untouched closed-form value, so
        # no pseudocount was folded in anywhere.
        x = math.log2(GENOME_SIZES["Pana"])
        self.assertAlmostEqual(
            float(row["median_log2_observed_over_expected"]),
            SENSITIVITY_INTERCEPT_INDEPENDENT - x, places=12)
        self.assertAlmostEqual(float(row["median_observed_over_expected"]),
                               2.0 ** (SENSITIVITY_INTERCEPT_INDEPENDENT - x),
                               places=12)

    def test_exclusion_is_paired_across_arms(self):
        fixture = self._corrupt(["0"])
        out = self.succeed(fixture)
        rows = {(r["scope"], r["arm"], r["target_label"]): r
                for r in read_tsv(out / "metaphlan_genome_size_target_level.tsv")}
        sensitivity = rows[("independent", "read_proportional_sensitivity", "Pana")]
        primary = rows[("independent", "genome_equivalent_primary", "Pana")]
        self.assertEqual(sensitivity["eligible_observations"],
                         primary["eligible_observations"])
        self.assertEqual(sensitivity["excluded_observations"],
                         primary["excluded_observations"])

    def test_entirely_unusable_target_fails(self):
        fixture = self.fixture()
        for row in fixture.rows:
            if row["target_feature"] == METAPHLAN_ALIAS["Dpne"]:
                row["primary_observed_over_expected"] = "0"
        self.fail(fixture, "entirely unusable")

    def test_exclusion_totals_reach_the_validation_table(self):
        fixture = self._corrupt(["0", "nan"])
        out = self.succeed(fixture)
        validation = {r["metric"]: r["value"]
                      for r in read_tsv(out / "metaphlan_genome_size_audit_validation.tsv")}
        self.assertEqual(validation["paired_excluded_observations"], "2")
        self.assertEqual(validation["pseudocount_applied"], "NONE")
        self.assertIn("nonpositive", validation["exclusion_reasons"])


# ------------------------------------------------------------------ case F
class PairingAndProvenanceTest(AuditTestCase):

    def test_missing_primary_value_fails(self):
        fixture = self.fixture()
        fixture.rows[0]["primary_observed_over_expected"] = ""
        self.fail(fixture, "missing primary_observed_over_expected")

    def test_missing_sensitivity_value_fails(self):
        fixture = self.fixture()
        fixture.rows[0]["sensitivity_observed_over_expected"] = ""
        self.fail(fixture, "missing sensitivity_observed_over_expected")

    def test_missing_selected_provenance_fields_fail(self):
        for column in ("primary_row_reference_type",
                       "sensitivity_row_reference_type",
                       "primary_selected_reference_type",
                       "sensitivity_selected_reference_type",
                       "primary_reference_scale",
                       "sensitivity_reference_scale"):
            fixture = self.fixture()
            fixture.rows[0][column] = ""
            self.fail(fixture, f"missing {column}", name=f"blank_{column}")

    def test_missing_arm_column_fails(self):
        fixture = self.fixture()
        comparison, validation, sizes, panel, aliases = fixture.paths()
        rows = read_tsv(comparison)
        fields = [f for f in COMPARISON_FIELDS
                  if f != "primary_observed_over_expected"]
        for row in rows:
            row.pop("primary_observed_over_expected")
        write_tsv(comparison, fields, rows)
        out = self.root / "nocol"
        done = subprocess.run(
            [sys.executable, str(SCRIPT), "--comparison", str(comparison),
             "--comparison-validation", str(validation),
             "--target-genome-sizes", str(sizes), "--spike-panel", str(panel),
             "--feature-aliases", str(aliases), "--outdir", str(out)],
            text=True, capture_output=True)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("lacks column(s): primary_observed_over_expected",
                      done.stdout + done.stderr)

    def test_duplicate_observation_key_fails(self):
        fixture = self.fixture()
        fixture.rows.append(dict(fixture.rows[0]))
        self.fail(fixture, "duplicate physical observation key")

    def test_missing_bracken_gate_fails(self):
        fixture = self.fixture()
        fixture.validation_rows = [r for r in fixture.validation_rows
                                   if r["metric"] != "bracken_identical"]
        self.fail(fixture, "no bracken_identical metric")

    def test_failed_bracken_gate_fails(self):
        fixture = self.fixture()
        for row in fixture.validation_rows:
            if row["metric"] == "bracken_identical":
                row["value"] = "FAIL"
        self.fail(fixture, "Bracken identity gate is 'FAIL'")

    def test_single_profiler_fails(self):
        fixture = self.fixture()
        fixture.rows = [r for r in fixture.rows
                        if r["profiler"] != "kraken2_bracken"]
        self.fail(fixture, "must contain both profilers")

    def test_bracken_rows_never_enter_the_regression(self):
        plain = self.fixture()
        moved = self.fixture()
        for row in moved.rows:
            if row["profiler"] == "kraken2_bracken":
                row["primary_observed_over_expected"] = "1234.5"
                row["sensitivity_observed_over_expected"] = "0.0001"
        left = read_tsv(self.succeed(plain, name="plain")
                        / "metaphlan_genome_size_regression.tsv")
        right = read_tsv(self.succeed(moved, name="moved")
                         / "metaphlan_genome_size_regression.tsv")
        self.assertEqual(left, right,
                         "Bracken ratios changed a MetaPhlAn-only regression")
        validation = {r["metric"]: r["value"] for r in read_tsv(
            self.root / "plain/metaphlan_genome_size_audit_validation.tsv")}
        self.assertEqual(validation["bracken_rows_in_regression"], "0")
        self.assertEqual(validation["bracken_observations"], "2")
        self.assertEqual(validation["metaphlan_observations"], "20")

    def test_checksums_and_seed_are_recorded(self):
        import hashlib

        fixture = self.fixture()
        out = self.succeed(fixture)
        comparison, validation_path, sizes, panel, aliases = fixture.paths()
        validation = {r["metric"]: r["value"]
                      for r in read_tsv(out / "metaphlan_genome_size_audit_validation.tsv")}
        for metric, path in (("comparison_sha256", comparison),
                             ("comparison_validation_sha256", validation_path),
                             ("target_genome_sizes_sha256", sizes),
                             ("spike_panel_sha256", panel),
                             ("feature_aliases_sha256", aliases)):
            self.assertEqual(validation[metric],
                             hashlib.sha256(path.read_bytes()).hexdigest(), metric)
        self.assertEqual(validation["bootstrap_seed"], "20260920")
        self.assertEqual(validation["statistical_unit"], "implanted_target")
        self.assertEqual(validation["bootstrap_unit"], "implanted_target")
        self.assertEqual(validation["fitted_genome_size_constant_used"], "NONE")
        self.assertEqual(validation["development_status"], "DEVELOPMENT_ONLY")
        self.assertEqual(validation["target_count"], "10")
        self.assertEqual(validation["expected_target_count"], "10")
        checksums = (out / "metaphlan_genome_size_audit.sha256").read_text()
        self.assertIn(str(comparison.resolve()), checksums)
        self.assertIn("metaphlan_genome_size_regression.tsv", checksums)
        development = (out / "DEVELOPMENT_ONLY.txt").read_text()
        for needle in ("yachida", "use_for_manuscript=NO", "feng_replication",
                       "zeller_replication", "three_cohort_validation"):
            self.assertIn(needle, development)

    def test_nonempty_outdir_fails(self):
        fixture = self.fixture()
        out = self.root / "occupied"
        out.mkdir()
        (out / "stale.tsv").write_text("x\n", encoding="utf-8")
        done = fixture.run(out)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("must be new or empty", done.stdout + done.stderr)


# ------------------------------------------------------------------ case G
class BootstrapTest(AuditTestCase):

    def test_repeated_runs_are_byte_identical(self):
        fixture = self.fixture()
        first = self.succeed(fixture, name="first", replicates=500)
        second = self.succeed(fixture, name="second", replicates=500)
        for name in ("metaphlan_genome_size_target_level.tsv",
                     "metaphlan_genome_size_regression.tsv",
                     "metaphlan_genome_size_paired_change.tsv"):
            self.assertEqual((first / name).read_bytes(),
                             (second / name).read_bytes(), name)

    def test_replicate_counts_are_internally_consistent(self):
        out = self.succeed(self.fixture(), replicates=500)
        for row in read_tsv(out / "metaphlan_genome_size_regression.tsv"):
            requested = int(row["bootstrap_replicates_requested"])
            valid = int(row["bootstrap_replicates_valid"])
            discarded = int(row["bootstrap_replicates_discarded"])
            self.assertEqual(requested, 500)
            self.assertEqual(valid + discarded, requested)
            self.assertGreaterEqual(valid, math.ceil(0.95 * requested))
            self.assertEqual(row["bootstrap_seed"], "20260920")

    def test_interval_contains_the_fitted_slope(self):
        out = self.succeed(self.fixture(residual={"Hhat": 0.4, "Dpne": -0.3}),
                           replicates=1000)
        for row in read_tsv(out / "metaphlan_genome_size_regression.tsv"):
            slope = float(row["slope"])
            lower = float(row["bootstrap_ci_lower"])
            upper = float(row["bootstrap_ci_upper"])
            self.assertLessEqual(lower, upper)
            self.assertLessEqual(lower, slope + 1e-12)
            self.assertGreaterEqual(upper, slope - 1e-12)

    def test_degenerate_bootstrap_request_is_rejected(self):
        fixture = self.fixture()
        out = self.root / "one"
        done = fixture.run(out, replicates=1)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("at least 2", done.stdout + done.stderr)


# --------------------------------------------------------- mutation sensitivity
class MutationSensitivityTest(AuditTestCase):
    """Show concretely which assertions catch the two defects of concern."""

    def test_swapping_the_arms_is_detected(self):
        """Caught by KnownSlopeTest.test_slopes_match_the_closed_form."""
        fixture = self.fixture()
        for row in fixture.rows:
            for suffix in ("observed_over_expected", "absolute_relative_error"):
                left, right = f"primary_{suffix}", f"sensitivity_{suffix}"
                row[left], row[right] = row[right], row[left]
        rows = self.regression(self.succeed(fixture, name="swapped"))
        sensitivity = float(rows[("pooled", "read_proportional_sensitivity")]["slope"])
        primary = float(rows[("pooled", "genome_equivalent_primary")]["slope"])
        # The arms now carry each other's data, so the very assertions in
        # KnownSlopeTest (-1 for sensitivity, 0 for primary) would fail.
        self.assertAlmostEqual(sensitivity, 0.0, places=12)
        self.assertAlmostEqual(primary, -1.0, places=12)
        self.assertNotAlmostEqual(sensitivity, -1.0, places=6)
        self.assertNotAlmostEqual(primary, 0.0, places=6)

    def test_observation_level_regression_would_differ(self):
        """Caught by ReplicationTest.test_unbalanced_replication_cannot_reweight."""
        residual = {"Hhat": 0.6}
        fixture = self.fixture(residual=residual, replication={"Hhat": 40},
                               community_replication={"Hhat": 40})
        produced = float(self.regression(self.succeed(fixture))
                         [("independent", "read_proportional_sensitivity")]["slope"])
        target_points = [(math.log2(GENOME_SIZES[label]),
                          SENSITIVITY_INTERCEPT_INDEPENDENT
                          - math.log2(GENOME_SIZES[label]) + residual.get(label, 0.0))
                         for label in LABELS]
        observation_points = []
        for x, y in target_points:
            count = 40 if abs(y + x - SENSITIVITY_INTERCEPT_INDEPENDENT - 0.6) < 1e-12 else 1
            observation_points.extend([(x, y)] * count)
        self.assertAlmostEqual(produced, ols(target_points), places=12)
        self.assertNotAlmostEqual(produced, ols(observation_points), places=3)


# ---------------------------------------------- defect 1: selected reference
class SelectedReferenceTest(AuditTestCase):
    """Row provenance and the selected estimand are different things."""

    def test_genome_equivalent_row_in_the_sensitivity_arm_is_read_proportional(self):
        """The exact case the old fixture hid: the MetaPhlAn sensitivity arm
        keeps genome_equivalent ROW provenance while measuring read-proportional."""
        fixture = self.fixture()
        for row in fixture.rows:
            if row["profiler"] == "metaphlan4":
                self.assertEqual(row["sensitivity_row_reference_type"],
                                 "genome_equivalent")
        out = self.succeed(fixture)
        rows = {(r["scope"], r["arm"], r["target_label"]): r
                for r in read_tsv(out / "metaphlan_genome_size_target_level.tsv")}
        for label in LABELS:
            sensitivity = rows[("pooled", "read_proportional_sensitivity", label)]
            primary = rows[("pooled", "genome_equivalent_primary", label)]
            # The reported estimand, not the source row's provenance.
            self.assertEqual(sensitivity["reference_type"], "read_proportional")
            self.assertEqual(primary["reference_type"], "genome_equivalent")
            # Provenance is preserved beside it, not lost.
            self.assertEqual(sensitivity["row_reference_type"], "genome_equivalent")
            self.assertEqual(sensitivity["selected_reference_type"],
                             "read_proportional")
            self.assertEqual(sensitivity["reference_scale"], "read_proportional")
            self.assertEqual(primary["row_reference_type"], "genome_equivalent")
            self.assertEqual(primary["selected_reference_type"],
                             "profiler_scale_primary")
            self.assertEqual(primary["reference_scale"], "profiler_scale")

    def test_sensitivity_selected_type_claiming_genome_equivalent_fails(self):
        fixture = self.fixture()
        for row in fixture.rows:
            row["sensitivity_selected_reference_type"] = "profiler_scale_primary"
        self.fail(fixture, "read_proportional' and 'read_proportional' are required")

    def test_sensitivity_scale_claiming_profiler_scale_fails(self):
        fixture = self.fixture()
        for row in fixture.rows:
            row["sensitivity_reference_scale"] = "profiler_scale"
        self.fail(fixture, "the read_proportional_sensitivity arm carries")

    def test_primary_selected_type_claiming_read_proportional_fails(self):
        fixture = self.fixture()
        for row in fixture.rows:
            row["primary_selected_reference_type"] = "read_proportional"
        self.fail(fixture, "profiler_scale' and 'profiler_scale_primary' are required")

    def test_primary_scale_claiming_read_proportional_fails(self):
        fixture = self.fixture()
        for row in fixture.rows:
            row["primary_reference_scale"] = "read_proportional"
        self.fail(fixture, "the genome_equivalent_primary arm carries")

    def test_metaphlan_row_provenance_must_be_genome_equivalent(self):
        fixture = self.fixture()
        for row in fixture.rows:
            if row["profiler"] == "metaphlan4":
                row["primary_row_reference_type"] = "read_proportional"
        self.fail(fixture, "MetaPhlAn source rows must be 'genome_equivalent'")

    def test_mutating_only_selected_metadata_fails(self):
        """No quantitative cell changes; only the selection labels are swapped."""
        baseline = self.fixture()
        numbers = [{k: v for k, v in row.items() if "reference" not in k
                    and "scale" not in k} for row in baseline.rows]
        mutated = self.fixture()
        for row in mutated.rows:
            row["primary_selected_reference_type"], \
                row["sensitivity_selected_reference_type"] = (
                    row["sensitivity_selected_reference_type"],
                    row["primary_selected_reference_type"])
            row["primary_reference_scale"], row["sensitivity_reference_scale"] = (
                row["sensitivity_reference_scale"], row["primary_reference_scale"])
        self.assertEqual(
            numbers, [{k: v for k, v in row.items() if "reference" not in k
                       and "scale" not in k} for row in mutated.rows],
            "the mutation must touch metadata only")
        self.succeed(baseline, name="clean")
        self.fail(mutated, "arm carries", name="metadata_only")

    def test_ambiguous_legacy_alias_is_not_used_by_the_audit(self):
        """Corrupting only the legacy alias must not change any result."""
        clean = read_tsv(self.succeed(self.fixture(), name="clean")
                         / "metaphlan_genome_size_target_level.tsv")
        legacy = self.fixture()
        for row in legacy.rows:
            row["primary_reference_type"] = "NONSENSE"
            row["sensitivity_reference_type"] = "NONSENSE"
        after = read_tsv(self.succeed(legacy, name="legacy")
                         / "metaphlan_genome_size_target_level.tsv")
        self.assertEqual(clean, after)


# --------------------------------------------------------- defect 2: cohort
class CohortScopeTest(AuditTestCase):

    def _set_cohort(self, fixture, cohort, predicate=lambda row: True):
        for row in fixture.rows:
            if predicate(row):
                row["cohort"] = cohort
        return fixture

    def test_yachida_only_passes_and_is_recorded(self):
        out = self.succeed(self.fixture())
        validation = {r["metric"]: r["value"] for r in read_tsv(
            out / "metaphlan_genome_size_audit_validation.tsv")}
        self.assertEqual(validation["observed_cohorts"], "yachida")
        self.assertEqual(validation["validated_cohort_count"], "1")
        self.assertEqual(validation["required_cohort"], "yachida")
        self.assertEqual(validation["cohort_scope"], "yachida")

    def test_blank_cohort_fails(self):
        message = self.fail(self._set_cohort(self.fixture(), ""), "blank cohort")
        self.assertIn("cohort", message)

    def test_feng_only_fails(self):
        self.fail(self._set_cohort(self.fixture(), "feng"),
                  "restricted to the 'yachida' cohort")

    def test_zeller_only_fails(self):
        self.fail(self._set_cohort(self.fixture(), "zeller"), "'zeller'")

    def test_mixed_cohorts_fail(self):
        fixture = self._set_cohort(
            self.fixture(), "feng",
            lambda row: row["target_feature"] == METAPHLAN_ALIAS["Bfrag"])
        message = self.fail(fixture, "restricted to the 'yachida' cohort")
        self.assertIn("'feng'", message)
        self.assertIn("'yachida'", message)

    def test_one_mutated_row_fails(self):
        fixture = self.fixture()
        fixture.rows[3]["cohort"] = "feng"
        self.fail(fixture, "restricted to the 'yachida' cohort")

    def test_a_non_yachida_bracken_row_also_fails(self):
        """Every row is inspected, including rows never regressed."""
        fixture = self.fixture()
        for row in fixture.rows:
            if row["profiler"] == "kraken2_bracken":
                row["cohort"] = "feng"
        self.fail(fixture, "restricted to the 'yachida' cohort")

    def test_non_yachida_rows_are_not_silently_dropped(self):
        clean = self.succeed(self.fixture(), name="clean")
        clean_rows = read_tsv(clean / "metaphlan_genome_size_target_level.tsv")
        fixture = self.fixture()
        fixture.rows.append(comparison_row(
            "Bfrag", "community", "feng_extra", 1.0, PRIMARY_LOG2))
        fixture.rows[-1]["cohort"] = "feng"
        self.fail(fixture, "restricted to the 'yachida' cohort")
        self.assertTrue(clean_rows)


# ------------------------------------------------- defect 3: identity fields
class ObservationIdentityTest(AuditTestCase):

    KEY_FIELDS = ("cohort", "study", "sample_id", "condition",
                  "analysis_population", "profiler", "baseline_id",
                  "target_feature", "dose_index", "nominal_target_fraction")

    def test_every_identity_field_must_be_nonblank(self):
        for index, field in enumerate(self.KEY_FIELDS):
            fixture = self.fixture()
            fixture.rows[0][field] = ""
            message = self.fail(fixture, f"blank {field}", name=f"blank{index}")
            self.assertIn("line 2", message)
            self.assertIn(field, message)

    def test_whitespace_only_identity_field_is_blank(self):
        fixture = self.fixture()
        fixture.rows[0]["baseline_id"] = "   "
        self.fail(fixture, "blank baseline_id")

    def test_blank_gate_fires_before_the_duplicate_gate(self):
        """Two distinct rows collapse onto one key only because of blanks."""
        fixture = self.fixture()
        first = dict(fixture.rows[0])
        second = dict(fixture.rows[0])
        first["sample_id"] = ""
        second["sample_id"] = ""
        second["primary_observed_over_expected"] = first[
            "primary_observed_over_expected"]
        fixture.rows = [first, second] + fixture.rows[1:]
        message = self.fail(fixture, "blank sample_id")
        self.assertNotIn("duplicate physical observation key", message)

    def test_complete_duplicate_keys_still_fail(self):
        fixture = self.fixture()
        fixture.rows.append(dict(fixture.rows[0]))
        self.fail(fixture, "duplicate physical observation key")


# --------------------------------------------- defect 4: error consistency
class AbsoluteRelativeErrorTest(AuditTestCase):

    def _first_metaphlan(self, fixture):
        for row in fixture.rows:
            if row["profiler"] == "metaphlan4":
                return row
        raise AssertionError("no MetaPhlAn row")

    def test_consistent_errors_pass(self):
        out = self.succeed(self.fixture())
        validation = {r["metric"]: r["value"] for r in read_tsv(
            out / "metaphlan_genome_size_audit_validation.tsv")}
        self.assertIn("abs(observed_over_expected - 1)",
                      validation["absolute_relative_error_tolerance"])
        self.assertIn("never", validation["absolute_relative_error_policy"])

    def test_negative_error_fails(self):
        fixture = self.fixture()
        self._first_metaphlan(fixture)["primary_absolute_relative_error"] = "-0.1"
        self.fail(fixture, "negative primary_absolute_relative_error")

    def test_inconsistent_finite_error_fails(self):
        fixture = self.fixture()
        self._first_metaphlan(fixture)["sensitivity_absolute_relative_error"] = "0.42"
        self.fail(fixture, "is inconsistent with abs(observed_over_expected - 1)")

    def test_round_trip_difference_within_tolerance_passes(self):
        fixture = self.fixture()
        row = self._first_metaphlan(fixture)
        exact = abs(float(row["primary_observed_over_expected"]) - 1.0)
        row["primary_absolute_relative_error"] = g(exact + 1e-10)
        self.succeed(fixture)

    def test_difference_outside_tolerance_fails(self):
        fixture = self.fixture()
        row = self._first_metaphlan(fixture)
        exact = abs(float(row["primary_observed_over_expected"]) - 1.0)
        row["primary_absolute_relative_error"] = g(exact + 1e-6)
        self.fail(fixture, "round-trip tolerance")

    def test_mutating_the_primary_error_alone_fails(self):
        fixture = self.fixture()
        row = self._first_metaphlan(fixture)
        row["primary_absolute_relative_error"] = g(
            float(row["primary_absolute_relative_error"]) * 1.5 + 0.01)
        self.fail(fixture, "primary_absolute_relative_error")

    def test_mutating_the_sensitivity_error_alone_fails(self):
        fixture = self.fixture()
        row = self._first_metaphlan(fixture)
        row["sensitivity_absolute_relative_error"] = g(
            float(row["sensitivity_absolute_relative_error"]) * 1.5 + 0.01)
        self.fail(fixture, "sensitivity_absolute_relative_error")

    def test_nonnumeric_error_fails(self):
        fixture = self.fixture()
        self._first_metaphlan(fixture)["primary_absolute_relative_error"] = "oops"
        self.fail(fixture, "nonnumeric primary_absolute_relative_error")

    def test_excluded_ratio_is_not_rescued_by_its_error(self):
        """A paired-excluded observation must neither fail nor be rehabilitated."""
        fixture = self.fixture(replication={"Pana": 2})
        changed = 0
        for row in fixture.rows:
            if (row["target_feature"] == METAPHLAN_ALIAS["Pana"]
                    and row["analysis_population"] == "independent"
                    and changed == 0):
                # Unusable ratio AND a wildly inconsistent error on both arms.
                row["sensitivity_observed_over_expected"] = "0"
                row["sensitivity_absolute_relative_error"] = "99"
                row["primary_absolute_relative_error"] = "77"
                changed = 1
        self.assertEqual(changed, 1)
        out = self.succeed(fixture)
        rows = {(r["scope"], r["arm"], r["target_label"]): r
                for r in read_tsv(out / "metaphlan_genome_size_target_level.tsv")}
        row = rows[("independent", "read_proportional_sensitivity", "Pana")]
        self.assertEqual(int(row["excluded_observations"]), 1)
        self.assertEqual(int(row["eligible_observations"]), 1)
        self.assertIn("nonpositive", row["exclusion_reasons"])
        # The surviving observation keeps its untouched closed-form value.
        x = math.log2(GENOME_SIZES["Pana"])
        self.assertAlmostEqual(
            float(row["median_log2_observed_over_expected"]),
            SENSITIVITY_INTERCEPT_INDEPENDENT - x, places=12)
        self.assertNotIn("99", row["median_absolute_relative_error"])


# ------------------------------------------------------- runner integration
class RunnerIntegrationTest(unittest.TestCase):
    """The audit must be a resumable, hard-gated stage of the development run."""

    STAGE7 = "# ---------- stage 7: MetaPhlAn genome-size residual audit"
    REPORTING = "# ---------- reporting"

    @classmethod
    def setUpClass(cls):
        cls.runner = (ROOT / "analysis_v2/run_geff_propagation_development.sh")
        cls.text = cls.runner.read_text()
        assert cls.text.count(cls.STAGE7) == 1
        cls.stage = cls.text[cls.text.index(cls.STAGE7):
                             cls.text.index(cls.REPORTING)]

    def test_stage_is_wired_with_the_expected_inputs(self):
        self.assertIn("audit_metaphlan_genome_size_residual.py", self.text)
        for flag, value in (
            ("--comparison",
             "$COMPARISON_DIR/target_recovery_reference_comparison.tsv"),
            ("--comparison-validation",
             "$COMPARISON_DIR/reference_comparison_validation.tsv"),
            ("--target-genome-sizes", "$TARGET_GENOME_SIZES"),
            ("--spike-panel", "$SPIKE_PANEL"),
            ("--feature-aliases", "$ALIASES"),
            ("--outdir", "$AUDIT_DIR"),
        ):
            self.assertIn(f'{flag} "{value}"', self.text, flag)
        self.assertIn('SPIKE_PANEL="spikes/spike_panel.tsv"', self.text)
        self.assertIn(
            'AUDIT_DIR="$RUN_ROOT/metaphlan_genome_size_residual_audit"', self.text)
        self.assertIn('COMPARISON_DIR="$RUN_ROOT/reference_comparison"', self.text)

    def test_stage_runs_after_the_comparison_and_is_hard_gated(self):
        comparison_gate = self.text.index('test -s "$COMPARISON_DIR/SUCCESS"')
        audit_call = self.text.index("audit_metaphlan_genome_size_residual.py")
        self.assertLess(comparison_gate, audit_call,
                        "the audit must run after the reference comparison gate")
        self.assertIn('test -s "$AUDIT_DIR/SUCCESS" || { echo "FAIL genome-size '
                      'residual audit"; exit 1; }', self.text)
        # Both stages also re-assert schema compatibility after regenerating.
        self.assertIn('stage_schema_reason reference_comparison "$COMPARISON_DIR"',
                      self.text)
        self.assertIn('stage_schema_reason metaphlan_genome_size_residual_audit '
                      '"$AUDIT_DIR"', self.text)

    def test_stage_is_resumable_and_quarantines_partial_output(self):
        self.assertIn("source analysis_v2/lib/stage_compatibility.sh", self.text)
        for stage, directory, label in (
            ("reference_comparison", "$COMPARISON_DIR", "reference comparison"),
            ("metaphlan_genome_size_residual_audit", "$AUDIT_DIR",
             "genome-size residual audit"),
        ):
            self.assertIn(
                f'reuse_or_quarantine {stage} "{directory}" \\\n'
                f'       "$RUN_ROOT" "$STAMP" "{label}"', self.text, stage)
        # The reuse decision and the quarantine both live in the shared gate.
        library = (ROOT / "analysis_v2/lib/stage_compatibility.sh").read_text()
        self.assertIn("[REUSE]", library)
        self.assertIn('mv "$dir" "$quarantine"', library)
        self.assertIn("failed_attempts", library)

    def test_stage_uses_the_container_python(self):
        invocations = [line.strip() for line in self.text.splitlines()
                       if "audit_metaphlan_genome_size_residual.py" in line]
        self.assertEqual(len(invocations), 1, invocations)
        self.assertTrue(
            invocations[0].startswith("analysis_python "),
            f"the audit must run through analysis_python: {invocations[0]!r}")

    def test_stage_does_not_rebuild_upstream_work(self):
        stage = self.stage
        for upstream in ("build_perturbation_response_input.py",
                         "summarize_target_recovery.py",
                         "analyze_perturbation_response.py",
                         "compare_target_recovery_references.py",
                         "derive_endpoints_with_references"):
            self.assertNotIn(upstream, stage,
                             f"the audit stage must not rerun {upstream}")

    def test_comparison_migration_never_touches_the_recovery_arms(self):
        """Quarantining a stale comparison must not remove its input arms."""
        comparison = self.text[self.text.index("# ---------- stage 6:"):
                               self.text.index(self.STAGE7)]
        for arm in ("$RUN_ROOT/target_recovery/target_recovery_observations.tsv",
                    "$RUN_ROOT/target_recovery_read_sensitivity/"
                    "target_recovery_observations.tsv"):
            self.assertIn(arm, comparison)
        # Inspect commands, not prose: "per-arm reference" contains "rm ".
        for line in comparison.splitlines():
            code = line.split("#", 1)[0].strip()
            self.assertNotIn(code.split(" ")[0] if code else "",
                             {"rm", "rmdir", "shred", "truncate"}, line)
        # Only the comparison directory is ever moved.
        self.assertNotIn('mv "$RUN_ROOT/target_recovery', comparison)


if __name__ == "__main__":
    unittest.main()
