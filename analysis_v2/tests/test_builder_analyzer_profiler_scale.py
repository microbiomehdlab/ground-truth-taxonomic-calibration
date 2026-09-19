#!/usr/bin/env python3
"""End-to-end: real builder output -> real analyzer, on the profiler scale.

This deliberately does NOT fabricate already-consistent analyzer input. It runs
`build_perturbation_response_input.py` and feeds its genuine Parquet output to
`analyze_perturbation_response.py`.

Every MetaPhlAn quantity is checked against an independent, hand-derived
reference implementation (`Reference` below), never against the analyzer's own
helpers:

    q_it = f_it * G_eff,i / G_t ; Q_i = sum_k q_ik ; D_i = (1-F_i) + Q_i
    retained_ps = (1-F)o / D_i ; implanted_ps = q_it / D_i
    expected_ps = retained_ps + implanted_ps ; signal_ps = a - retained_ps

The fixture is deliberately discriminating: two cohorts, both profilers,
matched independent and community perturbations, unequal target genome sizes,
sample-specific G_eff, and a non-implanted feature that is nonetheless reported
inside every perturbation, so the operator's regressor cannot be confused with
the row's own `target_fraction_for_feature`.
"""

from __future__ import annotations

import csv
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import duckdb

ROOT = Path(__file__).resolve().parents[2]
BUILDER = ROOT / "analysis_v2/scripts/build_perturbation_response_input.py"
ANALYZER = ROOT / "analysis_v2/scripts/analyze_perturbation_response.py"
sys.path.insert(0, str(ROOT / "analysis_v2/scripts"))

M_FIELDS = ["cohort", "study", "analysis_population", "sample_id", "condition",
            "target_label", "assembly_arm", "profiler", "profile_id",
            "baseline_profile_id", "spike_fraction_target", "dose_level",
            "source_profile", "target_feature"]
E_FIELDS = ["cohort", "study", "analysis_population", "sample_id", "condition",
            "target_label", "assembly_arm", "profiler", "profile_id",
            "baseline_profile_id", "spike_fraction_total", "spike_fraction_target",
            "source_baseline_profile", "source_profile"]

# --- frozen fixture constants ---------------------------------------------
G = {"A": 2_000_000.0, "B": 6_000_000.0}          # unequal implanted genomes
FEATURE = {"A": "Species A", "B": "Species B"}
BYSTANDER = "Bystander"
FEATURES = (BYSTANDER, FEATURE["A"], FEATURE["B"])
COHORTS = ("yachida", "feng")
PROFILERS = {"metaphlan4": "M", "kraken2_bracken": "K"}

# Community dose_06 matches independent dose_04 under the analyzer's ordinal
# superposition join (community dose_index - 2).
COMMUNITY_DOSE = "dose_06"
INDEPENDENT_DOSE = "dose_04"
COMMUNITY_FRACTIONS = {"A": 0.03, "B": 0.02}
COMMUNITY_TOTAL = 0.05
INDEPENDENT_FRACTION = dict(COMMUNITY_FRACTIONS)   # exact superposition components

# Sample-specific G_eff: two clearly different values, neither equal to any G_t.
GEFF = {"M1": 3_477_000.0, "K1": 3_477_000.0,
        "M2": 4_284_000.0, "K2": 4_284_000.0}
SAMPLE = {(cohort, profiler): f"{prefix}{index}"
          for index, cohort in enumerate(COHORTS, start=1)
          for profiler, prefix in PROFILERS.items()}

BASELINE = {"yachida": {BYSTANDER: 0.10}, "feng": {BYSTANDER: 0.08}}
OBSERVED_COMMUNITY = {
    "yachida": {BYSTANDER: 0.09, "Species A": 0.020, "Species B": 0.010},
    "feng": {BYSTANDER: 0.072, "Species A": 0.026, "Species B": 0.013},
}
# Every independent perturbation also reports BOTH implanted species, so each
# observation contains a non-implanted feature whose driver must still be the
# perturbation's own implanted dose.
OBSERVED_INDEPENDENT = {
    "yachida": {"A": {BYSTANDER: 0.095, "Species A": 0.025, "Species B": 0.0010},
                "B": {BYSTANDER: 0.092, "Species A": 0.0008, "Species B": 0.011}},
    "feng": {"A": {BYSTANDER: 0.076, "Species A": 0.031, "Species B": 0.0012},
             "B": {BYSTANDER: 0.074, "Species A": 0.0009, "Species B": 0.014}},
}


def write(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_fixture(root):
    m_rows, e_rows, abundance = [], [], []
    for cohort in COHORTS:
        for profiler in PROFILERS:
            sample = SAMPLE[(cohort, profiler)]
            base = f"/p/{cohort}/{profiler}/base"
            community = f"/p/{cohort}/{profiler}/comm"
            common = dict(cohort=cohort, study="S", condition="Control",
                          assembly_arm="original", profiler=profiler,
                          sample_id=sample, baseline_profile_id=sample)

            for label, fraction in COMMUNITY_FRACTIONS.items():
                row = {**common, "analysis_population": "community",
                       "target_label": label, "profile_id": f"{sample}_comm",
                       "spike_fraction_target": fraction,
                       "source_profile": community}
                m_rows.append({**row, "dose_level": COMMUNITY_DOSE,
                               "target_feature": FEATURE[label]})
                e_rows.append({**row, "spike_fraction_total": COMMUNITY_TOTAL,
                               "source_baseline_profile": base})

            for label, fraction in INDEPENDENT_FRACTION.items():
                profile = f"/p/{cohort}/{profiler}/ind_{label}"
                row = {**common, "analysis_population": "independent",
                       "target_label": label, "profile_id": f"{sample}_ind{label}",
                       "spike_fraction_target": fraction, "source_profile": profile}
                m_rows.append({**row, "dose_level": INDEPENDENT_DOSE,
                               "target_feature": FEATURE[label]})
                e_rows.append({**row, "spike_fraction_total": fraction,
                               "source_baseline_profile": base})

            for feature, value in BASELINE[cohort].items():
                abundance.append({"profiler": profiler, "source_profile": base,
                                  "feature": feature, "abundance_fraction": value})
            for feature, value in OBSERVED_COMMUNITY[cohort].items():
                abundance.append({"profiler": profiler, "source_profile": community,
                                  "feature": feature, "abundance_fraction": value})
            for label, observed in OBSERVED_INDEPENDENT[cohort].items():
                profile = f"/p/{cohort}/{profiler}/ind_{label}"
                for feature, value in observed.items():
                    abundance.append({"profiler": profiler, "source_profile": profile,
                                      "feature": feature, "abundance_fraction": value})

    manifest, endpoints, abund = (root / "manifest.tsv", root / "endpoints.tsv",
                                  root / "abundance.tsv")
    write(manifest, M_FIELDS, m_rows)
    write(endpoints, E_FIELDS, e_rows)
    write(abund, ["profiler", "source_profile", "feature", "abundance_fraction"],
          abundance)
    targets = root / "targets.tsv"
    write(targets, ["target_label", "genome_size_bp"],
          [{"target_label": k, "genome_size_bp": v} for k, v in G.items()])
    geff = root / "geff.tsv"
    write(geff, ["cohort", "sample_id", "effective_genome_size_bp"],
          [{"cohort": cohort, "sample_id": SAMPLE[(cohort, profiler)],
            "effective_genome_size_bp": GEFF[SAMPLE[(cohort, profiler)]]}
           for cohort in COHORTS for profiler in PROFILERS])
    return manifest, endpoints, abund, targets, geff


# --- independent hand-derived reference ------------------------------------
class Reference:
    """Closed-form expectations, written from the frozen formulas alone."""

    def __init__(self, scale):
        self.scale = scale

    def genome_equivalent(self, profiler):
        return self.scale == "profiler_scale" and profiler == "metaphlan4"

    def denominator(self, cohort, profiler, labels, total):
        """D_i = (1-F) + Q_i, or 1 on the read-proportional scale."""
        if not self.genome_equivalent(profiler):
            return 1.0
        sample = SAMPLE[(cohort, profiler)]
        fractions = (COMMUNITY_FRACTIONS if len(labels) > 1
                     else INDEPENDENT_FRACTION)
        q = sum(fractions[label] * GEFF[sample] / G[label] for label in labels)
        return (1.0 - total) + q

    def terms(self, cohort, profiler, labels, total, feature, baseline, observed):
        """(retained, implanted, expected, signal) on the selected scale."""
        denominator = self.denominator(cohort, profiler, labels, total)
        retained = (1.0 - total) * baseline / denominator
        implanted = 0.0
        for label in labels:
            if FEATURE[label] != feature:
                continue
            fractions = (COMMUNITY_FRACTIONS if len(labels) > 1
                         else INDEPENDENT_FRACTION)
            if self.genome_equivalent(profiler):
                sample = SAMPLE[(cohort, profiler)]
                implanted += fractions[label] * GEFF[sample] / G[label] / denominator
            else:
                implanted += fractions[label]
        return retained, implanted, retained + implanted, observed - retained

    def independent(self, cohort, profiler, label, feature):
        baseline = BASELINE[cohort].get(feature, 0.0)
        observed = OBSERVED_INDEPENDENT[cohort][label].get(feature, 0.0)
        return self.terms(cohort, profiler, (label,), INDEPENDENT_FRACTION[label],
                          feature, baseline, observed)

    def community(self, cohort, profiler, feature):
        baseline = BASELINE[cohort].get(feature, 0.0)
        observed = OBSERVED_COMMUNITY[cohort].get(feature, 0.0)
        return self.terms(cohort, profiler, tuple(COMMUNITY_FRACTIONS),
                          COMMUNITY_TOTAL, feature, baseline, observed)

    def driver(self, cohort, profiler, label):
        """Selected-scale implanted driver: q_it/D_i for MetaPhlAn, else f_it."""
        return self.independent(cohort, profiler, label, FEATURE[label])[1]

    def slope(self, cohort, profiler, label, feature):
        """Through-origin slope of one cohort's single independent observation."""
        _, _, _, signal = self.independent(cohort, profiler, label, feature)
        return signal / self.driver(cohort, profiler, label)


def run_builder(manifest, endpoints, abundance, outdir, extra):
    return subprocess.run(
        [sys.executable, str(BUILDER), "--profile-manifest", str(manifest),
         "--endpoints", str(endpoints), "--abundance", str(abundance),
         "--outdir", str(outdir), "--threads", "1", "--memory-limit", "1GB",
         "--expected-community-targets", "2", *map(str, extra)],
        text=True, capture_output=True)


def run_analyzer(responses, outdir, scale):
    return subprocess.run(
        [sys.executable, str(ANALYZER), "--responses", str(responses),
         "--outdir", str(outdir), "--reference-scale", scale],
        text=True, capture_output=True)


def load(parquet):
    con = duckdb.connect()
    names = [c[0] for c in con.execute(
        "DESCRIBE SELECT * FROM read_parquet(?)", [str(parquet)]).fetchall()]
    out = {}
    for row in con.execute("SELECT * FROM read_parquet(?)", [str(parquet)]).fetchall():
        record = dict(zip(names, row))
        out[(record["cohort"], record["profiler"], record["analysis_population"],
             record["implanted_targets"], record["feature"])] = record
    return out


class BuilderAnalyzerTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        cls.root = Path(cls._tmp.name)
        manifest, endpoints, abundance, targets, geff = build_fixture(cls.root)
        cls.built = cls.root / "response_input"
        done = run_builder(manifest, endpoints, abundance, cls.built,
                           ["--metaphlan-reference", "genome_equivalent",
                            "--target-genome-sizes", targets,
                            "--effective-genome-sizes", geff])
        assert done.returncode == 0, done.stdout + done.stderr
        cls.responses = cls.built / "paired_feature_responses.parquet"
        cls.rows = load(cls.responses)
        cls.out = {}
        for scale in ("profiler_scale", "read_proportional"):
            directory = cls.root / f"analysis_{scale}"
            result = run_analyzer(cls.responses, directory, scale)
            assert result.returncode == 0, result.stdout + result.stderr
            cls.out[scale] = directory

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def key(self, cohort, profiler, population, targets, feature):
        return self.rows[(cohort, profiler, population, targets, feature)]

    def independent_row(self, cohort, profiler, label, feature):
        return self.key(cohort, profiler, "independent", label, feature)

    def community_row(self, cohort, profiler, feature):
        return self.key(cohort, profiler, "community", "A;B", feature)

    def table(self, scale, name):
        return read_tsv(self.out[scale] / name)

    # ---------------------------------------------------------------- builder
    def test_builder_profiler_scale_matches_hand_calculation(self):
        reference = Reference("profiler_scale")
        for cohort in COHORTS:
            for feature in FEATURES:
                retained, implanted, expected, signal = reference.community(
                    cohort, "metaphlan4", feature)
                row = self.community_row(cohort, "metaphlan4", feature)
                self.assertAlmostEqual(row["retained_baseline_profiler_scale"],
                                       retained, places=14)
                self.assertAlmostEqual(row["implanted_signal_profiler_scale"],
                                       implanted, places=14)
                self.assertAlmostEqual(row["expected_abundance_profiler_scale"],
                                       expected, places=14)
                self.assertAlmostEqual(row["response_signal_profiler_scale"],
                                       signal, places=14)
                for label in INDEPENDENT_FRACTION:
                    retained, implanted, expected, signal = reference.independent(
                        cohort, "metaphlan4", label, feature)
                    row = self.independent_row(cohort, "metaphlan4", label, feature)
                    self.assertAlmostEqual(
                        row["retained_baseline_profiler_scale"], retained, places=14)
                    self.assertAlmostEqual(
                        row["implanted_signal_profiler_scale"], implanted, places=14)
                    self.assertAlmostEqual(
                        row["expected_abundance_profiler_scale"], expected, places=14)
                    self.assertAlmostEqual(
                        row["response_signal_profiler_scale"], signal, places=14)

    def test_fixture_is_discriminating(self):
        """(F) The two references must differ clearly for MetaPhlAn."""
        selected = Reference("profiler_scale")
        read = Reference("read_proportional")
        for label in INDEPENDENT_FRACTION:
            q = selected.driver("yachida", "metaphlan4", label)
            f = read.driver("yachida", "metaphlan4", label)
            self.assertGreater(abs(q - f) / f, 0.3,
                               f"target {label} does not discriminate the scales")
        # Unequal genome sizes must give unequal per-read scaling.
        self.assertGreater(
            selected.driver("yachida", "metaphlan4", "A") / COMMUNITY_FRACTIONS["A"],
            selected.driver("yachida", "metaphlan4", "B") / COMMUNITY_FRACTIONS["B"])
        # Sample-specific G_eff must actually move the reference.
        self.assertNotAlmostEqual(selected.driver("yachida", "metaphlan4", "A"),
                                  selected.driver("feng", "metaphlan4", "A"))

    # ------------------------------------------------------- (A) retained row
    def test_response_row_retained_baseline_is_the_supplied_value(self):
        """(A) ResponseRow.retained_baseline is the supplied selected column."""
        from analyze_perturbation_response import iter_rows

        reference = Reference("profiler_scale")
        seen = 0
        for row in iter_rows(self.responses, "profiler_scale"):
            if row.profiler != "metaphlan4":
                continue
            supplied = self.rows[(row.cohort, row.profiler, row.population,
                                  ";".join(row.targets), row.feature)]
            self.assertAlmostEqual(row.retained_baseline,
                                   supplied["retained_baseline_profiler_scale"],
                                   places=15)
            # It is emphatically NOT the read-proportional reconstruction.
            reconstructed = (1.0 - row.actual_total) * row.baseline
            if row.baseline > 0:
                self.assertNotAlmostEqual(row.retained_baseline, reconstructed,
                                          places=6)
            seen += 1
        self.assertEqual(seen, 18)
        # Cross-check one cell against the closed form.
        retained, _, _, _ = reference.community("yachida", "metaphlan4", BYSTANDER)
        self.assertAlmostEqual(
            (1 - COMMUNITY_TOTAL) * BASELINE["yachida"][BYSTANDER]
            / reference.denominator("yachida", "metaphlan4",
                                    tuple(COMMUNITY_FRACTIONS), COMMUNITY_TOTAL),
            retained, places=15)

    # ------------------------------------------------- (B) operator regressor
    def test_operator_x_uses_genome_equivalent_driver(self):
        """(B) The operator's x is q_it/D_i, not the raw read fraction f."""
        reference = Reference("profiler_scale")
        read = Reference("read_proportional")
        operators = {(r["cohort"], r["profiler"], r["target_label"], r["feature"]): r
                     for r in self.table("profiler_scale", "response_operator.tsv")}
        checked = 0
        relative_gaps = []
        for cohort in COHORTS:
            for label in INDEPENDENT_FRACTION:
                for feature in FEATURES:
                    row = operators[(cohort, "metaphlan4", label, feature)]
                    expected = reference.slope(cohort, "metaphlan4", label, feature)
                    wrong = read.slope(cohort, "metaphlan4", label, feature)
                    self.assertAlmostEqual(float(row["operator_slope"]), expected,
                                           delta=abs(expected) * 1e-9 + 1e-15,
                                           msg=f"{cohort}/{label}/{feature}")
                    if abs(expected) > 1e-12:
                        gap = abs(expected - wrong) / abs(expected)
                        self.assertGreater(
                            gap, 0.05,
                            f"{cohort}/{label}/{feature} barely discriminates")
                        relative_gaps.append(gap)
                    checked += 1
        self.assertEqual(checked, 12)
        # The two references must be grossly, not marginally, different.
        self.assertGreater(max(relative_gaps), 0.3)
        # The driver itself is what moves; check it directly on both targets.
        for label in INDEPENDENT_FRACTION:
            self.assertGreater(
                abs(reference.driver("yachida", "metaphlan4", label)
                    - read.driver("yachida", "metaphlan4", label))
                / read.driver("yachida", "metaphlan4", label), 0.3)
        # The non-implanted feature is driven by the perturbation's dose, not by
        # its own (structurally zero) target_fraction_for_feature.
        row = self.independent_row("yachida", "metaphlan4", "A", BYSTANDER)
        self.assertAlmostEqual(row["implanted_signal_profiler_scale"], 0.0)
        self.assertNotAlmostEqual(
            float(operators[("yachida", "metaphlan4", "A", BYSTANDER)]
                  ["operator_slope"]), 0.0)

    # ------------------------------------------------- (C) holdout prediction
    def test_heldout_operator_validation_matches_hand_calculation(self):
        """(C, E) Prediction = selected retained baseline + driver * slope."""
        for scale in ("profiler_scale", "read_proportional"):
            reference = Reference(scale)
            expected = {}
            for cohort in COHORTS:
                training = [c for c in COHORTS if c != cohort]
                for profiler in PROFILERS:
                    for label in INDEPENDENT_FRACTION:
                        driver = reference.driver(cohort, profiler, label)
                        for feature in FEATURES:
                            retained, _, expectation, _ = reference.independent(
                                cohort, profiler, label, feature)
                            # Cross-cohort through-origin slope pooled over the
                            # training cohorts (one observation each here).
                            numerator = sum(
                                reference.driver(t, profiler, label)
                                * reference.independent(t, profiler, label,
                                                        feature)[3]
                                for t in training)
                            denominator = sum(
                                reference.driver(t, profiler, label) ** 2
                                for t in training)
                            slope = numerator / denominator
                            prediction = retained + driver * slope
                            observed = OBSERVED_INDEPENDENT[cohort][label].get(
                                feature, 0.0)
                            role = ("implanted_target"
                                    if feature == FEATURE[label] else "bystander")
                            bucket = expected.setdefault(
                                (cohort, profiler, role),
                                {"n": 0, "composition": 0.0, "operator": 0.0,
                                 "square": 0.0})
                            bucket["n"] += 1
                            bucket["composition"] += abs(observed - expectation)
                            bucket["operator"] += abs(observed - prediction)
                            bucket["square"] += (observed - prediction) ** 2
            observed_table = {
                (r["holdout_cohort"], r["profiler"], r["feature_role"]): r
                for r in self.table(scale, "heldout_operator_validation.tsv")}
            self.assertEqual(set(observed_table), set(expected))
            for key, bucket in expected.items():
                row = observed_table[key]
                self.assertEqual(int(row["eligible_contexts"]), bucket["n"], key)
                self.assertAlmostEqual(
                    float(row["composition_only_mae"]),
                    bucket["composition"] / bucket["n"], places=12, msg=key)
                self.assertAlmostEqual(
                    float(row["operator_prediction_mae"]),
                    bucket["operator"] / bucket["n"], places=12, msg=key)
                self.assertAlmostEqual(
                    float(row["operator_prediction_rmse"]),
                    math.sqrt(bucket["square"] / bucket["n"]), places=12, msg=key)

    # --------------------------------------------------- (D) superposition
    def test_superposition_summary_matches_hand_calculation(self):
        """(D, E) predicted = community retained_ps + sum(independent signal_ps)."""
        for scale in ("profiler_scale", "read_proportional"):
            reference = Reference(scale)
            expected = {}
            for cohort in COHORTS:
                for profiler in PROFILERS:
                    for feature in FEATURES:
                        retained, _, expectation, _ = reference.community(
                            cohort, profiler, feature)
                        prediction = retained + sum(
                            reference.independent(cohort, profiler, label, feature)[3]
                            for label in INDEPENDENT_FRACTION)
                        observed = OBSERVED_COMMUNITY[cohort].get(feature, 0.0)
                        role = ("implanted_target"
                                if feature in FEATURE.values() else "bystander")
                        bucket = expected.setdefault(
                            (cohort, profiler, role),
                            {"n": 0, "composition": 0.0, "operator": 0.0,
                             "square": 0.0})
                        bucket["n"] += 1
                        bucket["composition"] += abs(observed - expectation)
                        bucket["operator"] += abs(observed - prediction)
                        bucket["square"] += (observed - prediction) ** 2
            observed_table = {(r["cohort"], r["profiler"], r["feature_role"]): r
                              for r in self.table(scale, "superposition_summary.tsv")}
            self.assertEqual(set(observed_table), set(expected))
            for key, bucket in expected.items():
                row = observed_table[key]
                self.assertEqual(row["dose_level"], COMMUNITY_DOSE)
                self.assertEqual(int(row["eligible_contexts"]), bucket["n"], key)
                self.assertAlmostEqual(
                    float(row["composition_only_mae"]),
                    bucket["composition"] / bucket["n"], places=12, msg=key)
                self.assertAlmostEqual(
                    float(row["operator_prediction_mae"]),
                    bucket["operator"] / bucket["n"], places=12, msg=key)
                self.assertAlmostEqual(
                    float(row["operator_prediction_rmse"]),
                    math.sqrt(bucket["square"] / bucket["n"]), places=12, msg=key)

    def test_superposition_never_reads_read_proportional_columns(self):
        """A profiler-scale superposition must fail if only the old columns exist."""
        from analyze_perturbation_response import exact_superposition_rows

        stripped = self.root / "stripped.parquet"
        con = duckdb.connect()
        con.execute(
            "CREATE TABLE t AS SELECT * EXCLUDE (response_signal_profiler_scale, "
            "expected_abundance_profiler_scale, retained_baseline_profiler_scale) "
            "FROM read_parquet(?)", [str(self.responses)])
        con.execute(f"COPY t TO '{stripped}' (FORMAT PARQUET)")
        with self.assertRaises(ValueError) as caught:
            exact_superposition_rows(stripped, "profiler_scale")
        message = str(caught.exception)
        for column in ("response_signal_profiler_scale",
                       "expected_abundance_profiler_scale",
                       "retained_baseline_profiler_scale"):
            self.assertIn(column, message)
        # The read-proportional selection still works on the same stripped table.
        self.assertTrue(exact_superposition_rows(stripped, "read_proportional"))

    # ------------------------------------------- (F) scales must differ
    def test_metaphlan_results_differ_between_scales(self):
        for name, key_fields, value_field in (
            ("superposition_summary.tsv", ("cohort", "profiler", "feature_role"),
             "operator_prediction_mae"),
            ("heldout_operator_validation.tsv",
             ("holdout_cohort", "profiler", "feature_role"),
             "operator_prediction_mae"),
            ("response_operator.tsv",
             ("cohort", "profiler", "target_label", "feature"), "operator_slope"),
        ):
            def index(scale):
                return {tuple(r[f] for f in key_fields): r
                        for r in self.table(scale, name)}
            selected, read = index("profiler_scale"), index("read_proportional")
            differing = [key for key in selected
                         if "metaphlan4" in key
                         and selected[key][value_field] != read[key][value_field]]
            self.assertTrue(differing,
                            f"{name} is identical under both scales: not discriminating")

    # ------------------------------------------- (G) Bracken is invariant
    def test_bracken_identical_between_scale_selections(self):
        for name, profiler_field in (
            ("superposition_summary.tsv", "profiler"),
            ("heldout_operator_validation.tsv", "profiler"),
            ("response_operator.tsv", "profiler"),
            ("quantitative_response_summary.tsv", "profiler"),
            ("reliability_certificates.tsv", "profiler"),
        ):
            selected = [r for r in self.table("profiler_scale", name)
                        if r[profiler_field] == "kraken2_bracken"]
            read = [r for r in self.table("read_proportional", name)
                    if r[profiler_field] == "kraken2_bracken"]
            self.assertTrue(selected)
            self.assertEqual(selected, read, f"Bracken changed in {name}")

    # ------------------------------------------- (H) detection is invariant
    def test_detection_identical_between_scale_selections(self):
        self.assertEqual(self.table("profiler_scale", "detection_response_summary.tsv"),
                         self.table("read_proportional", "detection_response_summary.tsv"))
        for key, row in self.rows.items():
            self.assertEqual(row["baseline_detected"],
                             int(row["baseline_abundance_fraction"] > 0), key)
            self.assertEqual(row["observed_detected"],
                             int(row["observed_abundance_fraction"] > 0), key)

    def test_identity_holds_for_every_row(self):
        for key, row in self.rows.items():
            self.assertAlmostEqual(
                row["expected_abundance_profiler_scale"],
                row["retained_baseline_profiler_scale"]
                + row["implanted_signal_profiler_scale"],
                places=12, msg=f"identity failed for {key}")

    def test_read_proportional_columns_unchanged_for_every_row(self):
        for key, row in self.rows.items():
            total = row["effective_total_fraction"]
            self.assertAlmostEqual(
                row["dilution_retained_baseline"],
                (1 - total) * row["baseline_abundance_fraction"], places=15, msg=key)


class FailClosedTest(unittest.TestCase):
    """(I) Corrupting the selected retained baseline or driver must fail."""

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        cls.root = Path(cls._tmp.name)
        manifest, endpoints, abundance, targets, geff = build_fixture(cls.root)
        done = run_builder(manifest, endpoints, abundance, cls.root / "built",
                           ["--metaphlan-reference", "genome_equivalent",
                            "--target-genome-sizes", targets,
                            "--effective-genome-sizes", geff])
        assert done.returncode == 0, done.stdout + done.stderr
        cls.responses = cls.root / "built/paired_feature_responses.parquet"

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def corrupt(self, name, expression, predicate="profiler='metaphlan4'"):
        path = self.root / f"corrupt_{name}.parquet"
        con = duckdb.connect()
        con.execute(
            f"CREATE TABLE t AS SELECT * EXCLUDE ({name}), "
            f"CASE WHEN {predicate} THEN {expression} ELSE {name} END AS {name} "
            "FROM read_parquet(?)", [str(self.responses)])
        con.execute(f"COPY t TO '{path}' (FORMAT PARQUET)")
        return run_analyzer(path, self.root / f"out_{name}", "profiler_scale")

    def test_corrupted_retained_baseline_fails_closed(self):
        done = self.corrupt("retained_baseline_profiler_scale",
                            f"{'retained_baseline_profiler_scale'} * 1.5")
        self.assertNotEqual(done.returncode, 0,
                            "analyzer accepted a corrupted retained baseline")
        self.assertIn("expected abundance", done.stdout + done.stderr)

    def test_corrupted_driver_fails_closed(self):
        done = self.corrupt("implanted_signal_profiler_scale",
                            "implanted_signal_profiler_scale * 1.5",
                            "profiler='metaphlan4' AND is_direct_target")
        self.assertNotEqual(done.returncode, 0,
                            "analyzer accepted a corrupted implanted driver")
        self.assertIn("expected abundance", done.stdout + done.stderr)

    def test_corrupted_bounded_error_fails_closed(self):
        done = self.corrupt("signed_bounded_error_profiler_scale", "0.42",
                            f"profiler='metaphlan4' AND feature='{BYSTANDER}'")
        self.assertNotEqual(done.returncode, 0,
                            "analyzer accepted a corrupted bounded error")
        self.assertIn("bounded error", done.stdout + done.stderr)

    def test_missing_direct_target_row_fails_closed(self):
        """Without its direct-target row an independent perturbation has no driver."""
        path = self.root / "no_direct.parquet"
        con = duckdb.connect()
        con.execute(
            "CREATE TABLE t AS SELECT * FROM read_parquet(?) "
            "WHERE NOT (profiler='metaphlan4' AND analysis_population='independent' "
            "           AND implanted_targets='A' AND is_direct_target)",
            [str(self.responses)])
        con.execute(f"COPY t TO '{path}' (FORMAT PARQUET)")
        done = run_analyzer(path, self.root / "out_no_direct", "profiler_scale")
        self.assertNotEqual(done.returncode, 0,
                            "analyzer accepted an observation with no driver")
        self.assertIn("direct-target row", done.stdout + done.stderr)

    def test_non_positive_driver_fails_closed(self):
        done = self.corrupt("implanted_signal_profiler_scale", "0.0",
                            "profiler='metaphlan4' AND is_direct_target "
                            "AND analysis_population='independent'")
        self.assertNotEqual(done.returncode, 0,
                            "analyzer accepted a non-positive driver")

    def test_duplicate_identical_direct_target_row_fails_closed(self):
        """Exactly one driver row is required, even when duplicates agree."""
        path = self.root / "duplicate_direct.parquet"
        con = duckdb.connect()
        con.execute(
            "CREATE TABLE source AS SELECT * FROM read_parquet(?)",
            [str(self.responses)])
        con.execute("CREATE TABLE t AS SELECT * FROM source")
        con.execute(
            "INSERT INTO t SELECT * FROM source "
            "WHERE profiler='metaphlan4' AND analysis_population='independent' "
            "AND implanted_targets='A' AND is_direct_target LIMIT 1")
        con.execute(f"COPY t TO '{path}' (FORMAT PARQUET)")
        done = run_analyzer(path, self.root / "out_duplicate_direct",
                            "profiler_scale")
        self.assertNotEqual(done.returncode, 0,
                            "analyzer accepted duplicate identical driver rows")
        self.assertIn("more than one direct-target driver row",
                      done.stdout + done.stderr)


if __name__ == "__main__":
    unittest.main()
