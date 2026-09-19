#!/usr/bin/env python3
"""Genome-equivalent MetaPhlAn reference in the perturbation-response input.

Bracken keeps e_read = (1-F)o + f. MetaPhlAn uses
  q_it = f_it * G_eff/G_t ,  Q_i = sum_k q_ik
  e_MP,it = ((1-F)o + q_it) / ((1-F) + Q_i)      implanted target
  e_MP,ij = ((1-F)o)        / ((1-F) + Q_i)      non-implanted feature

The non-implanted case is the point: under Bracken a non-implanted feature is
only diluted, but under MetaPhlAn the whole composition is renormalised.
"""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import duckdb

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis_v2/scripts/build_perturbation_response_input.py"

M_FIELDS = ["cohort", "study", "analysis_population", "sample_id", "condition",
            "target_label", "assembly_arm", "profiler", "profile_id",
            "baseline_profile_id", "spike_fraction_target", "dose_level",
            "source_profile", "target_feature"]
E_FIELDS = ["cohort", "study", "analysis_population", "sample_id", "condition",
            "target_label", "assembly_arm", "profiler", "profile_id",
            "baseline_profile_id", "spike_fraction_total", "spike_fraction_target",
            "source_baseline_profile", "source_profile"]
GEFF = 3_477_000.0


def write(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build(root, targets, profiler="metaphlan4", population="community",
          total=None, baseline_bystander=0.10, observed_bystander=0.095,
          sample="X", cohort="yachida"):
    total = sum(f for _, _, f in targets) if total is None else total
    m_rows, e_rows = [], []
    for label, feature, fraction in targets:
        common = dict(cohort=cohort, study="S", analysis_population=population,
                      sample_id=sample, condition="Control", target_label=label,
                      assembly_arm="original", profiler=profiler,
                      profile_id=f"{sample}_mix", baseline_profile_id=sample,
                      spike_fraction_target=fraction, source_profile="/profiles/mix")
        m_rows.append({**common, "dose_level": "dose_06", "target_feature": feature})
        e_rows.append({**common, "spike_fraction_total": total,
                       "source_baseline_profile": "/profiles/base"})
    manifest, endpoints = root / "manifest.tsv", root / "endpoints.tsv"
    abundance = root / "abundance.tsv"
    write(manifest, M_FIELDS, m_rows)
    write(endpoints, E_FIELDS, e_rows)
    rows = [{"profiler": profiler, "source_profile": "/profiles/base",
             "feature": "Bystander", "abundance_fraction": baseline_bystander},
            {"profiler": profiler, "source_profile": "/profiles/mix",
             "feature": "Bystander", "abundance_fraction": observed_bystander}]
    for _, feature, _ in targets:
        rows.append({"profiler": profiler, "source_profile": "/profiles/mix",
                     "feature": feature, "abundance_fraction": 0.02})
    write(abundance, ["profiler", "source_profile", "feature", "abundance_fraction"], rows)
    return manifest, endpoints, abundance


def sizes(root, mapping):
    path = root / "target_sizes.tsv"
    write(path, ["target_label", "genome_size_bp"],
          [{"target_label": k, "genome_size_bp": v} for k, v in mapping.items()])
    return path


def geff(root, value=GEFF, sample="X", cohort="yachida"):
    path = root / "geff.tsv"
    write(path, ["cohort", "sample_id", "effective_genome_size_bp"],
          [{"cohort": cohort, "sample_id": sample, "effective_genome_size_bp": value}])
    return path


def run(manifest, endpoints, abundance, out, extra=(), n_targets=2):
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--profile-manifest", str(manifest),
         "--endpoints", str(endpoints), "--abundance", str(abundance),
         "--outdir", str(out), "--threads", "1", "--memory-limit", "1GB",
         "--expected-community-targets", str(n_targets), *map(str, extra)],
        text=True, capture_output=True)


def responses(out):
    con = duckdb.connect()
    cols = ("feature, reference_type, expected_abundance_fraction, "
            "read_proportional_reference, expected_abundance_profiler_scale, "
            "retained_baseline_profiler_scale, effective_community_genome_size_bp, "
            "total_implanted_genome_equivalent_fraction, "
            "implanted_genome_equivalent_fraction, target_genome_size_bp, "
            "is_direct_target, observed_detected, baseline_detected")
    rows = con.execute(f"SELECT {cols} FROM read_parquet(?) ORDER BY feature",
                       [str(out / "paired_feature_responses.parquet")]).fetchall()
    keys = [c.strip() for c in cols.split(",")]
    return {r[0]: dict(zip(keys, r)) for r in rows}


class GenomeEquivalentTest(unittest.TestCase):

    def test_no_silent_default_for_metaphlan(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            m, e, a = build(root, [("A", "Species A", 0.03), ("B", "Species B", 0.02)])
            done = run(m, e, a, root / "out")
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("--metaphlan-reference", done.stdout + done.stderr)

    def test_bracken_unchanged_and_read_proportional(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            m, e, a = build(root, [("A", "Species A", 0.03), ("B", "Species B", 0.02)],
                            profiler="kraken2_bracken")
            out = root / "out"
            self.assertEqual(run(m, e, a, out).returncode, 0)
            got = responses(out)
            # Non-implanted Bracken feature: pure dilution, no renormalisation.
            self.assertAlmostEqual(got["Bystander"]["expected_abundance_fraction"],
                                   (1 - 0.05) * 0.10)
            self.assertEqual(got["Bystander"]["reference_type"], "read_proportional")
            for feature in got:
                self.assertAlmostEqual(
                    got[feature]["expected_abundance_profiler_scale"],
                    got[feature]["expected_abundance_fraction"],
                    msg="Bracken profiler scale must equal the read reference")

    def test_metaphlan_genome_equivalent_target_and_non_implanted(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            targets = [("A", "Species A", 0.03), ("B", "Species B", 0.02)]
            g = {"A": 2_000_000, "B": 6_000_000}      # unequal genome sizes
            m, e, a = build(root, targets)
            out = root / "out"
            done = run(m, e, a, out, ["--metaphlan-reference", "genome_equivalent",
                                      "--target-genome-sizes", sizes(root, g),
                                      "--effective-genome-sizes", geff(root)])
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            got = responses(out)
            F = 0.05
            q = {k: frac * GEFF / g[k] for k, (_, _, frac) in
                 zip(["A", "B"], targets)}
            Q = sum(q.values())
            denom = (1 - F) + Q
            for feature in got:
                self.assertEqual(got[feature]["reference_type"], "genome_equivalent")
                self.assertAlmostEqual(
                    got[feature]["total_implanted_genome_equivalent_fraction"], Q)
                self.assertAlmostEqual(
                    got[feature]["effective_community_genome_size_bp"], GEFF)
            # Non-implanted: e = (1-F)o / ((1-F)+Q), strictly below plain dilution.
            self.assertAlmostEqual(got["Bystander"]["expected_abundance_profiler_scale"],
                                   (1 - F) * 0.10 / denom)
            self.assertLess(got["Bystander"]["expected_abundance_profiler_scale"],
                            got["Bystander"]["read_proportional_reference"])
            # Implanted target: e = ((1-F)o + q) / ((1-F)+Q); baseline absent -> o=0.
            self.assertAlmostEqual(got["Species A"]["expected_abundance_profiler_scale"],
                                   q["A"] / denom)
            self.assertAlmostEqual(got["Species A"]["target_genome_size_bp"], g["A"])
            # Smaller genome -> larger genome-equivalent fraction at equal reads.
            self.assertGreater(got["Species A"]["implanted_genome_equivalent_fraction"]
                               / 0.03,
                               got["Species B"]["implanted_genome_equivalent_fraction"]
                               / 0.02)
            # Read-proportional fields retain their original meaning.
            self.assertAlmostEqual(got["Bystander"]["expected_abundance_fraction"],
                                   (1 - F) * 0.10)

    def test_equal_genome_sizes_reduce_to_scaled_read_fraction(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            targets = [("A", "Species A", 0.03), ("B", "Species B", 0.02)]
            m, e, a = build(root, targets)
            out = root / "out"
            self.assertEqual(run(m, e, a, out, [
                "--metaphlan-reference", "genome_equivalent",
                "--target-genome-sizes", sizes(root, {"A": GEFF, "B": GEFF}),
                "--effective-genome-sizes", geff(root)]).returncode, 0)
            got = responses(out)
            # G_t == G_eff => q == f => Q == F, denominator == 1.
            self.assertAlmostEqual(
                got["Bystander"]["expected_abundance_profiler_scale"],
                got["Bystander"]["read_proportional_reference"])

    def test_individual_perturbation(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            m, e, a = build(root, [("A", "Species A", 0.01)], population="independent")
            out = root / "out"
            done = run(m, e, a, out, [
                "--metaphlan-reference", "genome_equivalent",
                "--target-genome-sizes", sizes(root, {"A": 2_000_000}),
                "--effective-genome-sizes", geff(root)], n_targets=2)
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            got = responses(out)
            Q = 0.01 * GEFF / 2_000_000
            self.assertAlmostEqual(
                got["Bystander"]["expected_abundance_profiler_scale"],
                (1 - 0.01) * 0.10 / ((1 - 0.01) + Q))

    def test_zero_baseline_abundance(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            m, e, a = build(root, [("A", "Species A", 0.03), ("B", "Species B", 0.02)],
                            baseline_bystander=0.0, observed_bystander=0.0)
            out = root / "out"
            self.assertEqual(run(m, e, a, out, [
                "--metaphlan-reference", "genome_equivalent",
                "--target-genome-sizes", sizes(root, {"A": 2e6, "B": 6e6}),
                "--effective-genome-sizes", geff(root)]).returncode, 0)
            got = responses(out)
            self.assertAlmostEqual(
                got["Bystander"]["expected_abundance_profiler_scale"], 0.0)
            self.assertAlmostEqual(
                got["Bystander"]["retained_baseline_profiler_scale"], 0.0)

    def test_missing_target_genome_size_fails(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            m, e, a = build(root, [("A", "Species A", 0.03), ("B", "Species B", 0.02)])
            done = run(m, e, a, root / "out", [
                "--metaphlan-reference", "genome_equivalent",
                "--target-genome-sizes", sizes(root, {"A": 2_000_000}),
                "--effective-genome-sizes", geff(root)])
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("lack a genome size", done.stdout + done.stderr)

    def test_missing_geff_fails(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            m, e, a = build(root, [("A", "Species A", 0.03), ("B", "Species B", 0.02)])
            done = run(m, e, a, root / "out", [
                "--metaphlan-reference", "genome_equivalent",
                "--target-genome-sizes", sizes(root, {"A": 2e6, "B": 6e6}),
                "--effective-genome-sizes", geff(root, sample="OTHER")])
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("lack an effective genome size", done.stdout + done.stderr)

    def test_incomplete_community_membership_fails(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            # Two implanted members present but ten declared.
            m, e, a = build(root, [("A", "Species A", 0.03), ("B", "Species B", 0.02)])
            done = run(m, e, a, root / "out", [
                "--metaphlan-reference", "genome_equivalent",
                "--target-genome-sizes", sizes(root, {"A": 2e6, "B": 6e6}),
                "--effective-genome-sizes", geff(root)], n_targets=10)
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("malformed physical perturbation", done.stdout + done.stderr)

    def test_achieved_total_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            m, e, a = build(root, [("A", "Species A", 0.03), ("B", "Species B", 0.02)],
                            total=0.5)   # recorded total disagrees with sum(f)
            done = run(m, e, a, root / "out", [
                "--metaphlan-reference", "genome_equivalent",
                "--target-genome-sizes", sizes(root, {"A": 2e6, "B": 6e6}),
                "--effective-genome-sizes", geff(root)])
            self.assertNotEqual(done.returncode, 0)

    def test_detection_fields_invariant_to_reference(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            targets = [("A", "Species A", 0.03), ("B", "Species B", 0.02)]
            m, e, a = build(root, targets)
            read_out, ge_out = root / "read", root / "ge"
            self.assertEqual(run(m, e, a, read_out, [
                "--metaphlan-reference", "read_proportional"]).returncode, 0)
            self.assertEqual(run(m, e, a, ge_out, [
                "--metaphlan-reference", "genome_equivalent",
                "--target-genome-sizes", sizes(root, {"A": 2e6, "B": 6e6}),
                "--effective-genome-sizes", geff(root)]).returncode, 0)
            read_rows, ge_rows = responses(read_out), responses(ge_out)
            for feature in read_rows:
                for field in ("observed_detected", "baseline_detected"):
                    self.assertEqual(read_rows[feature][field], ge_rows[feature][field],
                                     f"detection changed for {feature}.{field}")
                self.assertEqual(read_rows[feature]["expected_abundance_fraction"],
                                 ge_rows[feature]["expected_abundance_fraction"],
                                 "read-proportional field must not change meaning")


if __name__ == "__main__":
    unittest.main()
