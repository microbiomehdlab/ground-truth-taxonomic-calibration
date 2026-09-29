#!/usr/bin/env python3
"""Fixture suite for the upstream evidence package builder.

Every case drives the real builder as a subprocess against synthetic
`production_seal_v2` seals shaped exactly like the ones
`seal_cohort_upstream.py` writes. No real cluster output is touched, and
nothing here reimplements the builder's logic.

The 15 cases required by UPSTREAM_EVIDENCE_PACKAGE_SPEC.md section 10 map onto
the `test_caseNN_*` methods; the remaining classes cover the additional gates.
"""
from __future__ import annotations

import csv
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILDER = ROOT / "analysis_v2/scripts/build_upstream_evidence_package.py"
PLOTTER = ROOT / "analysis_v2/scripts/plot_upstream_evidence.R"

SEAL_CONTRACT = "upstream_seal_v2"
COMMIT = "a" * 40
CONDITIONS = ("Control", "Adenoma", "CRC")
DESIGNS = ("baseline", "independent", "community")
COHORTS = {
    "yachida": {"study": "YachidaS_2019", "samples": 201, "community": 1407},
    "feng": {"study": "Public_study__FengQ_2015", "samples": 154, "community": 1078},
    "zeller": {"study": "Public_study__ZellerG_2014", "samples": 156, "community": 1092},
}
INDEPENDENT = 30
INDEPENDENT_PROFILES = 1800
SAMPLE_FLOW_FIELDS = [
    "sample_id", "study", "condition", "independent_subset",
    "batch_id",
    "expected_baseline_profiles", "observed_baseline_profiles",
    "expected_independent_profiles", "observed_independent_profiles",
    "expected_community_profiles", "observed_community_profiles",
    "expected_profiles", "observed_profiles", "retained_output_files",
    "verified_marker", "retained_output_receipt", "input_provenance",
    "input_provenance_mode", "provenance_error", "sample_success",
    "baseline_profile_count", "independent_profile_count",
    "community_profile_count",
    "completion_table", "manifest_independent_flag", "status",
    "failure_reasons", "receipt_error", "completion_error",
]
COVARIATE_FIELDS = ["field", "missing", "present", "distinct_nonmissing"]
CANONICAL = ["sample_id", "condition", "study", "independent_subset"]
SEAL_MEMBERS = ("sample_flow.tsv", "covariate_audit.tsv",
                "production_manifest.tsv",
                "production_manifest.independent.tsv",
                "source_seal_inventory.tsv", "SUCCESS",
                "production_seal.sha256")
PROVENANCE_FIELDS = ["category", "asset_id", "name", "version_or_release",
                     "sha256", "availability", "notes"]
DIGEST_CATEGORIES = ("source_code", "container_upstream",
                     "container_analysis", "database_kraken2",
                     "database_metaphlan", "host_reference")
PARAMETER_SENTINEL = "not_a_file"
REQUIRED_PARAMETERS = ("read_length", "bracken_threshold", "profiler_threads",
                       "spike_fractions", "pool_coverage", "seeds")
SPIKE_TARGETS = ("Bfrag", "Csym", "Dpne", "Fnuc", "Hhat",
                 "Pmic", "Pana", "Psto", "Porp", "Pint")
ANALYSIS_IMAGE_SHA256 = "c" * 64
LEDGER_FIELDS = ["cohort", "audit_job_id", "audit_date_utc", "state",
                 "exit_code", "samples", "seal_sha256", "notes"]
AUDIT_JOBS = {"yachida": "3097679", "feng": "3097680", "zeller": "3097681"}


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_tsv(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(fields)]
    for row in rows:
        lines.append("\t".join(str(row.get(name, "")) for name in fields))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_tsv(path):
    lines = Path(path).read_text(encoding="utf-8").rstrip("\n").split("\n")
    header = lines[0].split("\t")
    return header, [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def tree_digest(root):
    root = Path(root)
    return {str(p.relative_to(root)): sha256_file(p)
            for p in sorted(root.rglob("*")) if p.is_file()}


class Fixture:
    """Three synthetic production_seal_v2 seals plus the assembly experiment."""

    def __init__(self, root):
        self.root = Path(root)
        self.seals = {}
        self.manifests = {}
        self.independent_manifests = {}
        self.samples = {}
        self.conditions = {}
        for cohort in COHORTS:
            self.build_cohort(cohort)
        self.build_assembly()
        self.build_side_inputs()

    # ------------------------------------------------------------- cohorts
    def plan(self, cohort):
        """Sample ids, conditions and subset membership for one cohort."""
        total = COHORTS[cohort]["samples"]
        samples, conditions, independent = [], {}, []
        per_condition = {name: 0 for name in CONDITIONS}
        for index in range(total):
            sample = "%s_S%04d" % (cohort[:3].upper(), index)
            condition = CONDITIONS[index % 3]
            samples.append(sample)
            conditions[sample] = condition
            if per_condition[condition] < 10:
                per_condition[condition] += 1
                independent.append(sample)
        return samples, conditions, independent

    def build_cohort(self, cohort):
        samples, conditions, independent = self.plan(cohort)
        self.samples[cohort] = samples
        self.conditions[cohort] = conditions
        study = COHORTS[cohort]["study"]
        seal = self.root / "seals" / cohort
        seal.mkdir(parents=True, exist_ok=True)
        self.seals[cohort] = seal

        flow = []
        for sample in samples:
            in_subset = sample in independent
            independent_profiles = 60 if in_subset else 0
            total = 1 + independent_profiles + 7
            row = {name: "" for name in SAMPLE_FLOW_FIELDS}
            row.update({
                "sample_id": sample, "study": study,
                "condition": conditions[sample],
                "independent_subset": "1" if in_subset else "0",
                "batch_id": "", "expected_baseline_profiles": 1,
                "observed_baseline_profiles": 1,
                "expected_independent_profiles": independent_profiles,
                "observed_independent_profiles": independent_profiles,
                "expected_community_profiles": 7,
                "observed_community_profiles": 7,
                "expected_profiles": total, "observed_profiles": total,
                "retained_output_files": 2, "verified_marker": 1,
                "retained_output_receipt": 1, "input_provenance": 1,
                "input_provenance_mode": ("sealed_manifest_and_qc_receipt"
                                          if cohort == "yachida"
                                          else "state_file"),
                "provenance_error": "", "sample_success": 1,
                "baseline_profile_count": 1, "independent_profile_count": 1,
                "community_profile_count": 1, "completion_table": 1,
                "manifest_independent_flag": 1, "status": "PASS",
                "failure_reasons": "", "receipt_error": "",
                "completion_error": "",
            })
            flow.append(row)
        write_tsv(seal / "sample_flow.tsv", SAMPLE_FLOW_FIELDS, flow)

        write_tsv(seal / "covariate_audit.tsv", COVARIATE_FIELDS, [
            {"field": "condition", "missing": 0,
             "present": len(samples), "distinct_nonmissing": 3},
            {"field": "study", "missing": 0, "present": len(samples),
             "distinct_nonmissing": 1},
            {"field": "independent_subset", "missing": 0,
             "present": len(samples), "distinct_nonmissing": 2},
            {"field": "age", "missing": 2, "present": len(samples) - 2,
             "distinct_nonmissing": 25},
        ])

        manifest_rows = [
            {"sample_id": sample, "condition": conditions[sample],
             "study": study,
             "independent_subset": "1" if sample in independent else "0",
             "source_age": "55", "source_sex": "Female", "source_bmi": "23.4"}
            for sample in samples]
        fields = CANONICAL + ["source_age", "source_sex", "source_bmi"]
        write_tsv(seal / "production_manifest.tsv", fields, manifest_rows)
        write_tsv(seal / "production_manifest.independent.tsv", fields,
                  [row for row in manifest_rows
                   if row["sample_id"] in independent])
        write_tsv(seal / "source_seal_inventory.tsv",
                  ["member", "sha256", "bytes", "status", "source_audit_job"],
                  [{"member": "native_SUCCESS", "sha256": "0" * 64,
                    "bytes": 32, "status": "VERIFIED",
                    "source_audit_job": AUDIT_JOBS[cohort]}])

        (seal / "SUCCESS").write_text(
            "seal_contract\t%s\ncohort\t%s\nsamples\t%d\n"
            "independent_samples\t%d\nbaseline_profiles\t%d\n"
            "independent_profiles\t%d\ncommunity_profiles\t%d\nstatus\tPASS\n"
            % (SEAL_CONTRACT, cohort, len(samples), INDEPENDENT, len(samples),
               INDEPENDENT_PROFILES, COHORTS[cohort]["community"]),
            encoding="utf-8")
        self.reseal(cohort)

        manifests = self.root / "manifests"
        manifests.mkdir(parents=True, exist_ok=True)
        self.manifests[cohort] = manifests / ("%s_production.tsv" % cohort)
        self.independent_manifests[cohort] = (
            manifests / ("%s_independent.tsv" % cohort))
        write_tsv(self.manifests[cohort], ["sample_id", "condition"],
                  [{"sample_id": s, "condition": conditions[s]}
                   for s in samples])
        write_tsv(self.independent_manifests[cohort],
                  ["sample_id", "condition"],
                  [{"sample_id": s, "condition": conditions[s]}
                   for s in independent])

    def reseal(self, cohort):
        """Recompute the seal's own checksum manifest over its six members."""
        seal = self.seals[cohort]
        lines = []
        for member in SEAL_MEMBERS:
            if member == "production_seal.sha256":
                continue
            lines.append("%s  %s\n" % (sha256_file(seal / member), member))
        (seal / "production_seal.sha256").write_text("".join(lines),
                                                     encoding="utf-8")

    # -------------------------------------------------------- side inputs
    def build_assembly(self):
        seal = self.root / "seals" / "assembly"
        seal.mkdir(parents=True, exist_ok=True)
        self.assembly = seal
        (seal / "summary.tsv").write_text(
            "field\tvalue\narms\toriginal;clean\n", encoding="utf-8")
        # The real experiment seal's documented SUCCESS schema.
        (seal / "SUCCESS").write_text(
            "experiment\tyachida_assembly_choice_sensitivity\n"
            "samples\t30\nassembly_arms\t2\nfractions_per_arm\t6\n"
            "expected_profiles\t360\nobserved_profiles\t360\nstatus\tPASS\n",
            encoding="utf-8")
        sidecar_rows = []
        fractions = ("0.0001", "0.0005", "0.001", "0.005", "0.01", "0.05")
        for sample_index in range(30):
            sample = "assembly_%02d" % sample_index
            for label, clean_label in (("Pana", "Pana_clean_GCA_000381525.1"),
                                       ("Pint", "Pint_clean_GCA_001953955.1")):
                for fraction_index, fraction in enumerate(fractions):
                    seed = str(100000 + sample_index * 12 + fraction_index)
                    sidecar_rows.append({
                        "sample_id": sample, "study": "YachidaS_2019",
                        "original_label": label, "clean_label": clean_label,
                        "fraction": fraction, "original_seed": seed,
                        "clean_seed": seed, "status": "PASS",
                    })
        write_tsv(seal / "matched_seed_audit.tsv",
                  ["sample_id", "study", "original_label", "clean_label",
                   "fraction", "original_seed", "clean_seed", "status"],
                  sidecar_rows[:126])
        self.reseal_assembly()

    def reseal_assembly(self):
        lines = []
        for member in ("SUCCESS", "summary.tsv"):
            lines.append("%s  %s\n" % (sha256_file(self.assembly / member),
                                       member))
        (self.assembly / "experiment_inputs_and_summary.sha256").write_text(
            "".join(lines), encoding="utf-8")

    def build_side_inputs(self):
        side = self.root / "inputs"
        side.mkdir(parents=True, exist_ok=True)
        self.spike_panel = side / "spike_panel.tsv"
        write_tsv(self.spike_panel, ["label", "taxon_name", "assembly",
                                     "weight"],
                  [{"label": label, "taxon_name": "Species %s" % label,
                    "assembly": "GCF_%03d.1" % index, "weight": "1"}
                   for index, label in enumerate(SPIKE_TARGETS, start=1)])
        # Genuinely comma-delimited, exactly like the frozen alias file.
        self.taxon_aliases = side / "spike_taxon_aliases.csv"
        self.taxon_aliases.write_text(
            "canonical,alias,tool,spike_label\n"
            + "".join("Species %s,Species %s,metaphlan4,\n" % (label, label)
                      for label in SPIKE_TARGETS)
            + "".join("Species %s,Species %s alias,kraken2_bracken,\n" % (label, label)
                      for label in SPIKE_TARGETS), encoding="utf-8")
        self.provenance = side / "upstream_provenance_metadata.tsv"
        rows = []
        for category in DIGEST_CATEGORIES:
            rows.append({
                "category": category, "asset_id": "%s_id" % category,
                "name": ("%s name" % category if category != "source_code"
                         else "repository at %s" % COMMIT),
                "version_or_release": "v1",
                "sha256": (ANALYSIS_IMAGE_SHA256
                           if category == "container_analysis" else "b" * 64),
                "availability": "public", "notes": "fixture"})
        for label in SPIKE_TARGETS:
            rows.append({
                "category": "spike_reference", "asset_id": label,
                "name": "Species %s assembly" % label,
                "version_or_release": "GCF_v1", "sha256": "d" * 64,
                "availability": "public", "notes": "fixture"})
        for parameter in REQUIRED_PARAMETERS:
            rows.append({
                "category": "upstream_parameter", "asset_id": parameter,
                "name": parameter, "version_or_release": "frozen",
                "sha256": PARAMETER_SENTINEL, "availability": "documented",
                "notes": "fixture"})
        write_tsv(self.provenance, PROVENANCE_FIELDS, rows)
        self.ledger = side / "audit_ledger.tsv"
        write_tsv(self.ledger, LEDGER_FIELDS,
                  [{"cohort": cohort, "audit_job_id": AUDIT_JOBS[cohort],
                    "audit_date_utc": "2026-09-27", "state": "COMPLETED",
                    "exit_code": "0:0",
                    "samples": COHORTS[cohort]["samples"],
                    "seal_sha256": sha256_file(
                        self.seals[cohort] / "production_seal.sha256"),
                    "notes": "unified upstream_seal_v2 audit"}
                   for cohort in COHORTS])

    # --------------------------------------------------------------- run
    def command(self, outdir, extra=()):
        command = [sys.executable, str(BUILDER),
                   "--assembly-sensitivity-seal", str(self.assembly),
                   "--spike-panel", str(self.spike_panel),
                   "--taxon-aliases", str(self.taxon_aliases),
                   "--provenance-metadata", str(self.provenance),
                   "--audit-ledger", str(self.ledger),
                   "--source-commit", COMMIT,
                   "--analysis-image-sha256", ANALYSIS_IMAGE_SHA256,
                   "--plotter", str(PLOTTER),
                   "--outdir", str(outdir)]
        for cohort in COHORTS:
            command += ["--%s-seal" % cohort, str(self.seals[cohort]),
                        "--%s-manifest" % cohort, str(self.manifests[cohort]),
                        "--%s-independent-manifest" % cohort,
                        str(self.independent_manifests[cohort])]
        command.extend(extra)
        return command

    def run(self, outdir, extra=()):
        return subprocess.run(self.command(outdir, extra), text=True,
                              capture_output=True)


class PackageTestCase(unittest.TestCase):

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self._count = 0

    def tearDown(self):
        self._tmp.cleanup()

    def fixture(self):
        self._count += 1
        directory = self.root / ("fix%02d" % self._count)
        directory.mkdir()
        return Fixture(directory)

    def outdir(self, name="package"):
        return self.root / name

    def build_ok(self, fixture, name="package", extra=()):
        out = self.outdir(name)
        done = fixture.run(out, extra)
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertTrue((out / "SUCCESS").is_file())
        self.assertFalse(out.with_name(out.name + ".incomplete").exists(),
                         "the staging directory must not survive")
        return out

    def build_fails(self, fixture, needle, name="bad", extra=()):
        out = self.outdir(name)
        done = fixture.run(out, extra)
        self.assertNotEqual(done.returncode, 0,
                            "builder accepted bad input (%s)" % needle)
        message = done.stdout + done.stderr
        self.assertIn(needle, message)
        self.assertFalse((out / "SUCCESS").exists())
        self.assertFalse(out.with_name(out.name + ".incomplete").exists(),
                         "a failed build must clean up its staging directory")
        return message


# ------------------------------------------------------------- cases 1-15
class SpecificationTest(PackageTestCase):

    def test_case01_valid_three_cohort_fixture_succeeds(self):
        fixture = self.fixture()
        out = self.build_ok(fixture)
        success = dict(line.split("\t") for line in
                       (out / "SUCCESS").read_text().rstrip("\n").split("\n"))
        self.assertEqual(success["status"], "PASS")
        self.assertEqual(success["seal_contract"], SEAL_CONTRACT)
        self.assertEqual(success["source_commit"], COMMIT)
        self.assertEqual(success["yachida_samples"], "201")
        self.assertEqual(success["feng_samples"], "154")
        self.assertEqual(success["zeller_samples"], "156")
        self.assertEqual(success["assembly_sensitivity_profiles"], "360")
        self.assertGreater(int(success["files"]), 20)
        for relative in ("MANIFEST.tsv", "SHA256SUMS", "README.md",
                         "zenodo/cohort_registry.tsv",
                         "zenodo/cohort_sample_flow.tsv",
                         "zenodo/profile_completeness.tsv",
                         "zenodo/covariate_missingness.tsv",
                         "zenodo/assembly_sensitivity_completeness.tsv",
                         "zenodo/source_seal_inventory.tsv",
                         "zenodo/audit_ledger.tsv",
                         "zenodo/upstream_provenance_metadata.tsv",
                         "zenodo/spike_panel.tsv", "zenodo/taxon_aliases.tsv",
                         "manuscript_supplement/figure_source_data.tsv",
                         "manuscript_supplement/supplementary_sample_flow.tsv",
                         "provenance/build_parameters.tsv",
                         "provenance/input_checksums.tsv",
                         "provenance/repository_commit.txt",
                         "provenance/validation_report.tsv",
                         "source_seal_projections/yachida/SUCCESS",
                         "source_seal_projections/yachida/projection.sha256",
                         "source_seal_projections/yachida/README.md",
                         "source_seal_projections/yachida/"
                         "original_source_seal_inventory.tsv",
                         "source_seal_projections/yachida_assembly_sensitivity/"
                         "SUCCESS"):
            self.assertTrue((out / relative).is_file(), relative)

    def test_case02_missing_sample_fails(self):
        fixture = self.fixture()
        seal = fixture.seals["feng"]
        header, rows = read_tsv(seal / "sample_flow.tsv")
        write_tsv(seal / "sample_flow.tsv", header, rows[:-1])
        fixture.reseal("feng")
        self.build_fails(fixture, "cover different samples")

    def test_case03_wrong_condition_balance_fails(self):
        fixture = self.fixture()
        seal = fixture.seals["zeller"]
        header, rows = read_tsv(seal / "production_manifest.independent.tsv")
        rows[0]["condition"] = "CRC"
        write_tsv(seal / "production_manifest.independent.tsv", header, rows)
        fixture.reseal("zeller")
        self.build_fails(fixture, "independent subset is not balanced")

    def test_case04_mismatched_identity_with_right_count_fails(self):
        fixture = self.fixture()
        seal = fixture.seals["yachida"]
        for member in ("sample_flow.tsv", "production_manifest.tsv"):
            header, rows = read_tsv(seal / member)
            rows[5]["sample_id"] = "IMPOSTOR_0005"
            write_tsv(seal / member, header, rows)
        fixture.reseal("yachida")
        message = self.build_fails(fixture, "cover different samples")
        self.assertIn("IMPOSTOR_0005", message)

    def test_case05_wrong_profile_topology_fails(self):
        fixture = self.fixture()
        seal = fixture.seals["feng"]
        header, rows = read_tsv(seal / "sample_flow.tsv")
        rows[0]["observed_community_profiles"] = "6"
        write_tsv(seal / "sample_flow.tsv", header, rows)
        fixture.reseal("feng")
        self.build_fails(fixture, "observed 6 community profiles but expected 7")

    def test_case06_corrupted_seal_member_fails(self):
        fixture = self.fixture()
        member = fixture.seals["zeller"] / "covariate_audit.tsv"
        member.write_text(member.read_text(encoding="utf-8") + "tampered\t0\t0\t0\n",
                          encoding="utf-8")
        self.build_fails(fixture, "checksum mismatch")

    def test_case07_missing_or_stale_seal_fails(self):
        missing = self.fixture()
        (missing.seals["feng"] / "SUCCESS").unlink()
        self.build_fails(missing, "lacks member(s): SUCCESS")

        stale = self.fixture()
        (stale.seals["feng"] / "AUDIT_IN_PROGRESS").write_text(
            "status\tIN_PROGRESS\n", encoding="utf-8")
        self.build_fails(stale, "mid-audit", name="stale")

        absent = self.fixture()
        shutil.rmtree(str(absent.seals["zeller"]))
        self.build_fails(absent, "is not a directory", name="absent")

    def test_case08_duplicate_sample_across_cohorts_fails(self):
        fixture = self.fixture()
        shared = fixture.samples["yachida"][0]
        original = fixture.samples["feng"][0]
        seal = fixture.seals["feng"]
        # Rename consistently everywhere, so only the cross-cohort collision
        # is wrong and no earlier nesting or identity gate fires first.
        for member in ("sample_flow.tsv", "production_manifest.tsv",
                       "production_manifest.independent.tsv"):
            header, rows = read_tsv(seal / member)
            for row in rows:
                if row["sample_id"] == original:
                    row["sample_id"] = shared
            write_tsv(seal / member, header, rows)
        for path in (fixture.manifests["feng"],
                     fixture.independent_manifests["feng"]):
            header, rows = read_tsv(path)
            for row in rows:
                if row["sample_id"] == original:
                    row["sample_id"] = shared
            write_tsv(path, header, rows)
        fixture.reseal("feng")
        message = self.build_fails(fixture, "appears in both")
        self.assertIn(shared, message)

    def test_case09_incomplete_assembly_sensitivity_fails(self):
        for field, value in (("observed_profiles", "359"), ("samples", "29"),
                             ("assembly_arms", "1"), ("fractions_per_arm", "5")):
            fixture = self.fixture()
            text = (fixture.assembly / "SUCCESS").read_text(encoding="utf-8")
            lines = [line if not line.startswith(field + "\t")
                     else "%s\t%s" % (field, value)
                     for line in text.rstrip("\n").split("\n")]
            (fixture.assembly / "SUCCESS").write_text(
                "\n".join(lines) + "\n", encoding="utf-8")
            fixture.reseal_assembly()
            self.build_fails(fixture, "the frozen experiment requires",
                             name="assembly_%s" % field)

    def test_case10_blank_provenance_identity_fails(self):
        for field in ("asset_id", "name", "version_or_release", "availability"):
            fixture = self.fixture()
            header, rows = read_tsv(fixture.provenance)
            rows[0][field] = ""
            write_tsv(fixture.provenance, header, rows)
            self.build_fails(fixture, "has a blank %s" % field,
                             name="blank_%s" % field)

        missing = self.fixture()
        header, rows = read_tsv(missing.provenance)
        write_tsv(missing.provenance, header,
                  [row for row in rows if row["category"] != "host_reference"])
        self.build_fails(missing, "lacks required categor", name="nocat")

    def test_case11_no_absolute_paths_in_release_outputs(self):
        out = self.build_ok(self.fixture())
        offenders = []
        for path in sorted(out.rglob("*")):
            if not path.is_file() or path.suffix not in (".tsv", ".md", ".txt"):
                continue
            text = path.read_text(encoding="utf-8")
            for token in text.replace("\t", " ").split():
                if token.startswith("/") and token.count("/") > 1:
                    offenders.append("%s: %s" % (path.name, token))
        self.assertEqual(offenders, [], "absolute paths leaked into the package")
        # The build directory itself is never named in a release table.
        for name in ("zenodo/cohort_registry.tsv",
                     "provenance/input_checksums.tsv"):
            self.assertNotIn(str(self.root), (out / name).read_text())

    def test_case11b_cluster_path_in_supplied_metadata_is_refused(self):
        """The realistic leak vector: an operator pastes a path into a note.

        This failure also happens after the staging directory exists, so it
        proves the staging cleanup as well.
        """
        for table, field in (("provenance", "notes"), ("ledger", "notes")):
            fixture = self.fixture()
            path = getattr(fixture, table)
            header, rows = read_tsv(path)
            rows[0][field] = ("copied from "
                              "/mnt/nfs/microbiomehd/crc-lab/work/seal")
            write_tsv(path, header, rows)
            out = self.outdir("leak_%s" % table)
            done = fixture.run(out)
            self.assertNotEqual(done.returncode, 0,
                                "a cluster path reached a release table")
            message = done.stdout + done.stderr
            self.assertIn("absolute paths", message)
            self.assertIn("/mnt/nfs/", message)
            self.assertFalse(out.exists())
            self.assertFalse(out.with_name(out.name + ".incomplete").exists(),
                             "a post-staging failure must clean up")

    def test_case12_attempted_overwrite_fails(self):
        fixture = self.fixture()
        out = self.build_ok(fixture)
        before = tree_digest(out)
        done = fixture.run(out)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("already exists and is not empty", done.stdout + done.stderr)
        self.assertEqual(tree_digest(out), before,
                         "a refused rebuild must not disturb the package")

    def test_case13_source_seals_unchanged_on_success_and_failure(self):
        fixture = self.fixture()
        before = tree_digest(self.root / "fix01" / "seals")
        self.assertTrue(before)
        self.build_ok(fixture)
        self.assertEqual(tree_digest(self.root / "fix01" / "seals"), before,
                         "a successful build modified a source seal")
        header, rows = read_tsv(fixture.provenance)
        rows[0]["availability"] = ""
        write_tsv(fixture.provenance, header, rows)
        self.build_fails(fixture, "has a blank availability", name="p2")
        self.assertEqual(tree_digest(self.root / "fix01" / "seals"), before,
                         "a failed build modified a source seal")

    def test_case14_tampering_breaks_sha256sum_c(self):
        out = self.build_ok(self.fixture())
        if shutil.which("sha256sum") is None:
            self.skipTest("sha256sum is unavailable")
        clean = subprocess.run(["sha256sum", "-c", "--quiet", "SHA256SUMS"],
                               cwd=str(out), capture_output=True, text=True)
        self.assertEqual(clean.returncode, 0, clean.stdout + clean.stderr)
        target = out / "zenodo" / "cohort_registry.tsv"
        target.write_text(target.read_text(encoding="utf-8") + "tampered\n",
                          encoding="utf-8")
        dirty = subprocess.run(["sha256sum", "-c", "--quiet", "SHA256SUMS"],
                               cwd=str(out), capture_output=True, text=True)
        self.assertNotEqual(dirty.returncode, 0)
        self.assertIn("cohort_registry.tsv", dirty.stdout + dirty.stderr)

        out = self.build_ok(self.fixture(), name="tamper_success")
        success = out / "SUCCESS"
        success.write_text(success.read_text(encoding="utf-8") + "tampered\n",
                           encoding="utf-8")
        dirty = subprocess.run(["sha256sum", "-c", "--quiet", "SHA256SUMS"],
                               cwd=str(out), capture_output=True, text=True)
        self.assertNotEqual(dirty.returncode, 0)
        self.assertIn("SUCCESS", dirty.stdout + dirty.stderr)

    def test_case15_figure_source_data_matches_the_tables(self):
        out = self.build_ok(self.fixture())
        _, source = read_tsv(out / "manuscript_supplement/figure_source_data.tsv")
        indexed = {(r["panel"], r["cohort"], r["group"], r["measure"]): r["value"]
                   for r in source}
        self.assertEqual(len(indexed), len(source), "duplicate source rows")

        _, flow = read_tsv(out / "zenodo/cohort_sample_flow.tsv")
        for row in flow:
            if row["condition"] == "TOTAL":
                continue
            self.assertEqual(
                indexed[("A", row["cohort"], row["condition"], "samples")],
                row["manifest_samples"])
        _, completeness = read_tsv(out / "zenodo/profile_completeness.tsv")
        for row in completeness:
            for measure in ("expected_profiles", "observed_profiles"):
                self.assertEqual(
                    indexed[("B", row["cohort"], row["design"], measure)],
                    row[measure])
        _, balance = read_tsv(out / "zenodo/independent_subset_balance.tsv")
        for row in balance:
            for measure in ("expected_samples", "observed_samples"):
                self.assertEqual(
                    indexed[("C", row["cohort"], row["condition"], measure)],
                    row[measure])
        _, covariates = read_tsv(out / "zenodo/covariate_missingness.tsv")
        for row in covariates:
            self.assertEqual(
                indexed[("D", row["cohort"], row["field"], "missing")],
                row["missing"])
        # Every source row is accounted for by one of the four panels.
        self.assertEqual({row["panel"] for row in source}, {"A", "B", "C", "D"})


# --------------------------------------------------------- additional gates
class AdditionalGateTest(PackageTestCase):

    def test_wrong_cohort_identity_fails(self):
        fixture = self.fixture()
        seal = fixture.seals["feng"]
        text = (seal / "SUCCESS").read_text(encoding="utf-8")
        (seal / "SUCCESS").write_text(
            text.replace("cohort\tfeng", "cohort\tzeller"), encoding="utf-8")
        fixture.reseal("feng")
        self.build_fails(fixture, "SUCCESS declares cohort 'zeller'")

    def test_wrong_seal_contract_fails(self):
        fixture = self.fixture()
        seal = fixture.seals["yachida"]
        text = (seal / "SUCCESS").read_text(encoding="utf-8")
        (seal / "SUCCESS").write_text(
            text.replace(SEAL_CONTRACT, "legacy_seal_v1"), encoding="utf-8")
        fixture.reseal("yachida")
        self.build_fails(fixture, "declares contract")

    def test_non_pass_sample_row_fails(self):
        fixture = self.fixture()
        seal = fixture.seals["zeller"]
        header, rows = read_tsv(seal / "sample_flow.tsv")
        rows[3]["status"] = "FAIL"
        write_tsv(seal / "sample_flow.tsv", header, rows)
        fixture.reseal("zeller")
        message = self.build_fails(fixture, "is not PASS")
        self.assertIn("ZEL_S0003", message)

        # A recorded failure reason is refused even if status still says PASS.
        reasoned = self.fixture()
        header, rows = read_tsv(reasoned.seals["zeller"] / "sample_flow.tsv")
        rows[3]["failure_reasons"] = "completion_table"
        write_tsv(reasoned.seals["zeller"] / "sample_flow.tsv", header, rows)
        reasoned.reseal("zeller")
        self.build_fails(reasoned, "nonblank failure_reasons", name="reasoned")

        # So is a cleared boolean audit field.
        cleared = self.fixture()
        header, rows = read_tsv(cleared.seals["zeller"] / "sample_flow.tsv")
        rows[2]["input_provenance"] = "0"
        write_tsv(cleared.seals["zeller"] / "sample_flow.tsv", header, rows)
        cleared.reseal("zeller")
        self.build_fails(cleared, "has input_provenance='0', not 1",
                         name="cleared")

    def test_sample_flow_schema_drift_fails(self):
        fixture = self.fixture()
        seal = fixture.seals["feng"]
        header, rows = read_tsv(seal / "sample_flow.tsv")
        header = [name for name in header if name != "provenance_error"]
        write_tsv(seal / "sample_flow.tsv", header, rows)
        fixture.reseal("feng")
        self.build_fails(fixture, "sample_flow.tsv schema drifted")

    def test_covariate_schema_drift_fails(self):
        fixture = self.fixture()
        seal = fixture.seals["feng"]
        header, rows = read_tsv(seal / "covariate_audit.tsv")
        write_tsv(seal / "covariate_audit.tsv", header + ["extra"], rows)
        fixture.reseal("feng")
        self.build_fails(fixture, "covariate_audit.tsv schema drifted")

    def test_canonical_columns_must_lead_in_order(self):
        for column in CANONICAL[1:]:
            fixture = self.fixture()
            seal = fixture.seals["yachida"]
            header, rows = read_tsv(seal / "production_manifest.tsv")
            write_tsv(seal / "production_manifest.tsv",
                      [name for name in header if name != column], rows)
            fixture.reseal("yachida")
            self.build_fails(fixture, "must begin with sample_id, condition, "
                                      "study, independent_subset",
                             name="canon_%s" % column)

        # Present but out of order is also refused.
        reordered = self.fixture()
        seal = reordered.seals["yachida"]
        header, rows = read_tsv(seal / "production_manifest.tsv")
        swapped = ["sample_id", "study", "condition", "independent_subset"]
        write_tsv(seal / "production_manifest.tsv",
                  swapped + [n for n in header if n not in swapped], rows)
        reordered.reseal("yachida")
        self.build_fails(reordered, "must begin with", name="reordered")

    def test_unsafe_checksum_member_fails(self):
        for member in ("/etc/passwd", "../escape.tsv"):
            fixture = self.fixture()
            seal = fixture.seals["yachida"]
            checksum = seal / "production_seal.sha256"
            checksum.write_text(
                checksum.read_text(encoding="utf-8")
                + "%s  %s\n" % ("0" * 64, member), encoding="utf-8")
            self.build_fails(fixture, "unsafe member path",
                             name="unsafe_%s" % member.strip("/.").replace("/", "_"))

    def test_duplicate_checksum_member_fails(self):
        fixture = self.fixture()
        checksum = fixture.seals["feng"] / "production_seal.sha256"
        first = checksum.read_text(encoding="utf-8").splitlines()[0]
        checksum.write_text(
            checksum.read_text(encoding="utf-8") + first + "\n", encoding="utf-8")
        self.build_fails(fixture, "twice")

    def test_flow_and_manifest_condition_must_agree(self):
        fixture = self.fixture()
        seal = fixture.seals["feng"]
        header, rows = read_tsv(seal / "production_manifest.tsv")
        rows[0]["condition"] = "Healthy"
        write_tsv(seal / "production_manifest.tsv", header, rows)
        fixture.reseal("feng")
        self.build_fails(fixture, "in sample_flow.tsv but 'Healthy' in the "
                                  "canonical manifest")

    def test_unexpected_condition_label_fails(self):
        """Changed consistently in both, so the label gate is what refuses it."""
        fixture = self.fixture()
        seal = fixture.seals["feng"]
        for member in ("production_manifest.tsv", "sample_flow.tsv"):
            header, rows = read_tsv(seal / member)
            rows[0]["condition"] = "Healthy"
            write_tsv(seal / member, header, rows)
        fixture.reseal("feng")
        self.build_fails(fixture, "unexpected condition label")

    def test_short_commit_is_refused(self):
        fixture = self.fixture()
        out = self.outdir("shortcommit")
        command = fixture.command(out)
        command[command.index("--source-commit") + 1] = "abc123"
        done = subprocess.run(command, text=True, capture_output=True)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("40-character", done.stdout + done.stderr)

    def test_credential_in_provenance_is_refused(self):
        fixture = self.fixture()
        header, rows = read_tsv(fixture.provenance)
        rows[0]["notes"] = "password: hunter2"
        write_tsv(fixture.provenance, header, rows)
        self.build_fails(fixture, "credential or private key")

    def test_audit_ledger_must_cover_every_cohort(self):
        fixture = self.fixture()
        header, rows = read_tsv(fixture.ledger)
        write_tsv(fixture.ledger, header,
                  [row for row in rows if row["cohort"] != "zeller"])
        self.build_fails(fixture, "audit ledger lacks row(s) for: zeller")

    def test_registry_is_release_safe_and_complete(self):
        out = self.build_ok(self.fixture())
        header, rows = read_tsv(out / "zenodo/cohort_registry.tsv")
        self.assertEqual(header, [
            "cohort", "study", "samples", "independent_samples",
            "baseline_profiles", "community_profiles", "independent_profiles",
            "seal_contract", "seal_status", "input_provenance_mode",
            "source_audit_job"])
        registry = {row["cohort"]: row for row in rows}
        self.assertEqual(sorted(registry), ["feng", "yachida", "zeller"])
        for cohort, expected in COHORTS.items():
            row = registry[cohort]
            self.assertEqual(row["study"], expected["study"])
            self.assertEqual(row["samples"], str(expected["samples"]))
            self.assertEqual(row["independent_samples"], "30")
            self.assertEqual(row["baseline_profiles"], str(expected["samples"]))
            self.assertEqual(row["community_profiles"],
                             str(expected["community"]))
            self.assertEqual(row["independent_profiles"], "1800")
            self.assertEqual(row["seal_contract"], SEAL_CONTRACT)
            self.assertEqual(row["seal_status"], "PASS")
            self.assertEqual(row["source_audit_job"], AUDIT_JOBS[cohort])
        self.assertEqual(registry["yachida"]["input_provenance_mode"],
                         "sealed_manifest_and_qc_receipt")
        self.assertEqual(registry["feng"]["input_provenance_mode"], "state_file")

    def test_individual_covariates_are_withheld_by_default(self):
        out = self.build_ok(self.fixture())
        for path in sorted(out.rglob("*.tsv")):
            header, _ = read_tsv(path)
            for forbidden in ("source_age", "source_sex", "source_bmi",
                              "age", "sex", "bmi"):
                if path.name == "covariate_missingness.tsv":
                    continue  # aggregate rows name the field, not its values
                self.assertNotIn(forbidden, header,
                                 "%s exposes %s" % (path.name, forbidden))
        header, _ = read_tsv(
            out / "source_seal_projections/yachida/production_manifest.canonical.tsv")
        self.assertEqual(header, CANONICAL)
        self.assertFalse(
            (out / "source_seal_projections/yachida/production_manifest.tsv").exists())
        parameters = dict(
            (row["parameter"], row["value"]) for row in
            read_tsv(out / "provenance/build_parameters.tsv")[1])
        self.assertEqual(parameters["individual_covariates_released"], "no")

    def test_individual_covariates_require_a_review_reference(self):
        fixture = self.fixture()
        self.build_fails(fixture, "requires --redistribution-review",
                         extra=["--include-individual-covariates"])

    def test_manifest_and_sums_cover_the_package(self):
        out = self.build_ok(self.fixture())
        header, entries = read_tsv(out / "MANIFEST.tsv")
        self.assertEqual(header, ["relative_path", "bytes", "sha256",
                                  "media_type", "release_role"])
        listed = {row["relative_path"] for row in entries}
        self.assertNotIn("MANIFEST.tsv", listed)
        self.assertNotIn("SHA256SUMS", listed)
        self.assertIn("SUCCESS", listed)
        for row in entries:
            target = out / row["relative_path"]
            self.assertTrue(target.is_file(), row["relative_path"])
            self.assertEqual(row["sha256"], sha256_file(target))
            self.assertIn(row["release_role"],
                          {"manuscript_supplement", "zenodo_release",
                           "provenance_source_seal", "source_seal_projection",
                           "build_provenance",
                           "package_root"})
        summed = set()
        for line in (out / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
            recorded, relative = line.split("  ", 1)
            summed.add(relative)
            self.assertEqual(recorded, sha256_file(out / relative), relative)
        self.assertIn("MANIFEST.tsv", summed)
        self.assertIn("SUCCESS", summed)
        self.assertEqual(listed - summed, set())

    def test_staging_directory_is_removed_after_failure(self):
        fixture = self.fixture()
        out = self.outdir("atomic")
        header, rows = read_tsv(fixture.provenance)
        rows[0]["name"] = ""
        write_tsv(fixture.provenance, header, rows)
        done = fixture.run(out)
        self.assertNotEqual(done.returncode, 0)
        self.assertFalse(out.exists(), "no partial package may be promoted")
        self.assertFalse(out.with_name(out.name + ".incomplete").exists())

    def test_preexisting_staging_directory_is_refused(self):
        fixture = self.fixture()
        out = self.outdir("staged")
        staging = out.with_name(out.name + ".incomplete")
        staging.mkdir(parents=True)
        (staging / "leftover.tsv").write_text("x\n", encoding="utf-8")
        done = fixture.run(out)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("staging directory is still present",
                      done.stdout + done.stderr)
        self.assertTrue((staging / "leftover.tsv").is_file(),
                        "an operator's staging directory must be preserved")


# ---------------------------------------------------------------- figure
class FigureTest(PackageTestCase):

    def test_figure_is_rendered_or_explicitly_skipped(self):
        out = self.build_ok(self.fixture())
        parameters = dict(
            (row["parameter"], row["value"]) for row in
            read_tsv(out / "provenance/build_parameters.tsv")[1])
        validation = {row["check"]: row for row in
                      read_tsv(out / "provenance/validation_report.tsv")[1]}
        if shutil.which("Rscript") is None:
            self.assertEqual(parameters["figure_rendered"], "no")
            self.assertEqual(validation["figure"]["status"], "SKIPPED")
            self.skipTest("Rscript unavailable; structural checks still ran")
        if parameters["figure_rendered"] != "yes":
            self.assertEqual(validation["figure"]["status"], "SKIPPED")
            self.skipTest("ggplot2/patchwork unavailable: %s"
                          % parameters["figure_note"])
        for name in ("upstream_completeness.pdf", "upstream_completeness.svg",
                     "upstream_completeness.png"):
            target = out / "manuscript_supplement" / name
            self.assertTrue(target.is_file(), name)
            self.assertGreater(target.stat().st_size, 0)
        self.assertEqual(validation["figure"]["status"], "PASS")

    def test_require_figure_fails_when_r_is_unavailable(self):
        if shutil.which("Rscript") is not None:
            self.skipTest("Rscript is available; the skip path cannot be forced")
        self.build_fails(self.fixture(), "Rscript is unavailable",
                         extra=["--require-figure"])


# --------------------------------------------- Codex release-integrity gates
class PackagedSealCopyTest(PackageTestCase):
    """Item 1: a packaged directory must verify against its own manifest."""

    def test_every_projection_verifies_independently(self):
        out = self.build_ok(self.fixture())
        root = out / "source_seal_projections"
        self.assertTrue(root.is_dir())
        self.assertFalse((out / "source_seals").exists(),
                         "projection and exact-seal semantics must not mix")
        components = sorted(item.name for item in root.iterdir())
        self.assertEqual(components, ["feng", "yachida",
                                      "yachida_assembly_sensitivity", "zeller"])
        for component in components:
            directory = root / component
            manifest = directory / "projection.sha256"
            self.assertTrue(manifest.is_file(), component)
            # No original manifest may sit where it could not verify.
            self.assertFalse((directory / "production_seal.sha256").exists())
            self.assertFalse(
                (directory / "experiment_inputs_and_summary.sha256").exists())
            listed = {}
            for line in manifest.read_text(encoding="utf-8").splitlines():
                recorded, name = line.split("  ", 1)
                listed[name] = recorded
                self.assertEqual(sha256_file(directory / name), recorded,
                                 "%s/%s" % (component, name))
            present = sorted(item.name for item in directory.iterdir()
                             if item.name != "projection.sha256")
            self.assertEqual(sorted(listed), present,
                             "%s: manifest and contents disagree" % component)
            if shutil.which("sha256sum") is not None:
                done = subprocess.run(
                    ["sha256sum", "-c", "--quiet", "projection.sha256"],
                    cwd=str(directory), capture_output=True, text=True)
                self.assertEqual(done.returncode, 0,
                                 done.stdout + done.stderr)
            readme = (directory / "README.md").read_text(encoding="utf-8")
            self.assertIn("NOT the original authoritative seal", readme)

    def test_projection_records_every_original_member(self):
        fixture = self.fixture()
        out = self.build_ok(fixture)
        for cohort in COHORTS:
            header, rows = read_tsv(
                out / "source_seal_projections" / cohort
                / "original_source_seal_inventory.tsv")
            self.assertEqual(header, ["original_member", "sha256", "bytes",
                                      "copied_into_projection"])
            recorded = {row["original_member"]: row for row in rows}
            self.assertEqual(sorted(recorded), sorted(SEAL_MEMBERS))
            for member, row in sorted(recorded.items()):
                self.assertEqual(row["sha256"],
                                 sha256_file(fixture.seals[cohort] / member),
                                 member)
                self.assertEqual(
                    int(row["bytes"]),
                    (fixture.seals[cohort] / member).stat().st_size)
            self.assertEqual(recorded["production_manifest.tsv"]
                             ["copied_into_projection"], "no")
            self.assertEqual(recorded["SUCCESS"]["copied_into_projection"],
                             "yes")

    def test_approved_exact_copy_keeps_the_original_manifest_valid(self):
        out = self.build_ok(self.fixture(), name="exact", extra=[
            "--include-individual-covariates",
            "--redistribution-review", "REVIEW-2026-09-28"])
        self.assertFalse((out / "source_seal_projections").exists())
        root = out / "source_seals"
        note = (root / "README.md").read_text(encoding="utf-8")
        self.assertIn("byte-identical copy", note)
        for cohort in COHORTS:
            directory = root / cohort
            present = sorted(item.name for item in directory.iterdir())
            self.assertEqual(present, sorted(SEAL_MEMBERS),
                             "an exact copy must be exactly the seal")
            self.assertTrue((directory / "production_manifest.tsv").is_file())
            if shutil.which("sha256sum") is not None:
                done = subprocess.run(
                    ["sha256sum", "-c", "--quiet", "production_seal.sha256"],
                    cwd=str(directory), capture_output=True, text=True)
                self.assertEqual(done.returncode, 0,
                                 "the copied original manifest must verify")
        parameters = dict(
            (row["parameter"], row["value"]) for row in
            read_tsv(out / "provenance/build_parameters.tsv")[1])
        self.assertEqual(parameters["packaged_seal_mode"], "exact_copy")
        self.assertEqual(parameters["individual_covariates_released"], "yes")


class ExactSealMembershipTest(PackageTestCase):
    """Item 2: exactly the documented members, no more and no fewer."""

    def test_unexpected_file_in_a_seal_fails(self):
        fixture = self.fixture()
        (fixture.seals["feng"] / "notes.txt").write_text("x\n", encoding="utf-8")
        self.build_fails(fixture, "holds unexpected file(s): notes.txt")

    def test_unexpected_subdirectory_in_a_seal_fails(self):
        fixture = self.fixture()
        (fixture.seals["zeller"] / "extra").mkdir()
        self.build_fails(fixture, "unexpected subdirector")

    def test_checksum_manifest_must_cover_exactly_six_members(self):
        omitted = self.fixture()
        checksum = omitted.seals["yachida"] / "production_seal.sha256"
        kept = [line for line in checksum.read_text().splitlines()
                if "covariate_audit.tsv" not in line]
        checksum.write_text("\n".join(kept) + "\n", encoding="utf-8")
        self.build_fails(omitted, "checksum manifest omits covariate_audit.tsv")

        extra = self.fixture()
        seal = extra.seals["yachida"]
        (seal / "production_seal.sha256").write_text(
            (seal / "production_seal.sha256").read_text()
            + "%s  SUCCESS\n" % sha256_file(seal / "SUCCESS"),
            encoding="utf-8")
        self.build_fails(extra, "lists %s twice" % "SUCCESS", name="extra")

    def test_assembly_member_absent_from_its_manifest_fails(self):
        fixture = self.fixture()
        (fixture.assembly / "stray.tsv").write_text("x\n", encoding="utf-8")
        self.build_fails(fixture, "observed unlisted members")

    def test_historical_matched_seed_sidecar_is_required_and_validated(self):
        missing = self.fixture()
        (missing.assembly / "matched_seed_audit.tsv").unlink()
        self.build_fails(missing, "historical unsealed sidecar",
                         name="missing_seed_sidecar")

        mismatched = self.fixture()
        path = mismatched.assembly / "matched_seed_audit.tsv"
        header, rows = read_tsv(path)
        rows[0]["clean_seed"] = str(int(rows[0]["original_seed"]) + 1)
        write_tsv(path, header, rows)
        self.build_fails(mismatched, "does not preserve the seed",
                         name="mismatched_seed_sidecar")

    def test_assembly_success_must_be_checksummed(self):
        fixture = self.fixture()
        checksum = fixture.assembly / "experiment_inputs_and_summary.sha256"
        kept = [line for line in checksum.read_text().splitlines()
                if not line.endswith("  SUCCESS")]
        checksum.write_text("\n".join(kept) + "\n", encoding="utf-8")
        self.build_fails(fixture, "does not checksum its own SUCCESS")

    def test_assembly_success_schema_is_required(self):
        for field in ("experiment", "samples", "assembly_arms",
                      "fractions_per_arm", "expected_profiles",
                      "observed_profiles", "status"):
            fixture = self.fixture()
            text = (fixture.assembly / "SUCCESS").read_text(encoding="utf-8")
            kept = [line for line in text.rstrip("\n").split("\n")
                    if not line.startswith(field + "\t")]
            (fixture.assembly / "SUCCESS").write_text(
                "\n".join(kept) + "\n", encoding="utf-8")
            fixture.reseal_assembly()
            self.build_fails(fixture, "SUCCESS", name="asm_%s" % field)


class PerSampleTopologyTest(PackageTestCase):
    """Item 3: an aggregate-preserving per-sample error must still fail."""

    def test_offset_preserving_community_corruption_fails(self):
        fixture = self.fixture()
        seal = fixture.seals["feng"]
        header, rows = read_tsv(seal / "sample_flow.tsv")
        rows[0]["observed_community_profiles"] = "6"
        rows[0]["observed_profiles"] = str(int(rows[0]["observed_profiles"]) - 1)
        rows[1]["observed_community_profiles"] = "8"
        rows[1]["observed_profiles"] = str(int(rows[1]["observed_profiles"]) + 1)
        write_tsv(seal / "sample_flow.tsv", header, rows)
        fixture.reseal("feng")
        # The cohort total is still 1078, so only per-sample checking catches it.
        self.assertEqual(
            sum(int(r["observed_community_profiles"]) for r in rows),
            COHORTS["feng"]["community"])
        self.build_fails(fixture, "observed 6 community profiles but expected 7")

    def test_offset_preserving_independent_corruption_fails(self):
        fixture = self.fixture()
        seal = fixture.seals["yachida"]
        header, rows = read_tsv(seal / "sample_flow.tsv")
        subset = [row for row in rows if row["independent_subset"] == "1"]
        subset[0]["observed_independent_profiles"] = "59"
        subset[0]["observed_profiles"] = str(
            int(subset[0]["observed_profiles"]) - 1)
        subset[1]["observed_independent_profiles"] = "61"
        subset[1]["observed_profiles"] = str(
            int(subset[1]["observed_profiles"]) + 1)
        write_tsv(seal / "sample_flow.tsv", header, rows)
        fixture.reseal("yachida")
        self.assertEqual(
            sum(int(r["observed_independent_profiles"]) for r in rows),
            INDEPENDENT_PROFILES)
        self.build_fails(fixture, "observed 59 independent profiles")

    def test_offset_preserving_baseline_corruption_fails(self):
        fixture = self.fixture()
        seal = fixture.seals["zeller"]
        header, rows = read_tsv(seal / "sample_flow.tsv")
        rows[0]["observed_baseline_profiles"] = "0"
        rows[0]["observed_profiles"] = str(int(rows[0]["observed_profiles"]) - 1)
        rows[1]["observed_baseline_profiles"] = "2"
        rows[1]["observed_profiles"] = str(int(rows[1]["observed_profiles"]) + 1)
        write_tsv(seal / "sample_flow.tsv", header, rows)
        fixture.reseal("zeller")
        self.build_fails(fixture, "baseline profiles but expected")

    def test_component_sum_mismatch_fails(self):
        fixture = self.fixture()
        seal = fixture.seals["feng"]
        header, rows = read_tsv(seal / "sample_flow.tsv")
        rows[0]["observed_profiles"] = str(int(rows[0]["observed_profiles"]) + 1)
        write_tsv(seal / "sample_flow.tsv", header, rows)
        fixture.reseal("feng")
        self.build_fails(fixture, "but its components sum to")

    def test_independent_manifest_must_equal_the_flagged_rows(self):
        fixture = self.fixture()
        seal = fixture.seals["feng"]
        header, rows = read_tsv(seal / "production_manifest.tsv")
        flagged = [row for row in rows if row["independent_subset"] == "1"]
        unflagged = [row for row in rows if row["independent_subset"] == "0"]
        # Swap one membership flag in both the manifest and the flow ledger.
        demoted, promoted = flagged[0]["sample_id"], unflagged[0]["sample_id"]
        flagged[0]["independent_subset"] = "0"
        unflagged[0]["independent_subset"] = "1"
        write_tsv(seal / "production_manifest.tsv", header, rows)
        # Keep the ledger fully self-consistent with the swapped flags, so the
        # per-sample topology gate passes and only the independent manifest
        # still disagrees.
        flow_header, flow = read_tsv(seal / "sample_flow.tsv")
        by_id = {row["sample_id"]: row for row in flow}
        for sample, subset in ((demoted, "0"), (promoted, "1")):
            row = by_id[sample]
            count = 60 if subset == "1" else 0
            row["independent_subset"] = subset
            row["expected_independent_profiles"] = str(count)
            row["observed_independent_profiles"] = str(count)
            row["expected_profiles"] = str(1 + count + 7)
            row["observed_profiles"] = str(1 + count + 7)
        write_tsv(seal / "sample_flow.tsv", flow_header, flow)
        fixture.reseal("feng")
        self.build_fails(fixture, "does not match the rows flagged")

    def test_wrong_input_provenance_mode_fails(self):
        fixture = self.fixture()
        seal = fixture.seals["yachida"]
        header, rows = read_tsv(seal / "sample_flow.tsv")
        for row in rows:
            row["input_provenance_mode"] = "state_file"
        write_tsv(seal / "sample_flow.tsv", header, rows)
        fixture.reseal("yachida")
        self.build_fails(fixture, "reports input-provenance mode 'state_file'")

    def test_inventory_schema_and_values_are_validated(self):
        for field, value, needle in (
                ("sha256", "nothex", "invalid sha256"),
                ("bytes", "0", "non-positive byte count"),
                ("status", "UNVERIFIED", "not VERIFIED")):
            fixture = self.fixture()
            seal = fixture.seals["feng"]
            header, rows = read_tsv(seal / "source_seal_inventory.tsv")
            rows[0][field] = value
            write_tsv(seal / "source_seal_inventory.tsv", header, rows)
            fixture.reseal("feng")
            self.build_fails(fixture, needle, name="inv_%s" % field)

        drifted = self.fixture()
        seal = drifted.seals["zeller"]
        header, rows = read_tsv(seal / "source_seal_inventory.tsv")
        write_tsv(seal / "source_seal_inventory.tsv", header[:-1], rows)
        drifted.reseal("zeller")
        self.build_fails(drifted, "source_seal_inventory.tsv schema must be",
                         name="inv_schema")

    def test_inventory_member_and_native_audit_job_are_valid(self):
        duplicated = self.fixture()
        seal = duplicated.seals["feng"]
        header, rows = read_tsv(seal / "source_seal_inventory.tsv")
        rows.append(dict(rows[0]))
        write_tsv(seal / "source_seal_inventory.tsv", header, rows)
        duplicated.reseal("feng")
        self.build_fails(duplicated, "repeats member", name="inv_duplicate")

        invalid = self.fixture()
        seal = invalid.seals["feng"]
        header, rows = read_tsv(seal / "source_seal_inventory.tsv")
        rows[0]["source_audit_job"] = "job-3097680"
        write_tsv(seal / "source_seal_inventory.tsv", header, rows)
        invalid.reseal("feng")
        self.build_fails(invalid, "invalid source_audit_job", name="inv_job")

        legacy = self.fixture()
        seal = legacy.seals["yachida"]
        header, rows = read_tsv(seal / "source_seal_inventory.tsv")
        rows[0]["source_audit_job"] = ""
        write_tsv(seal / "source_seal_inventory.tsv", header, rows)
        legacy.reseal("yachida")
        header, rows = read_tsv(legacy.ledger)
        for row in rows:
            if row["cohort"] == "yachida":
                row["seal_sha256"] = sha256_file(
                    seal / "production_seal.sha256")
        write_tsv(legacy.ledger, header, rows)
        self.build_ok(legacy, name="legacy_blank_native_job")

        mixed = self.fixture()
        seal = mixed.seals["feng"]
        header, rows = read_tsv(seal / "source_seal_inventory.tsv")
        second = dict(rows[0])
        second["member"] = "another_source_member.tsv"
        second["source_audit_job"] = "999999"
        rows.append(second)
        write_tsv(seal / "source_seal_inventory.tsv", header, rows)
        mixed.reseal("feng")
        self.build_fails(mixed, "reports multiple audit jobs",
                         name="mixed_native_jobs")


class AuditLedgerSemanticsTest(PackageTestCase):
    """Item 4: the ledger must describe these audits of these seals."""

    def corrupt(self, fixture, field, value, cohort="feng"):
        header, rows = read_tsv(fixture.ledger)
        for row in rows:
            if row["cohort"] == cohort:
                row[field] = value
        write_tsv(fixture.ledger, header, rows)
        return fixture

    def test_wrong_seal_checksum_fails(self):
        fixture = self.corrupt(self.fixture(), "seal_sha256", "e" * 64)
        self.build_fails(fixture, "seal_sha256 for feng does not match")

    def test_wrong_job_identifier_fails(self):
        fixture = self.corrupt(self.fixture(), "audit_job_id", "9999999")
        self.build_fails(fixture, "records job 9999999 for feng")

    def test_non_numeric_job_identifier_fails(self):
        fixture = self.corrupt(self.fixture(), "audit_job_id", "job-1")
        self.build_fails(fixture, "non-numeric audit_job_id")

    def test_expected_job_can_be_overridden_explicitly(self):
        fixture = self.corrupt(self.fixture(), "audit_job_id", "4200001")
        self.build_ok(fixture, extra=["--feng-audit-job", "4200001"])

    def test_non_completed_state_fails(self):
        self.build_fails(self.corrupt(self.fixture(), "state", "FAILED"),
                         "has state 'FAILED', not COMPLETED")

    def test_nonzero_exit_code_fails(self):
        self.build_fails(self.corrupt(self.fixture(), "exit_code", "1:0"),
                         "has exit_code '1:0', not 0:0")

    def test_wrong_sample_count_fails(self):
        self.build_fails(self.corrupt(self.fixture(), "samples", "153"),
                         "records 153 samples for feng")

    def test_invalid_utc_stamp_fails(self):
        for stamp in ("27/09/2026", "2026-09", "yesterday"):
            fixture = self.corrupt(self.fixture(), "audit_date_utc", stamp)
            self.build_fails(fixture, "invalid UTC audit_date_utc",
                             name="stamp_%s" % abs(hash(stamp)))

    def test_unknown_or_duplicate_cohort_fails(self):
        unknown = self.fixture()
        header, rows = read_tsv(unknown.ledger)
        rows.append(dict(rows[0], cohort="gutbiome"))
        write_tsv(unknown.ledger, header, rows)
        self.build_fails(unknown, "names the unknown cohort 'gutbiome'")

        duplicate = self.fixture()
        header, rows = read_tsv(duplicate.ledger)
        rows.append(dict(rows[0]))
        write_tsv(duplicate.ledger, header, rows)
        self.build_fails(duplicate, "more than one row for", name="dupe")

    def test_registry_job_comes_from_the_validated_ledger(self):
        out = self.build_ok(self.fixture())
        registry = {row["cohort"]: row for row in
                    read_tsv(out / "zenodo/cohort_registry.tsv")[1]}
        for cohort, job in AUDIT_JOBS.items():
            self.assertEqual(registry[cohort]["source_audit_job"], job)


class ProvenanceMetadataTest(PackageTestCase):
    """Item 5: asset identities must be real, unique and complete."""

    def edit(self, fixture, predicate, field, value):
        header, rows = read_tsv(fixture.provenance)
        for row in rows:
            if predicate(row):
                row[field] = value
        write_tsv(fixture.provenance, header, rows)
        return fixture

    def test_short_or_blank_digest_fails(self):
        for value in ("", "abc", "z" * 64):
            fixture = self.edit(self.fixture(),
                                lambda row: row["category"] == "database_kraken2",
                                "sha256", value)
            self.build_fails(fixture, "needs a 64-character sha256",
                             name="digest_%s" % (value[:4] or "blank"))

    def test_parameter_rows_need_the_documented_sentinel(self):
        fixture = self.edit(self.fixture(),
                            lambda row: row["category"] == "upstream_parameter",
                            "sha256", "f" * 64)
        self.build_fails(fixture, "must record sha256 'not_a_file'")

    def test_duplicate_asset_id_fails(self):
        fixture = self.fixture()
        header, rows = read_tsv(fixture.provenance)
        rows[1]["asset_id"] = rows[0]["asset_id"]
        write_tsv(fixture.provenance, header, rows)
        self.build_fails(fixture, "reuses asset_id")

    def test_unknown_category_fails(self):
        fixture = self.fixture()
        header, rows = read_tsv(fixture.provenance)
        rows.append(dict(rows[0], category="mystery", asset_id="mystery_id"))
        write_tsv(fixture.provenance, header, rows)
        self.build_fails(fixture, "undocumented category 'mystery'")

    def test_source_code_must_name_the_commit(self):
        fixture = self.edit(self.fixture(),
                            lambda row: row["category"] == "source_code",
                            "name", "repository at some other revision")
        self.build_fails(fixture, "no source_code row identifies the supplied")

    def test_every_spike_target_needs_a_reference(self):
        fixture = self.fixture()
        header, rows = read_tsv(fixture.provenance)
        write_tsv(fixture.provenance, header,
                  [row for row in rows if row["asset_id"] != "Hhat"])
        self.build_fails(fixture, "exactly one spike_reference for: Hhat")

    def test_every_frozen_parameter_is_required(self):
        for parameter in REQUIRED_PARAMETERS:
            fixture = self.fixture()
            header, rows = read_tsv(fixture.provenance)
            write_tsv(fixture.provenance, header,
                      [row for row in rows if row["asset_id"] != parameter])
            self.build_fails(fixture, "exactly one frozen upstream_parameter "
                                      "row for: %s" % parameter,
                             name="param_%s" % parameter)

    def test_spike_references_and_parameters_are_one_to_one(self):
        duplicated = self.fixture()
        header, rows = read_tsv(duplicated.provenance)
        extra = dict(next(row for row in rows
                          if row["category"] == "spike_reference"))
        extra["asset_id"] = "extra_spike_reference"
        extra["name"] = "Bfrag"
        rows.append(extra)
        write_tsv(duplicated.provenance, header, rows)
        self.build_fails(duplicated, "exactly one spike_reference for: Bfrag",
                         name="duplicate_spike_reference")

        duplicated = self.fixture()
        header, rows = read_tsv(duplicated.provenance)
        extra = dict(next(row for row in rows
                          if row["category"] == "upstream_parameter"))
        extra["asset_id"] = "extra_parameter"
        extra["name"] = "read_length"
        rows.append(extra)
        write_tsv(duplicated.provenance, header, rows)
        self.build_fails(duplicated,
                         "exactly one frozen upstream_parameter row for: "
                         "read_length", name="duplicate_parameter")

    def test_analysis_image_checksum_is_required(self):
        fixture = self.fixture()
        command = fixture.command(self.outdir("missing_image_sha"))
        index = command.index("--analysis-image-sha256")
        del command[index:index + 2]
        done = subprocess.run(command, text=True, capture_output=True)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("--analysis-image-sha256 is required",
                      done.stdout + done.stderr)

    def test_analysis_image_checksum_must_agree(self):
        fixture = self.fixture()
        out = self.outdir("image")
        command = fixture.command(out)
        command[command.index("--analysis-image-sha256") + 1] = "a" * 64
        done = subprocess.run(command, text=True, capture_output=True)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("container_analysis sha256 does not match",
                      done.stdout + done.stderr)

    def test_analysis_image_checksum_is_recorded(self):
        out = self.build_ok(self.fixture())
        parameters = dict(
            (row["parameter"], row["value"]) for row in
            read_tsv(out / "provenance/build_parameters.tsv")[1])
        self.assertEqual(parameters["analysis_image_sha256"],
                         ANALYSIS_IMAGE_SHA256)


class PanelAndAliasTest(PackageTestCase):
    """Item 6: parse by the real delimiter and validate both inputs."""

    def test_alias_output_is_genuinely_tab_delimited(self):
        fixture = self.fixture()
        # The input really is comma-delimited.
        raw = fixture.taxon_aliases.read_text(encoding="utf-8")
        self.assertIn(",", raw.splitlines()[0])
        self.assertNotIn("\t", raw)
        out = self.build_ok(fixture)
        written = (out / "zenodo/taxon_aliases.tsv").read_text(encoding="utf-8")
        header = written.splitlines()[0]
        self.assertEqual(header, "canonical\talias\ttool")
        self.assertNotIn(",", written)
        for line in written.splitlines()[1:]:
            self.assertEqual(len(line.split("\t")), 3, line)
        _, rows = read_tsv(out / "zenodo/taxon_aliases.tsv")
        self.assertEqual(len(rows), 2 * len(SPIKE_TARGETS))

    def test_alias_schema_is_validated(self):
        fixture = self.fixture()
        fixture.taxon_aliases.write_text(
            "canonical,alias\nSpecies Fnuc,Species Fnuc\n", encoding="utf-8")
        self.build_fails(fixture, "taxon alias table lacks column(s): tool")

        blank = self.fixture()
        blank.taxon_aliases.write_text(
            "canonical,alias,tool\nSpecies Fnuc,,metaphlan4\n",
            encoding="utf-8")
        self.build_fails(blank, "taxon alias line 2 has a blank alias",
                         name="blank_alias")

    def test_spike_panel_schema_and_counts_are_validated(self):
        missing_weight = self.fixture()
        header, rows = read_tsv(missing_weight.spike_panel)
        header.remove("weight")
        write_tsv(missing_weight.spike_panel, header, rows)
        self.build_fails(missing_weight, "lacks column(s): weight",
                         name="missing_weight")

        short = self.fixture()
        header, rows = read_tsv(short.spike_panel)
        write_tsv(short.spike_panel, header, rows[:-1])
        self.build_fails(short, "spike panel has 9 targets")

        duplicated = self.fixture()
        header, rows = read_tsv(duplicated.spike_panel)
        rows[1]["label"] = rows[0]["label"]
        write_tsv(duplicated.spike_panel, header, rows)
        self.build_fails(duplicated, "spike panel repeats label",
                         name="dupe_label")

        blank = self.fixture()
        header, rows = read_tsv(blank.spike_panel)
        rows[0]["assembly"] = ""
        write_tsv(blank.spike_panel, header, rows)
        self.build_fails(blank, "has a blank assembly", name="blank_assembly")

        weighted = self.fixture()
        header, rows = read_tsv(weighted.spike_panel)
        rows[0]["weight"] = "0"
        write_tsv(weighted.spike_panel, header, rows)
        self.build_fails(weighted, "non-positive weight", name="weight")

        nonnumeric = self.fixture()
        header, rows = read_tsv(nonnumeric.spike_panel)
        rows[0]["weight"] = "heavy"
        write_tsv(nonnumeric.spike_panel, header, rows)
        self.build_fails(nonnumeric, "nonnumeric weight", name="weight2")


class MidBuildMutationTest(PackageTestCase):
    """Item 8: the builder itself must detect a source seal changing."""

    def mutating_build(self, mutate):
        """Run the builder while a background thread mutates a source seal.

        The mutation is applied when the builder's own staging directory
        appears, which is after the seals have been read and long before
        SUCCESS, so the revalidation is what must catch it.
        """
        import threading

        fixture = self.fixture()
        out = self.outdir("midbuild")
        staging = out.with_name(out.name + ".incomplete")
        stop = threading.Event()
        applied = []

        def watcher():
            while not stop.is_set():
                if staging.exists():
                    mutate(fixture)
                    applied.append(True)
                    return
                stop.wait(0.01)

        thread = threading.Thread(target=watcher)
        thread.start()
        try:
            done = fixture.run(out)
        finally:
            stop.set()
            thread.join(timeout=30)
        return fixture, out, done, applied

    def test_a_file_added_during_the_build_prevents_success(self):
        def mutate(fixture):
            (fixture.seals["feng"] / "sneaked_in.tsv").write_text(
                "added\n", encoding="utf-8")

        fixture, out, done, applied = self.mutating_build(mutate)
        self.assertTrue(applied, "the mutation never ran")
        self.assertNotEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertIn("changed during the build", done.stdout + done.stderr)
        self.assertIn("added sneaked_in.tsv", done.stdout + done.stderr)
        self.assertFalse((out / "SUCCESS").exists())
        self.assertFalse(out.exists())

    def test_a_file_changed_during_the_build_prevents_success(self):
        def mutate(fixture):
            target = fixture.seals["zeller"] / "covariate_audit.tsv"
            target.write_text(target.read_text(encoding="utf-8") + "x\t0\t0\t0\n",
                              encoding="utf-8")

        fixture, out, done, applied = self.mutating_build(mutate)
        self.assertTrue(applied)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("changed during the build", done.stdout + done.stderr)
        self.assertFalse(out.exists())

    def test_audit_in_progress_during_the_build_prevents_success(self):
        def mutate(fixture):
            (fixture.seals["yachida"] / "AUDIT_IN_PROGRESS").write_text(
                "status\tIN_PROGRESS\n", encoding="utf-8")

        fixture, out, done, applied = self.mutating_build(mutate)
        self.assertTrue(applied)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("mid-audit during the build", done.stdout + done.stderr)
        self.assertFalse(out.exists())


class TemplateTest(unittest.TestCase):
    """Item 9: documented templates, with placeholders rather than hashes."""

    ROOT = ROOT / "analysis_v2/templates/upstream_evidence"

    def test_templates_exist_outside_the_untracked_inputs_directory(self):
        self.assertTrue(self.ROOT.is_dir())
        for name in ("upstream_provenance_metadata.tsv", "audit_ledger.tsv",
                     "README.md"):
            self.assertTrue((self.ROOT / name).is_file(), name)
        self.assertNotIn("inputs", self.ROOT.parts)

    def test_template_schemas_match_the_builder(self):
        header, rows = read_tsv(self.ROOT / "upstream_provenance_metadata.tsv")
        self.assertEqual(header, PROVENANCE_FIELDS)
        categories = {row["category"] for row in rows}
        for required in DIGEST_CATEGORIES + ("spike_reference",
                                             "upstream_parameter"):
            self.assertIn(required, categories)
        parameters = {row["asset_id"] for row in rows
                      if row["category"] == "upstream_parameter"}
        self.assertEqual(parameters, set(REQUIRED_PARAMETERS))
        spikes = {row["asset_id"] for row in rows
                  if row["category"] == "spike_reference"}
        self.assertEqual(spikes, set(SPIKE_TARGETS))

        header, rows = read_tsv(self.ROOT / "audit_ledger.tsv")
        self.assertEqual(header, LEDGER_FIELDS)
        self.assertEqual({row["cohort"] for row in rows}, set(COHORTS))
        for row in rows:
            self.assertEqual(row["audit_job_id"], AUDIT_JOBS[row["cohort"]])
            self.assertEqual(row["state"], "COMPLETED")
            self.assertEqual(row["exit_code"], "0:0")
            self.assertEqual(int(row["samples"]),
                             COHORTS[row["cohort"]]["samples"])

    def test_templates_contain_placeholders_not_invented_hashes(self):
        import re

        for name in ("upstream_provenance_metadata.tsv", "audit_ledger.tsv"):
            text = (self.ROOT / name).read_text(encoding="utf-8")
            self.assertIn("<", text, "%s has no placeholders" % name)
            for token in re.findall(r"\b[0-9a-f]{64}\b", text):
                self.fail("%s contains an invented digest: %s" % (name, token))
            # The parameter sentinel is a documented value, not a hash.
            if name.startswith("upstream"):
                self.assertIn(PARAMETER_SENTINEL, text)


if __name__ == "__main__":
    unittest.main(verbosity=2)
