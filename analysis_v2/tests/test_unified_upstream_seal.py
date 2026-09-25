#!/usr/bin/env python3
"""Fixture suite for the unified upstream cohort seal.

Every case drives the real auditor as a subprocess against a synthetic cohort
tree. Nothing here reimplements the auditor's logic, and no real cluster output
is touched.

The 18 cases required by UNIFIED_UPSTREAM_SEAL_SPEC.md section 7 map onto the
`test_caseNN_*` methods below.
"""
from __future__ import annotations

import hashlib
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
AUDITOR = ROOT / "analysis_v2/scripts/seal_cohort_upstream.py"

SEAL_MEMBERS = [
    "sample_flow.tsv", "covariate_audit.tsv", "production_manifest.tsv",
    "production_manifest.independent.tsv", "source_seal_inventory.tsv",
    "SUCCESS", "production_seal.sha256",
]
CONDITIONS = ("Control", "Adenoma", "CRC")
SELECTION_TRIPLET = ("selection_rank", "selection_hash", "selection_seed")
COMMUNITY_PER_SAMPLE = 7
INDEPENDENT_PER_SUBSET_SAMPLE = 60

# Boundary settings per cohort, mirroring the production manifests.
COHORTS = {
    "yachida": {
        "study": "YachidaS_2019",
        "condition_column": "Target_Condition",
        "study_column": None,
        "batch_column": "batch_id",
        # The sealed 29-column production manifest carries all four.
        "covariates": ("age", "sex", "bmi", "batch_id"),
        # Batch provenance is production-only; the sealed 26-column
        # independent manifest has none of it.
        "independent_drop": ("batch_id",),
        # And its selection triplet appears twice: pilot, then independent.
        "historical_duplicate_selection": True,
        "completion_has_study": False,
        "seal_manifest_member": "pilot_batched.tsv",
        "seal_independent_member": "independent_10_per_condition.tsv",
        "seal_extra": ("dataset_completion.tsv", "dataset_completion.tsv.sha256"),
        "success_schema": "yachida_dataset",
    },
    "feng": {
        "study": "FengQ_2015",
        "condition_column": "condition",
        "study_column": "study",
        "batch_column": None,
        "covariates": ("age", "sex", "bmi"),
        "independent_drop": (),
        "historical_duplicate_selection": False,
        "completion_has_study": True,
        "seal_manifest_member": "production_manifest.tsv",
        "seal_independent_member": "production_manifest.independent.tsv",
        "seal_extra": ("sample_flow.tsv", "covariate_audit.tsv"),
        "success_schema": "cohort_profiles",
    },
    "zeller": {
        "study": "ZellerG_2014",
        "condition_column": "condition",
        "study_column": "study",
        "batch_column": None,
        "covariates": ("age", "sex", "bmi"),
        "independent_drop": (),
        "historical_duplicate_selection": False,
        "completion_has_study": True,
        "seal_manifest_member": "production_manifest.tsv",
        "seal_independent_member": "production_manifest.independent.tsv",
        "seal_extra": ("sample_flow.tsv", "covariate_audit.tsv"),
        "success_schema": "cohort_profiles",
    },
}


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(fields)]
    for row in rows:
        lines.append("\t".join(str(row.get(name, "")) for name in fields))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_tsv(path):
    lines = path.read_text(encoding="utf-8").rstrip("\n").split("\n")
    header = lines[0].split("\t")
    return header, [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def tree_digest(root):
    """Every file under `root`, keyed by relative path, for byte-identity proof."""
    return {
        str(path.relative_to(root)): sha256_file(path)
        for path in sorted(root.rglob("*")) if path.is_file()
    }


class Fixture:
    """A synthetic but structurally complete upstream cohort."""

    def __init__(self, root, cohort="feng", samples=6, independent=3):
        self.root = Path(root)
        self.cohort = cohort
        self.config = COHORTS[cohort]
        self.study = self.config["study"]
        self.samples = ["S%02d" % index for index in range(1, samples + 1)]
        self.independent = self.samples[:independent]
        self.condition = {}
        for index, sample in enumerate(self.samples):
            self.condition[sample] = CONDITIONS[index % len(CONDITIONS)]
        self.state = self.root / "state"
        self.results = self.root / "results"
        self.qc = self.root / "qc"
        self.scratch = self.root / "scratch"
        self.seal = self.root / "native_seal"
        self.outdir = self.root / "production_seal_v2"
        self.manifest = self.root / "manifests" / "production.tsv"
        self.independent_manifest = self.root / "manifests" / "independent.tsv"
        for directory in (self.state / "samples", self.results, self.qc,
                          self.scratch, self.seal):
            directory.mkdir(parents=True, exist_ok=True)
        self.build()

    # ------------------------------------------------------------ manifests
    def manifest_fields(self):
        fields = ["sample_id", self.config["condition_column"]]
        if self.config["study_column"]:
            fields.append(self.config["study_column"])
        if self.config["batch_column"]:
            fields.append(self.config["batch_column"])
        for name in self.config["covariates"]:
            if name not in fields:
                fields.append(name)
        if self.cohort != "yachida":
            fields.append("independent_subset")
        else:
            # Pilot selection provenance, as the real manifest carries it.
            fields.extend(SELECTION_TRIPLET)
        return fields

    def manifest_row(self, sample):
        row = {"sample_id": sample,
               self.config["condition_column"]: self.condition[sample]}
        if self.config["study_column"]:
            row[self.config["study_column"]] = self.study
        if self.config["batch_column"]:
            row[self.config["batch_column"]] = "batch_%s" % (
                "a" if sample in self.independent else "b")
        number = self.samples.index(sample)
        for name in self.config["covariates"]:
            if name == "batch_id":
                continue
            # One sample has no recorded covariates, so missingness is nonzero.
            if number == 0:
                row[name] = ""
            elif name == "age":
                row[name] = str(50 + number)
            elif name == "sex":
                row[name] = "Female" if number % 2 else "Male"
            else:
                row[name] = "%.1f" % (20.0 + number)
        if self.cohort != "yachida":
            row["independent_subset"] = "1" if sample in self.independent else "0"
        else:
            for name in SELECTION_TRIPLET:
                row[name] = "pilot_%s_%d" % (name, number)
        return row

    def independent_selection_value(self, name, position):
        return "independent_%s_%d" % (name, position)

    def write_manifests(self, independent_extra=None, independent_drop=(),
                        independent_overrides=None, duplicate_selection=None,
                        duplicate_fields=SELECTION_TRIPLET, occurrences=2,
                        extra_shared=None):
        """Write both manifests, reproducing the real Yachida asymmetry.

        The independent manifest is written positionally so the documented
        historical duplicated header can be reproduced exactly.
        """
        fields = self.manifest_fields()
        rows = [self.manifest_row(sample) for sample in self.samples]
        if extra_shared:
            for name, values in extra_shared.items():
                fields = fields + [name]
                for position, row in enumerate(rows):
                    row[name] = values[position]
        write_tsv(self.manifest, fields, rows)

        drop = set(independent_drop) | set(self.config["independent_drop"])
        subset = [row for row in rows if row["sample_id"] in self.independent]
        columns = [(name, [row[name] for row in subset])
                   for name in fields if name not in drop]
        if duplicate_selection is None:
            duplicate_selection = self.config["historical_duplicate_selection"]
        if duplicate_selection:
            # Reproduce the real historical selector: it overwrote the
            # inherited values in the row dict before writing, and its
            # fieldnames listed the triplet twice, so csv.DictWriter emitted
            # the SAME independent value into both occurrences. The pilot
            # values survive only in the production manifest.
            for _ in range(occurrences - 1):
                for name in duplicate_fields:
                    columns.append((name, [
                        self.independent_selection_value(name, position)
                        for position in range(len(subset))]))
            for name in duplicate_fields:
                if name not in SELECTION_TRIPLET:
                    continue
                for index, (column, values) in enumerate(columns):
                    if column == name:
                        columns[index] = (column, [
                            self.independent_selection_value(name, position)
                            for position in range(len(subset))])
        if independent_extra:
            for name in independent_extra:
                columns.append((name, ["%s_%d" % (name, position)
                                       for position in range(len(subset))]))
        if independent_overrides:
            for name, value in independent_overrides.items():
                for index, (column, values) in enumerate(columns):
                    if column == name:
                        values[0] = value
                        break
                else:
                    raise AssertionError("no such column: %s" % name)
        header = [name for name, _ in columns]
        lines = ["\t".join(header)]
        for position in range(len(subset)):
            lines.append("\t".join(values[position] for _, values in columns))
        self.independent_manifest.parent.mkdir(parents=True, exist_ok=True)
        self.independent_manifest.write_text(
            "\n".join(lines) + "\n", encoding="utf-8")

    # -------------------------------------------------------------- cohort
    def sample_root(self, sample):
        return self.results / self.study / sample

    def qc_root(self, sample):
        return self.qc / self.study / sample

    def receipt_path(self, sample):
        return self.state / "samples" / ("%s.retained_outputs.tsv" % sample)

    def write_receipt(self, sample, payloads=None, mutate=None):
        if payloads is None:
            payloads = [self.sample_root(sample) / "sample_completion.tsv",
                        self.qc_root(sample) / "qc_report.txt"]
        rows = []
        for payload in payloads:
            rows.append({"path": str(payload), "sha256": sha256_file(payload),
                         "bytes": payload.stat().st_size})
        if mutate is not None:
            rows = mutate(rows)
        write_tsv(self.receipt_path(sample), ["path", "sha256", "bytes"], rows)

    def write_completion(self, sample, overrides=None):
        in_independent = sample in self.independent
        expected_independent = (
            INDEPENDENT_PER_SUBSET_SAMPLE if in_independent else 0)
        total = 1 + expected_independent + COMMUNITY_PER_SAMPLE
        values = [("sample_id", sample)]
        if self.config["completion_has_study"]:
            values.append(("study", self.study))
        values.extend([
            ("condition", self.condition[sample]),
            ("independent_subset", "1" if in_independent else "0"),
            ("expected_profiles", str(total)),
            ("observed_profiles", str(total)),
            ("community_design_rows", str(COMMUNITY_PER_SAMPLE)),
            ("independent_design_rows", str(expected_independent)),
        ])
        if overrides:
            values = [(name, overrides.get(name, value)) for name, value in values]
        text = "field\tvalue\n" + "".join(
            "%s\t%s\n" % (name, value) for name, value in values)
        path = self.sample_root(sample) / "sample_completion.tsv"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")

    def build(self):
        self.write_manifests()
        for sample in self.samples:
            sample_root = self.sample_root(sample)
            sample_root.mkdir(parents=True, exist_ok=True)
            (sample_root / "SUCCESS").touch()
            baseline = sample_root / "profiles/baseline" / sample / "SUCCESS"
            baseline.parent.mkdir(parents=True, exist_ok=True)
            baseline.touch()
            for index in range(COMMUNITY_PER_SAMPLE):
                marker = sample_root / "profiles/community" / str(index) / "SUCCESS"
                marker.parent.mkdir(parents=True, exist_ok=True)
                marker.touch()
            if sample in self.independent:
                for index in range(INDEPENDENT_PER_SUBSET_SAMPLE):
                    marker = (sample_root / "profiles/independent"
                              / str(index) / "SUCCESS")
                    marker.parent.mkdir(parents=True, exist_ok=True)
                    marker.touch()
            qc_root = self.qc_root(sample)
            qc_root.mkdir(parents=True, exist_ok=True)
            (qc_root / "qc_report.txt").write_text(
                "sample\t%s\nstatus\tPASS\n" % sample, encoding="utf-8")
            (self.scratch / sample).mkdir(parents=True, exist_ok=True)
            for suffix in ("verified", "input_provenance.tsv"):
                (self.state / "samples" / ("%s.%s" % (sample, suffix))).write_text(
                    "status\tPASS\n", encoding="utf-8")
            self.write_completion(sample)
            self.write_receipt(sample)
        self.write_source_seal()

    # ---------------------------------------------------------- native seal
    def native_success(self, **overrides):
        """The documented native SUCCESS for this cohort's schema."""
        samples = len(self.samples)
        independent = len(self.independent)
        values = {
            "samples": samples,
            "independent": independent,
            "baseline": samples,
            "independent_profiles": independent * INDEPENDENT_PER_SUBSET_SAMPLE,
            "community": samples * COMMUNITY_PER_SAMPLE,
            "identity": (self.study if self.config["success_schema"]
                         == "yachida_dataset" else self.cohort),
            "status": "PASS",
        }
        values.update(overrides)
        counts = {}
        for sample in self.samples:
            counts[self.condition[sample]] = counts.get(
                self.condition[sample], 0) + 1
        if self.config["success_schema"] == "yachida_dataset":
            return (
                "dataset\t%(identity)s\n"
                "samples\t%(samples)s\n"
                "batches\t2\n"
                "independent_subset\t%(independent)s\n"
                "profiles\tbaseline=%(baseline)s;independent="
                "%(independent_profiles)s;community=%(community)s\n"
                "conditions\t" + ";".join(
                    "%s=%d" % item for item in sorted(counts.items()))
                + "\nstatus\t%(status)s\n") % values
        return (
            "cohort\t%(identity)s\n"
            "samples\t%(samples)s\n"
            "independent_samples\t%(independent)s\n"
            "baseline_profiles\t%(baseline)s\n"
            "independent_profiles\t%(independent_profiles)s\n"
            "community_profiles\t%(community)s\n"
            "status\t%(status)s\n") % values

    def write_source_seal(self, corrupt_member=None, drop_member=None,
                          in_progress=False, success=True,
                          success_overrides=None):
        for item in list(self.seal.iterdir()):
            if item.is_file():
                item.unlink()
        members = []
        for name in self.config["seal_extra"]:
            (self.seal / name).write_text(
                "field\tvalue\nnative\t%s\n" % name, encoding="utf-8")
            members.append(name)
        shutil.copyfile(str(self.manifest),
                        str(self.seal / self.config["seal_manifest_member"]))
        shutil.copyfile(str(self.independent_manifest),
                        str(self.seal / self.config["seal_independent_member"]))
        members.append(self.config["seal_manifest_member"])
        members.append(self.config["seal_independent_member"])
        if success:
            (self.seal / "SUCCESS").write_text(
                self.native_success(**(success_overrides or {})),
                encoding="utf-8")
            members.append("SUCCESS")
        members = [name for name in members if name != drop_member]
        lines = []
        for name in members:
            lines.append("%s  %s\n" % (sha256_file(self.seal / name), name))
        (self.seal / "production_seal.sha256").write_text(
            "".join(lines), encoding="utf-8")
        if corrupt_member:
            path = self.seal / corrupt_member
            path.write_text(path.read_text(encoding="utf-8") + "tampered\n",
                            encoding="utf-8")
        if in_progress:
            (self.seal / "AUDIT_IN_PROGRESS").write_text(
                "status\tIN_PROGRESS\n", encoding="utf-8")

    # ---------------------------------------------------------------- run
    def command(self, extra=()):
        command = [
            sys.executable, str(AUDITOR),
            "--cohort", self.cohort,
            "--manifest", str(self.manifest),
            "--independent-manifest", str(self.independent_manifest),
            "--source-seal", str(self.seal),
            "--state-dir", str(self.state),
            "--results-root", str(self.results),
            "--qc-root", str(self.qc),
            "--scratch-root", str(self.scratch),
            "--outdir", str(self.outdir),
            "--expected-samples", str(len(self.samples)),
            "--expected-independent", str(len(self.independent)),
            "--expected-conditions", self.expected_conditions(),
        ]
        command.extend(extra)
        return command

    def expected_conditions(self):
        counts = {}
        for sample in self.samples:
            counts[self.condition[sample]] = counts.get(
                self.condition[sample], 0) + 1
        return ",".join("%s=%d" % item for item in sorted(counts.items()))

    def run(self, extra=()):
        return subprocess.run(self.command(extra), text=True, capture_output=True)


class SealFixtureTestCase(unittest.TestCase):
    """Shared fixture plumbing; carries no cases of its own."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self._count = 0

    def tearDown(self):
        self._tmp.cleanup()

    def fixture(self, cohort="feng", samples=6, independent=3):
        self._count += 1
        directory = self.root / ("fix%02d" % self._count)
        directory.mkdir()
        return Fixture(directory, cohort=cohort, samples=samples,
                       independent=independent)

    def seal_ok(self, fixture, extra=()):
        done = fixture.run(extra)
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertTrue((fixture.outdir / "SUCCESS").is_file())
        self.assertTrue((fixture.outdir / "production_seal.sha256").is_file())
        self.assertFalse((fixture.outdir / "AUDIT_IN_PROGRESS").exists())
        return fixture.outdir

    def seal_fails(self, fixture, needle, extra=()):
        done = fixture.run(extra)
        self.assertNotEqual(done.returncode, 0,
                            "auditor accepted bad input (%s)" % needle)
        message = done.stdout + done.stderr
        self.assertIn(needle, message)
        self.assertFalse((fixture.outdir / "SUCCESS").exists())
        self.assertFalse((fixture.outdir / "production_seal.sha256").exists())
        return message


class UnifiedSealTestCase(SealFixtureTestCase):
    """The 18 fixture cases required by the frozen specification."""

    # ------------------------------------------------- cases 1, 2, 3
    def test_case01_yachida_configuration_seals(self):
        fixture = self.fixture(cohort="yachida")
        outdir = self.seal_ok(fixture)
        success = dict(
            line.split("\t") for line in
            (outdir / "SUCCESS").read_text().rstrip("\n").split("\n"))
        self.assertEqual(success["seal_contract"], "upstream_seal_v2")
        self.assertEqual(success["cohort"], "yachida")
        self.assertEqual(success["samples"], "6")
        self.assertEqual(success["independent_samples"], "3")
        self.assertEqual(success["baseline_profiles"], "6")
        self.assertEqual(success["independent_profiles"], "180")
        self.assertEqual(success["community_profiles"], "42")
        self.assertEqual(success["status"], "PASS")
        header, rows = read_tsv(outdir / "sample_flow.tsv")
        self.assertTrue(all(row["batch_id"].startswith("batch_") for row in rows))

    def test_case02_feng_configuration_seals(self):
        outdir = self.seal_ok(self.fixture(cohort="feng"))
        self.assertIn("cohort\tfeng", (outdir / "SUCCESS").read_text())

    def test_case03_zeller_configuration_seals(self):
        outdir = self.seal_ok(self.fixture(cohort="zeller"))
        self.assertIn("cohort\tzeller", (outdir / "SUCCESS").read_text())

    # ------------------------------------------------------------- case 4
    def test_case04_schemas_identical_across_modes(self):
        seals = {}
        for cohort in ("yachida", "feng", "zeller"):
            seals[cohort] = self.seal_ok(self.fixture(cohort=cohort))
        reference = None
        for cohort, outdir in sorted(seals.items()):
            members = sorted(item.name for item in outdir.iterdir())
            flow_header, _ = read_tsv(outdir / "sample_flow.tsv")
            covariate_header, _ = read_tsv(outdir / "covariate_audit.tsv")
            success_fields = [
                line.split("\t")[0] for line in
                (outdir / "SUCCESS").read_text().rstrip("\n").split("\n")]
            shape = (members, flow_header, covariate_header, success_fields)
            if reference is None:
                reference = shape
                self.assertEqual(members, sorted(SEAL_MEMBERS))
                self.assertIn("batch_id", flow_header)
            else:
                self.assertEqual(shape, reference, cohort)
        # Cohorts without upstream batches keep the column and leave it blank.
        _, feng_rows = read_tsv(seals["feng"] / "sample_flow.tsv")
        self.assertTrue(all(row["batch_id"] == "" for row in feng_rows))

    # ------------------------------------------------------------- case 5
    def test_case05_condition_column_normalization(self):
        yachida = self.seal_ok(self.fixture(cohort="yachida"))
        feng = self.seal_ok(self.fixture(cohort="feng"))
        for outdir in (yachida, feng):
            header, rows = read_tsv(outdir / "production_manifest.tsv")
            self.assertEqual(header[:4],
                             ["sample_id", "condition", "study", "independent_subset"])
            self.assertEqual([row["condition"] for row in rows][:3],
                             ["Control", "Adenoma", "CRC"])
        # Original columns survive under collision-safe names.
        yachida_header, _ = read_tsv(yachida / "production_manifest.tsv")
        self.assertIn("Target_Condition", yachida_header)
        feng_header, feng_rows = read_tsv(feng / "production_manifest.tsv")
        self.assertIn("source_condition", feng_header)
        self.assertIn("source_study", feng_header)
        self.assertIn("source_independent_subset", feng_header)
        self.assertEqual(feng_rows[0]["condition"], feng_rows[0]["source_condition"])

    # ------------------------------------------------------------- case 6
    def test_case06_study_value_versus_study_column(self):
        yachida = self.seal_ok(self.fixture(cohort="yachida"))
        feng = self.seal_ok(self.fixture(cohort="feng"))
        _, yachida_rows = read_tsv(yachida / "production_manifest.tsv")
        _, feng_rows = read_tsv(feng / "production_manifest.tsv")
        self.assertTrue(all(row["study"] == "YachidaS_2019" for row in yachida_rows))
        self.assertTrue(all(row["study"] == "FengQ_2015" for row in feng_rows))
        # Supplying both, or neither, is ambiguous and refused.
        fixture = self.fixture(cohort="feng")
        self.seal_fails(fixture, "ambiguous study configuration",
                        extra=["--study-column", "study",
                               "--study-value", "FengQ_2015"])

    # ------------------------------------------------------------- case 7
    def test_case07_condition_counts_and_subset_balance(self):
        fixture = self.fixture(cohort="feng", samples=30, independent=30)
        self.seal_ok(fixture)
        success = dict(
            line.split("\t") for line in
            (fixture.outdir / "SUCCESS").read_text().rstrip("\n").split("\n"))
        self.assertEqual(success["independent_samples"], "30")
        self.assertEqual(success["independent_profiles"], "1800")
        # Exact condition counts are enforced.
        self.seal_fails(fixture, "unexpected condition counts",
                        extra=["--expected-conditions",
                               "Control=11,Adenoma=9,CRC=10"])
        # The frozen 30-sample subset must be balanced 10/10/10.
        unbalanced = self.fixture(cohort="feng", samples=30, independent=30)
        unbalanced.condition[unbalanced.samples[0]] = "CRC"
        unbalanced.write_manifests()
        unbalanced.write_source_seal()
        self.seal_fails(unbalanced, "independent subset is not balanced")

    # ------------------------------------------------------------- case 8
    def test_case08_missing_baseline_fails(self):
        fixture = self.fixture()
        sample = fixture.samples[0]
        (fixture.sample_root(sample) / "profiles/baseline" / sample
         / "SUCCESS").unlink()
        self.seal_fails(fixture, "sample(s) are incomplete")
        _, rows = read_tsv(fixture.outdir / "sample_flow.tsv")
        failing = [row for row in rows if row["sample_id"] == sample][0]
        self.assertEqual(failing["observed_baseline_profiles"], "0")
        self.assertIn("baseline_profile_count", failing["failure_reasons"])

    # ------------------------------------------------------------- case 9
    def test_case09_community_baseline_offset_fails(self):
        """The total is unchanged, so only per-component checks can catch it."""
        fixture = self.fixture()
        sample = fixture.samples[0]
        root = fixture.sample_root(sample)
        (root / "profiles/community/0/SUCCESS").unlink()
        extra = root / "profiles/baseline/extra/SUCCESS"
        extra.parent.mkdir(parents=True)
        extra.touch()
        self.seal_fails(fixture, "sample(s) are incomplete")
        _, rows = read_tsv(fixture.outdir / "sample_flow.tsv")
        failing = [row for row in rows if row["sample_id"] == sample][0]
        self.assertEqual(failing["observed_profiles"], failing["expected_profiles"])
        self.assertEqual(failing["observed_baseline_profiles"], "2")
        self.assertEqual(failing["observed_community_profiles"], "6")
        self.assertEqual(failing["failure_reasons"],
                         "baseline_profile_count;community_profile_count")

    # ------------------------------------------------------------ case 10
    def test_case10_independent_profile_count_failures(self):
        missing = self.fixture()
        sample = missing.independent[0]
        (missing.sample_root(sample) / "profiles/independent/0/SUCCESS").unlink()
        self.seal_fails(missing, "sample(s) are incomplete")
        _, rows = read_tsv(missing.outdir / "sample_flow.tsv")
        row = [item for item in rows if item["sample_id"] == sample][0]
        self.assertEqual(row["observed_independent_profiles"], "59")
        self.assertIn("independent_profile_count", row["failure_reasons"])

        extra = self.fixture()
        sample = extra.independent[0]
        marker = extra.sample_root(sample) / "profiles/independent/extra/SUCCESS"
        marker.parent.mkdir(parents=True)
        marker.touch()
        self.seal_fails(extra, "sample(s) are incomplete")
        _, rows = read_tsv(extra.outdir / "sample_flow.tsv")
        row = [item for item in rows if item["sample_id"] == sample][0]
        self.assertEqual(row["observed_independent_profiles"], "61")

        # A non-subset sample must have exactly zero independent profiles.
        stray = self.fixture()
        sample = [s for s in stray.samples if s not in stray.independent][0]
        marker = stray.sample_root(sample) / "profiles/independent/0/SUCCESS"
        marker.parent.mkdir(parents=True)
        marker.touch()
        self.seal_fails(stray, "sample(s) are incomplete")

    # ------------------------------------------------------------ case 11
    def test_case11_completion_mismatch_fails(self):
        for field, value in (("condition", "WRONG"),
                             ("independent_subset", "0"),
                             ("community_design_rows", "6"),
                             ("observed_profiles", "1")):
            fixture = self.fixture()
            sample = fixture.independent[0]
            fixture.write_completion(sample, overrides={field: value})
            fixture.write_receipt(sample)
            self.seal_fails(fixture, "sample(s) are incomplete")
            _, rows = read_tsv(fixture.outdir / "sample_flow.tsv")
            row = [item for item in rows if item["sample_id"] == sample][0]
            self.assertIn("completion_table", row["failure_reasons"], field)
            self.assertIn(field, row["completion_error"])

        missing = self.fixture()
        sample = missing.independent[0]
        (missing.sample_root(sample) / "sample_completion.tsv").unlink()
        self.seal_fails(missing, "sample(s) are incomplete")
        _, rows = read_tsv(missing.outdir / "sample_flow.tsv")
        row = [item for item in rows if item["sample_id"] == sample][0]
        self.assertEqual(row["completion_error"], "missing_or_empty_completion_table")

        # Cohorts whose completion tables carry a study must have it checked.
        drifted = self.fixture(cohort="feng")
        sample = drifted.samples[0]
        drifted.write_completion(sample, overrides={"study": "WrongStudy"})
        drifted.write_receipt(sample)
        self.seal_fails(drifted, "sample(s) are incomplete")
        _, rows = read_tsv(drifted.outdir / "sample_flow.tsv")
        row = [item for item in rows if item["sample_id"] == sample][0]
        self.assertIn("study:WrongStudy!=FengQ_2015", row["completion_error"])

        absent = self.fixture(cohort="feng")
        sample = absent.samples[0]
        path = absent.sample_root(sample) / "sample_completion.tsv"
        path.write_text(
            "".join(line + "\n" for line
                    in path.read_text(encoding="utf-8").rstrip("\n").split("\n")
                    if not line.startswith("study\t")), encoding="utf-8")
        absent.write_receipt(sample)
        self.seal_fails(absent, "sample(s) are incomplete")
        _, rows = read_tsv(absent.outdir / "sample_flow.tsv")
        row = [item for item in rows if item["sample_id"] == sample][0]
        self.assertIn("study:<missing>", row["completion_error"])

        # Yachida completion tables have no study field, and that is correct.
        yachida = self.fixture(cohort="yachida")
        self.assertNotIn(
            "study\t",
            (yachida.sample_root(yachida.samples[0])
             / "sample_completion.tsv").read_text(encoding="utf-8"))
        self.seal_ok(yachida)

    # ------------------------------------------------------------ case 12
    def test_case12_missing_input_provenance_fails(self):
        fixture = self.fixture()
        sample = fixture.samples[0]
        (fixture.state / "samples" / ("%s.input_provenance.tsv" % sample)).unlink()
        self.seal_fails(fixture, "sample(s) are incomplete")
        _, rows = read_tsv(fixture.outdir / "sample_flow.tsv")
        row = [item for item in rows if item["sample_id"] == sample][0]
        self.assertIn("input_provenance", row["failure_reasons"])

        empty = self.fixture()
        sample = empty.samples[0]
        (empty.state / "samples" / ("%s.input_provenance.tsv" % sample)).write_text("")
        self.seal_fails(empty, "sample(s) are incomplete")

        unverified = self.fixture()
        sample = unverified.samples[0]
        (unverified.state / "samples" / ("%s.verified" % sample)).unlink()
        self.seal_fails(unverified, "sample(s) are incomplete")
        _, rows = read_tsv(unverified.outdir / "sample_flow.tsv")
        row = [item for item in rows if item["sample_id"] == sample][0]
        self.assertIn("verified_marker", row["failure_reasons"])

    # ------------------------------------------------------------ case 13
    def test_case13_receipt_failures(self):
        def receipt_error(mutate=None, payloads=None, delete=False):
            fixture = self.fixture()
            sample = fixture.samples[0]
            if delete:
                fixture.receipt_path(sample).unlink()
            else:
                fixture.write_receipt(sample, payloads=payloads, mutate=mutate)
            self.seal_fails(fixture, "sample(s) are incomplete")
            _, rows = read_tsv(fixture.outdir / "sample_flow.tsv")
            row = [item for item in rows if item["sample_id"] == sample][0]
            self.assertIn("retained_output_receipt", row["failure_reasons"])
            return fixture, row["receipt_error"]

        # A missing retained output.
        fixture = self.fixture()
        sample = fixture.samples[0]
        (fixture.qc_root(sample) / "qc_report.txt").unlink()
        self.seal_fails(fixture, "sample(s) are incomplete")
        _, rows = read_tsv(fixture.outdir / "sample_flow.tsv")
        self.assertIn("missing_or_empty_output",
                      [item for item in rows
                       if item["sample_id"] == sample][0]["receipt_error"])

        def wrong_size(rows):
            rows[0]["bytes"] = int(rows[0]["bytes"]) + 1
            return rows

        def wrong_digest(rows):
            rows[0]["sha256"] = "0" * 64
            return rows

        def duplicated(rows):
            return [rows[0], dict(rows[0])]

        _, error = receipt_error(mutate=wrong_size)
        self.assertIn("size_mismatch", error)
        _, error = receipt_error(mutate=wrong_digest)
        self.assertIn("sha256_mismatch", error)
        _, error = receipt_error(mutate=duplicated)
        self.assertIn("duplicate_receipt_path", error)
        _, error = receipt_error(delete=True)
        self.assertEqual(error, "missing_or_empty_receipt")

        # A retained output inside disposable scratch.
        scratch_fixture = self.fixture()
        sample = scratch_fixture.samples[0]
        disposable = scratch_fixture.scratch / sample / "tmp.txt"
        disposable.parent.mkdir(parents=True, exist_ok=True)
        disposable.write_text("scratch\n", encoding="utf-8")
        scratch_fixture.write_receipt(sample, payloads=[disposable])
        self.seal_fails(scratch_fixture, "sample(s) are incomplete")
        _, rows = read_tsv(scratch_fixture.outdir / "sample_flow.tsv")
        self.assertIn("retained_output_inside_scratch",
                      [item for item in rows
                       if item["sample_id"] == sample][0]["receipt_error"])

        # A retained output outside the sample's permitted roots.
        outside = self.fixture()
        sample = outside.samples[0]
        stray = outside.root / "elsewhere.txt"
        stray.write_text("stray\n", encoding="utf-8")
        outside.write_receipt(sample, payloads=[stray])
        self.seal_fails(outside, "sample(s) are incomplete")
        _, rows = read_tsv(outside.outdir / "sample_flow.tsv")
        self.assertIn("retained_output_outside_sample_roots",
                      [item for item in rows
                       if item["sample_id"] == sample][0]["receipt_error"])

        # Another sample's results directory is not a permitted root either.
        neighbour = self.fixture()
        sample, other = neighbour.samples[0], neighbour.samples[1]
        neighbour.write_receipt(
            sample,
            payloads=[neighbour.sample_root(other) / "sample_completion.tsv"])
        self.seal_fails(neighbour, "sample(s) are incomplete")
        _, rows = read_tsv(neighbour.outdir / "sample_flow.tsv")
        self.assertIn("retained_output_outside_sample_roots",
                      [item for item in rows
                       if item["sample_id"] == sample][0]["receipt_error"])

    # ------------------------------------------------------------ case 14
    def test_case14_source_seal_failures(self):
        corrupt = self.fixture()
        corrupt.write_source_seal(
            corrupt_member=corrupt.config["seal_manifest_member"])
        self.seal_fails(corrupt, "source-seal checksum mismatch")

        stale = self.fixture()
        stale.write_source_seal(in_progress=True)
        self.seal_fails(stale, "mid-audit")

        incomplete = self.fixture()
        incomplete.write_source_seal(success=False)
        self.seal_fails(incomplete, "no nonempty SUCCESS")

        dropped = self.fixture()
        dropped.write_source_seal(
            drop_member=dropped.config["seal_independent_member"])
        self.seal_fails(dropped, "does not contain the independent manifest copy")

        absent = self.fixture()
        for item in absent.seal.iterdir():
            item.unlink()
        self.seal_fails(absent, "no nonempty SUCCESS")

        no_checksum = self.fixture()
        (no_checksum.seal / "production_seal.sha256").unlink()
        self.seal_fails(no_checksum, "no nonempty production_seal.sha256")

        missing_member = self.fixture()
        (missing_member.seal
         / missing_member.config["seal_manifest_member"]).unlink()
        self.seal_fails(missing_member, "source-seal member missing or empty")

    # ------------------------------------------------------------ case 15
    def test_case15_manifest_seal_identity_mismatch_fails(self):
        drifted = self.fixture()
        rows = drifted.manifest.read_text(encoding="utf-8")
        drifted.manifest.write_text(rows.replace("S01", "S99"), encoding="utf-8")
        self.seal_fails(drifted,
                        "not byte-identical to the native source-seal copy")

        # A whitespace-only difference is still a difference.
        whitespace = self.fixture()
        whitespace.manifest.write_text(
            whitespace.manifest.read_text(encoding="utf-8") + "\n",
            encoding="utf-8")
        self.seal_fails(whitespace,
                        "not byte-identical to the native source-seal copy")

        # The independent manifest must stay an exact row-consistent subset.
        inconsistent = self.fixture()
        header, rows = read_tsv(inconsistent.independent_manifest)
        rows[0][inconsistent.config["condition_column"]] = "WRONG"
        write_tsv(inconsistent.independent_manifest, header, rows)
        inconsistent.write_source_seal()
        self.seal_fails(inconsistent, "but the production manifest has")

        # A sample that is not in the production manifest cannot be in the subset.
        unnested = self.fixture()
        header, rows = read_tsv(unnested.independent_manifest)
        rows[0]["sample_id"] = "GHOST"
        write_tsv(unnested.independent_manifest, header, rows)
        unnested.write_source_seal()
        self.seal_fails(unnested, "not nested in the production manifest")

    # ------------------------------------------------------------ case 16
    def test_case16_stale_authority_invalidated(self):
        fixture = self.fixture()
        outdir = self.seal_ok(fixture)
        first_success = (outdir / "SUCCESS").read_bytes()
        first_checksum = (outdir / "production_seal.sha256").read_bytes()
        self.assertTrue(first_success and first_checksum)

        # Break the cohort and rerun into the same directory.
        sample = fixture.samples[0]
        (fixture.sample_root(sample) / "profiles/baseline" / sample
         / "SUCCESS").unlink()
        done = fixture.run()
        self.assertNotEqual(done.returncode, 0)
        self.assertFalse((outdir / "SUCCESS").exists(),
                         "a failed rerun must not leave the older SUCCESS current")
        self.assertFalse((outdir / "production_seal.sha256").exists())
        self.assertTrue((outdir / "AUDIT_IN_PROGRESS").is_file())
        self.assertIn("IN_PROGRESS",
                      (outdir / "AUDIT_IN_PROGRESS").read_text(encoding="utf-8"))
        # The diagnostic ledger is still refreshed so the failure is auditable.
        _, rows = read_tsv(outdir / "sample_flow.tsv")
        self.assertIn("baseline_profile_count",
                      [item for item in rows
                       if item["sample_id"] == sample][0]["failure_reasons"])

    # ------------------------------------------------------------ case 17
    def test_case17_source_seal_unchanged(self):
        fixture = self.fixture()
        before = tree_digest(fixture.seal)
        self.assertTrue(before)
        self.seal_ok(fixture)
        self.assertEqual(tree_digest(fixture.seal), before,
                         "a successful audit modified the native source seal")
        sample = fixture.samples[0]
        (fixture.sample_root(sample) / "profiles/baseline" / sample
         / "SUCCESS").unlink()
        done = fixture.run()
        self.assertNotEqual(done.returncode, 0)
        self.assertEqual(tree_digest(fixture.seal), before,
                         "a failed audit modified the native source seal")

        # Also unchanged when the seal itself is what fails verification.
        corrupt = self.fixture()
        corrupt.write_source_seal(
            corrupt_member=corrupt.config["seal_manifest_member"])
        snapshot = tree_digest(corrupt.seal)
        self.seal_fails(corrupt, "source-seal checksum mismatch")
        self.assertEqual(tree_digest(corrupt.seal), snapshot)

    # ------------------------------------------------------------ case 18
    def test_case18_checksum_manifest_verifies(self):
        for cohort in ("yachida", "feng", "zeller"):
            outdir = self.seal_ok(self.fixture(cohort=cohort))
            listed = {}
            for line in (outdir / "production_seal.sha256").read_text(
                    encoding="utf-8").splitlines():
                recorded, name = line.split()
                listed[name] = recorded
            self.assertEqual(
                sorted(listed), sorted(set(SEAL_MEMBERS) - {"production_seal.sha256"}),
                cohort)
            for name, recorded in sorted(listed.items()):
                member = outdir / name
                self.assertTrue(member.is_file(), name)
                self.assertEqual(sha256_file(member), recorded, name)
            # Every top-level member except the checksum itself is covered.
            present = sorted(item.name for item in outdir.iterdir())
            self.assertEqual(present, sorted(SEAL_MEMBERS), cohort)
            inventory_header, inventory = read_tsv(
                outdir / "source_seal_inventory.tsv")
            self.assertEqual(
                inventory_header,
                ["member", "sha256", "bytes", "status", "source_audit_job"])
            self.assertTrue(all(row["status"] == "VERIFIED" for row in inventory))
            self.assertFalse(any("/" in row["member"] for row in inventory))
            expected_job = {"yachida": "", "feng": "3097587",
                            "zeller": "3097588"}[cohort]
            self.assertTrue(
                all(row["source_audit_job"] == expected_job for row in inventory),
                cohort)


class ConfigurationRefusalTest(SealFixtureTestCase):
    """Boundary configuration must be explicit and is refused when ambiguous."""

    def test_unknown_cohort_is_refused(self):
        fixture = self.fixture()
        command = fixture.command()
        command[command.index("--cohort") + 1] = "mystery"
        done = subprocess.run(command, text=True, capture_output=True)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("unknown cohort", done.stdout + done.stderr)

    def test_missing_required_column_is_refused(self):
        fixture = self.fixture()
        self.seal_fails(fixture, "lacks required column(s): absent",
                        extra=["--covariate-column", "absent"])

    def test_duplicate_manifest_column_is_refused(self):
        fixture = self.fixture()
        text = fixture.manifest.read_text(encoding="utf-8").split("\n")
        text[0] = text[0] + "\tcondition"
        fixture.manifest.write_text(
            "\n".join([text[0]] + [row + "\tX" if row else row
                                   for row in text[1:]]), encoding="utf-8")
        fixture.write_source_seal()
        self.seal_fails(fixture, "duplicate column(s): condition")

    def test_outdir_must_not_overlap_the_native_seal(self):
        fixture = self.fixture()
        command = fixture.command()
        command[command.index("--outdir") + 1] = str(fixture.seal)
        done = subprocess.run(command, text=True, capture_output=True)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("must not overlap the native source seal",
                      done.stdout + done.stderr)

    def test_relative_paths_are_refused(self):
        fixture = self.fixture()
        command = fixture.command()
        command[command.index("--outdir") + 1] = "relative/seal"
        done = subprocess.run(command, text=True, capture_output=True)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("must be an absolute path", done.stdout + done.stderr)

    def test_unexpected_outdir_member_is_refused(self):
        fixture = self.fixture()
        fixture.outdir.mkdir(parents=True)
        (fixture.outdir / "leftover.tsv").write_text("x\n", encoding="utf-8")
        self.seal_fails(fixture, "unexpected member(s): leftover.tsv")

    def test_expected_sample_count_is_enforced(self):
        # The native SUCCESS is checked against the expectation first.
        fixture = self.fixture()
        command = fixture.command()
        command[command.index("--expected-samples") + 1] = "5"
        done = subprocess.run(command, text=True, capture_output=True)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("native source SUCCESS reports samples=6",
                      done.stdout + done.stderr)
        # With a seal that agrees with the wrong expectation, the manifest
        # row-count gate is what refuses it.
        consistent = self.fixture()
        consistent.write_source_seal(
            success_overrides={"samples": 5, "baseline": 5,
                               "community": 5 * COMMUNITY_PER_SAMPLE})
        command = consistent.command()
        command[command.index("--expected-samples") + 1] = "5"
        done = subprocess.run(command, text=True, capture_output=True)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("expected 5 samples; found 6", done.stdout + done.stderr)


# ------------------------------------------- Codex review corrections
class SourceSealOverlapTest(SealFixtureTestCase):
    """Issue 1: overlap is refused before anything in --outdir is touched."""

    def run_with_outdir(self, fixture, outdir):
        command = fixture.command()
        command[command.index("--outdir") + 1] = str(outdir)
        return subprocess.run(command, text=True, capture_output=True)

    def assert_refused(self, fixture, outdir, created_check=None):
        before = tree_digest(fixture.seal)
        listing = sorted(item.name for item in fixture.seal.iterdir())
        done = self.run_with_outdir(fixture, outdir)
        self.assertNotEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertIn("must not overlap the native source seal",
                      done.stdout + done.stderr)
        self.assertEqual(tree_digest(fixture.seal), before,
                         "the native seal was modified by a refused run")
        self.assertEqual(sorted(item.name for item in fixture.seal.iterdir()),
                         listing, "a refused run added or removed seal members")
        if created_check is not None:
            self.assertFalse(created_check.exists(),
                             "a refused run created %s" % created_check)

    def test_outdir_equal_to_source_seal_is_refused(self):
        fixture = self.fixture()
        self.assert_refused(fixture, fixture.seal)

    def test_outdir_inside_source_seal_is_refused(self):
        fixture = self.fixture()
        nested = fixture.seal / "production_seal_v2"
        self.assert_refused(fixture, nested, created_check=nested)

    def test_source_seal_inside_outdir_is_refused(self):
        fixture = self.fixture()
        parent = fixture.seal.parent / "enclosing"
        parent.mkdir()
        seal = parent / "native_seal"
        shutil.copytree(str(fixture.seal), str(seal))
        fixture.seal = seal
        self.assert_refused(fixture, parent)

    def test_refusal_does_not_clear_existing_authority(self):
        """Stale authority elsewhere must survive a path-refused invocation."""
        fixture = self.fixture()
        outdir = self.seal_ok(fixture)
        before = tree_digest(outdir)
        nested = fixture.seal / "inner"
        done = self.run_with_outdir(fixture, nested)
        self.assertNotEqual(done.returncode, 0)
        self.assertFalse(nested.exists())
        self.assertEqual(tree_digest(outdir), before)


class CovariateAuditTest(SealFixtureTestCase):
    """Issue 2: covariates are read through the original-to-canonical mapping."""

    def audit_rows(self, outdir):
        header, rows = read_tsv(outdir / "covariate_audit.tsv")
        self.assertEqual(header, ["field", "missing", "present",
                                  "distinct_nonmissing"])
        return {row["field"]: row for row in rows}

    def test_feng_covariates_report_real_values(self):
        fixture = self.fixture(cohort="feng")
        rows = self.audit_rows(self.seal_ok(fixture))
        # Six samples; the first has blank age/sex/bmi, the rest are populated.
        for field, distinct in (("age", 5), ("sex", 2), ("bmi", 5)):
            self.assertEqual(rows[field]["missing"], "1", field)
            self.assertEqual(rows[field]["present"], "5", field)
            self.assertEqual(rows[field]["distinct_nonmissing"], str(distinct),
                             field)
        self.assertEqual(rows["condition"]["missing"], "0")
        self.assertEqual(rows["condition"]["present"], "6")
        self.assertEqual(rows["condition"]["distinct_nonmissing"], "3")
        self.assertEqual(rows["study"]["present"], "6")
        self.assertEqual(rows["study"]["distinct_nonmissing"], "1")
        self.assertEqual(rows["independent_subset"]["present"], "6")
        self.assertEqual(rows["independent_subset"]["distinct_nonmissing"], "2")

    def test_yachida_audits_the_fields_its_manifest_really_has(self):
        """The sealed production manifest carries age, sex, bmi and batch_id."""
        fixture = self.fixture(cohort="yachida")
        rows = self.audit_rows(self.seal_ok(fixture))
        self.assertEqual(
            sorted(rows),
            ["age", "batch_id", "bmi", "condition", "independent_subset",
             "sex", "study"])
        self.assertEqual(rows["batch_id"]["missing"], "0")
        self.assertEqual(rows["batch_id"]["present"], "6")
        self.assertEqual(rows["batch_id"]["distinct_nonmissing"], "2")
        for field, distinct in (("age", 5), ("sex", 2), ("bmi", 5)):
            self.assertEqual(rows[field]["missing"], "1", field)
            self.assertEqual(rows[field]["present"], "5", field)
            self.assertEqual(rows[field]["distinct_nonmissing"], str(distinct),
                             field)

    def test_a_colliding_covariate_audits_the_original_column(self):
        """A manifest column named `condition` is renamed `source_condition` in
        the canonical copy; auditing it must still read the original values."""
        fixture = self.fixture(cohort="yachida")
        fixture.write_manifests(extra_shared={
            "condition": ["" if index < 2 else "recorded"
                          for index in range(len(fixture.samples))]})
        fixture.write_source_seal()
        outdir = self.seal_ok(fixture,
                              extra=["--covariate-column", "condition"])
        rows = self.audit_rows(outdir)
        # Reported under the canonical manifest's own column name, with the
        # ORIGINAL values: two blanks, one distinct value.
        self.assertEqual(rows["source_condition"]["missing"], "2")
        self.assertEqual(rows["source_condition"]["present"], "4")
        self.assertEqual(rows["source_condition"]["distinct_nonmissing"], "1")
        # The normalized condition is still audited separately and is complete.
        self.assertEqual(rows["condition"]["missing"], "0")
        self.assertEqual(rows["condition"]["distinct_nonmissing"], "3")
        # Not the canonical condition, which is complete with three values.
        manifest_header, manifest_rows = read_tsv(
            outdir / "production_manifest.tsv")
        self.assertIn("source_condition", manifest_header)
        self.assertEqual(manifest_rows[0]["condition"], "Control")
        self.assertEqual(manifest_rows[0]["source_condition"], "")

    def test_a_wholly_missing_covariate_is_reported_as_missing(self):
        fixture = self.fixture(cohort="feng")
        header, rows = read_tsv(fixture.manifest)
        for row in rows:
            row["bmi"] = ""
        write_tsv(fixture.manifest, header, rows)
        header, rows = read_tsv(fixture.independent_manifest)
        for row in rows:
            row["bmi"] = ""
        write_tsv(fixture.independent_manifest, header, rows)
        fixture.write_source_seal()
        audited = self.audit_rows(self.seal_ok(fixture))
        self.assertEqual(audited["bmi"]["missing"], "6")
        self.assertEqual(audited["bmi"]["present"], "0")
        self.assertEqual(audited["bmi"]["distinct_nonmissing"], "0")


class IndependentSelectionProvenanceTest(SealFixtureTestCase):
    """Issue 3: independent-only selection columns are legitimate."""

    def test_independent_only_columns_are_accepted_and_preserved(self):
        fixture = self.fixture(cohort="yachida")
        fixture.write_manifests(
            independent_extra=("selection_note", "selection_stage"))
        fixture.write_source_seal()
        outdir = self.seal_ok(fixture)
        header, rows = read_tsv(outdir / "production_manifest.independent.tsv")
        self.assertEqual(header[:4],
                         ["sample_id", "condition", "study", "independent_subset"])
        self.assertIn("selection_note", header)
        self.assertIn("selection_stage", header)
        self.assertEqual(rows[0]["selection_note"], "selection_note_0")
        # The production manifest keeps its own, narrower header.
        production_header, _ = read_tsv(outdir / "production_manifest.tsv")
        self.assertNotIn("selection_note", production_header)

    def test_selection_column_colliding_with_canonical_is_renamed(self):
        fixture = self.fixture(cohort="yachida")
        fixture.write_manifests(independent_extra=("study",))
        fixture.write_source_seal()
        outdir = self.seal_ok(fixture)
        header, rows = read_tsv(outdir / "production_manifest.independent.tsv")
        self.assertIn("source_study", header)
        self.assertEqual(rows[0]["study"], "YachidaS_2019")
        self.assertEqual(rows[0]["source_study"], "study_0")

    def test_altered_shared_value_fails(self):
        # batch_id is production-only, so an altered shared value is tested on
        # a column the independent manifest genuinely carries.
        for column, value in (("Target_Condition", "WRONG"), ("age", "999")):
            fixture = self.fixture(cohort="yachida")
            fixture.write_manifests(
                independent_extra=("selection_note",),
                independent_overrides={column: value})
            fixture.write_source_seal()
            self.seal_fails(fixture, "but the production manifest has")

        covariate = self.fixture(cohort="feng")
        covariate.write_manifests(independent_overrides={"bmi": "99.9"})
        covariate.write_source_seal()
        self.seal_fails(covariate, "but the production manifest has")

    def test_missing_shared_required_column_fails(self):
        fixture = self.fixture(cohort="feng")
        fixture.write_manifests(independent_drop=("study",))
        fixture.write_source_seal()
        self.seal_fails(fixture, "lacks required column(s): study")

        condition = self.fixture(cohort="yachida")
        condition.write_manifests(independent_drop=("Target_Condition",))
        condition.write_source_seal()
        self.seal_fails(condition, "lacks required column(s): Target_Condition")

    def test_duplicate_independent_header_still_fails(self):
        """Even under the Yachida adapter, any other duplicate is refused."""
        fixture = self.fixture(cohort="yachida")
        fixture.write_manifests(duplicate_fields=("age",))
        fixture.write_source_seal()
        self.seal_fails(
            fixture, "outside the documented historical selection triplet: age")

        feng = self.fixture(cohort="feng")
        feng.write_manifests(duplicate_selection=True, duplicate_fields=("age",))
        feng.write_source_seal()
        self.seal_fails(feng, "duplicate column(s): age")

    def test_extra_columns_do_not_weaken_identity_or_subset_checks(self):
        fixture = self.fixture(cohort="yachida")
        fixture.write_manifests(independent_extra=("selection_note",),
                                independent_overrides={"sample_id": "GHOST"})
        fixture.write_source_seal()
        self.seal_fails(fixture, "not nested in the production manifest")


class InterruptedTemporaryTest(SealFixtureTestCase):
    """Issue 4: a recognized interrupted temporary is recoverable."""

    def test_recognized_temporaries_are_cleared_and_the_rerun_succeeds(self):
        fixture = self.fixture()
        fixture.outdir.mkdir(parents=True)
        leftovers = [fixture.outdir / "sample_flow.tsv.tmp",
                     fixture.outdir / "production_seal.sha256.tmp",
                     fixture.outdir / "SUCCESS.tmp"]
        for path in leftovers:
            path.write_text("partial\n", encoding="utf-8")
        outdir = self.seal_ok(fixture)
        for path in leftovers:
            self.assertFalse(path.exists(), path.name)
        self.assertEqual(sorted(item.name for item in outdir.iterdir()),
                         sorted(SEAL_MEMBERS))

    def test_interrupted_run_leaves_a_recoverable_directory(self):
        """Stale authority plus a stray temporary still reruns cleanly."""
        fixture = self.fixture()
        self.seal_ok(fixture)
        (fixture.outdir / "covariate_audit.tsv.tmp").write_text(
            "partial\n", encoding="utf-8")
        outdir = self.seal_ok(fixture)
        self.assertFalse((outdir / "covariate_audit.tsv.tmp").exists())
        self.assertEqual(sorted(item.name for item in outdir.iterdir()),
                         sorted(SEAL_MEMBERS))

    def test_arbitrary_unexpected_files_are_still_refused(self):
        fixture = self.fixture()
        fixture.outdir.mkdir(parents=True)
        (fixture.outdir / "notes.txt.tmp").write_text("x\n", encoding="utf-8")
        self.seal_fails(fixture, "unexpected member(s): notes.txt.tmp")

    def test_stale_authority_is_dropped_even_when_temporaries_exist(self):
        """A failing run never writes SUCCESS, so only the explicit cleanup can
        remove a leftover SUCCESS.tmp."""
        fixture = self.fixture()
        outdir = self.seal_ok(fixture)
        for name in ("SUCCESS.tmp", "production_seal.sha256.tmp",
                     "sample_flow.tsv.tmp"):
            (outdir / name).write_text("partial\n", encoding="utf-8")
        sample = fixture.samples[0]
        (fixture.sample_root(sample) / "profiles/baseline" / sample
         / "SUCCESS").unlink()
        done = fixture.run()
        self.assertNotEqual(done.returncode, 0)
        self.assertFalse((outdir / "SUCCESS").exists())
        self.assertFalse((outdir / "production_seal.sha256").exists())
        for name in ("SUCCESS.tmp", "production_seal.sha256.tmp",
                     "sample_flow.tsv.tmp"):
            self.assertFalse((outdir / name).exists(), name)
        self.assertTrue((outdir / "AUDIT_IN_PROGRESS").is_file())
        # And the directory is still recoverable by a clean rerun.
        marker = (fixture.sample_root(sample) / "profiles/baseline" / sample
                  / "SUCCESS")
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.touch()
        self.assertEqual(sorted(item.name for item in self.seal_ok(
            fixture).iterdir()), sorted(SEAL_MEMBERS))


class NativeSuccessSemanticsTest(SealFixtureTestCase):
    """Issue 5: a checksummed SUCCESS must also describe this cohort."""

    def rewrite_success(self, fixture, **overrides):
        fixture.write_source_seal(success_overrides=overrides)

    def test_wrong_cohort_identity_fails(self):
        fixture = self.fixture(cohort="feng")
        self.rewrite_success(fixture, identity="zeller")
        self.seal_fails(fixture, "native source SUCCESS cohort is 'zeller'")

    def test_wrong_dataset_identity_fails(self):
        fixture = self.fixture(cohort="yachida")
        self.rewrite_success(fixture, identity="SomeOther_2019")
        self.seal_fails(fixture, "native source SUCCESS dataset is 'SomeOther_2019'")

    def test_wrong_sample_count_fails(self):
        for cohort in ("yachida", "feng"):
            fixture = self.fixture(cohort=cohort)
            self.rewrite_success(fixture, samples=99)
            self.seal_fails(fixture, "reports samples=99")

    def test_wrong_independent_count_fails(self):
        fixture = self.fixture(cohort="feng")
        self.rewrite_success(fixture, independent=7)
        self.seal_fails(fixture, "reports independent_samples=7")

    def test_wrong_profile_topology_fails(self):
        feng = self.fixture(cohort="feng")
        self.rewrite_success(feng, community=999)
        self.seal_fails(feng, "reports community_profiles=999")

        yachida = self.fixture(cohort="yachida")
        self.rewrite_success(yachida, independent_profiles=12345)
        self.seal_fails(yachida, "reports independent profiles=12345")

    def test_non_pass_status_fails(self):
        for cohort in ("yachida", "feng"):
            fixture = self.fixture(cohort=cohort)
            self.rewrite_success(fixture, status="FAIL")
            self.seal_fails(fixture, "status is 'FAIL', not PASS")

    def test_condition_counts_in_the_native_success_must_agree(self):
        fixture = self.fixture(cohort="yachida")
        text = (fixture.seal / "SUCCESS").read_text(encoding="utf-8")
        fixture.write_source_seal()
        replaced = text.replace("Adenoma=2", "Adenoma=5")
        self.assertNotEqual(replaced, text)
        (fixture.seal / "SUCCESS").write_text(replaced, encoding="utf-8")
        # Re-checksum so only the semantics are wrong, not the integrity.
        lines = []
        for line in (fixture.seal / "production_seal.sha256").read_text(
                encoding="utf-8").splitlines():
            _, member = line.split()
            lines.append("%s  %s\n" % (sha256_file(fixture.seal / member), member))
        (fixture.seal / "production_seal.sha256").write_text(
            "".join(lines), encoding="utf-8")
        self.seal_fails(fixture, "conditions")

    def test_missing_required_success_field_fails(self):
        fixture = self.fixture(cohort="feng")
        (fixture.seal / "SUCCESS").write_text(
            "cohort\tfeng\nstatus\tPASS\n", encoding="utf-8")
        fixture.write_source_seal()
        (fixture.seal / "SUCCESS").write_text(
            "cohort\tfeng\nstatus\tPASS\n", encoding="utf-8")
        lines = []
        for line in (fixture.seal / "production_seal.sha256").read_text(
                encoding="utf-8").splitlines():
            _, member = line.split()
            lines.append("%s  %s\n" % (sha256_file(fixture.seal / member), member))
        (fixture.seal / "production_seal.sha256").write_text(
            "".join(lines), encoding="utf-8")
        self.seal_fails(fixture, "lacks 'samples' for the cohort_profiles schema")


class SourceSealRevalidationTest(SealFixtureTestCase):
    """Issue 6: the seal is reverified before authority is granted."""

    def module(self):
        import importlib.util

        spec = importlib.util.spec_from_file_location(
            "seal_cohort_upstream_under_test", str(AUDITOR))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def test_revalidation_detects_a_seal_that_changed(self):
        seal = self.module()
        fixture = self.fixture()
        adapter = dict(seal.ADAPTERS[fixture.cohort])
        baseline = seal.source_seal_inventory(fixture.seal, adapter)
        seal.reverify_source_seal(fixture.seal, adapter, fixture.manifest,
                                  fixture.independent_manifest, baseline)

        member = fixture.seal / fixture.config["seal_extra"][0]
        member.write_text(member.read_text(encoding="utf-8") + "drift\n",
                          encoding="utf-8")
        with self.assertRaises(seal.SealError) as caught:
            seal.reverify_source_seal(fixture.seal, adapter, fixture.manifest,
                                      fixture.independent_manifest, baseline)
        self.assertIn("source-seal checksum mismatch", str(caught.exception))

    def test_revalidation_detects_a_new_audit_in_progress(self):
        seal = self.module()
        fixture = self.fixture()
        adapter = dict(seal.ADAPTERS[fixture.cohort])
        baseline = seal.source_seal_inventory(fixture.seal, adapter)
        (fixture.seal / "AUDIT_IN_PROGRESS").write_text(
            "status\tIN_PROGRESS\n", encoding="utf-8")
        with self.assertRaises(seal.SealError) as caught:
            seal.reverify_source_seal(fixture.seal, adapter, fixture.manifest,
                                      fixture.independent_manifest, baseline)
        self.assertIn("became mid-audit during this audit", str(caught.exception))

    def test_mutation_during_the_audit_blocks_authority(self):
        """End to end: the seal is swapped while the audit is running.

        The revalidation hook is wrapped rather than stubbed, so removing the
        call from the auditor makes this test fail instead of silently pass.
        """
        seal = self.module()
        fixture = self.fixture()
        mutated = []
        original = seal.reverify_source_seal

        def mutate_then_reverify(seal_root, adapter, manifest, independent,
                                 baseline):
            if not mutated:
                mutated.append(True)
                member = seal_root / fixture.config["seal_extra"][0]
                member.write_text(
                    member.read_text(encoding="utf-8") + "mid-audit drift\n",
                    encoding="utf-8")
            return original(seal_root, adapter, manifest, independent, baseline)

        seal.reverify_source_seal = mutate_then_reverify
        argv = fixture.command()[1:]
        saved = sys.argv
        try:
            sys.argv = argv
            with self.assertRaises(SystemExit) as caught:
                seal.main()
        finally:
            sys.argv = saved
        self.assertTrue(mutated, "the auditor never reverified the source seal")
        self.assertIn("changed during this audit", str(caught.exception))
        self.assertFalse((fixture.outdir / "SUCCESS").exists())
        self.assertFalse((fixture.outdir / "production_seal.sha256").exists())
        self.assertTrue((fixture.outdir / "AUDIT_IN_PROGRESS").is_file())


# ------------------------------------ observed real Yachida header evidence
class HistoricalYachidaHeaderTest(SealFixtureTestCase):
    """Issue 3: the sealed duplicated header is translated, never rewritten."""

    def independent_header(self, fixture):
        with fixture.independent_manifest.open(encoding="utf-8") as handle:
            return handle.readline().rstrip("\n").split("\t")

    def test_observed_header_shape_is_reproduced_by_the_fixture(self):
        fixture = self.fixture(cohort="yachida")
        production = self.independent_header(fixture)
        with fixture.manifest.open(encoding="utf-8") as handle:
            production_header = handle.readline().rstrip("\n").split("\t")
        # Production: one selection triplet, and batch provenance.
        for field in SELECTION_TRIPLET:
            self.assertEqual(production_header.count(field), 1, field)
        self.assertIn("batch_id", production_header)
        for field in ("age", "sex", "bmi"):
            self.assertIn(field, production_header)
        # Independent: two triplets, no batch provenance.
        for field in SELECTION_TRIPLET:
            self.assertEqual(production.count(field), 2, field)
        self.assertNotIn("batch_id", production)
        for field in ("age", "sex", "bmi"):
            self.assertIn(field, production)

    def test_both_historical_copies_hold_the_independent_value(self):
        """The old selector overwrote the inherited values before writing."""
        fixture = self.fixture(cohort="yachida")
        header = self.independent_header(fixture)
        with fixture.independent_manifest.open(encoding="utf-8") as handle:
            handle.readline()
            row = handle.readline().rstrip("\n").split("\t")
        for field in SELECTION_TRIPLET:
            first = header.index(field)
            second = len(header) - 1 - header[::-1].index(field)
            self.assertNotEqual(first, second)
            self.assertEqual(row[first], row[second],
                             "both historical copies must be identical")
            self.assertTrue(row[first].startswith("independent_"), field)
        # The pilot value exists only in the production manifest.
        production_header, production_rows = read_tsv(fixture.manifest)
        self.assertEqual(production_rows[0]["selection_hash"],
                         "pilot_selection_hash_0")
        self.assertNotIn("pilot_selection_hash_0", row)

    def test_canonical_output_reconstructs_pilot_provenance(self):
        fixture = self.fixture(cohort="yachida")
        outdir = self.seal_ok(fixture)
        header, rows = read_tsv(outdir / "production_manifest.independent.tsv")
        for field in SELECTION_TRIPLET:
            self.assertNotIn(field, header, field)
            self.assertIn("pilot_%s" % field, header)
            self.assertIn("independent_%s" % field, header)
            self.assertEqual(header.count("pilot_%s" % field), 1)
            self.assertEqual(header.count("independent_%s" % field), 1)
        self.assertEqual(len(header), len(set(header)),
                         "the canonical header must be unambiguous")
        # Pilot provenance is reconstructed from the production manifest row;
        # the independent triplet is the historical duplicated value.
        _, production_rows = read_tsv(fixture.manifest)
        production = {row["sample_id"]: row for row in production_rows}
        for row in rows:
            source = production[row["sample_id"]]
            for field in SELECTION_TRIPLET:
                self.assertEqual(row["pilot_%s" % field], source[field], field)
                self.assertTrue(
                    row["independent_%s" % field].startswith("independent_"),
                    field)
                self.assertNotEqual(row["pilot_%s" % field],
                                    row["independent_%s" % field], field)
        self.assertEqual(rows[0]["pilot_selection_rank"],
                         "pilot_selection_rank_0")
        self.assertEqual(rows[0]["independent_selection_rank"],
                         "independent_selection_rank_0")
        # The sealed file itself is untouched and still duplicated.
        self.assertEqual(
            self.independent_header(fixture).count("selection_rank"), 2)

    def corrupt_occurrence(self, fixture, field, which):
        lines = fixture.independent_manifest.read_text(
            encoding="utf-8").rstrip("\n").split("\n")
        header = lines[0].split("\t")
        if which == "first":
            position = header.index(field)
        else:
            position = len(header) - 1 - header[::-1].index(field)
        row = lines[1].split("\t")
        row[position] = "tampered"
        lines[1] = "\t".join(row)
        fixture.independent_manifest.write_text(
            "\n".join(lines) + "\n", encoding="utf-8")
        fixture.write_source_seal()

    def test_changing_the_first_historical_copy_fails(self):
        """Neither raw copy is pilot provenance; they must simply agree."""
        fixture = self.fixture(cohort="yachida")
        self.corrupt_occurrence(fixture, "selection_hash", "first")
        message = self.seal_fails(
            fixture, "the two historical selection_hash columns differ")
        self.assertIn("tampered", message)

    def test_changing_the_second_historical_copy_fails(self):
        fixture = self.fixture(cohort="yachida")
        self.corrupt_occurrence(fixture, "selection_seed", "second")
        message = self.seal_fails(
            fixture, "the two historical selection_seed columns differ")
        self.assertIn("tampered", message)

    def test_every_historical_field_is_cross_checked(self):
        for field in SELECTION_TRIPLET:
            for which in ("first", "second"):
                fixture = self.fixture(cohort="yachida")
                self.corrupt_occurrence(fixture, field, which)
                self.seal_fails(
                    fixture, "the two historical %s columns differ" % field)

    def test_production_must_still_carry_the_pilot_triplet(self):
        """Format A's pilot provenance has nowhere else to come from."""
        fixture = self.fixture(cohort="yachida")
        header, rows = read_tsv(fixture.manifest)
        header = [name for name in header if name != "selection_hash"]
        write_tsv(fixture.manifest, header, rows)
        fixture.write_source_seal()
        self.seal_fails(fixture, "lacks the pilot selection column(s)")

    PREFIXED = tuple("pilot_%s" % field for field in SELECTION_TRIPLET) + tuple(
        "independent_%s" % field for field in SELECTION_TRIPLET)

    def prefixed_fixture(self, extra=PREFIXED, overrides=None):
        """A format-B manifest: prefixed provenance, no bare triplet."""
        fixture = self.fixture(cohort="yachida")
        fixture.write_manifests(duplicate_selection=False,
                                independent_drop=SELECTION_TRIPLET,
                                independent_extra=extra,
                                independent_overrides=overrides)
        fixture.write_source_seal()
        return fixture

    def test_single_bare_occurrence_is_refused(self):
        """A lone bare triplet is neither the historical nor the corrected
        format, so its provenance is ambiguous and it is refused."""
        fixture = self.fixture(cohort="yachida")
        fixture.write_manifests(duplicate_selection=False)
        fixture.write_source_seal()
        message = self.seal_fails(fixture, "must appear exactly twice each")
        self.assertIn("selection_rank=1", message)
        self.assertIn("pilot_selection_rank=0", message)

    def test_partial_pilot_triplet_is_refused(self):
        fixture = self.prefixed_fixture(
            extra=("pilot_selection_rank", "pilot_selection_hash")
            + tuple("independent_%s" % field for field in SELECTION_TRIPLET))
        message = self.seal_fails(fixture, "complete pilot/independent pair")
        self.assertIn("pilot_selection_seed=0", message)

    def test_partial_independent_triplet_is_refused(self):
        fixture = self.prefixed_fixture(
            extra=tuple("pilot_%s" % field for field in SELECTION_TRIPLET)
            + ("independent_selection_rank",))
        message = self.seal_fails(fixture, "complete pilot/independent pair")
        self.assertIn("independent_selection_seed=0", message)

    def test_missing_pilot_provenance_is_refused(self):
        fixture = self.prefixed_fixture(
            extra=tuple("independent_%s" % field for field in SELECTION_TRIPLET))
        self.seal_fails(fixture, "complete pilot/independent pair")

    def test_missing_independent_provenance_is_refused(self):
        fixture = self.prefixed_fixture(
            extra=tuple("pilot_%s" % field for field in SELECTION_TRIPLET))
        self.seal_fails(fixture, "complete pilot/independent pair")

    def test_no_selection_provenance_at_all_is_refused(self):
        fixture = self.fixture(cohort="yachida")
        fixture.write_manifests(duplicate_selection=False,
                                independent_drop=SELECTION_TRIPLET)
        fixture.write_source_seal()
        self.seal_fails(fixture, "must appear exactly twice each")

    def test_bare_and_prefixed_mixture_is_refused(self):
        """One bare triplet beside a complete prefixed pair is still ambiguous."""
        fixture = self.fixture(cohort="yachida")
        fixture.write_manifests(duplicate_selection=False,
                                independent_extra=self.PREFIXED)
        fixture.write_source_seal()
        message = self.seal_fails(fixture, "complete pilot/independent pair")
        self.assertIn("selection_rank=1", message)

    def test_three_occurrences_fail(self):
        fixture = self.fixture(cohort="yachida")
        fixture.write_manifests(occurrences=3)
        fixture.write_source_seal()
        self.seal_fails(fixture, "must appear exactly twice each")

    def test_incomplete_duplicated_triplet_fails(self):
        for partial in (("selection_rank",),
                        ("selection_rank", "selection_hash")):
            fixture = self.fixture(cohort="yachida")
            fixture.write_manifests(duplicate_fields=partial)
            self.assertTrue(fixture.independent_manifest.is_file())
            fixture.write_source_seal()
            self.seal_fails(fixture, "must appear exactly twice each")

    def test_historical_triplet_beside_a_prefixed_field_is_refused(self):
        """The duplicated triplet cannot be translated over an existing name."""
        fixture = self.fixture(cohort="yachida")
        fixture.write_manifests(independent_extra=("pilot_selection_rank",))
        fixture.write_source_seal()
        message = self.seal_fails(
            fixture, "mixes the historical duplicated selection triplet")
        self.assertIn("pilot_selection_rank=1", message)

    def test_feng_independent_duplicates_are_never_translated(self):
        """The occurrence adapter is Yachida-only; Feng gets no translation."""
        fixture = self.fixture(cohort="feng")
        extra = {name: ["%s_%d" % (name, index)
                        for index in range(len(fixture.samples))]
                 for name in SELECTION_TRIPLET}
        fixture.write_manifests(extra_shared=extra, duplicate_selection=True)
        fixture.write_source_seal()
        message = self.seal_fails(fixture, "duplicate column(s)")
        self.assertNotIn("pilot_", message)
        self.assertNotIn("historical", message)

    def test_production_manifest_duplicates_are_always_refused(self):
        fixture = self.fixture(cohort="yachida")
        text = fixture.manifest.read_text(encoding="utf-8").split("\n")
        text[0] = text[0] + "\tselection_rank"
        fixture.manifest.write_text(
            "\n".join([text[0]] + [row + "\tX" if row else row
                                   for row in text[1:]]), encoding="utf-8")
        fixture.write_source_seal()
        self.seal_fails(fixture, "duplicate column(s): selection_rank")

    def test_future_prefixed_manifest_is_accepted(self):
        """The corrected selector's unambiguous output also seals."""
        outdir = self.seal_ok(self.prefixed_fixture())
        header, rows = read_tsv(outdir / "production_manifest.independent.tsv")
        for field in SELECTION_TRIPLET:
            self.assertIn("pilot_%s" % field, header)
            self.assertIn("independent_%s" % field, header)
        self.assertEqual(len(header), len(set(header)))

    def test_prefixed_pilot_values_are_still_verified(self):
        fixture = self.prefixed_fixture(
            overrides={"pilot_selection_rank": "wrong"})
        message = self.seal_fails(fixture, "pilot_selection_rank='wrong'")
        self.assertIn("the production manifest has selection_rank=", message)

    def test_prefixed_independent_values_are_not_compared(self):
        fixture = self.prefixed_fixture(
            overrides={"independent_selection_seed": "a-new-seed"})
        outdir = self.seal_ok(fixture)
        _, rows = read_tsv(outdir / "production_manifest.independent.tsv")
        self.assertEqual(rows[0]["independent_selection_seed"], "a-new-seed")


class ManifestAsymmetryTest(SealFixtureTestCase):
    """Issue 2: production-only batch fields are not required downstream."""

    def test_batch_provenance_absent_from_the_independent_manifest_is_fine(self):
        fixture = self.fixture(cohort="yachida")
        with fixture.independent_manifest.open(encoding="utf-8") as handle:
            header = handle.readline().rstrip("\n").split("\t")
        self.assertNotIn("batch_id", header)
        outdir = self.seal_ok(fixture)
        _, rows = read_tsv(outdir / "sample_flow.tsv")
        self.assertTrue(all(row["batch_id"].startswith("batch_") for row in rows))

    def test_other_production_only_batch_fields_are_not_required(self):
        fixture = self.fixture(cohort="yachida")
        extra = {name: ["%s_%d" % (name, index)
                        for index in range(len(fixture.samples))]
                 for name in ("batch_hash", "batch_position", "batch_size",
                              "batch_seed", "processing_order")}
        fixture.write_manifests(extra_shared=extra,
                                independent_drop=tuple(extra))
        fixture.write_source_seal()
        self.seal_ok(fixture)

    def test_batch_column_missing_from_production_still_fails(self):
        fixture = self.fixture(cohort="yachida")
        header, rows = read_tsv(fixture.manifest)
        position = header.index("batch_id")
        header.pop(position)
        write_tsv(fixture.manifest, header, rows)
        fixture.write_source_seal()
        self.seal_fails(fixture, "lacks required column(s): batch_id")

    def test_shared_covariate_missing_from_independent_still_fails(self):
        fixture = self.fixture(cohort="yachida")
        fixture.write_manifests(independent_drop=("bmi",))
        fixture.write_source_seal()
        self.seal_fails(fixture, "lacks required column(s): bmi")


class StaleAuthorityBeforeInspectionTest(SealFixtureTestCase):
    """Issue 5: authority dies before the directory contents are judged."""

    def test_unexpected_member_in_a_valid_seal_revokes_authority(self):
        fixture = self.fixture()
        outdir = self.seal_ok(fixture)
        seal_before = tree_digest(fixture.seal)
        intruder = outdir / "leftover_notes.txt"
        intruder.write_text("operator note\n", encoding="utf-8")
        done = fixture.run()
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("unexpected member(s): leftover_notes.txt",
                      done.stdout + done.stderr)
        self.assertFalse((outdir / "SUCCESS").exists(),
                         "stale SUCCESS survived a refused rerun")
        self.assertFalse((outdir / "production_seal.sha256").exists())
        self.assertTrue(intruder.is_file(),
                        "the unexpected file must be preserved for the operator")
        self.assertEqual(intruder.read_text(encoding="utf-8"), "operator note\n")
        self.assertEqual(tree_digest(fixture.seal), seal_before)
        # No AUDIT_IN_PROGRESS either: the run stopped before claiming the run.
        self.assertFalse((outdir / "AUDIT_IN_PROGRESS").exists())


class NativeSuccessCompletenessTest(SealFixtureTestCase):
    """Issue 6: topology fields are required, not optionally inspected."""

    def reseal_with(self, fixture, text):
        (fixture.seal / "SUCCESS").write_text(text, encoding="utf-8")
        fixture.write_source_seal()
        (fixture.seal / "SUCCESS").write_text(text, encoding="utf-8")
        lines = []
        for line in (fixture.seal / "production_seal.sha256").read_text(
                encoding="utf-8").splitlines():
            _, member = line.split()
            lines.append("%s  %s\n" % (sha256_file(fixture.seal / member), member))
        (fixture.seal / "production_seal.sha256").write_text(
            "".join(lines), encoding="utf-8")

    def drop_line(self, fixture, field):
        text = (fixture.seal / "SUCCESS").read_text(encoding="utf-8")
        kept = [line for line in text.rstrip("\n").split("\n")
                if not line.startswith(field + "\t")]
        self.assertEqual(len(kept), len(text.rstrip("\n").split("\n")) - 1, field)
        self.reseal_with(fixture, "\n".join(kept) + "\n")

    def test_every_cohort_profiles_topology_field_is_required(self):
        for field in ("cohort", "samples", "independent_samples",
                      "baseline_profiles", "independent_profiles",
                      "community_profiles", "status"):
            fixture = self.fixture(cohort="feng")
            self.drop_line(fixture, field)
            self.seal_fails(fixture, "native source SUCCESS")

    def test_every_yachida_topology_field_is_required(self):
        for field in ("dataset", "samples", "independent_subset", "profiles",
                      "conditions", "status"):
            fixture = self.fixture(cohort="yachida")
            self.drop_line(fixture, field)
            self.seal_fails(fixture, "native source SUCCESS")

    def test_incomplete_yachida_profile_summary_fails(self):
        fixture = self.fixture(cohort="yachida")
        text = (fixture.seal / "SUCCESS").read_text(encoding="utf-8")
        lines = []
        for line in text.rstrip("\n").split("\n"):
            if line.startswith("profiles\t"):
                line = "profiles\tbaseline=6;community=42"
            lines.append(line)
        self.reseal_with(fixture, "\n".join(lines) + "\n")
        self.seal_fails(fixture, "profiles lacks independent")

    def test_incomplete_yachida_condition_summary_fails(self):
        fixture = self.fixture(cohort="yachida")
        text = (fixture.seal / "SUCCESS").read_text(encoding="utf-8")
        lines = []
        for line in text.rstrip("\n").split("\n"):
            if line.startswith("conditions\t"):
                line = "conditions\tControl=2;CRC=2"
            lines.append(line)
        self.reseal_with(fixture, "\n".join(lines) + "\n")
        self.seal_fails(fixture, "conditions lacks Adenoma")


if __name__ == "__main__":
    unittest.main(verbosity=2)
