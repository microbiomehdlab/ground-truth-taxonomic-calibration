#!/usr/bin/env python3
"""Build the immutable, checksummed upstream evidence package.

The package proves which upstream inputs were complete before definitive
downstream analysis began. It is an evidence *index*: it copies compact sealed
tables and records checksums and stable identifiers for large external assets.
It never copies raw reads, host-cleaned reads, reference databases, container
images, scratch files or unrestricted scheduler logs.

Its only scientific inputs are the three `production_seal_v2` directories
produced by `seal_cohort_upstream.py` plus the sealed Yachida
assembly-sensitivity experiment. The three native cohort seals are provenance
sources already indexed by each v2 seal and are never parsed here.

Every authoritative path is an explicit named argument: nothing is discovered
by searching `work/`, taking the newest directory, or inferring a cohort from a
filename. Standard library only, with postponed annotation evaluation so the
cluster's older Python can run it.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime, timedelta, tzinfo
from pathlib import Path

PACKAGE_VERSION = "upstream_evidence_v1"
SEAL_CONTRACT = "upstream_seal_v2"

# --- the unified seal contract this package consumes ----------------------
SEAL_MEMBERS = (
    "sample_flow.tsv", "covariate_audit.tsv", "production_manifest.tsv",
    "production_manifest.independent.tsv", "source_seal_inventory.tsv",
    "SUCCESS", "production_seal.sha256",
)
SUCCESS_NAME = "SUCCESS"
CHECKSUM_NAME = "production_seal.sha256"
IN_PROGRESS_NAME = "AUDIT_IN_PROGRESS"
SEAL_SUCCESS_FIELDS = (
    "seal_contract", "cohort", "samples", "independent_samples",
    "baseline_profiles", "independent_profiles", "community_profiles", "status",
)
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
COVARIATE_AUDIT_FIELDS = ["field", "missing", "present", "distinct_nonmissing"]
CANONICAL_MANIFEST_COLUMNS = ("sample_id", "condition", "study",
                              "independent_subset")
CONDITIONS = ("Control", "Adenoma", "CRC")
DESIGNS = ("baseline", "independent", "community")
SAMPLE_BOOLEAN_FIELDS = (
    "verified_marker", "retained_output_receipt", "input_provenance",
    "sample_success", "baseline_profile_count", "independent_profile_count",
    "community_profile_count", "completion_table", "manifest_independent_flag",
)
SAMPLE_BLANK_FIELDS = ("failure_reasons", "receipt_error", "completion_error",
                       "provenance_error")
EXPECTED_PROVENANCE_MODE = {"yachida": "sealed_manifest_and_qc_receipt",
                            "feng": "state_file", "zeller": "state_file"}

# --- frozen scientific expectations ---------------------------------------
BASELINE_PER_SAMPLE = 1
COMMUNITY_PER_SAMPLE = 7
INDEPENDENT_PER_SUBSET_SAMPLE = 60
FROZEN_INDEPENDENT_SAMPLES = 30
FROZEN_SUBSET_BALANCE = {"Control": 10, "Adenoma": 10, "CRC": 10}
FROZEN_COHORTS = {
    "yachida": {"study": "YachidaS_2019", "samples": 201,
                "baseline": 201, "community": 1407, "independent": 1800},
    "feng": {"study": "Public_study__FengQ_2015", "samples": 154,
             "baseline": 154, "community": 1078, "independent": 1800},
    "zeller": {"study": "Public_study__ZellerG_2014", "samples": 156,
               "baseline": 156, "community": 1092, "independent": 1800},
}
COHORT_ORDER = ("yachida", "feng", "zeller")

# --- assembly-sensitivity expectations ------------------------------------
ASSEMBLY_CHECKSUM_NAME = "experiment_inputs_and_summary.sha256"
ASSEMBLY_EXPECTED = {"samples": 30, "assembly_arms": 2, "fractions_per_arm": 6,
                     "expected_profiles": 360, "observed_profiles": 360}
# The real experiment seal's SUCCESS schema, confirmed against the sealed
# 2026-09-01 run. No key override is needed for it; the flags exist only so a
# differently labelled historical seal is adapted explicitly rather than
# guessed at.
ASSEMBLY_DEFAULT_KEYS = {"samples": "samples", "assembly_arms": "assembly_arms",
                         "fractions_per_arm": "fractions_per_arm",
                         "expected_profiles": "expected_profiles",
                         "observed_profiles": "observed_profiles"}
ASSEMBLY_SUCCESS_FIELDS = ("experiment", "samples", "assembly_arms",
                           "fractions_per_arm", "expected_profiles",
                           "observed_profiles", "status")
ASSEMBLY_SIDECAR_NAME = "matched_seed_audit.tsv"
ASSEMBLY_SIDECAR_FIELDS = [
    "sample_id", "study", "original_label", "clean_label", "fraction",
    "original_seed", "clean_seed", "status",
]
ASSEMBLY_CLEAN_LABELS = {
    "Pana": "Pana_clean_GCA_000381525.1",
    "Pint": "Pint_clean_GCA_001953955.1",
}
ASSEMBLY_LABELS = set(ASSEMBLY_CLEAN_LABELS)
ASSEMBLY_FRACTIONS = {"0.0001", "0.0005", "0.001", "0.005", "0.01", "0.05"}

# --- provenance metadata contract -----------------------------------------
PROVENANCE_FIELDS = ["category", "asset_id", "name", "version_or_release",
                     "sha256", "availability", "notes"]
REQUIRED_PROVENANCE_CATEGORIES = (
    "source_code", "container_upstream", "container_analysis",
    "database_kraken2", "database_metaphlan", "host_reference",
    "spike_reference", "upstream_parameter",
)
PROVENANCE_REQUIRED_VALUES = ("asset_id", "name", "version_or_release",
                              "availability")
AUDIT_LEDGER_FIELDS = ["cohort", "audit_job_id", "audit_date_utc", "state",
                       "exit_code", "samples", "seal_sha256", "notes"]
LEDGER_REQUIRED_STATE = "COMPLETED"
LEDGER_REQUIRED_EXIT = "0:0"
# The unified upstream_seal_v2 audits for this definitive package. Frozen here
# and overridable per cohort, so an arbitrary job identifier is never silently
# accepted as the provenance of a released package.
EXPECTED_AUDIT_JOBS = {"yachida": "3097679", "feng": "3097680",
                       "zeller": "3097681"}

SEAL_INVENTORY_FIELDS = ["member", "sha256", "bytes", "status",
                         "source_audit_job"]
SHA256_PATTERN = re.compile(r"\A[0-9a-fA-F]{64}\Z")
JOB_PATTERN = re.compile(r"\A[0-9]+\Z")
UTC_STAMP_PATTERNS = (
    re.compile(r"\A\d{4}-\d{2}-\d{2}\Z"),
    re.compile(r"\A\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}(:\d{2})?Z\Z"),
)

# Provenance categories whose asset must carry a real content digest. An
# `upstream_parameter` row records a frozen numeric setting rather than a file,
# so its sha256 must be the documented sentinel instead of a blank cell.
PROVENANCE_DIGEST_CATEGORIES = (
    "source_code", "container_upstream", "container_analysis",
    "database_kraken2", "database_metaphlan", "host_reference",
    "spike_reference",
)
PARAMETER_SHA256_SENTINEL = "not_a_file"
REQUIRED_PARAMETERS = ("read_length", "bracken_threshold", "profiler_threads",
                       "spike_fractions", "pool_coverage", "seeds")
SPIKE_PANEL_COLUMNS = ("label", "taxon_name", "assembly", "weight")
EXPECTED_SPIKE_TARGETS = 10
ALIAS_COLUMNS = ("canonical", "alias", "tool")
ALIAS_OUTPUT_COLUMNS = ["canonical", "alias", "tool"]

PROJECTION_ROOT = "source_seal_projections"
PROJECTION_CHECKSUM_NAME = "projection.sha256"
PROJECTION_INVENTORY_NAME = "original_source_seal_inventory.tsv"
PROJECTION_README_NAME = "README.md"
EXACT_SEAL_ROOT = "source_seals"

# --- privacy screening -----------------------------------------------------
URL_PATTERN = re.compile(r"[A-Za-z][A-Za-z0-9+.\-]*://\S+")
ABSOLUTE_PATH_PATTERN = re.compile(r"(?:^|[\s,;=\"'])(/(?:[\w.\-]+/)+[\w.\-]*)")
CREDENTIAL_PATTERNS = (
    re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----"),
    re.compile(r"(?i)\b(?:password|passwd|secret_key|api[_-]?key|access[_-]?token"
               r"|aws_secret_access_key)\b\s*[:=]"),
    re.compile(r"(?i)\bssh-rsa\s+AAAA"),
)

MEDIA_TYPES = {".tsv": "text/tab-separated-values", ".md": "text/markdown",
               ".txt": "text/plain", ".pdf": "application/pdf",
               ".svg": "image/svg+xml", ".png": "image/png",
               ".sha256": "text/plain", "": "text/plain"}


class _Utc(tzinfo):
    """`datetime.UTC` is too new for the cluster interpreter."""

    def utcoffset(self, when):
        return timedelta(0)

    def tzname(self, when):
        return "UTC"

    def dst(self, when):
        return timedelta(0)


UTC = _Utc()


class PackageError(Exception):
    """A fail-closed contract violation."""


# ------------------------------------------------------------------- io
def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def nonempty_file(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def require_absolute(path, label: str) -> Path:
    expanded = Path(os.path.expanduser(str(path)))
    if not expanded.is_absolute():
        raise PackageError("%s must be an absolute path: %s" % (label, path))
    if ".." in expanded.parts:
        raise PackageError("%s must not contain '..': %s" % (label, path))
    return expanded.resolve()


def read_header(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        line = handle.readline()
    if not line.strip():
        raise PackageError("%s has no header" % path.name)
    header = line.rstrip("\r\n").split("\t")
    duplicates = sorted({name for name in header if header.count(name) > 1})
    if duplicates:
        raise PackageError(
            "%s has duplicate column(s): %s" % (path.name, ", ".join(duplicates)))
    return header


def read_table(path: Path, label: str):
    if not nonempty_file(path):
        raise PackageError("missing or empty %s: %s" % (label, path))
    header = read_header(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise PackageError("%s has no data rows: %s" % (label, path))
    return header, rows


def read_key_value(path: Path, label: str):
    if not nonempty_file(path):
        raise PackageError("missing or empty %s: %s" % (label, path))
    values = {}
    for number, line in enumerate(
            path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != 2 or not parts[0].strip():
            raise PackageError("malformed %s line %d: %r" % (label, number, line))
        if parts[0] in values:
            raise PackageError("duplicate %s field: %s" % (label, parts[0]))
        values[parts[0]] = parts[1]
    if not values:
        raise PackageError("%s is empty: %s" % (label, path))
    return values


def write_tsv(path: Path, fields, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(fields)]
    for row in rows:
        lines.append("\t".join(str(row.get(name, "")) for name in fields))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def integer(value, label: str) -> int:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        raise PackageError("%s is not an integer: %r" % (label, value))


# --------------------------------------------------------- checksum manifests
def verify_checksum_manifest(root: Path, checksum: Path, label: str):
    """Verify every member a seal's checksum manifest lists.

    Absolute and parent-traversal member names are rejected before any file is
    opened, so a malformed manifest cannot reach outside its own seal.
    """
    if not nonempty_file(checksum):
        raise PackageError("%s has no nonempty %s" % (label, checksum.name))
    verified = {}
    for number, line in enumerate(
            checksum.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) != 2:
            raise PackageError(
                "%s checksum line %d is malformed: %r" % (label, number, line))
        recorded, member = parts
        if member in verified:
            raise PackageError("%s lists %s twice" % (label, member))
        candidate = Path(member)
        if candidate.is_absolute() or ".." in candidate.parts:
            raise PackageError("%s has an unsafe member path: %s" % (label, member))
        target = root / member
        if not nonempty_file(target):
            raise PackageError("%s member missing or empty: %s" % (label, member))
        actual = digest(target)
        if actual != recorded:
            raise PackageError("%s checksum mismatch for %s" % (label, member))
        verified[member] = (actual, target.stat().st_size)
    if not verified:
        raise PackageError("%s checksum manifest lists nothing" % label)
    return verified


# --------------------------------------------------------------- seal reading
def load_cohort_seal(cohort: str, seal_root: Path, manifest: Path,
                     independent_manifest: Path):
    """Verify one production_seal_v2 directory and return its evidence."""
    label = "%s production_seal_v2" % cohort
    if not seal_root.is_dir():
        raise PackageError("%s is not a directory: %s" % (label, seal_root))
    if (seal_root / IN_PROGRESS_NAME).exists():
        raise PackageError(
            "%s is mid-audit (%s present)" % (label, IN_PROGRESS_NAME))
    present = sorted(item.name for item in seal_root.iterdir() if item.is_file())
    missing = [name for name in SEAL_MEMBERS if name not in present]
    if missing:
        raise PackageError("%s lacks member(s): %s" % (label, ", ".join(missing)))
    # Exactly the seven documented members, no more.
    unexpected = [name for name in present if name not in SEAL_MEMBERS]
    if unexpected:
        raise PackageError(
            "%s holds unexpected file(s): %s" % (label, ", ".join(unexpected)))
    directories = sorted(item.name for item in seal_root.iterdir()
                         if not item.is_file())
    if directories:
        raise PackageError(
            "%s holds unexpected subdirector(ies): %s"
            % (label, ", ".join(directories)))

    verified = verify_checksum_manifest(
        seal_root, seal_root / CHECKSUM_NAME, label)
    covered = set(SEAL_MEMBERS) - {CHECKSUM_NAME}
    # The checksum manifest must cover exactly the other six members.
    unlisted = sorted(covered - set(verified))
    if unlisted:
        raise PackageError(
            "%s checksum manifest omits %s" % (label, ", ".join(unlisted)))
    extra = sorted(set(verified) - covered)
    if extra:
        raise PackageError(
            "%s checksum manifest lists unexpected member(s): %s"
            % (label, ", ".join(extra)))

    success = read_key_value(seal_root / SUCCESS_NAME, "%s SUCCESS" % label)
    missing_keys = [name for name in SEAL_SUCCESS_FIELDS if name not in success]
    if missing_keys:
        raise PackageError(
            "%s SUCCESS lacks %s" % (label, ", ".join(missing_keys)))
    if success["seal_contract"] != SEAL_CONTRACT:
        raise PackageError(
            "%s SUCCESS declares contract %r, not %r"
            % (label, success["seal_contract"], SEAL_CONTRACT))
    if success["cohort"] != cohort:
        raise PackageError(
            "%s SUCCESS declares cohort %r" % (label, success["cohort"]))
    if success["status"] != "PASS":
        raise PackageError(
            "%s SUCCESS status is %r, not PASS" % (label, success["status"]))

    frozen = FROZEN_COHORTS[cohort]
    for key, expected in (("samples", frozen["samples"]),
                          ("independent_samples", FROZEN_INDEPENDENT_SAMPLES),
                          ("baseline_profiles", frozen["baseline"]),
                          ("independent_profiles", frozen["independent"]),
                          ("community_profiles", frozen["community"])):
        observed = integer(success[key], "%s SUCCESS %s" % (label, key))
        if observed != expected:
            raise PackageError(
                "%s SUCCESS reports %s=%d; the frozen design requires %d"
                % (label, key, observed, expected))

    flow_header, flow = read_table(seal_root / "sample_flow.tsv",
                                   "%s sample_flow.tsv" % label)
    if flow_header != SAMPLE_FLOW_FIELDS:
        raise PackageError(
            "%s sample_flow.tsv schema drifted from the unified contract"
            % label)
    covariate_header, covariates = read_table(
        seal_root / "covariate_audit.tsv", "%s covariate_audit.tsv" % label)
    if covariate_header != COVARIATE_AUDIT_FIELDS:
        raise PackageError(
            "%s covariate_audit.tsv schema drifted from the unified contract"
            % label)

    canonical_header, canonical = read_table(
        seal_root / "production_manifest.tsv", "%s production manifest" % label)
    independent_header, independent_rows = read_table(
        seal_root / "production_manifest.independent.tsv",
        "%s independent manifest" % label)
    # The canonical four must lead, in order, in both manifests.
    for name, header in (("production_manifest.tsv", canonical_header),
                         ("production_manifest.independent.tsv",
                          independent_header)):
        if tuple(header[:4]) != CANONICAL_MANIFEST_COLUMNS:
            raise PackageError(
                "%s %s must begin with %s; observed %s"
                % (label, name, ", ".join(CANONICAL_MANIFEST_COLUMNS),
                   ", ".join(header[:4])))

    # Sample identity, compared as sets rather than counts.
    flow_ids = [row["sample_id"] for row in flow]
    if len(set(flow_ids)) != len(flow_ids):
        raise PackageError("%s sample_flow.tsv has duplicate sample ids" % label)
    canonical_ids = {row["sample_id"] for row in canonical}
    if canonical_ids != set(flow_ids):
        raise PackageError(
            "%s sealed manifest and sample flow cover different samples" % label)
    frozen_ids = frozen_manifest_ids(manifest, "%s frozen manifest" % label)
    if frozen_ids != set(flow_ids):
        only_manifest = sorted(frozen_ids - set(flow_ids))[:3]
        only_seal = sorted(set(flow_ids) - frozen_ids)[:3]
        raise PackageError(
            "%s frozen manifest and seal cover different samples "
            "(manifest only: %s; seal only: %s)"
            % (label, only_manifest or "none", only_seal or "none"))
    if len(flow_ids) != frozen["samples"]:
        raise PackageError(
            "%s has %d samples; the frozen design requires %d"
            % (label, len(flow_ids), frozen["samples"]))

    # Per-sample cross-check between the ledger and the canonical manifest,
    # then the full frozen per-sample topology. Aggregate totals cannot mask a
    # sample-level error, because every component is checked on every row.
    canonical_by_id = {row["sample_id"]: row for row in canonical}
    declared_independent = {row["sample_id"] for row in canonical
                            if row["independent_subset"] == "1"}
    for row in flow:
        sample = row["sample_id"]
        source = canonical_by_id[sample]
        for column in CANONICAL_MANIFEST_COLUMNS:
            if row[column] != source[column]:
                raise PackageError(
                    "%s sample %s has %s=%r in sample_flow.tsv but %r in the "
                    "canonical manifest"
                    % (label, sample, column, row[column], source[column]))
        in_subset = sample in declared_independent
        expected = {
            "baseline": BASELINE_PER_SAMPLE,
            "community": COMMUNITY_PER_SAMPLE,
            "independent": (INDEPENDENT_PER_SUBSET_SAMPLE if in_subset else 0),
        }
        for design in DESIGNS:
            declared = integer(row["expected_%s_profiles" % design],
                               "%s %s expected_%s_profiles"
                               % (label, sample, design))
            observed = integer(row["observed_%s_profiles" % design],
                               "%s %s observed_%s_profiles"
                               % (label, sample, design))
            if declared != expected[design]:
                raise PackageError(
                    "%s sample %s expects %d %s profiles; the frozen design "
                    "requires %d" % (label, sample, declared, design,
                                     expected[design]))
            if observed != declared:
                raise PackageError(
                    "%s sample %s observed %d %s profiles but expected %d"
                    % (label, sample, observed, design, declared))
        component_total = sum(expected.values())
        for column in ("expected_profiles", "observed_profiles"):
            total = integer(row[column], "%s %s %s" % (label, sample, column))
            if total != component_total:
                raise PackageError(
                    "%s sample %s has %s=%d but its components sum to %d"
                    % (label, sample, column, total, component_total))
        for column in SAMPLE_BOOLEAN_FIELDS:
            if row[column] != "1":
                raise PackageError(
                    "%s sample %s has %s=%r, not 1"
                    % (label, sample, column, row[column]))
        for column in SAMPLE_BLANK_FIELDS:
            if row[column].strip():
                raise PackageError(
                    "%s sample %s has a nonblank %s: %r"
                    % (label, sample, column, row[column]))
        if row["status"] != "PASS":
            raise PackageError(
                "%s sample %s is not PASS (%s)"
                % (label, sample, row["failure_reasons"] or "no reason"))
    conditions = Counter(row["condition"] for row in canonical)
    unexpected = sorted(set(conditions) - set(CONDITIONS))
    if unexpected:
        raise PackageError(
            "%s manifest uses unexpected condition label(s): %s"
            % (label, ", ".join(unexpected)))

    # Independent subset: exact size, exact nesting, exact balance.
    independent_ids = {row["sample_id"] for row in independent_rows}
    if len(independent_ids) != len(independent_rows):
        raise PackageError("%s independent manifest has duplicate ids" % label)
    if len(independent_ids) != FROZEN_INDEPENDENT_SAMPLES:
        raise PackageError(
            "%s independent subset has %d samples; %d are required"
            % (label, len(independent_ids), FROZEN_INDEPENDENT_SAMPLES))
    if not independent_ids <= set(flow_ids):
        raise PackageError(
            "%s independent subset is not nested in the cohort" % label)
    # The subset must be exactly the rows the manifest flags, not merely a
    # nested set of the right size.
    if independent_ids != declared_independent:
        only_subset = sorted(independent_ids - declared_independent)[:3]
        only_flagged = sorted(declared_independent - independent_ids)[:3]
        raise PackageError(
            "%s independent manifest does not match the rows flagged "
            "independent_subset=1 (subset only: %s; flagged only: %s)"
            % (label, only_subset or "none", only_flagged or "none"))
    frozen_independent = frozen_manifest_ids(
        independent_manifest, "%s frozen independent manifest" % label)
    if frozen_independent != independent_ids:
        raise PackageError(
            "%s frozen independent manifest and seal subset differ" % label)
    balance = Counter(row["condition"] for row in independent_rows)
    if dict(balance) != FROZEN_SUBSET_BALANCE:
        raise PackageError(
            "%s independent subset is not balanced %s; observed %s"
            % (label, FROZEN_SUBSET_BALANCE, dict(balance)))

    # Topology, summed from the per-sample ledger and cross-checked.
    totals = {"baseline": 0, "independent": 0, "community": 0}
    expected_totals = {"baseline": 0, "independent": 0, "community": 0}
    for row in flow:
        for design in DESIGNS:
            expected_totals[design] += integer(
                row["expected_%s_profiles" % design],
                "%s expected_%s_profiles" % (label, design))
            totals[design] += integer(
                row["observed_%s_profiles" % design],
                "%s observed_%s_profiles" % (label, design))
    for design in DESIGNS:
        if totals[design] != frozen[design]:
            raise PackageError(
                "%s observed %d %s profiles; the frozen design requires %d"
                % (label, totals[design], design, frozen[design]))
        if expected_totals[design] != frozen[design]:
            raise PackageError(
                "%s expected %d %s profiles; the frozen design requires %d"
                % (label, expected_totals[design], design, frozen[design]))

    studies = {row["study"] for row in canonical}
    if studies != {frozen["study"]}:
        raise PackageError(
            "%s manifest study labels are %s; %r is required"
            % (label, sorted(studies), frozen["study"]))

    inventory_header, inventory = read_table(
        seal_root / "source_seal_inventory.tsv",
        "%s source_seal_inventory.tsv" % label)
    if inventory_header != SEAL_INVENTORY_FIELDS:
        raise PackageError(
            "%s source_seal_inventory.tsv schema must be exactly: %s"
            % (label, ", ".join(SEAL_INVENTORY_FIELDS)))
    inventory_members = set()
    inventory_jobs = set()
    for number, row in enumerate(inventory, start=2):
        if not (row["member"] or "").strip():
            raise PackageError(
                "%s source_seal_inventory.tsv line %d has no member"
                % (label, number))
        if not SHA256_PATTERN.match((row["sha256"] or "").strip()):
            raise PackageError(
                "%s source_seal_inventory.tsv line %d has an invalid sha256"
                % (label, number))
        size = integer(row["bytes"], "%s inventory line %d bytes"
                       % (label, number))
        if size <= 0:
            raise PackageError(
                "%s source_seal_inventory.tsv line %d has a non-positive byte "
                "count" % (label, number))
        if row["status"] != "VERIFIED":
            raise PackageError(
                "%s source_seal_inventory.tsv line %d is %r, not VERIFIED"
                % (label, number, row["status"]))
        member = row["member"].strip()
        if member in inventory_members:
            raise PackageError(
                "%s source_seal_inventory.tsv repeats member %r"
                % (label, member))
        inventory_members.add(member)
        job = (row["source_audit_job"] or "").strip()
        # This field belongs to the native/source seal, not to the later
        # unified-v2 audit recorded in audit_ledger.tsv. Historical Yachida
        # has no native Slurm audit identifier, so a consistently blank value
        # is legitimate; otherwise the value must be numeric.
        if job and not JOB_PATTERN.match(job):
            raise PackageError(
                "%s source_seal_inventory.tsv line %d has an invalid "
                "source_audit_job %r" % (label, number, job))
        inventory_jobs.add(job)
    if len(inventory_jobs) != 1:
        raise PackageError(
            "%s source_seal_inventory.tsv reports multiple audit jobs: %s"
            % (label, ", ".join(sorted(inventory_jobs))))

    modes = {row["input_provenance_mode"] for row in flow}
    if len(modes) != 1:
        raise PackageError(
            "%s reports multiple input-provenance modes: %s"
            % (label, ", ".join(sorted(modes))))
    expected_mode = EXPECTED_PROVENANCE_MODE[cohort]
    if sorted(modes)[0] != expected_mode:
        raise PackageError(
            "%s reports input-provenance mode %r; %r is required for this "
            "cohort" % (label, sorted(modes)[0], expected_mode))

    return {
        "cohort": cohort, "root": seal_root, "success": success, "flow": flow,
        "covariates": covariates, "canonical": canonical,
        "independent": independent_rows, "inventory": inventory,
        "verified": verified, "conditions": conditions,
        "balance": balance, "totals": totals,
        "expected_totals": expected_totals,
        "input_provenance_mode": sorted(modes)[0],
        "native_source_audit_job": next(iter(inventory_jobs)),
        "seal_sha256": digest(seal_root / CHECKSUM_NAME),
    }


def frozen_manifest_ids(path: Path, label: str):
    header, rows = read_table(path, label)
    if "sample_id" not in header:
        raise PackageError("%s lacks a sample_id column" % label)
    identifiers = [row["sample_id"].strip() for row in rows]
    if any(not value for value in identifiers):
        raise PackageError("%s has a blank sample_id" % label)
    if len(set(identifiers)) != len(identifiers):
        raise PackageError("%s has duplicate sample identifiers" % label)
    return set(identifiers)


def validate_assembly_sidecar(path: Path):
    """Validate the historical matched-seed audit not covered by the seal manifest."""
    header, rows = read_table(path, "assembly matched-seed audit")
    if header != ASSEMBLY_SIDECAR_FIELDS:
        raise PackageError(
            "assembly matched-seed audit schema must be exactly: %s"
            % ", ".join(ASSEMBLY_SIDECAR_FIELDS))
    if len(rows) != 126:
        raise PackageError(
            "assembly matched-seed audit has %d rows; the frozen audit has 126"
            % len(rows))
    samples = set()
    observed = set()
    for number, row in enumerate(rows, start=2):
        sample = row["sample_id"].strip()
        label = row["original_label"].strip()
        fraction = row["fraction"].strip()
        key = (sample, label, fraction)
        if not sample or key in observed:
            raise PackageError(
                "assembly matched-seed audit line %d has a blank sample or duplicate design key"
                % number)
        observed.add(key)
        samples.add(sample)
        if row["study"].strip() != "YachidaS_2019":
            raise PackageError(
                "assembly matched-seed audit line %d has the wrong study" % number)
        if label not in ASSEMBLY_LABELS or fraction not in ASSEMBLY_FRACTIONS:
            raise PackageError(
                "assembly matched-seed audit line %d has an unexpected label or fraction"
                % number)
        if row["clean_label"].strip() != ASSEMBLY_CLEAN_LABELS.get(label):
            raise PackageError(
                "assembly matched-seed audit line %d has the wrong clean label"
                % number)
        original_seed = row["original_seed"].strip()
        clean_seed = row["clean_seed"].strip()
        if not original_seed.isdigit() or original_seed != clean_seed:
            raise PackageError(
                "assembly matched-seed audit line %d does not preserve the seed"
                % number)
        if row["status"].strip() != "PASS":
            raise PackageError(
                "assembly matched-seed audit line %d is not PASS" % number)
    return digest(path), path.stat().st_size


def load_assembly_sensitivity(seal_root: Path, keys):
    label = "assembly-sensitivity experiment seal"
    if not seal_root.is_dir():
        raise PackageError("%s is not a directory: %s" % (label, seal_root))
    checksum = seal_root / ASSEMBLY_CHECKSUM_NAME
    if not nonempty_file(checksum):
        raise PackageError("%s lacks %s" % (label, ASSEMBLY_CHECKSUM_NAME))
    if (seal_root / IN_PROGRESS_NAME).exists():
        raise PackageError("%s is mid-audit (%s present)"
                           % (label, IN_PROGRESS_NAME))
    verified = verify_checksum_manifest(seal_root, checksum, label)
    if SUCCESS_NAME not in verified:
        raise PackageError(
            "%s does not checksum its own SUCCESS" % label)
    # The historical directory contains one post-seal matched-seed audit. It is
    # not silently treated as checksummed: require that exact sidecar, validate
    # its full design and assertions, and reject every other unlisted member.
    present = sorted(item.name for item in seal_root.iterdir() if item.is_file())
    unlisted = [name for name in present
                if name != ASSEMBLY_CHECKSUM_NAME and name not in verified]
    if unlisted != [ASSEMBLY_SIDECAR_NAME]:
        raise PackageError(
            "%s must have exactly the historical unsealed sidecar %s; observed unlisted members: %s"
            % (label, ASSEMBLY_SIDECAR_NAME, ", ".join(unlisted) or "none"))
    sidecar = validate_assembly_sidecar(seal_root / ASSEMBLY_SIDECAR_NAME)
    success = read_key_value(seal_root / SUCCESS_NAME, "%s SUCCESS" % label)
    missing_fields = [name for name in ASSEMBLY_SUCCESS_FIELDS
                      if name not in success]
    if missing_fields:
        raise PackageError(
            "%s SUCCESS lacks %s" % (label, ", ".join(missing_fields)))
    if success.get("status") != "PASS":
        raise PackageError(
            "%s SUCCESS status is %r, not PASS" % (label, success.get("status")))
    observed = {}
    for name, expected in sorted(ASSEMBLY_EXPECTED.items()):
        key = keys[name]
        if key not in success:
            raise PackageError(
                "%s SUCCESS lacks the %r field required for %s"
                % (label, key, name))
        value = integer(success[key], "%s %s" % (label, key))
        if value != expected:
            raise PackageError(
                "%s reports %s=%d; the frozen experiment requires %d"
                % (label, name, value, expected))
        observed[name] = value
    return {"root": seal_root, "success": success, "verified": verified,
            "sidecar": {ASSEMBLY_SIDECAR_NAME: sidecar},
            "observed": observed,
            "seal_sha256": digest(checksum)}


# ------------------------------------------------------- provenance metadata
def load_provenance_metadata(path: Path, commit: str, spike_targets,
                             analysis_image_sha256):
    """Validate the non-sensitive provenance metadata contract.

    File-backed assets must carry a real 64-character digest. An
    `upstream_parameter` row records a frozen numeric setting rather than a
    file, so the documented policy is that its `sha256` is the sentinel
    `not_a_file` -- explicit, not a silently accepted blank.
    """
    header, rows = read_table(path, "provenance metadata")
    if header != PROVENANCE_FIELDS:
        raise PackageError(
            "provenance metadata schema must be exactly: %s"
            % ", ".join(PROVENANCE_FIELDS))
    categories = set()
    asset_ids = {}
    spike_counts = Counter()
    parameter_counts = Counter()
    analysis_rows = []
    source_rows = []
    for number, row in enumerate(rows, start=2):
        category = (row.get("category") or "").strip()
        if not category:
            raise PackageError("provenance metadata line %d has no category" % number)
        if category not in REQUIRED_PROVENANCE_CATEGORIES:
            raise PackageError(
                "provenance metadata line %d has the undocumented category %r; "
                "document it in the specification before use"
                % (number, category))
        categories.add(category)
        for field in PROVENANCE_REQUIRED_VALUES:
            if not (row.get(field) or "").strip():
                raise PackageError(
                    "provenance metadata line %d (%s) has a blank %s"
                    % (number, category, field))
        asset_id = row["asset_id"].strip()
        if asset_id in asset_ids:
            raise PackageError(
                "provenance metadata reuses asset_id %r on lines %d and %d"
                % (asset_id, asset_ids[asset_id], number))
        asset_ids[asset_id] = number
        checksum = (row.get("sha256") or "").strip()
        if category in PROVENANCE_DIGEST_CATEGORIES:
            if not SHA256_PATTERN.match(checksum):
                raise PackageError(
                    "provenance metadata line %d (%s) needs a 64-character "
                    "sha256; observed %r" % (number, category, checksum))
        elif checksum != PARAMETER_SHA256_SENTINEL:
            raise PackageError(
                "provenance metadata line %d (%s) must record sha256 %r "
                "because it is not a file; observed %r"
                % (number, category, PARAMETER_SHA256_SENTINEL, checksum))
        if category == "spike_reference":
            matches = [target for target in spike_targets
                       if target in (asset_id, row["name"].strip())]
            if len(matches) != 1:
                raise PackageError(
                    "provenance metadata line %d must identify exactly one "
                    "implanted target; observed asset_id=%r name=%r"
                    % (number, asset_id, row["name"]))
            spike_counts[matches[0]] += 1
        if category == "upstream_parameter":
            matches = [name for name in REQUIRED_PARAMETERS
                       if name in (asset_id, row["name"].strip())]
            if len(matches) != 1:
                raise PackageError(
                    "provenance metadata line %d must identify exactly one "
                    "frozen upstream parameter; observed asset_id=%r name=%r"
                    % (number, asset_id, row["name"]))
            parameter_counts[matches[0]] += 1
        if category == "container_analysis":
            analysis_rows.append((number, checksum))
        if category == "source_code":
            source_rows.append((number, row))

    missing = [name for name in REQUIRED_PROVENANCE_CATEGORIES
               if name not in categories]
    if missing:
        raise PackageError(
            "provenance metadata lacks required categor(ies): %s"
            % ", ".join(missing))

    # source_code must identify the commit this package was built from.
    if not any(commit in "\t".join(row.values()) for _, row in source_rows):
        raise PackageError(
            "no source_code row identifies the supplied source commit %s"
            % commit)

    # One spike_reference per implanted target in the validated panel.
    wrong_spikes = [target for target in spike_targets
                    if spike_counts[target] != 1]
    if wrong_spikes:
        raise PackageError(
            "provenance metadata must contain exactly one spike_reference "
            "for: %s" % ", ".join(sorted(wrong_spikes)))

    wrong_parameters = [name for name in REQUIRED_PARAMETERS
                        if parameter_counts[name] != 1]
    if wrong_parameters:
        raise PackageError(
            "provenance metadata must contain exactly one frozen "
            "upstream_parameter row for: %s"
            % ", ".join(wrong_parameters))

    if analysis_image_sha256:
        mismatched = [number for number, checksum in analysis_rows
                      if checksum.lower() != analysis_image_sha256.lower()]
        if mismatched or not analysis_rows:
            raise PackageError(
                "the container_analysis sha256 does not match the analysis "
                "image checksum %s supplied by the runner"
                % analysis_image_sha256)
    return rows


def load_audit_ledger(path: Path, cohorts, expected_jobs, seal_checksums):
    """Validate the manually supplied operational evidence semantically.

    Exactly one row per cohort, each recording a completed audit of *this*
    cohort's actual seal. The expected job identifiers are frozen, so an
    arbitrary identifier can never become the released provenance.
    """
    header, rows = read_table(path, "audit ledger")
    if header != AUDIT_LEDGER_FIELDS:
        raise PackageError(
            "audit ledger schema must be exactly: %s"
            % ", ".join(AUDIT_LEDGER_FIELDS))
    listed = {}
    for number, row in enumerate(rows, start=2):
        cohort = (row.get("cohort") or "").strip()
        if not cohort:
            raise PackageError("audit ledger line %d has no cohort" % number)
        if cohort not in cohorts:
            raise PackageError(
                "audit ledger line %d names the unknown cohort %r"
                % (number, cohort))
        if cohort in listed:
            raise PackageError(
                "audit ledger has more than one row for %s" % cohort)
        for field in ("audit_job_id", "audit_date_utc", "state", "exit_code",
                      "samples", "seal_sha256"):
            if not (row.get(field) or "").strip():
                raise PackageError(
                    "audit ledger line %d (%s) has a blank %s"
                    % (number, cohort, field))
        job = row["audit_job_id"].strip()
        if not JOB_PATTERN.match(job):
            raise PackageError(
                "audit ledger line %d (%s) has a non-numeric audit_job_id %r"
                % (number, cohort, job))
        if job != expected_jobs[cohort]:
            raise PackageError(
                "audit ledger records job %s for %s; this package expects %s"
                % (job, cohort, expected_jobs[cohort]))
        stamp = row["audit_date_utc"].strip()
        if not any(pattern.match(stamp) for pattern in UTC_STAMP_PATTERNS):
            raise PackageError(
                "audit ledger line %d (%s) has an invalid UTC audit_date_utc "
                "%r; use YYYY-MM-DD or YYYY-MM-DDThh:mm:ssZ"
                % (number, cohort, stamp))
        if row["state"].strip() != LEDGER_REQUIRED_STATE:
            raise PackageError(
                "audit ledger line %d (%s) has state %r, not %s"
                % (number, cohort, row["state"], LEDGER_REQUIRED_STATE))
        if row["exit_code"].strip() != LEDGER_REQUIRED_EXIT:
            raise PackageError(
                "audit ledger line %d (%s) has exit_code %r, not %s"
                % (number, cohort, row["exit_code"], LEDGER_REQUIRED_EXIT))
        samples = integer(row["samples"], "audit ledger %s samples" % cohort)
        if samples != FROZEN_COHORTS[cohort]["samples"]:
            raise PackageError(
                "audit ledger records %d samples for %s; the frozen design has "
                "%d" % (samples, cohort, FROZEN_COHORTS[cohort]["samples"]))
        recorded = row["seal_sha256"].strip()
        if not SHA256_PATTERN.match(recorded):
            raise PackageError(
                "audit ledger line %d (%s) has an invalid seal_sha256"
                % (number, cohort))
        if recorded.lower() != seal_checksums[cohort].lower():
            raise PackageError(
                "audit ledger seal_sha256 for %s does not match that cohort's "
                "production_seal.sha256" % cohort)
        listed[cohort] = row
    missing = [cohort for cohort in cohorts if cohort not in listed]
    if missing:
        raise PackageError(
            "audit ledger lacks row(s) for: %s" % ", ".join(sorted(missing)))
    return rows


def load_spike_panel(path: Path):
    """Validate the frozen implanted panel and return its target labels."""
    header, rows = read_table(path, "spike panel")
    missing = [name for name in SPIKE_PANEL_COLUMNS if name not in header]
    if missing:
        raise PackageError(
            "spike panel lacks column(s): %s" % ", ".join(missing))
    labels = []
    for number, row in enumerate(rows, start=2):
        for field in SPIKE_PANEL_COLUMNS:
            if not (row.get(field) or "").strip():
                raise PackageError(
                    "spike panel line %d has a blank %s" % (number, field))
        label = row["label"].strip()
        if label in labels:
            raise PackageError("spike panel repeats label %r" % label)
        labels.append(label)
        raw = (row.get("weight") or "").strip()
        try:
            weight = float(raw)
        except ValueError:
            raise PackageError(
                "spike panel line %d has a nonnumeric weight %r"
                % (number, raw))
        if weight <= 0:
            raise PackageError(
                "spike panel line %d has a non-positive weight" % number)
    if len(labels) != EXPECTED_SPIKE_TARGETS:
        raise PackageError(
            "spike panel has %d targets; the frozen panel has %d"
            % (len(labels), EXPECTED_SPIKE_TARGETS))
    return labels, rows


def load_taxon_aliases(path: Path):
    """Read the alias table by its real delimiter and normalise it to TSV.

    The frozen alias file is comma-delimited, so copying it to a `.tsv` name
    would produce a falsely labelled file.
    """
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        header = list(reader.fieldnames or [])
        rows = [row for row in reader]
    if not header:
        raise PackageError("taxon alias table has no header: %s" % path)
    if not rows:
        raise PackageError("taxon alias table has no data rows: %s" % path)
    missing = [name for name in ALIAS_COLUMNS if name not in header]
    if missing:
        raise PackageError(
            "taxon alias table lacks column(s): %s" % ", ".join(missing))
    seen = set()
    for number, row in enumerate(rows, start=2):
        for field in ALIAS_COLUMNS:
            if not (row.get(field) or "").strip():
                raise PackageError(
                    "taxon alias line %d has a blank %s" % (number, field))
        key = (row["canonical"].strip(), row["tool"].strip())
        if key in seen:
            raise PackageError(
                "taxon alias table repeats %s for tool %s" % key)
        seen.add(key)
    return [dict((name, (row.get(name) or "").strip())
                 for name in ALIAS_OUTPUT_COLUMNS) for row in rows]


def scan_for_credentials(rows, label: str) -> None:
    for number, row in enumerate(rows, start=2):
        joined = "\t".join(str(value) for value in row.values())
        for pattern in CREDENTIAL_PATTERNS:
            if pattern.search(joined):
                raise PackageError(
                    "%s line %d looks like a credential or private key" % (label, number))


def absolute_paths_in(text: str):
    stripped = URL_PATTERN.sub(" ", text)
    return [match.group(1) for match in ABSOLUTE_PATH_PATTERN.finditer(stripped)]


def assert_release_safe(root: Path) -> None:
    """No release-facing table may carry a cluster-absolute path."""
    offenders = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.suffix not in (".tsv", ".md", ".txt"):
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        found = absolute_paths_in(text)
        if found:
            offenders.append("%s (%s)" % (path.name, found[0]))
    if offenders:
        raise PackageError(
            "release-facing output(s) contain absolute paths: "
            + "; ".join(sorted(offenders)))


# ----------------------------------------------------------------- git commit
def seal_inventory_snapshot(roots):
    """Full path -> digest for every file under every authoritative seal."""
    snapshot = {}
    for root in roots:
        for path in sorted(Path(root).rglob("*")):
            if path.is_file():
                snapshot[str(path)] = digest(path)
    return snapshot


def assert_sources_unchanged(roots, before) -> None:
    """Prove no authoritative seal changed while the package was built."""
    for root in roots:
        if (Path(root) / IN_PROGRESS_NAME).exists():
            raise PackageError(
                "%s became mid-audit during the build" % Path(root).name)
    after = seal_inventory_snapshot(roots)
    added = sorted(set(after) - set(before))
    removed = sorted(set(before) - set(after))
    changed = sorted(name for name in set(before) & set(after)
                     if before[name] != after[name])
    if added or removed or changed:
        detail = []
        if added:
            detail.append("added %s" % Path(added[0]).name)
        if removed:
            detail.append("removed %s" % Path(removed[0]).name)
        if changed:
            detail.append("changed %s" % Path(changed[0]).name)
        raise PackageError(
            "a source seal changed during the build: %s" % "; ".join(detail))


def resolve_commit(repository_root, supplied):
    if supplied:
        commit = supplied.strip()
    else:
        if repository_root is None:
            raise PackageError(
                "give either --source-commit or --repository-root")
        try:
            commit = subprocess.check_output(
                ["git", "-C", str(repository_root), "rev-parse", "HEAD"],
                stderr=subprocess.STDOUT).decode("utf-8").strip()
        except (OSError, subprocess.CalledProcessError) as error:
            raise PackageError("cannot resolve the source commit: %s" % error)
    if not re.match(r"\A[0-9a-f]{40}\Z", commit):
        raise PackageError(
            "the source commit must be a full 40-character hex sha: %r" % commit)
    return commit


# --------------------------------------------------------------- table builds
def build_tables(cohorts, assembly, provenance_rows, ledger_rows):
    """Every release table, derived once from the verified seals."""
    sample_flow = []
    condition_counts = []
    completeness = []
    balance_rows = []
    covariate_rows = []
    for cohort in COHORT_ORDER:
        data = cohorts[cohort]
        frozen = FROZEN_COHORTS[cohort]
        independent_by_condition = Counter(
            row["condition"] for row in data["independent"])
        total_independent = 0
        for condition in CONDITIONS:
            manifest_samples = data["conditions"].get(condition, 0)
            independent_samples = independent_by_condition.get(condition, 0)
            total_independent += independent_samples
            sample_flow.append({
                "cohort": cohort, "condition": condition,
                "manifest_samples": manifest_samples,
                "sealed_samples": manifest_samples,
                "independent_samples": independent_samples,
                "status": "PASS",
            })
            condition_counts.append({
                "cohort": cohort, "condition": condition,
                "samples": manifest_samples,
                "fraction_of_cohort": render_fraction(
                    manifest_samples, frozen["samples"]),
            })
            balance_rows.append({
                "cohort": cohort, "condition": condition,
                "expected_samples": FROZEN_SUBSET_BALANCE[condition],
                "observed_samples": independent_samples,
                "status": ("PASS" if independent_samples
                           == FROZEN_SUBSET_BALANCE[condition] else "FAIL"),
            })
        sample_flow.append({
            "cohort": cohort, "condition": "TOTAL",
            "manifest_samples": frozen["samples"],
            "sealed_samples": len(data["flow"]),
            "independent_samples": total_independent, "status": "PASS",
        })
        for design in DESIGNS:
            expected = data["expected_totals"][design]
            observed = data["totals"][design]
            completeness.append({
                "cohort": cohort, "design": design,
                "expected_profiles": expected, "observed_profiles": observed,
                "completion_fraction": render_fraction(observed, expected),
                "status": "PASS" if observed == expected else "FAIL",
            })
        for row in data["covariates"]:
            covariate_rows.append({
                "cohort": cohort, "field": row["field"],
                "missing": row["missing"], "present": row["present"],
                "distinct_nonmissing": row["distinct_nonmissing"],
            })

    assembly_rows = [{
        "experiment": "yachida_assembly_choice_sensitivity",
        "samples": assembly["observed"]["samples"],
        "assembly_arms": assembly["observed"]["assembly_arms"],
        "fractions_per_arm": assembly["observed"]["fractions_per_arm"],
        "expected_profiles": assembly["observed"]["expected_profiles"],
        "observed_profiles": assembly["observed"]["observed_profiles"],
        "status": "PASS",
    }]

    registry = []
    for cohort in COHORT_ORDER:
        data = cohorts[cohort]
        frozen = FROZEN_COHORTS[cohort]
        jobs = sorted({row["audit_job_id"] for row in ledger_rows
                       if row["cohort"] == cohort})
        registry.append({
            "cohort": cohort, "study": frozen["study"],
            "samples": frozen["samples"],
            "independent_samples": FROZEN_INDEPENDENT_SAMPLES,
            "baseline_profiles": data["totals"]["baseline"],
            "community_profiles": data["totals"]["community"],
            "independent_profiles": data["totals"]["independent"],
            "seal_contract": SEAL_CONTRACT,
            "seal_status": data["success"]["status"],
            "input_provenance_mode": data["input_provenance_mode"],
            "source_audit_job": ";".join(jobs),
        })

    inventory = []
    for cohort in COHORT_ORDER:
        data = cohorts[cohort]
        for member, (sha, size) in sorted(data["verified"].items()):
            inventory.append({
                "component": cohort, "logical_file": member, "sha256": sha,
                "bytes": size, "source_seal_checksum": data["seal_sha256"],
                "status": "VERIFIED",
            })
        for row in data["inventory"]:
            inventory.append({
                "component": "%s_native_source_seal" % cohort,
                "logical_file": row["member"], "sha256": row["sha256"],
                "bytes": row["bytes"],
                "source_seal_checksum": data["seal_sha256"],
                "status": row["status"],
            })
    for member, (sha, size) in sorted(assembly["verified"].items()):
        inventory.append({
            "component": "yachida_assembly_sensitivity",
            "logical_file": member, "sha256": sha, "bytes": size,
            "source_seal_checksum": assembly["seal_sha256"],
            "status": "VERIFIED",
        })
    for member, (sha, size) in sorted(assembly["sidecar"].items()):
        inventory.append({
            "component": "yachida_assembly_sensitivity",
            "logical_file": member, "sha256": sha, "bytes": size,
            "source_seal_checksum": "not_covered_by_historical_manifest",
            "status": "AUDITED_UNSEALED",
        })

    figure_source = []
    for row in sample_flow:
        if row["condition"] == "TOTAL":
            continue
        figure_source.append({
            "panel": "A", "cohort": row["cohort"], "group": row["condition"],
            "measure": "samples", "value": row["manifest_samples"],
        })
    for row in completeness:
        for measure in ("expected_profiles", "observed_profiles"):
            figure_source.append({
                "panel": "B", "cohort": row["cohort"], "group": row["design"],
                "measure": measure, "value": row[measure],
            })
    for row in balance_rows:
        for measure in ("expected_samples", "observed_samples"):
            figure_source.append({
                "panel": "C", "cohort": row["cohort"],
                "group": row["condition"], "measure": measure,
                "value": row[measure],
            })
    for row in covariate_rows:
        figure_source.append({
            "panel": "D", "cohort": row["cohort"], "group": row["field"],
            "measure": "missing", "value": row["missing"],
        })

    return {
        "cohort_sample_flow": sample_flow,
        "cohort_condition_counts": condition_counts,
        "profile_completeness": completeness,
        "independent_subset_balance": balance_rows,
        "covariate_missingness": covariate_rows,
        "assembly_sensitivity_completeness": assembly_rows,
        "cohort_registry": registry,
        "source_seal_inventory": inventory,
        "figure_source_data": figure_source,
    }


def render_fraction(numerator, denominator) -> str:
    if not denominator:
        return ""
    return format(float(numerator) / float(denominator), ".6f")


# --------------------------------------------------------------- package build
def media_type(path: Path) -> str:
    return MEDIA_TYPES.get(path.suffix, "application/octet-stream")


def release_role(relative: str) -> str:
    if relative.startswith("manuscript_supplement/"):
        return "manuscript_supplement"
    if relative.startswith("zenodo/"):
        return "zenodo_release"
    if relative.startswith(PROJECTION_ROOT + "/"):
        return "source_seal_projection"
    if relative.startswith(EXACT_SEAL_ROOT + "/"):
        return "provenance_source_seal"
    if relative.startswith("provenance/"):
        return "build_provenance"
    return "package_root"


PROJECTION_README = """# Source-seal projection: {component}

**This is NOT the original authoritative seal.** It is a release-safe
projection of it, produced by the upstream evidence package builder.

{explanation}

- `{checksum}` covers exactly the files in this directory and nothing else.
  Verify it from here with `sha256sum -c {checksum}`.
- `{inventory}` records every member of the **original** seal with its
  original SHA-256 and byte count. Those hashes are also in
  `zenodo/source_seal_inventory.tsv` and `provenance/input_checksums.tsv`.
- The original seal's own checksum manifest is deliberately **not** copied
  here, because it indexes members this projection does not contain and would
  therefore fail to verify.
"""

EXACT_SEAL_README = """# Source seal: {component} (exact copy)

This is a byte-identical copy of the original authoritative seal, made because
individual covariate release was explicitly approved. Every original member is
present and the original `{checksum}` verifies unchanged from this directory.
"""


def write_projection(destination: Path, component: str, explanation: str,
                     copied, projected, original_members):
    """Write one release-safe projection with its own verifiable checksums.

    `copied` maps output name -> source path for byte-identical members.
    `projected` maps output name -> (fields, rows) for derived members.
    `original_members` maps original member name -> (sha256, bytes).
    """
    destination.mkdir(parents=True, exist_ok=True)
    for name, source in sorted(copied.items()):
        shutil.copyfile(str(source), str(destination / name))
    for name, (fields, rows) in sorted(projected.items()):
        write_tsv(destination / name, fields, rows)
    write_tsv(destination / PROJECTION_INVENTORY_NAME,
              ["original_member", "sha256", "bytes", "copied_into_projection"],
              [{"original_member": member, "sha256": sha, "bytes": size,
                "copied_into_projection":
                    "yes" if member in copied else "no"}
               for member, (sha, size) in sorted(original_members.items())])
    (destination / PROJECTION_README_NAME).write_text(
        PROJECTION_README.format(
            component=component, explanation=explanation,
            checksum=PROJECTION_CHECKSUM_NAME,
            inventory=PROJECTION_INVENTORY_NAME), encoding="utf-8")
    # The projection's own manifest covers exactly what is here.
    lines = []
    for path in sorted(destination.iterdir()):
        if not path.is_file() or path.name == PROJECTION_CHECKSUM_NAME:
            continue
        lines.append("%s  %s\n" % (digest(path), path.name))
    (destination / PROJECTION_CHECKSUM_NAME).write_text(
        "".join(lines), encoding="utf-8")


def copy_exact_seal(destination: Path, component: str, source_root: Path,
                    members, checksum_name: str):
    """Copy a complete seal byte-for-byte, so its own manifest still verifies.

    Nothing else is placed in the directory: an added README would be a file
    the original manifest cannot cover, which is precisely the defect this
    design avoids. The explanatory note is written beside the components.
    """
    destination.mkdir(parents=True, exist_ok=True)
    for member in sorted(members):
        shutil.copyfile(str(source_root / member), str(destination / member))
    for member in sorted(members):
        if digest(destination / member) != digest(source_root / member):
            raise PackageError("copy of %s/%s changed during packaging"
                               % (component, member))
    shutil.copyfile(str(source_root / checksum_name),
                    str(destination / checksum_name))
    # Prove here, not only upstream, that the copied manifest verifies.
    verify_checksum_manifest(destination, destination / checksum_name,
                             "%s exact copy" % component)


def copy_source_seals(staging: Path, cohorts, assembly, include_individual):
    """Copy release-safe evidence from every source seal.

    By default nothing claims to be an authoritative seal: each component gets
    an explicitly named projection carrying only release-safe members, its own
    checksum manifest over exactly those files, and an inventory of the
    original members' hashes. The original checksum manifest is never copied
    into a directory where it could not verify.

    When individual covariate release is approved, a complete byte-identical
    seal is copied instead, with every original member and its own verifying
    manifest. The two semantics are never mixed.
    """
    root = staging / (EXACT_SEAL_ROOT if include_individual else PROJECTION_ROOT)
    root.mkdir(parents=True, exist_ok=True)
    if include_individual:
        (root / PROJECTION_README_NAME).write_text(
            EXACT_SEAL_README.format(component="all components",
                                     checksum=CHECKSUM_NAME),
            encoding="utf-8")
    for cohort in COHORT_ORDER:
        data = cohorts[cohort]
        destination = root / cohort
        members = dict(
            (member, (sha, size))
            for member, (sha, size) in data["verified"].items())
        members[CHECKSUM_NAME] = (data["seal_sha256"],
                                  (data["root"] / CHECKSUM_NAME).stat().st_size)
        if include_individual:
            copy_exact_seal(destination, cohort, data["root"],
                            [name for name in SEAL_MEMBERS
                             if name != CHECKSUM_NAME], CHECKSUM_NAME)
            continue
        copied = dict(
            (member, data["root"] / member)
            for member in ("SUCCESS", "source_seal_inventory.tsv",
                           "sample_flow.tsv", "covariate_audit.tsv"))
        projected = {
            "production_manifest.canonical.tsv": (
                list(CANONICAL_MANIFEST_COLUMNS), data["canonical"]),
            "production_manifest.independent.canonical.tsv": (
                list(CANONICAL_MANIFEST_COLUMNS), data["independent"]),
        }
        write_projection(
            destination, cohort,
            "The two manifest members are projected to "
            "`sample_id, condition, study, independent_subset`, because the "
            "sealed originals also carry individual covariates that are not "
            "released by default.",
            copied, projected, members)

    destination = root / "yachida_assembly_sensitivity"
    members = dict(assembly["verified"])
    members.update(assembly["sidecar"])
    members[ASSEMBLY_CHECKSUM_NAME] = (
        assembly["seal_sha256"],
        (assembly["root"] / ASSEMBLY_CHECKSUM_NAME).stat().st_size)
    if include_individual:
        copy_exact_seal(
            destination, "yachida_assembly_sensitivity", assembly["root"],
            sorted(set(assembly["verified"]) | set(assembly["sidecar"])),
            ASSEMBLY_CHECKSUM_NAME)
    else:
        write_projection(
            destination, "yachida_assembly_sensitivity",
            "Only the experiment's SUCCESS sentinel is copied; the remaining "
            "sealed members stay in the original seal.",
            {"SUCCESS": assembly["root"] / SUCCESS_NAME}, {}, members)
    return root.name


def verify_packaged_seal_copies(staging: Path) -> int:
    """Independently verify every packaged projection or exact seal copy."""
    checked = 0
    for parent in (staging / PROJECTION_ROOT, staging / EXACT_SEAL_ROOT):
        if not parent.is_dir():
            continue
        for component in sorted(parent.iterdir()):
            if not component.is_dir():
                continue
            manifests = [name for name in
                         (PROJECTION_CHECKSUM_NAME, CHECKSUM_NAME,
                          ASSEMBLY_CHECKSUM_NAME)
                         if (component / name).is_file()]
            if len(manifests) != 1:
                raise PackageError(
                    "%s/%s must carry exactly one checksum manifest; found %s"
                    % (parent.name, component.name,
                       ", ".join(manifests) or "none"))
            verified = verify_checksum_manifest(
                component, component / manifests[0],
                "packaged %s/%s" % (parent.name, component.name))
            present = sorted(item.name for item in component.iterdir()
                             if item.is_file() and item.name != manifests[0])
            unlisted = [name for name in present if name not in verified]
            allowed_unlisted = ({ASSEMBLY_SIDECAR_NAME}
                                if parent.name == EXACT_SEAL_ROOT and
                                component.name == "yachida_assembly_sensitivity"
                                else set())
            unlisted = [name for name in unlisted if name not in allowed_unlisted]
            if unlisted:
                raise PackageError(
                    "packaged %s/%s holds file(s) its checksum manifest does "
                    "not cover: %s"
                    % (parent.name, component.name, ", ".join(unlisted)))
            checked += len(verified)
    if not checked:
        raise PackageError("no packaged seal copy was verified")
    return checked


def write_manifest_and_sums(staging: Path):
    skip = {"MANIFEST.tsv", "SHA256SUMS"}
    entries = []
    for path in sorted(staging.rglob("*")):
        if not path.is_file():
            continue
        relative = str(path.relative_to(staging))
        if relative in skip:
            continue
        entries.append({
            "relative_path": relative, "bytes": path.stat().st_size,
            "sha256": digest(path), "media_type": media_type(path),
            "release_role": release_role(relative),
        })
    write_tsv(staging / "MANIFEST.tsv",
              ["relative_path", "bytes", "sha256", "media_type",
               "release_role"], entries)
    sums = []
    for path in sorted(staging.rglob("*")):
        if not path.is_file():
            continue
        relative = str(path.relative_to(staging))
        if relative == "SHA256SUMS":
            continue
        sums.append("%s  %s\n" % (digest(path), relative))
    (staging / "SHA256SUMS").write_text("".join(sums), encoding="utf-8")
    return entries


def recheck_sums(staging: Path) -> int:
    checked = 0
    for line in (staging / "SHA256SUMS").read_text(
            encoding="utf-8").splitlines():
        if not line.strip():
            continue
        recorded, relative = line.split("  ", 1)
        target = staging / relative
        if not target.is_file():
            raise PackageError("SHA256SUMS lists a missing file: %s" % relative)
        if digest(target) != recorded:
            raise PackageError("SHA256SUMS mismatch for %s" % relative)
        checked += 1
    if not checked:
        raise PackageError("SHA256SUMS is empty")
    return checked


def render_figure(staging: Path, plotter: Path, require_plot: bool):
    """Draw the supplementary figure from the exported source data."""
    source = staging / "manuscript_supplement" / "figure_source_data.tsv"
    outdir = staging / "manuscript_supplement"
    if shutil.which("Rscript") is None:
        if require_plot:
            raise PackageError(
                "Rscript is unavailable and --require-figure was given")
        return False, "Rscript unavailable; figure skipped"
    done = subprocess.run(
        ["Rscript", str(plotter), "--source-data", str(source),
         "--outdir", str(outdir)], stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT)
    if done.returncode != 0:
        detail = done.stdout.decode("utf-8", "replace").strip().splitlines()
        message = detail[-1] if detail else "no output"
        if require_plot:
            raise PackageError("figure rendering failed: %s" % message)
        return False, "figure rendering failed: %s" % message
    missing = [name for name in ("upstream_completeness.pdf",
                                 "upstream_completeness.svg",
                                 "upstream_completeness.png")
               if not nonempty_file(outdir / name)]
    if missing:
        if require_plot:
            raise PackageError(
                "figure rendering produced no %s" % ", ".join(missing))
        return False, "figure rendering produced no %s" % ", ".join(missing)
    return True, "rendered"


README_TEMPLATE = """# Upstream evidence package

**Package version:** {version}
**Source commit:** {commit}
**Built (UTC):** {stamp}
**Status:** {status}

This package is an immutable, checksummed index of the upstream evidence that
was complete before definitive downstream analysis began. It is **not** a copy
of raw reads, reference databases, container images or native profiles: those
are recorded as logical asset identifiers with checksums in
`zenodo/upstream_provenance_metadata.tsv`.

Completeness is an analytical-input property. **It does not validate biological
accuracy**, cellular abundance or biomarker truth.

## Cohorts

| Cohort | Study | Samples | Independent | Baseline | Community | Independent profiles |
|---|---|---:|---:|---:|---:|---:|
{registry_rows}

Assembly-sensitivity experiment: {assembly_samples} samples, {assembly_arms}
arms, {assembly_fractions} fractions per arm,
{assembly_observed}/{assembly_expected} profiles.

## Contents

- `manuscript_supplement/` — the supplementary figure, its exported source
  data, and the two supplementary tables.
- `zenodo/` — release-facing machine-readable tables, including
  `cohort_registry.tsv`.
- `{seal_directory}/` — {seal_explanation}
- `provenance/` — build parameters, input checksums, source commit and the
  validation report.
- `MANIFEST.tsv`, `SHA256SUMS`, `SUCCESS`.

Verify with `sha256sum -c SHA256SUMS` from the package root.

## Privacy

Sample identifiers are the pseudonymous accessions already present in the
sealed manifests. Individual age, sex and BMI are **not** released by default:
`zenodo/covariate_missingness.tsv` carries aggregate counts only, and the
copied manifests under `{seal_directory}/` are projected to
`sample_id, condition, study, independent_subset`. Release-facing tables carry
no cluster-absolute paths. Individual covariates: {covariate_policy}

## Draft supplementary caption

**Upstream cohort and profile completeness.** (A) Samples per cohort and
clinical condition. (B) Expected and observed profiles per cohort and design
(baseline, independent, community). (C) Independent-subset samples per cohort
and condition; the frozen expectation of 10 per condition is marked. (D)
Aggregate covariate missingness per cohort and field. All values are exact
counts taken from `manuscript_supplement/figure_source_data.tsv`, which is
exported from the three verified `production_seal_v2` seals. Completeness of
analytical inputs does not by itself establish biological accuracy.
"""

ZENODO_README = """# Upstream evidence — machine-readable release tables

Every table here derives from three verified `production_seal_v2` cohort seals
and the sealed Yachida assembly-sensitivity experiment. No cluster-absolute
path, raw read, reference database or container image is included. Sample
identifiers are pseudonymous accessions; individual covariates are not
released by default.

| Table | Contents |
|---|---|
| `cohort_registry.tsv` | one release-safe row per cohort |
| `cohort_sample_flow.tsv` | samples per cohort and condition, plus totals |
| `cohort_condition_counts.tsv` | condition counts and cohort fractions |
| `profile_completeness.tsv` | expected versus observed profiles per design |
| `independent_subset_balance.tsv` | the frozen 10/10/10 subset balance |
| `covariate_missingness.tsv` | aggregate missingness only |
| `assembly_sensitivity_completeness.tsv` | the 360-profile arm experiment |
| `spike_panel.tsv` | the frozen implanted panel |
| `taxon_aliases.tsv` | the frozen profiler taxon-alias policy |
| `upstream_provenance_metadata.tsv` | logical asset identities and checksums |
| `audit_ledger.tsv` | operational audit evidence, manually supplied |
| `source_seal_inventory.tsv` | verified members of every original seal |
"""


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for cohort in COHORT_ORDER:
        parser.add_argument("--%s-seal" % cohort, required=True, type=Path,
                            help="%s production_seal_v2 directory" % cohort)
        parser.add_argument("--%s-manifest" % cohort, required=True, type=Path)
        parser.add_argument("--%s-independent-manifest" % cohort,
                            required=True, type=Path)
    parser.add_argument("--assembly-sensitivity-seal", required=True, type=Path)
    parser.add_argument("--spike-panel", required=True, type=Path)
    parser.add_argument("--taxon-aliases", required=True, type=Path)
    parser.add_argument("--provenance-metadata", required=True, type=Path)
    parser.add_argument("--audit-ledger", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--repository-root", type=Path)
    parser.add_argument("--source-commit")
    parser.add_argument("--package-version", default=PACKAGE_VERSION)
    parser.add_argument("--plotter", type=Path,
                        default=Path(__file__).resolve().parent
                        / "plot_upstream_evidence.R")
    parser.add_argument("--analysis-image-sha256",
                        help="SHA-256 of the pinned analysis container, "
                             "cross-checked against container_analysis.")
    for cohort in COHORT_ORDER:
        parser.add_argument("--%s-audit-job" % cohort,
                            default=EXPECTED_AUDIT_JOBS[cohort],
                            help="Expected unified audit job for %s" % cohort)
    parser.add_argument("--require-figure", action="store_true",
                        help="Fail instead of skipping when R is unavailable.")
    parser.add_argument("--include-individual-covariates", action="store_true",
                        help="Only after a redistribution review; also requires "
                             "--redistribution-review.")
    parser.add_argument("--redistribution-review",
                        help="Reference for the review that approved releasing "
                             "individual covariates.")
    for name in sorted(ASSEMBLY_DEFAULT_KEYS):
        parser.add_argument("--assembly-%s-key" % name.replace("_", "-"),
                            default=ASSEMBLY_DEFAULT_KEYS[name],
                            help="SUCCESS field holding %s" % name)
    return parser.parse_args(argv)


def build(args) -> Path:
    if args.include_individual_covariates and not args.redistribution_review:
        raise PackageError(
            "--include-individual-covariates requires --redistribution-review; "
            "individual covariates are not released by default")

    outdir = require_absolute(args.outdir, "--outdir")
    if outdir.exists() and any(outdir.iterdir()):
        raise PackageError(
            "--outdir already exists and is not empty: %s" % outdir)
    staging = outdir.with_name(outdir.name + ".incomplete")
    if staging.exists():
        raise PackageError(
            "a previous staging directory is still present: %s" % staging)

    inputs = {"spike_panel": require_absolute(args.spike_panel, "--spike-panel"),
              "taxon_aliases": require_absolute(args.taxon_aliases,
                                                "--taxon-aliases"),
              "provenance_metadata": require_absolute(
                  args.provenance_metadata, "--provenance-metadata"),
              "audit_ledger": require_absolute(args.audit_ledger,
                                               "--audit-ledger")}
    for label, path in sorted(inputs.items()):
        if not nonempty_file(path):
            raise PackageError("missing or empty %s: %s" % (label, path))

    repository_root = (require_absolute(args.repository_root, "--repository-root")
                       if args.repository_root else None)
    commit = resolve_commit(repository_root, args.source_commit)

    cohort_paths = {}
    for cohort in COHORT_ORDER:
        seal = require_absolute(getattr(args, "%s_seal" % cohort),
                                "--%s-seal" % cohort)
        manifest = require_absolute(getattr(args, "%s_manifest" % cohort),
                                    "--%s-manifest" % cohort)
        independent = require_absolute(
            getattr(args, "%s_independent_manifest" % cohort),
            "--%s-independent-manifest" % cohort)
        inputs["%s_seal" % cohort] = seal / CHECKSUM_NAME
        inputs["%s_manifest" % cohort] = manifest
        inputs["%s_independent_manifest" % cohort] = independent
        cohort_paths[cohort] = (seal, manifest, independent)

    assembly_root = require_absolute(args.assembly_sensitivity_seal,
                                     "--assembly-sensitivity-seal")
    inputs["assembly_sensitivity_seal"] = assembly_root / ASSEMBLY_CHECKSUM_NAME

    # Snapshot before the first source seal is parsed. A mutation during
    # initial validation must not become the unnoticed comparison baseline.
    source_roots = [cohort_paths[cohort][0] for cohort in COHORT_ORDER]
    source_roots.append(assembly_root)
    source_snapshot = seal_inventory_snapshot(source_roots)

    cohorts = {}
    seen_samples = {}
    for cohort in COHORT_ORDER:
        seal, manifest, independent = cohort_paths[cohort]
        cohorts[cohort] = load_cohort_seal(cohort, seal, manifest, independent)
        for row in cohorts[cohort]["flow"]:
            other = seen_samples.get(row["sample_id"])
            if other is not None:
                raise PackageError(
                    "sample %s appears in both %s and %s"
                    % (row["sample_id"], other, cohort))
            seen_samples[row["sample_id"]] = cohort

    assembly_keys = dict(
        (name, getattr(args, "assembly_%s_key" % name))
        for name in ASSEMBLY_DEFAULT_KEYS)
    assembly = load_assembly_sensitivity(assembly_root, assembly_keys)

    if not args.analysis_image_sha256:
        raise PackageError("--analysis-image-sha256 is required")
    if not SHA256_PATTERN.match(args.analysis_image_sha256.strip()):
        raise PackageError(
            "--analysis-image-sha256 must be 64 hexadecimal characters")
    spike_targets, _ = load_spike_panel(inputs["spike_panel"])
    alias_rows = load_taxon_aliases(inputs["taxon_aliases"])
    provenance_rows = load_provenance_metadata(
        inputs["provenance_metadata"], commit, spike_targets,
        (args.analysis_image_sha256 or "").strip())
    scan_for_credentials(provenance_rows, "provenance metadata")
    expected_jobs = dict(
        (cohort, str(getattr(args, "%s_audit_job" % cohort)).strip())
        for cohort in COHORT_ORDER)
    seal_checksums = dict((cohort, cohorts[cohort]["seal_sha256"])
                          for cohort in COHORT_ORDER)
    ledger_rows = load_audit_ledger(inputs["audit_ledger"], COHORT_ORDER,
                                    expected_jobs, seal_checksums)
    scan_for_credentials(ledger_rows, "audit ledger")

    tables = build_tables(cohorts, assembly, provenance_rows, ledger_rows)
    stamp = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")

    try:
        staging.mkdir(parents=True)
        zenodo = staging / "zenodo"
        supplement = staging / "manuscript_supplement"
        provenance = staging / "provenance"
        for directory in (zenodo, supplement, provenance):
            directory.mkdir(parents=True)

        write_tsv(zenodo / "cohort_sample_flow.tsv",
                  ["cohort", "condition", "manifest_samples", "sealed_samples",
                   "independent_samples", "status"],
                  tables["cohort_sample_flow"])
        write_tsv(zenodo / "cohort_condition_counts.tsv",
                  ["cohort", "condition", "samples", "fraction_of_cohort"],
                  tables["cohort_condition_counts"])
        write_tsv(zenodo / "profile_completeness.tsv",
                  ["cohort", "design", "expected_profiles", "observed_profiles",
                   "completion_fraction", "status"],
                  tables["profile_completeness"])
        write_tsv(zenodo / "independent_subset_balance.tsv",
                  ["cohort", "condition", "expected_samples",
                   "observed_samples", "status"],
                  tables["independent_subset_balance"])
        write_tsv(zenodo / "covariate_missingness.tsv",
                  ["cohort", "field", "missing", "present",
                   "distinct_nonmissing"], tables["covariate_missingness"])
        write_tsv(zenodo / "assembly_sensitivity_completeness.tsv",
                  ["experiment", "samples", "assembly_arms",
                   "fractions_per_arm", "expected_profiles",
                   "observed_profiles", "status"],
                  tables["assembly_sensitivity_completeness"])
        write_tsv(zenodo / "cohort_registry.tsv",
                  ["cohort", "study", "samples", "independent_samples",
                   "baseline_profiles", "community_profiles",
                   "independent_profiles", "seal_contract", "seal_status",
                   "input_provenance_mode", "source_audit_job"],
                  tables["cohort_registry"])
        write_tsv(zenodo / "source_seal_inventory.tsv",
                  ["component", "logical_file", "sha256", "bytes",
                   "source_seal_checksum", "status"],
                  tables["source_seal_inventory"])
        write_tsv(zenodo / "upstream_provenance_metadata.tsv",
                  PROVENANCE_FIELDS, provenance_rows)
        write_tsv(zenodo / "audit_ledger.tsv", AUDIT_LEDGER_FIELDS, ledger_rows)
        shutil.copyfile(str(inputs["spike_panel"]), str(zenodo / "spike_panel.tsv"))
        # The frozen alias file is comma-delimited; it is parsed and rewritten
        # rather than copied under a .tsv name it would not honour.
        write_tsv(zenodo / "taxon_aliases.tsv", ALIAS_OUTPUT_COLUMNS,
                  alias_rows)
        (zenodo / "README.md").write_text(ZENODO_README, encoding="utf-8")

        write_tsv(supplement / "figure_source_data.tsv",
                  ["panel", "cohort", "group", "measure", "value"],
                  tables["figure_source_data"])
        supplementary_flow = []
        for cohort in COHORT_ORDER:
            for row in cohorts[cohort]["flow"]:
                supplementary_flow.append({
                    "cohort": cohort, "sample_id": row["sample_id"],
                    "condition": row["condition"],
                    "independent_subset": row["independent_subset"],
                    "expected_profiles": row["expected_profiles"],
                    "observed_profiles": row["observed_profiles"],
                    "status": row["status"],
                })
        write_tsv(supplement / "supplementary_sample_flow.tsv",
                  ["cohort", "sample_id", "condition", "independent_subset",
                   "expected_profiles", "observed_profiles", "status"],
                  supplementary_flow)
        write_tsv(supplement / "supplementary_profile_completeness.tsv",
                  ["cohort", "design", "expected_profiles", "observed_profiles",
                   "completion_fraction", "status"],
                  tables["profile_completeness"])

        seal_directory = copy_source_seals(
            staging, cohorts, assembly, args.include_individual_covariates)

        figure_ok, figure_note = render_figure(
            staging, require_absolute(args.plotter, "--plotter"),
            args.require_figure)

        write_tsv(provenance / "build_parameters.tsv", ["parameter", "value"], [
            {"parameter": "package_version", "value": args.package_version},
            {"parameter": "seal_contract", "value": SEAL_CONTRACT},
            {"parameter": "built_utc", "value": stamp},
            {"parameter": "source_commit", "value": commit},
            {"parameter": "analysis_image_sha256",
             "value": (args.analysis_image_sha256 or "").strip()
                      or "not_supplied"},
            {"parameter": "expected_audit_jobs",
             "value": ";".join("%s=%s" % (cohort, expected_jobs[cohort])
                               for cohort in COHORT_ORDER)},
            {"parameter": "packaged_seal_mode",
             "value": "exact_copy" if args.include_individual_covariates
                      else "release_projection"},
            {"parameter": "individual_covariates_released",
             "value": "yes" if args.include_individual_covariates else "no"},
            {"parameter": "redistribution_review",
             "value": args.redistribution_review or "not_applicable"},
            {"parameter": "figure_rendered", "value": "yes" if figure_ok else "no"},
            {"parameter": "figure_note", "value": figure_note},
        ])
        write_tsv(provenance / "input_checksums.tsv",
                  ["logical_input", "basename", "sha256", "bytes"],
                  [{"logical_input": name, "basename": path.name,
                    "sha256": digest(path), "bytes": path.stat().st_size}
                   for name, path in sorted(inputs.items())])
        (provenance / "repository_commit.txt").write_text(
            commit + "\n", encoding="utf-8")

        validation = [
            {"check": "source_seal_checksums", "scope": "all", "status": "PASS",
             "detail": "every member of every source seal verified"},
            {"check": "cohort_identity_and_size", "scope": "all",
             "status": "PASS", "detail": "201/154/156 with 30 independent each"},
            {"check": "sample_set_agreement", "scope": "all", "status": "PASS",
             "detail": "frozen manifests and seals cover identical samples"},
            {"check": "profile_topology", "scope": "all", "status": "PASS",
             "detail": "baseline, community and independent totals exact"},
            {"check": "sample_status", "scope": "all", "status": "PASS",
             "detail": "every sealed sample row is PASS"},
            {"check": "independent_subset_balance", "scope": "all",
             "status": "PASS", "detail": "10 Control / 10 Adenoma / 10 CRC"},
            {"check": "assembly_sensitivity", "scope": "yachida",
             "status": "PASS", "detail": "30 samples, 2 arms, 6 fractions, 360/360"},
            {"check": "provenance_metadata", "scope": "all", "status": "PASS",
             "detail": "required categories present, no blank identities"},
            {"check": "credential_scan", "scope": "all", "status": "PASS",
             "detail": "no credential or private-key pattern found"},
            {"check": "figure", "scope": "manuscript_supplement",
             "status": "PASS" if figure_ok else "SKIPPED", "detail": figure_note},
            {"check": "audit_ledger", "scope": "all", "status": "PASS",
             "detail": "one COMPLETED 0:0 row per cohort, job and seal digest "
                       "cross-checked"},
            {"check": "packaged_seal_copies", "scope": "all", "status": "PASS",
             "detail": "every packaged projection or exact copy verifies "
                       "against its own checksum manifest"},
            {"check": "source_seal_immutability", "scope": "all",
             "status": "PASS",
             "detail": "no source seal file was added, removed or changed "
                       "during the build"},
        ]
        write_tsv(provenance / "validation_report.tsv",
                  ["check", "scope", "status", "detail"], validation)

        registry_rows = "\n".join(
            "| %s | %s | %d | %d | %d | %d | %d |"
            % (row["cohort"], row["study"], row["samples"],
               row["independent_samples"], row["baseline_profiles"],
               row["community_profiles"], row["independent_profiles"])
            for row in tables["cohort_registry"])
        covariate_policy = (
            "released under review %s" % args.redistribution_review
            if args.include_individual_covariates
            else "withheld pending a separate redistribution review.")
        seal_explanation = (
            "byte-identical copies of each authoritative seal, each verifying "
            "against its own original checksum manifest."
            if args.include_individual_covariates else
            "release-safe projections of each source seal. Each carries only "
            "release-safe members, its own `projection.sha256` covering "
            "exactly the files present, and an "
            "`original_source_seal_inventory.tsv` of every original member's "
            "SHA-256 and size. **These are not the original authoritative "
            "seals**; an original manifest is never copied where it could not "
            "verify.")
        (staging / "README.md").write_text(README_TEMPLATE.format(
            version=args.package_version, commit=commit, stamp=stamp,
            status="PASS", registry_rows=registry_rows,
            seal_directory=seal_directory, seal_explanation=seal_explanation,
            assembly_samples=assembly["observed"]["samples"],
            assembly_arms=assembly["observed"]["assembly_arms"],
            assembly_fractions=assembly["observed"]["fractions_per_arm"],
            assembly_observed=assembly["observed"]["observed_profiles"],
            assembly_expected=assembly["observed"]["expected_profiles"],
            covariate_policy=covariate_policy), encoding="utf-8")

        # No product of this package may carry a cluster-absolute path, so
        # the whole staged tree is screened, not only the zenodo tables.
        assert_release_safe(staging)

        # Every packaged projection or exact copy must verify on its own.
        packaged = verify_packaged_seal_copies(staging)
        # And no authoritative seal may have changed while we worked.
        assert_sources_unchanged(source_roots, source_snapshot)

        # SUCCESS is evidence, not a checksum file. Create it first so both
        # MANIFEST.tsv and SHA256SUMS protect it. MANIFEST.tsv is then covered
        # by SHA256SUMS; only SHA256SUMS itself is necessarily self-exempt.
        archived_files = sum(
            1 for path in staging.rglob("*") if path.is_file()) + 1
        checksummed_files = archived_files + 1
        (staging / SUCCESS_NAME).write_text(
            "package_version\t%s\n"
            "seal_contract\t%s\n"
            "source_commit\t%s\n"
            "yachida_samples\t%d\n"
            "feng_samples\t%d\n"
            "zeller_samples\t%d\n"
            "assembly_sensitivity_profiles\t%d\n"
            "files\t%d\n"
            "checksums_verified\t%d\n"
            "status\tPASS\n"
            % (args.package_version, SEAL_CONTRACT, commit,
               FROZEN_COHORTS["yachida"]["samples"],
               FROZEN_COHORTS["feng"]["samples"],
               FROZEN_COHORTS["zeller"]["samples"],
               assembly["observed"]["observed_profiles"],
               archived_files, checksummed_files), encoding="utf-8")
        entries = write_manifest_and_sums(staging)
        checked = recheck_sums(staging)
        if len(entries) != archived_files or checked != checksummed_files:
            raise PackageError(
                "internal package file-count mismatch: SUCCESS records %d/%d, "
                "observed %d/%d"
                % (archived_files, checksummed_files, len(entries), checked))
    except BaseException:
        # An incomplete package is never left behind for someone to mistake
        # for a finished one.
        shutil.rmtree(str(staging), ignore_errors=True)
        raise
    os.replace(str(staging), str(outdir))
    return outdir


def main() -> None:
    try:
        args = parse_args()
        outdir = build(args)
    except (PackageError, OSError, KeyError, ValueError) as error:
        raise SystemExit("[ERROR] %s" % error)
    print("[PASS] Upstream evidence package: %s" % outdir)


if __name__ == "__main__":
    main()
