#!/usr/bin/env python3
"""One configurable upstream auditor and sealer for Yachida, Feng and Zeller.

Every paper cohort is validated against the same scientific and integrity
contract and emits the same `production_seal_v2` layout. Cohort differences are
confined to an adapter at the input boundary (which column carries the
condition, whether the study is a manifest column or a fixed value, whether an
upstream batch exists, which logical members and which SUCCESS schema the
native seal uses). No check is relaxed for any cohort.

This is a validation and interface migration. It never recomputes a profile,
rewrites a retained output, or modifies a native seal: the native seal is an
immutable provenance input whose checksums are verified, whose SUCCESS is
semantically checked against the selected cohort, and whose manifest copies
must be byte-identical to the authoritative manifests. That seal is verified
again immediately before publication authority is granted, so a seal that
changes during the audit cannot be sealed over.

Standard library only, and `from __future__ import annotations` keeps the
annotations unevaluated so the cluster's older Python can run it.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import sys
from collections import Counter
from pathlib import Path

SEAL_CONTRACT = "upstream_seal_v2"

# Frozen scientific design. These are not configurable: a cohort that does not
# meet them is not sealable.
BASELINE_PER_SAMPLE = 1
COMMUNITY_PER_SAMPLE = 7
INDEPENDENT_PER_SUBSET_SAMPLE = 60
FROZEN_INDEPENDENT_SAMPLES = 30
FROZEN_SUBSET_BALANCE = {"Control": 10, "Adenoma": 10, "CRC": 10}

# Frozen production sizes and condition counts. Cross-checked only when the
# audit is actually run at production scale, so scaled fixtures stay usable
# while a real cohort can never be sealed against drifted expectations.
FROZEN_COHORTS = {
    "yachida": (201, {"Control": 67, "Adenoma": 67, "CRC": 67}),
    "feng": (154, {"Control": 61, "Adenoma": 47, "CRC": 46}),
    "zeller": (156, {"Control": 61, "Adenoma": 42, "CRC": 53}),
}

CANONICAL_COLUMNS = ("sample_id", "condition", "study", "independent_subset")

# One sample-flow schema for every cohort. `batch_id` is always present and is
# blank for cohorts without upstream batches.
SAMPLE_FLOW_FIELDS = [
    "sample_id", "study", "condition", "independent_subset",
    "batch_id",
    "expected_baseline_profiles", "observed_baseline_profiles",
    "expected_independent_profiles", "observed_independent_profiles",
    "expected_community_profiles", "observed_community_profiles",
    "expected_profiles", "observed_profiles", "retained_output_files",
    "verified_marker", "retained_output_receipt", "input_provenance",
    "sample_success",
    "baseline_profile_count", "independent_profile_count",
    "community_profile_count",
    "completion_table", "manifest_independent_flag", "status",
    "failure_reasons", "receipt_error", "completion_error",
]
SAMPLE_CHECKS = [
    "verified_marker", "retained_output_receipt", "input_provenance",
    "sample_success", "baseline_profile_count", "independent_profile_count",
    "community_profile_count", "completion_table", "manifest_independent_flag",
]
COVARIATE_FIELDS = ["field", "missing", "present", "distinct_nonmissing"]
INVENTORY_FIELDS = ["member", "sha256", "bytes", "status", "source_audit_job"]

SEAL_MEMBERS = (
    "sample_flow.tsv", "covariate_audit.tsv", "production_manifest.tsv",
    "production_manifest.independent.tsv", "source_seal_inventory.tsv",
)
SUCCESS_NAME = "SUCCESS"
CHECKSUM_NAME = "production_seal.sha256"
IN_PROGRESS_NAME = "AUDIT_IN_PROGRESS"
OUTDIR_MEMBERS = SEAL_MEMBERS + (SUCCESS_NAME, CHECKSUM_NAME, IN_PROGRESS_NAME)
# Only the temporary siblings this program itself can create are recoverable;
# anything else in the output directory is unexpected and fails closed.
RECOVERABLE_TEMPORARIES = frozenset(name + ".tmp" for name in OUTDIR_MEMBERS)
PERMITTED_OUTDIR_MEMBERS = frozenset(OUTDIR_MEMBERS) | RECOVERABLE_TEMPORARIES

# Boundary translation only. The Feng audit ran as Slurm job 3097587 and the
# Zeller audit as 3097588; no Yachida job identifier is on record, and one is
# never inferred, so it stays blank.
ADAPTERS = {
    "yachida": {
        "condition_column": "Target_Condition",
        "study_column": None,
        "study_value": "YachidaS_2019",
        "batch_column": "batch_id",
        "completion_has_study": False,
        # Audited covariates must exist in the frozen manifest. The sealed
        # 29-column production manifest carries all four of these; add others
        # with --covariate-column only if the frozen manifest gains them.
        "covariates": ("age", "sex", "bmi", "batch_id"),
        # Batch provenance exists only in the production manifest; the sealed
        # 26-column independent manifest carries none of it.
        "production_only_columns": ("batch_id", "batch_hash", "batch_position",
                                    "batch_size", "batch_seed",
                                    "processing_order"),
        "independent_header_schema": "yachida_historical_duplicate_selection",
        "seal_manifest_member": "pilot_batched.tsv",
        "seal_independent_member": "independent_10_per_condition.tsv",
        "source_success_schema": "yachida_dataset",
        "source_audit_job": "",
    },
    "feng": {
        "condition_column": "condition",
        "study_column": "study",
        "study_value": None,
        "batch_column": None,
        "completion_has_study": True,
        "covariates": ("age", "sex", "bmi"),
        "production_only_columns": (),
        "independent_header_schema": "unique",
        "seal_manifest_member": "production_manifest.tsv",
        "seal_independent_member": "production_manifest.independent.tsv",
        "source_success_schema": "cohort_profiles",
        "source_audit_job": "3097587",
    },
    "zeller": {
        "condition_column": "condition",
        "study_column": "study",
        "study_value": None,
        "batch_column": None,
        "completion_has_study": True,
        "covariates": ("age", "sex", "bmi"),
        "production_only_columns": (),
        "independent_header_schema": "unique",
        "seal_manifest_member": "production_manifest.tsv",
        "seal_independent_member": "production_manifest.independent.tsv",
        "source_success_schema": "cohort_profiles",
        "source_audit_job": "3097588",
    },
}
SOURCE_SUCCESS_SCHEMAS = ("yachida_dataset", "cohort_profiles")

# Independent-manifest header handling.
#   unique                                 -- duplicates are always rejected
#   yachida_historical_duplicate_selection -- accepts exactly two formats and
#       nothing else:
#         A. the historical sealed header, carrying selection_rank,
#            selection_hash and selection_seed exactly twice each. The
#            checksummed file is never rewritten, so it is read positionally and
#            translated by occurrence: first -> pilot_selection_*, second ->
#            independent_selection_*.
#         B. the corrected header the migrated selector now emits, carrying
#            exactly one pilot_selection_* triplet and exactly one
#            independent_selection_* triplet and no bare selection_* column.
#       A lone bare triplet, a partial or missing pilot/independent triplet, a
#       mixture of bare and prefixed fields, and any duplicate outside the exact
#       historical triplet are all refused. This is a narrow, explicitly
#       selected translation of one known format, not tolerance of duplicate
#       headers in general.
INDEPENDENT_HEADER_SCHEMAS = ("unique", "yachida_historical_duplicate_selection")
HISTORICAL_DUPLICATE_FIELDS = ("selection_rank", "selection_hash",
                               "selection_seed")
HISTORICAL_FIRST_PREFIX = "pilot"
HISTORICAL_SECOND_PREFIX = "independent"


class SealError(Exception):
    """A fail-closed contract violation."""


# --------------------------------------------------------------------- io
def digest_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def nonempty_file(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def marker_exists(path: Path) -> bool:
    """Production sentinels are created with touch, so zero bytes is valid."""
    return path.is_file()


def remove_if_present(path: Path) -> None:
    """`Path.unlink(missing_ok=)` is too new for the cluster interpreter."""
    try:
        path.unlink()
    except OSError:
        if path.exists():
            raise


def write_atomic(path: Path, payload: bytes) -> None:
    """Replace one member in a single step; readers never see a partial file.

    This is per-member atomicity, not whole-directory atomicity. Directory-level
    safety comes from the authority ordering in `main`: stale authority is
    dropped first, the checksum manifest is written before the SUCCESS it
    covers, and an interrupted run leaves AUDIT_IN_PROGRESS behind.
    """
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(str(temporary), str(path))


def tsv_bytes(fields, rows) -> bytes:
    lines = ["\t".join(fields)]
    for row in rows:
        lines.append("\t".join(str(row.get(name, "")) for name in fields))
    return ("\n".join(lines) + "\n").encode("utf-8")


def read_header(path: Path):
    """Return the raw header fields, rejecting duplicates DictReader would hide."""
    with path.open(newline="", encoding="utf-8") as handle:
        line = handle.readline()
    if not line.strip():
        raise SealError("%s has no header" % path.name)
    header = line.rstrip("\r\n").split("\t")
    duplicates = sorted({name for name in header if header.count(name) > 1})
    if duplicates:
        raise SealError(
            "%s has duplicate column(s): %s" % (path.name, ", ".join(duplicates)))
    if any(not name.strip() for name in header):
        raise SealError("%s has a blank column name" % path.name)
    return header


def read_table(path: Path):
    if not nonempty_file(path):
        raise SealError("missing or empty table: %s" % path)
    header = read_header(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SealError("%s has no data rows" % path.name)
    return header, rows


def raw_header(path: Path):
    """The header exactly as written, occurrences preserved."""
    with path.open(newline="", encoding="utf-8") as handle:
        line = handle.readline()
    if not line.strip():
        raise SealError("%s has no header" % path.name)
    header = line.rstrip("\r\n").split("\t")
    if any(not name.strip() for name in header):
        raise SealError("%s has a blank column name" % path.name)
    return header


def translate_historical_independent_header(path: Path, header):
    """Translate the one documented duplicated Yachida header by occurrence.

    The sealed independent manifest was produced by running the deterministic
    selector over a manifest that already carried pilot selection provenance,
    so it holds `selection_rank`, `selection_hash` and `selection_seed` twice:
    the first occurrence is the inherited pilot value, the second is the
    independent-subset value. The file is checksummed and is never rewritten,
    so it is read positionally and renamed into unambiguous canonical names.

    Returns the translated header plus the production column each translated
    name must be compared against.
    """
    counts = {}
    for name in header:
        counts[name] = counts.get(name, 0) + 1
    repeated = sorted(name for name, total in counts.items() if total > 1)
    unexpected = [name for name in repeated
                  if name not in HISTORICAL_DUPLICATE_FIELDS]
    if unexpected:
        raise SealError(
            "%s has duplicate column(s) outside the documented historical "
            "selection triplet: %s" % (path.name, ", ".join(unexpected)))

    def observed(prefix):
        return dict(
            (field, counts.get(
                "%s_%s" % (prefix, field) if prefix else field, 0))
            for field in HISTORICAL_DUPLICATE_FIELDS)

    bare = observed("")
    pilot = observed(HISTORICAL_FIRST_PREFIX)
    independent = observed(HISTORICAL_SECOND_PREFIX)
    summary = ", ".join(
        "%s=%d" % (name, total) for name, total in sorted(
            list(bare.items())
            + [("%s_%s" % (HISTORICAL_FIRST_PREFIX, name), total)
               for name, total in pilot.items()]
            + [("%s_%s" % (HISTORICAL_SECOND_PREFIX, name), total)
               for name, total in independent.items()]))
    prefixed_present = any(pilot.values()) or any(independent.values())

    if all(total == 1 for total in pilot.values()) and all(
            total == 1 for total in independent.values()) and not any(
            bare.values()):
        # Format B: already unambiguous, as the corrected selector emits. A
        # `pilot_<field>` column carries the same inherited provenance the
        # historical first occurrence does, so it is verified the same way.
        provenance = dict(
            ("%s_%s" % (HISTORICAL_FIRST_PREFIX, field), field)
            for field in HISTORICAL_DUPLICATE_FIELDS)
        return list(header), provenance

    if not all(total == 2 for total in bare.values()):
        if prefixed_present:
            raise SealError(
                "%s carries neither the historical duplicated selection "
                "triplet nor a complete pilot/independent pair; observed %s"
                % (path.name, summary))
        raise SealError(
            "%s does not match the documented historical header: the "
            "selection triplet must appear exactly twice each, observed %s"
            % (path.name, summary))
    if prefixed_present:
        raise SealError(
            "%s mixes the historical duplicated selection triplet with "
            "prefixed selection fields; observed %s" % (path.name, summary))

    translated = []
    provenance = {}
    seen = {}
    for name in header:
        if name in HISTORICAL_DUPLICATE_FIELDS:
            occurrence = seen.get(name, 0)
            seen[name] = occurrence + 1
            prefix = (HISTORICAL_FIRST_PREFIX if occurrence == 0
                      else HISTORICAL_SECOND_PREFIX)
            renamed = "%s_%s" % (prefix, name)
            if renamed in header:
                raise SealError(
                    "%s already contains %s, so the historical selection "
                    "triplet cannot be translated unambiguously"
                    % (path.name, renamed))
            translated.append(renamed)
            # The pilot occurrence is the production manifest's own value.
            if occurrence == 0:
                provenance[renamed] = name
        else:
            translated.append(name)
    duplicates = sorted({name for name in translated
                         if translated.count(name) > 1})
    if duplicates:
        raise SealError(
            "%s cannot be translated to a unique header; duplicate(s): %s"
            % (path.name, ", ".join(duplicates)))
    return translated, provenance


def read_independent_table(path: Path, adapter):
    """Read the independent manifest under the configured header schema."""
    if not nonempty_file(path):
        raise SealError("missing or empty table: %s" % path)
    header = raw_header(path)
    schema = adapter["independent_header_schema"]
    if schema == "unique":
        duplicates = sorted({name for name in header if header.count(name) > 1})
        if duplicates:
            raise SealError(
                "%s has duplicate column(s): %s"
                % (path.name, ", ".join(duplicates)))
        translated, provenance = list(header), {}
    elif schema == "yachida_historical_duplicate_selection":
        translated, provenance = translate_historical_independent_header(
            path, header)
    else:
        raise SealError("unknown independent header schema %r" % schema)
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        for number, values in enumerate(reader, start=2):
            if not values:
                continue
            if len(values) != len(translated):
                raise SealError(
                    "%s line %d has %d field(s) but the header has %d"
                    % (path.name, number, len(values), len(translated)))
            rows.append(dict(zip(translated, values)))
    if not rows:
        raise SealError("%s has no data rows" % path.name)
    return translated, rows, provenance


def read_key_value(path: Path, label: str):
    """Parse a two-column `field<TAB>value` table such as a native SUCCESS."""
    values = {}
    for number, line in enumerate(
            path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != 2 or not parts[0].strip():
            raise SealError("malformed %s line %d: %r" % (label, number, line))
        if parts[0] in values:
            raise SealError("duplicate %s field: %s" % (label, parts[0]))
        values[parts[0]] = parts[1]
    if not values:
        raise SealError("%s is empty" % label)
    return values


# ----------------------------------------------------------------- paths
def require_absolute(path, label: str) -> Path:
    expanded = Path(os.path.expanduser(str(path)))
    if not expanded.is_absolute():
        raise SealError("%s must be an absolute path: %s" % (label, path))
    if ".." in expanded.parts:
        raise SealError("%s must not contain '..': %s" % (label, path))
    return expanded.resolve()


def contains(parent: Path, child: Path) -> bool:
    return parent == child or parent in child.parents


def resolve_paths(args):
    """Resolve and safety-check every path BEFORE anything is created or removed.

    Overlap between the output directory and the immutable native seal is
    rejected here, so a bad invocation can never reach `mkdir`, `unlink` or
    `write_atomic`.
    """
    paths = {
        "manifest": require_absolute(args.manifest, "--manifest"),
        "independent_manifest": require_absolute(
            args.independent_manifest, "--independent-manifest"),
        "source_seal": require_absolute(args.source_seal, "--source-seal"),
        "state_dir": require_absolute(args.state_dir, "--state-dir"),
        "results_root": require_absolute(args.results_root, "--results-root"),
        "qc_root": require_absolute(args.qc_root, "--qc-root"),
        "scratch_root": require_absolute(args.scratch_root, "--scratch-root"),
        "outdir": require_absolute(args.outdir, "--outdir"),
    }
    seal_root = paths["source_seal"]
    outdir = paths["outdir"]
    if contains(seal_root, outdir) or contains(outdir, seal_root):
        raise SealError(
            "--outdir must not overlap the native source seal (outdir=%s, "
            "source-seal=%s); the native seal is never modified"
            % (outdir, seal_root))
    if contains(paths["scratch_root"], outdir):
        raise SealError("--outdir must not be inside disposable scratch")
    return paths


# --------------------------------------------------- native source seal
def source_seal_inventory(seal_root: Path, adapter):
    """Verify the native seal's checksum manifest and inventory its members."""
    if not seal_root.is_dir():
        raise SealError("native source seal is not a directory: %s" % seal_root)
    if (seal_root / IN_PROGRESS_NAME).exists():
        raise SealError(
            "native source seal is mid-audit (%s present): %s"
            % (IN_PROGRESS_NAME, seal_root))
    success = seal_root / SUCCESS_NAME
    if not nonempty_file(success):
        raise SealError("native source seal has no nonempty SUCCESS: %s" % seal_root)
    checksum = seal_root / CHECKSUM_NAME
    if not nonempty_file(checksum):
        raise SealError(
            "native source seal has no nonempty %s: %s" % (CHECKSUM_NAME, seal_root))

    inventory = []
    listed = []
    for number, line in enumerate(
            checksum.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) != 2:
            raise SealError(
                "malformed source-seal checksum line %d in %s" % (number, seal_root))
        recorded, member = parts
        if member in listed:
            raise SealError("duplicate source-seal member: %s" % member)
        if Path(member).is_absolute() or ".." in Path(member).parts:
            raise SealError("unsafe source-seal member path: %s" % member)
        listed.append(member)
        target = seal_root / member
        if not nonempty_file(target):
            raise SealError("source-seal member missing or empty: %s" % member)
        actual = digest(target)
        if actual != recorded:
            raise SealError(
                "source-seal checksum mismatch for %s in %s" % (member, seal_root))
        inventory.append({
            "member": member, "sha256": actual,
            "bytes": target.stat().st_size, "status": "VERIFIED",
            "source_audit_job": adapter["source_audit_job"],
        })
    if SUCCESS_NAME not in listed:
        raise SealError(
            "native source seal does not checksum its own SUCCESS: %s" % seal_root)
    for member, label in (
        (adapter["seal_manifest_member"], "production manifest"),
        (adapter["seal_independent_member"], "independent manifest"),
    ):
        if member not in listed:
            raise SealError(
                "native source seal does not contain the %s copy %s"
                % (label, member))
    inventory.append({
        "member": CHECKSUM_NAME, "sha256": digest(checksum),
        "bytes": checksum.stat().st_size, "status": "VERIFIED",
        "source_audit_job": adapter["source_audit_job"],
    })
    return inventory


def check_manifest_identity(seal_root: Path, adapter, manifest: Path,
                            independent_manifest: Path) -> None:
    for member, authoritative, label in (
        (adapter["seal_manifest_member"], manifest, "production manifest"),
        (adapter["seal_independent_member"], independent_manifest,
         "independent manifest"),
    ):
        if (seal_root / member).read_bytes() != authoritative.read_bytes():
            raise SealError(
                "%s is not byte-identical to the native source-seal copy %s"
                % (label, member))


def parse_profile_field(raw: str, label: str):
    """Parse Yachida's `baseline=N;independent=M;community=K` profile summary."""
    totals = {}
    for item in raw.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise SealError("malformed %s entry: %r" % (label, item))
        name, value = item.split("=", 1)
        try:
            totals[name.strip()] = int(value)
        except ValueError:
            raise SealError("non-integer %s entry: %r" % (label, item))
    return totals


def check_source_success(seal_root: Path, adapter, cohort: str,
                         expected_samples: int, expected_independent: int,
                         expected_conditions) -> None:
    """Semantically validate the native SUCCESS against the selected cohort.

    Checksum verification only proves the file is the one that was sealed. It
    must also describe *this* cohort at *these* expected counts, otherwise a
    perfectly intact seal for a different cohort or size would be accepted.
    The two documented native schemas are selected explicitly; nothing is
    inferred from the file's shape.
    """
    schema = adapter["source_success_schema"]
    if schema not in SOURCE_SUCCESS_SCHEMAS:
        raise SealError("unknown native SUCCESS schema %r" % schema)
    values = read_key_value(seal_root / SUCCESS_NAME, "native source SUCCESS")
    expected_totals = {
        "baseline": expected_samples * BASELINE_PER_SAMPLE,
        "independent": expected_independent * INDEPENDENT_PER_SUBSET_SAMPLE,
        "community": expected_samples * COMMUNITY_PER_SAMPLE,
    }

    def require(field):
        if field not in values:
            raise SealError(
                "native source SUCCESS lacks %r for the %s schema"
                % (field, schema))
        return values[field]

    def require_int(field, expected):
        raw = require(field)
        try:
            observed = int(raw)
        except ValueError:
            raise SealError(
                "native source SUCCESS has a non-integer %s: %r" % (field, raw))
        if observed != expected:
            raise SealError(
                "native source SUCCESS reports %s=%d; this audit expects %d"
                % (field, observed, expected))

    if require("status") != "PASS":
        raise SealError(
            "native source SUCCESS status is %r, not PASS" % values.get("status"))

    if schema == "yachida_dataset":
        identity = require("dataset")
        expected_identity = adapter["study_value"] or cohort
        if identity != expected_identity:
            raise SealError(
                "native source SUCCESS dataset is %r; this audit is configured "
                "for %r" % (identity, expected_identity))
        require_int("samples", expected_samples)
        require_int("independent_subset", expected_independent)
        totals = parse_profile_field(require("profiles"), "native profiles")
        missing = [name for name in sorted(expected_totals) if name not in totals]
        if missing:
            raise SealError(
                "native source SUCCESS profiles lacks %s" % ", ".join(missing))
        for name, expected in sorted(expected_totals.items()):
            if totals[name] != expected:
                raise SealError(
                    "native source SUCCESS reports %s profiles=%d; this "
                    "audit expects %d" % (name, totals[name], expected))
        observed = parse_profile_field(require("conditions"), "native conditions")
        missing = [name for name in sorted(FROZEN_SUBSET_BALANCE)
                   if name not in observed]
        if missing:
            raise SealError(
                "native source SUCCESS conditions lacks %s" % ", ".join(missing))
        if observed != expected_conditions:
            raise SealError(
                "native source SUCCESS conditions %s disagree with the "
                "expected %s" % (observed, expected_conditions))
    else:
        identity = require("cohort")
        if identity != cohort:
            raise SealError(
                "native source SUCCESS cohort is %r; this audit is configured "
                "for %r" % (identity, cohort))
        require_int("samples", expected_samples)
        require_int("independent_samples", expected_independent)
        for field, name in (("baseline_profiles", "baseline"),
                            ("independent_profiles", "independent"),
                            ("community_profiles", "community")):
            require_int(field, expected_totals[name])


def reverify_source_seal(seal_root: Path, adapter, manifest: Path,
                         independent_manifest: Path, baseline_inventory) -> None:
    """Prove the native seal did not change while this audit ran.

    Called immediately before publication authority is granted. A seal that was
    edited, re-audited or replaced mid-run must not be sealed over.
    """
    if (seal_root / IN_PROGRESS_NAME).exists():
        raise SealError(
            "native source seal became mid-audit during this audit: %s" % seal_root)
    try:
        current = source_seal_inventory(seal_root, adapter)
    except SealError as error:
        # It verified at the start, so any verification failure now means the
        # seal changed while this audit was running.
        raise SealError(
            "native source seal changed during this audit: %s" % error)
    if current != baseline_inventory:
        before = {row["member"]: row["sha256"] for row in baseline_inventory}
        after = {row["member"]: row["sha256"] for row in current}
        changed = sorted(
            name for name in set(before) | set(after)
            if before.get(name) != after.get(name))
        raise SealError(
            "native source seal changed during this audit: %s"
            % ", ".join(changed or ["member set differs"]))
    try:
        check_manifest_identity(seal_root, adapter, manifest, independent_manifest)
    except SealError as error:
        raise SealError(
            "native source seal changed during this audit: %s" % error)


# ------------------------------------------------------ retained outputs
def verify_receipt(receipt: Path, scratch: Path, permitted):
    """Rehash every retained output and enforce the permitted roots."""
    if not nonempty_file(receipt):
        return False, 0, "missing_or_empty_receipt"
    try:
        with receipt.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["path", "sha256", "bytes"]:
                return False, 0, "invalid_receipt_header"
            rows = list(reader)
        if not rows:
            return False, 0, "empty_receipt"
        seen = set()
        for index, row in enumerate(rows, start=2):
            raw = row.get("path", "")
            if not raw or not row.get("sha256") or not row.get("bytes"):
                return False, len(rows), "malformed_receipt_row_%d" % index
            path = Path(raw).resolve()
            if path in seen:
                return False, len(rows), "duplicate_receipt_path_%d" % index
            seen.add(path)
            if contains(scratch, path):
                return False, len(rows), "retained_output_inside_scratch_%d" % index
            if not any(contains(root, path) for root in permitted):
                return False, len(rows), (
                    "retained_output_outside_sample_roots_%d" % index)
            if not nonempty_file(path):
                return False, len(rows), "missing_or_empty_output_%d" % index
            try:
                expected_bytes = int(row["bytes"])
            except ValueError:
                return False, len(rows), "invalid_bytes_%d" % index
            if expected_bytes <= 0 or path.stat().st_size != expected_bytes:
                return False, len(rows), "size_mismatch_%d" % index
            if digest(path) != row["sha256"]:
                return False, len(rows), "sha256_mismatch_%d" % index
        return True, len(rows), ""
    except (OSError, csv.Error) as error:
        detail = str(error).replace("\t", " ").replace("\n", " ")
        return False, 0, "receipt_read_error:%s" % detail


def read_completion(path: Path):
    """Read the runner's two-column sample-completion contract."""
    if not nonempty_file(path):
        return {}, "missing_or_empty_completion_table"
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        if not rows or rows[0] != ["field", "value"]:
            return {}, "invalid_completion_header"
        values = {}
        for index, row in enumerate(rows[1:], start=2):
            if len(row) != 2 or not row[0] or row[0] in values:
                return {}, "malformed_or_duplicate_completion_row_%d" % index
            values[row[0]] = row[1]
        return values, ""
    except (OSError, csv.Error) as error:
        detail = str(error).replace("\t", " ").replace("\n", " ")
        return {}, "completion_read_error:%s" % detail


def count_success(root: Path) -> int:
    if not root.is_dir():
        return 0
    return sum(1 for path in root.rglob(SUCCESS_NAME) if marker_exists(path))


# ------------------------------------------------------- normalization
def canonical_names(header):
    """Map original columns to collision-safe names beside the canonical four."""
    fields = list(CANONICAL_COLUMNS)
    mapping = {}
    for name in header:
        target = name
        suffix = 0
        while target in fields:
            suffix += 1
            target = "source_%s" % name if suffix == 1 else "source_%s_%d" % (
                name, suffix - 1)
        fields.append(target)
        mapping[name] = target
    return fields, mapping


def canonical_rows(header, rows, study_of, condition_of, independent):
    fields, mapping = canonical_names(header)
    output = []
    for row in rows:
        sample = row["sample_id"]
        record = {
            "sample_id": sample,
            "condition": condition_of(row),
            "study": study_of(row),
            "independent_subset": "1" if sample in independent else "0",
        }
        for name in header:
            record[mapping[name]] = row.get(name, "")
        output.append(record)
    return fields, mapping, output


# ----------------------------------------------------------------- audit
def parse_expected_conditions(raw):
    """Explicit expected condition counts; there is no implicit skip."""
    if not raw or not raw.strip():
        raise SealError("--expected-conditions is required and must not be empty")
    counts = {}
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise SealError(
                "--expected-conditions entries must be NAME=COUNT: %r" % item)
        name, value = item.split("=", 1)
        name = name.strip()
        if not name or name in counts:
            raise SealError("blank or duplicate expected condition: %r" % item)
        try:
            counts[name] = int(value)
        except ValueError:
            raise SealError("non-integer expected condition count: %r" % item)
        if counts[name] < 0:
            raise SealError("negative expected condition count: %r" % item)
    if not counts:
        raise SealError("--expected-conditions is empty")
    return counts


def build_adapter(args):
    if args.cohort not in ADAPTERS:
        raise SealError(
            "unknown cohort %r; expected one of %s"
            % (args.cohort, ", ".join(sorted(ADAPTERS))))
    adapter = dict(ADAPTERS[args.cohort])
    if args.condition_column:
        adapter["condition_column"] = args.condition_column
    if args.study_column or args.study_value:
        adapter["study_column"] = args.study_column
        adapter["study_value"] = args.study_value
    if args.batch_column is not None:
        adapter["batch_column"] = args.batch_column or None
    if args.completion_has_study is not None:
        adapter["completion_has_study"] = args.completion_has_study
    if args.covariate_column is not None:
        adapter["covariates"] = tuple(args.covariate_column)
    if args.source_seal_manifest_member:
        adapter["seal_manifest_member"] = args.source_seal_manifest_member
    if args.source_seal_independent_member:
        adapter["seal_independent_member"] = args.source_seal_independent_member
    if args.source_success_schema:
        adapter["source_success_schema"] = args.source_success_schema
    if args.independent_header_schema:
        adapter["independent_header_schema"] = args.independent_header_schema
    if args.production_only_column is not None:
        adapter["production_only_columns"] = tuple(args.production_only_column)
    if args.source_audit_job is not None:
        adapter["source_audit_job"] = args.source_audit_job
    if bool(adapter["study_column"]) == bool(adapter["study_value"]):
        raise SealError(
            "ambiguous study configuration: give exactly one of a study column "
            "or a fixed study value")
    if not adapter["condition_column"]:
        raise SealError("a condition column is required")
    if adapter["source_success_schema"] not in SOURCE_SUCCESS_SCHEMAS:
        raise SealError(
            "unknown --source-success-schema %r; expected one of %s"
            % (adapter["source_success_schema"], ", ".join(SOURCE_SUCCESS_SCHEMAS)))
    if adapter["independent_header_schema"] not in INDEPENDENT_HEADER_SCHEMAS:
        raise SealError(
            "unknown --independent-header-schema %r; expected one of %s"
            % (adapter["independent_header_schema"],
               ", ".join(INDEPENDENT_HEADER_SCHEMAS)))
    return adapter


def prepare_outdir(outdir: Path) -> Path:
    """Make the output directory reusable without ever trusting stale authority.

    An interrupted `write_atomic` can leave `<member>.tmp` behind. Those are
    recognized and removed; anything else unexpected fails closed. Stale
    authority is dropped before validation so a failed rerun cannot leave an
    older seal looking current.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    # Authority dies first. Whatever else is wrong with this directory, an
    # older SUCCESS must not survive a run that is about to refuse it.
    for stale in (outdir / SUCCESS_NAME, outdir / CHECKSUM_NAME):
        remove_if_present(stale)
    present = sorted(item.name for item in outdir.iterdir())
    for name in present:
        if name in RECOVERABLE_TEMPORARIES:
            remove_if_present(outdir / name)
    unexpected = [name for name in present
                  if name not in PERMITTED_OUTDIR_MEMBERS]
    if unexpected:
        raise SealError(
            "--outdir holds unexpected member(s): %s" % ", ".join(unexpected))
    return outdir / IN_PROGRESS_NAME


def audit(args, adapter, paths, outdir: Path):
    manifest_path = paths["manifest"]
    independent_path = paths["independent_manifest"]
    seal_root = paths["source_seal"]
    state_dir = paths["state_dir"]
    results_root = paths["results_root"]
    qc_root = paths["qc_root"]
    scratch_root = paths["scratch_root"]

    expected_conditions = parse_expected_conditions(args.expected_conditions)
    frozen = FROZEN_COHORTS.get(args.cohort)
    if frozen and args.expected_samples == frozen[0]:
        if expected_conditions != frozen[1]:
            raise SealError(
                "frozen %s condition counts are %s; this run was given %s"
                % (args.cohort, frozen[1], expected_conditions))

    inventory = source_seal_inventory(seal_root, adapter)
    check_manifest_identity(seal_root, adapter, manifest_path, independent_path)
    check_source_success(seal_root, adapter, args.cohort, args.expected_samples,
                         args.expected_independent, expected_conditions)

    header, rows = read_table(manifest_path)
    independent_header, independent_manifest_rows, independent_provenance = (
        read_independent_table(independent_path, adapter))
    condition_column = adapter["condition_column"]
    study_column = adapter["study_column"]
    production_required = ["sample_id", condition_column]
    if study_column:
        production_required.append(study_column)
    if adapter["batch_column"]:
        production_required.append(adapter["batch_column"])
    for name in adapter["covariates"]:
        if name not in production_required:
            production_required.append(name)
    # Batch provenance is production-only: the sealed Yachida independent
    # manifest genuinely does not carry it, and requiring it there would be a
    # false failure rather than a stronger check.
    production_only = set(adapter["production_only_columns"])
    independent_required = [name for name in production_required
                            if name not in production_only]
    for table_name, table_header, needed in (
            (manifest_path.name, header, production_required),
            (independent_path.name, independent_header, independent_required)):
        missing = [name for name in needed if name not in table_header]
        if missing:
            raise SealError(
                "%s lacks required column(s): %s"
                % (table_name, ", ".join(missing)))

    def study_of(row):
        return row[study_column] if study_column else adapter["study_value"]

    def condition_of(row):
        return row[condition_column]

    identifiers = [row["sample_id"] for row in rows]
    if any(not value.strip() for value in identifiers):
        raise SealError("production manifest has a blank sample_id")
    if len(set(identifiers)) != len(identifiers):
        raise SealError("production manifest has duplicate sample identifiers")
    if len(rows) != args.expected_samples:
        raise SealError(
            "expected %d samples; found %d" % (args.expected_samples, len(rows)))

    observed_conditions = Counter(condition_of(row) for row in rows)
    if dict(observed_conditions) != expected_conditions:
        raise SealError(
            "unexpected condition counts: observed %s; expected %s"
            % (dict(observed_conditions), expected_conditions))

    independent_ids = [row["sample_id"] for row in independent_manifest_rows]
    independent = set(independent_ids)
    if len(independent) != len(independent_ids):
        raise SealError("independent manifest has duplicate sample identifiers")
    if len(independent) != args.expected_independent:
        raise SealError(
            "expected %d independent-subset samples; found %d"
            % (args.expected_independent, len(independent)))
    if not independent.issubset(set(identifiers)):
        raise SealError("independent manifest is not nested in the production manifest")
    if args.expected_independent == FROZEN_INDEPENDENT_SAMPLES:
        balance = Counter(condition_of(row) for row in independent_manifest_rows)
        if dict(balance) != FROZEN_SUBSET_BALANCE:
            raise SealError(
                "independent subset is not balanced %s; observed %s"
                % (FROZEN_SUBSET_BALANCE, dict(balance)))

    # The independent manifest may legitimately carry deterministic-selection
    # provenance columns the production manifest does not have. Every shared
    # production column must still agree exactly, row by row.
    production_by_id = {row["sample_id"]: row for row in rows}
    # A translated pilot column is compared against the production column it
    # was inherited from, so the historical first occurrence is still verified.
    comparable = []
    for name in independent_header:
        origin = independent_provenance.get(name, name)
        if origin in header:
            comparable.append((name, origin))
    for row in independent_manifest_rows:
        source = production_by_id[row["sample_id"]]
        for name, origin in comparable:
            if row.get(name) != source.get(origin):
                raise SealError(
                    "independent manifest row for %s has %s=%r but the "
                    "production manifest has %s=%r"
                    % (row["sample_id"], name, row.get(name), origin,
                       source.get(origin)))

    manifest_fields, manifest_mapping, manifest_canonical = canonical_rows(
        header, rows, study_of, condition_of, independent)
    independent_fields, _, independent_canonical = canonical_rows(
        independent_header, independent_manifest_rows, study_of, condition_of,
        independent)

    flow = []
    for row in rows:
        sample = row["sample_id"]
        study = study_of(row)
        in_independent = sample in independent
        sample_root = results_root / study / sample
        expected_independent_profiles = (
            INDEPENDENT_PER_SUBSET_SAMPLE if in_independent else 0)
        expected_profiles = (BASELINE_PER_SAMPLE + expected_independent_profiles
                             + COMMUNITY_PER_SAMPLE)
        profiles = sample_root / "profiles"
        observed_baseline = count_success(profiles / "baseline")
        observed_independent = count_success(profiles / "independent")
        observed_community = count_success(profiles / "community")
        observed_profiles = (
            observed_baseline + observed_independent + observed_community)

        completion, completion_read_error = read_completion(
            sample_root / "sample_completion.tsv")
        expected_completion = {
            "sample_id": sample,
            "condition": condition_of(row),
            "independent_subset": "1" if in_independent else "0",
            "expected_profiles": str(expected_profiles),
            "observed_profiles": str(expected_profiles),
            "community_design_rows": str(COMMUNITY_PER_SAMPLE),
            "independent_design_rows": str(expected_independent_profiles),
        }
        if adapter["completion_has_study"]:
            expected_completion["study"] = study
        mismatches = [
            "%s:%s!=%s" % (field, completion.get(field, "<missing>"), expected)
            for field, expected in sorted(expected_completion.items())
            if completion.get(field) != expected
        ]
        completion_error = completion_read_error or ";".join(mismatches)

        receipt_ok, receipt_files, receipt_error = verify_receipt(
            state_dir / "samples" / ("%s.retained_outputs.tsv" % sample),
            scratch_root / sample,
            (sample_root, qc_root / study / sample),
        )
        checks = {
            "verified_marker": nonempty_file(
                state_dir / "samples" / ("%s.verified" % sample)),
            "retained_output_receipt": receipt_ok,
            "input_provenance": nonempty_file(
                state_dir / "samples" / ("%s.input_provenance.tsv" % sample)),
            "sample_success": marker_exists(sample_root / SUCCESS_NAME),
            "baseline_profile_count": observed_baseline == BASELINE_PER_SAMPLE,
            "independent_profile_count": (
                observed_independent == expected_independent_profiles),
            "community_profile_count": observed_community == COMMUNITY_PER_SAMPLE,
            "completion_table": not completion_error,
            "manifest_independent_flag": True,
        }
        declared = row.get("independent_subset")
        if declared is not None:
            checks["manifest_independent_flag"] = (
                declared == ("1" if in_independent else "0"))
        record = {
            "sample_id": sample, "study": study, "condition": condition_of(row),
            "independent_subset": int(in_independent),
            "batch_id": row[adapter["batch_column"]] if adapter["batch_column"] else "",
            "expected_baseline_profiles": BASELINE_PER_SAMPLE,
            "observed_baseline_profiles": observed_baseline,
            "expected_independent_profiles": expected_independent_profiles,
            "observed_independent_profiles": observed_independent,
            "expected_community_profiles": COMMUNITY_PER_SAMPLE,
            "observed_community_profiles": observed_community,
            "expected_profiles": expected_profiles,
            "observed_profiles": observed_profiles,
            "retained_output_files": receipt_files,
            "status": "PASS" if all(checks.values()) else "FAIL",
            "failure_reasons": ";".join(
                name for name in SAMPLE_CHECKS if not checks[name]),
            "receipt_error": receipt_error, "completion_error": completion_error,
        }
        for name in SAMPLE_CHECKS:
            record[name] = int(checks[name])
        flow.append(record)

    # Covariates are read through the explicit original-to-canonical mapping.
    # Every audited field is named for a real column of the canonical manifest,
    # so a manifest covariate that collides with a canonical name is reported
    # from its own `source_<name>` column instead of being shadowed by the
    # normalized one.
    covariates = []
    audited = list(CANONICAL_COLUMNS[1:])
    for field in adapter["covariates"]:
        column = manifest_mapping[field]
        if column not in audited:
            audited.append(column)
    for column in audited:
        values = [str(row.get(column, "")).strip() for row in manifest_canonical]
        covariates.append({
            "field": column,
            "missing": sum(1 for value in values if not value),
            "present": sum(1 for value in values if value),
            "distinct_nonmissing": len({value for value in values if value}),
        })

    payloads = {
        "sample_flow.tsv": tsv_bytes(SAMPLE_FLOW_FIELDS, flow),
        "covariate_audit.tsv": tsv_bytes(COVARIATE_FIELDS, covariates),
        "production_manifest.tsv": tsv_bytes(manifest_fields, manifest_canonical),
        "production_manifest.independent.tsv": tsv_bytes(
            independent_fields, independent_canonical),
        "source_seal_inventory.tsv": tsv_bytes(INVENTORY_FIELDS, inventory),
    }
    for name in SEAL_MEMBERS:
        write_atomic(outdir / name, payloads[name])

    failed = [row["sample_id"] for row in flow if row["status"] != "PASS"]
    if failed:
        raise SealError(
            "%d sample(s) are incomplete; see %s (first: %s)"
            % (len(failed), outdir / "sample_flow.tsv", failed[0]))

    totals = {
        "baseline": sum(row["observed_baseline_profiles"] for row in flow),
        "independent": sum(row["observed_independent_profiles"] for row in flow),
        "community": sum(row["observed_community_profiles"] for row in flow),
    }
    expected_totals = {
        "baseline": len(rows) * BASELINE_PER_SAMPLE,
        "independent": len(independent) * INDEPENDENT_PER_SUBSET_SAMPLE,
        "community": len(rows) * COMMUNITY_PER_SAMPLE,
    }
    if totals != expected_totals:
        raise SealError(
            "aggregate profile topology mismatch: observed %s; expected %s"
            % (totals, expected_totals))
    return len(rows), len(independent), totals, inventory


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--independent-manifest", required=True, type=Path)
    parser.add_argument("--source-seal", required=True, type=Path,
                        help="Native production seal, read-only provenance.")
    parser.add_argument("--state-dir", required=True, type=Path)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--qc-root", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--expected-samples", required=True, type=int)
    parser.add_argument("--expected-independent", type=int,
                        default=FROZEN_INDEPENDENT_SAMPLES)
    parser.add_argument("--expected-conditions", required=True,
                        help="NAME=COUNT,... exact production condition counts.")
    parser.add_argument("--condition-column")
    parser.add_argument("--study-column")
    parser.add_argument("--study-value")
    parser.add_argument("--batch-column")
    parser.add_argument("--covariate-column", action="append")
    parser.add_argument("--source-seal-manifest-member")
    parser.add_argument("--source-seal-independent-member")
    parser.add_argument("--source-success-schema",
                        choices=list(SOURCE_SUCCESS_SCHEMAS))
    parser.add_argument("--independent-header-schema",
                        choices=list(INDEPENDENT_HEADER_SCHEMAS),
                        help="unique (default for Feng/Zeller) or the strict "
                             "yachida_historical_duplicate_selection adapter.")
    parser.add_argument("--production-only-column", action="append",
                        help="Column required in the production manifest but "
                             "not in the independent manifest; repeatable.")
    parser.add_argument("--source-audit-job")
    completion = parser.add_mutually_exclusive_group()
    completion.add_argument("--completion-has-study", dest="completion_has_study",
                            action="store_true", default=None)
    completion.add_argument("--no-completion-study", dest="completion_has_study",
                            action="store_false")
    return parser.parse_args(argv)


def run(argv=None) -> None:
    args = parse_args(argv)
    if args.expected_samples <= 0 or args.expected_independent <= 0:
        raise SealError("expected sample counts must be positive")
    adapter = build_adapter(args)
    # Nothing below this line may run before the paths are proven safe.
    paths = resolve_paths(args)
    outdir = paths["outdir"]
    in_progress = prepare_outdir(outdir)
    write_atomic(in_progress, (
        "seal_contract\t%s\ncohort\t%s\nstatus\tIN_PROGRESS\nmanifest\t%s\n"
        % (SEAL_CONTRACT, args.cohort, paths["manifest"])).encode("utf-8"))

    samples, independent_samples, totals, inventory = audit(
        args, adapter, paths, outdir)

    # The seal must still be the seal that was verified at the start.
    reverify_source_seal(paths["source_seal"], adapter, paths["manifest"],
                         paths["independent_manifest"], inventory)

    success = (
        "seal_contract\t%s\n"
        "cohort\t%s\n"
        "samples\t%d\n"
        "independent_samples\t%d\n"
        "baseline_profiles\t%d\n"
        "independent_profiles\t%d\n"
        "community_profiles\t%d\n"
        "status\tPASS\n"
    ) % (SEAL_CONTRACT, args.cohort, samples, independent_samples,
         totals["baseline"], totals["independent"], totals["community"])
    success_bytes = success.encode("utf-8")
    lines = []
    for name in SEAL_MEMBERS:
        lines.append("%s  %s\n" % (digest(outdir / name), name))
    lines.append("%s  %s\n" % (digest_bytes(success_bytes), SUCCESS_NAME))
    # The checksum lands first and SUCCESS last, so publication authority is
    # never asserted before the manifest that covers it exists.
    write_atomic(outdir / CHECKSUM_NAME, "".join(lines).encode("utf-8"))
    write_atomic(outdir / SUCCESS_NAME, success_bytes)
    remove_if_present(in_progress)
    print("[PASS] Sealed %s upstream (%s): %d samples, %d independent"
          % (args.cohort, SEAL_CONTRACT, samples, independent_samples))
    print("[INFO] Seal: %s" % outdir)


def main() -> None:
    try:
        run()
    except (SealError, OSError, KeyError, ValueError) as error:
        raise SystemExit("[ERROR] %s" % error)


if __name__ == "__main__":
    main()
