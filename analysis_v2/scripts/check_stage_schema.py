#!/usr/bin/env python3
"""Decide whether a completed pipeline stage directory may be reused.

A stage that finished before its producer gained new columns still carries a
valid `SUCCESS` marker, so `SUCCESS` alone is not evidence that the directory is
usable by the current code. This checker states, per stage, the exact files,
TSV header columns and validation metrics a reusable directory must contain.

Headers are parsed as tab-separated fields and compared exactly. A substring or
unquoted `grep` test would accept `primary_reference_type` where
`primary_row_reference_type` is required, which is precisely the confusion this
gate exists to catch.

Exit status: 0 reusable; 1 not reusable (the reason is printed); 2 usage error.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

COMPARISON = "reference_comparison"
AUDIT = "metaphlan_genome_size_residual_audit"

# The explicit per-arm reference provenance that compare_target_recovery_
# references.py began emitting on 20 September 2026. An older comparison has
# only the ambiguous `*_reference_type` alias and cannot feed the audit.
COMPARISON_COLUMNS = [
    "primary_row_reference_type", "sensitivity_row_reference_type",
    "primary_selected_reference_type", "sensitivity_selected_reference_type",
    "primary_reference_scale", "sensitivity_reference_scale",
    "primary_selected_estimand", "sensitivity_selected_estimand",
]
COMPARISON_METRICS = [
    "primary_selection", "sensitivity_selection", "selected_estimand_contract",
]

AUDIT_COLUMNS = [
    "reference_type", "row_reference_type", "selected_reference_type",
    "reference_scale",
]
AUDIT_METRICS = [
    "selected_reference_sensitivity", "selected_reference_primary",
    "observed_cohorts", "absolute_relative_error_policy", "audit_status",
]

STAGES = {
    COMPARISON: {
        "files": ["SUCCESS", "target_recovery_reference_comparison.tsv",
                  "reference_comparison_validation.tsv"],
        "headers": {"target_recovery_reference_comparison.tsv": COMPARISON_COLUMNS},
        "metrics": {"reference_comparison_validation.tsv": COMPARISON_METRICS},
    },
    AUDIT: {
        "files": ["SUCCESS", "metaphlan_genome_size_target_level.tsv",
                  "metaphlan_genome_size_regression.tsv",
                  "metaphlan_genome_size_paired_change.tsv",
                  "metaphlan_genome_size_audit_validation.tsv",
                  "metaphlan_genome_size_audit.sha256", "DEVELOPMENT_ONLY.txt"],
        "headers": {"metaphlan_genome_size_target_level.tsv": AUDIT_COLUMNS},
        "metrics": {"metaphlan_genome_size_audit_validation.tsv": AUDIT_METRICS},
    },
}


def header_fields(path: Path) -> list[str]:
    """Exact tab-separated header fields; never a substring test."""
    with path.open(newline="", encoding="utf-8") as handle:
        first = handle.readline()
    if not first.strip():
        raise ValueError(f"{path.name} has no header line")
    return first.rstrip("\r\n").split("\t")


def metric_names(path: Path) -> set[str]:
    """Metric names with a nonempty value, from a two-column metric/value TSV."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if "metric" not in fields or "value" not in fields:
            raise ValueError(
                f"{path.name} is not a metric/value table (header: "
                + ", ".join(fields) + ")")
        return {(row.get("metric") or "").strip() for row in reader
                if (row.get("metric") or "").strip()
                and (row.get("value") or "").strip()}


def incompatibility(stage: str, directory: Path) -> str | None:
    """Return the first reason the directory may not be reused, else None."""
    spec = STAGES[stage]
    if not directory.is_dir():
        return f"{directory} is not a directory"
    for name in spec["files"]:
        path = directory / name
        if not path.is_file():
            return f"missing {name}"
        if not path.stat().st_size:
            return f"empty {name}"
    for name, required in spec["headers"].items():
        try:
            present = set(header_fields(directory / name))
        except (OSError, ValueError) as error:
            return str(error)
        missing = [column for column in required if column not in present]
        if missing:
            return (f"{name} lacks column(s) " + ", ".join(missing)
                    + "; it predates the explicit reference-provenance schema")
    for name, required in spec["metrics"].items():
        try:
            present = metric_names(directory / name)
        except (OSError, ValueError) as error:
            return str(error)
        missing = [metric for metric in required if metric not in present]
        if missing:
            return (f"{name} lacks metric(s) " + ", ".join(missing)
                    + "; it predates the current selection contract")
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=sorted(STAGES), required=True)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--quiet", action="store_true",
                        help="Suppress the success line; reasons still print.")
    args = parser.parse_args()
    reason = incompatibility(args.stage, args.directory)
    if reason is None:
        if not args.quiet:
            print(f"[OK] {args.stage} directory is schema-compatible: "
                  f"{args.directory}")
        raise SystemExit(0)
    print(f"{reason}", file=sys.stdout)
    raise SystemExit(1)


if __name__ == "__main__":
    main()
