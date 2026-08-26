#!/usr/bin/env python3
"""Freeze metadata eligibility for a public CRC shotgun cohort."""

from __future__ import annotations

import argparse
import csv
import hashlib
import pathlib


ALLOWED_CONDITIONS = ("Control", "Adenoma", "CRC")


def checksum(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_table(path: pathlib.Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=pathlib.Path)
    parser.add_argument("--output", required=True, type=pathlib.Path)
    parser.add_argument("--study-name", required=True)
    parser.add_argument("--eligibility-version", default="crc-public-cohort-metadata-eligibility-v1")
    args = parser.parse_args()

    with args.metadata.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    required = {
        "Name", "Sample Source ID", "Organism", "Body Site", "Study name",
        "Study condition",
    }
    if missing := required.difference(fields):
        raise SystemExit(f"[ERROR] metadata lacks fields: {sorted(missing)}")

    name_counts: dict[str, int] = {}
    source_counts: dict[str, int] = {}
    for row in rows:
        name_counts[row["Name"]] = name_counts.get(row["Name"], 0) + 1
        source_counts[row["Sample Source ID"]] = source_counts.get(row["Sample Source ID"], 0) + 1

    eligible: list[dict[str, str]] = []
    excluded: list[dict[str, str]] = []
    output_fields = fields + ["metadata_eligibility", "metadata_exclusion_reason", "eligibility_version"]
    for row in rows:
        reasons = []
        name = row["Name"].strip()
        source_id = row["Sample Source ID"].strip()
        if not name:
            reasons.append("blank_sample_id")
        elif name_counts[row["Name"]] != 1:
            reasons.append("duplicate_sample_id")
        if not source_id:
            reasons.append("blank_source_id")
        elif source_counts[row["Sample Source ID"]] != 1:
            reasons.append("duplicate_source_id")
        if row["Organism"].strip().casefold() != "homo sapiens":
            reasons.append("not_human")
        if row["Body Site"].strip().casefold() not in {"feces", "faeces", "stool"}:
            reasons.append("not_fecal")
        if row["Study name"].strip() != args.study_name:
            reasons.append("wrong_study")
        if row["Study condition"].strip() not in ALLOWED_CONDITIONS:
            reasons.append("unsupported_condition")

        annotated = {
            **row,
            "metadata_eligibility": "EXCLUDED" if reasons else "ELIGIBLE",
            "metadata_exclusion_reason": ";".join(reasons),
            "eligibility_version": args.eligibility_version,
        }
        (excluded if reasons else eligible).append(annotated)

    eligible.sort(key=lambda row: (row["Study condition"], row["Name"]))
    excluded.sort(key=lambda row: row["Name"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    excluded_path = args.output.with_name(f"{args.output.stem}.excluded.tsv")
    provenance_path = args.output.with_name(f"{args.output.stem}.provenance.tsv")
    write_table(args.output, output_fields, eligible)
    write_table(excluded_path, output_fields, excluded)

    counts = {condition: 0 for condition in ALLOWED_CONDITIONS}
    for row in eligible:
        counts[row["Study condition"]] += 1
    provenance_rows = [
        {"key": "eligibility_version", "value": args.eligibility_version},
        {"key": "source_metadata", "value": str(args.metadata)},
        {"key": "source_metadata_sha256", "value": checksum(args.metadata)},
        {"key": "required_study", "value": args.study_name},
        {"key": "allowed_organism", "value": "Homo sapiens"},
        {"key": "allowed_body_sites", "value": "feces,faeces,stool"},
        {"key": "allowed_conditions", "value": ",".join(ALLOWED_CONDITIONS)},
        {"key": "eligible_samples", "value": str(len(eligible))},
        {"key": "excluded_samples", "value": str(len(excluded))},
        *({"key": f"eligible_{condition}", "value": str(counts[condition])} for condition in ALLOWED_CONDITIONS),
        {"key": "technical_failures", "value": "not part of metadata eligibility; retry and track separately"},
        {"key": "outcome_use", "value": "none"},
    ]
    write_table(provenance_path, ["key", "value"], provenance_rows)
    for path in (args.output, excluded_path, provenance_path):
        path.with_suffix(path.suffix + ".sha256").write_text(
            f"{checksum(path)}  {path.name}\n", encoding="utf-8"
        )

    print(f"[PASS] metadata eligibility frozen: {len(eligible)} eligible, {len(excluded)} excluded")
    print(f"[INFO] Counts: {counts}")
    print(f"[INFO] Eligible manifest: {args.output}")


if __name__ == "__main__":
    main()
