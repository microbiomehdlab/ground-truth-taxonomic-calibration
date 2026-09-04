#!/usr/bin/env python3
"""Validate and seal a canonical paired dose-response v2 TSV."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from collections import Counter
from pathlib import Path


SCHEMA_VERSION = "paired-dose-response-v2.1"
REQUIRED = [
    "schema_version", "cohort", "study", "sample_id", "condition",
    "analysis_population", "target_label", "target_taxon", "assembly_arm",
    "profiler", "profile_id", "baseline_profile_id", "spike_fraction_total",
    "spike_fraction_target", "implanted_read_pairs_target", "native_abundance",
    "native_unit", "abundance_fraction", "detected_native_nonzero",
    "source_profile", "source_design", "include", "exclusion_reason",
]
PROFILER_UNITS = {
    "kraken2_bracken": "fraction_total_reads",
    "metaphlan4": "relative_abundance_pct",
}
POPULATIONS = {"community", "independent"}
ARMS = {"original", "clean", "not_applicable"}


def numeric(row: dict[str, str], field: str, line: int) -> float:
    try:
        value = float(row[field])
    except ValueError as error:
        raise ValueError(f"line {line}: {field} is not numeric") from error
    if not math.isfinite(value):
        raise ValueError(f"line {line}: {field} must be finite")
    return value


def integer(row: dict[str, str], field: str, line: int) -> int:
    value = numeric(row, field, line)
    if value < 0 or not value.is_integer():
        raise ValueError(f"line {line}: {field} must be a non-negative integer")
    return int(value)


def key(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row[name] for name in (
        "cohort", "sample_id", "analysis_population", "target_label",
        "profiler", "assembly_arm",
    ))


def validate_row(row: dict[str, str], line: int) -> tuple[float, bool]:
    blank = [field for field in REQUIRED if field != "exclusion_reason" and not row[field].strip()]
    if blank:
        raise ValueError(f"line {line}: blank required fields: {', '.join(blank)}")
    if row["schema_version"] != SCHEMA_VERSION:
        raise ValueError(f"line {line}: unsupported schema_version")
    if row["analysis_population"] not in POPULATIONS:
        raise ValueError(f"line {line}: invalid analysis_population")
    if row["assembly_arm"] not in ARMS:
        raise ValueError(f"line {line}: invalid assembly_arm")
    expected_unit = PROFILER_UNITS.get(row["profiler"])
    if expected_unit is None or row["native_unit"] != expected_unit:
        raise ValueError(f"line {line}: profiler/native_unit mismatch")
    if row["include"] not in {"0", "1"}:
        raise ValueError(f"line {line}: include must be 0 or 1")
    if (row["include"] == "0") != bool(row["exclusion_reason"].strip()):
        raise ValueError(f"line {line}: exclusion reason/include mismatch")

    total = numeric(row, "spike_fraction_total", line)
    target = numeric(row, "spike_fraction_target", line)
    pairs = integer(row, "implanted_read_pairs_target", line)
    native = numeric(row, "native_abundance", line)
    fraction = numeric(row, "abundance_fraction", line)
    if not (0 <= target <= total < 1):
        raise ValueError(f"line {line}: invalid spike fractions")
    if native < 0 or not (0 <= fraction <= 1.000001):
        raise ValueError(f"line {line}: invalid abundance")
    expected_fraction = native if expected_unit == "fraction_total_reads" else native / 100
    if not math.isclose(fraction, expected_fraction, rel_tol=1e-9, abs_tol=1e-12):
        raise ValueError(f"line {line}: abundance_fraction is not a pure unit conversion")
    detected = fraction > 0
    if row["detected_native_nonzero"] != ("1" if detected else "0"):
        raise ValueError(f"line {line}: detected_native_nonzero disagrees with abundance")

    baseline = target == 0
    if baseline:
        if total != 0 or pairs != 0 or row["profile_id"] != row["baseline_profile_id"]:
            raise ValueError(f"line {line}: malformed zero-dose row")
    elif pairs == 0 or row["profile_id"] == row["baseline_profile_id"]:
        raise ValueError(f"line {line}: malformed positive-dose row")
    return target, baseline


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if not args.input.is_file() or args.input.stat().st_size == 0:
        raise SystemExit(f"[ERROR] missing or empty input: {args.input}")

    try:
        with args.input.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != REQUIRED:
                raise ValueError("header does not exactly match the canonical contract")
            rows = list(reader)
        if not rows:
            raise ValueError("canonical table contains no data rows")

        analytical_ids: set[tuple[str, ...]] = set()
        baselines: Counter[tuple[str, ...]] = Counter()
        positive: Counter[tuple[str, ...]] = Counter()
        included = excluded = 0
        for line, row in enumerate(rows, start=2):
            _, is_baseline = validate_row(row, line)
            analytical_id = key(row) + (row["profile_id"],)
            if analytical_id in analytical_ids:
                raise ValueError(f"line {line}: duplicate analytical row")
            analytical_ids.add(analytical_id)
            if row["include"] == "1":
                included += 1
                (baselines if is_baseline else positive)[key(row)] += 1
            else:
                excluded += 1
        for paired_key, count in baselines.items():
            if count != 1:
                raise ValueError(f"expected one baseline for key {paired_key}, found {count}")
        missing = [paired_key for paired_key in positive if baselines[paired_key] != 1]
        if missing:
            raise ValueError(f"positive-dose row lacks unique baseline: {missing[0]}")
    except ValueError as error:
        raise SystemExit(f"[ERROR] {error}") from error

    args.outdir.mkdir(parents=True, exist_ok=True)
    report = args.outdir / "canonical_input_validation.tsv"
    report.write_text(
        "metric\tvalue\n"
        f"schema_version\t{SCHEMA_VERSION}\n"
        f"rows\t{len(rows)}\n"
        f"included_rows\t{included}\n"
        f"excluded_rows\t{excluded}\n"
        f"baseline_rows\t{sum(baselines.values())}\n"
        f"positive_dose_rows\t{sum(positive.values())}\n"
        f"biological_samples\t{len({(r['cohort'], r['sample_id']) for r in rows})}\n"
        "status\tPASS\n",
        encoding="utf-8",
    )
    checksum = args.outdir / "canonical_input.sha256"
    checksum.write_text(
        f"{sha256(args.input)}  {args.input.resolve()}\n"
        f"{sha256(report)}  {report.resolve()}\n",
        encoding="utf-8",
    )
    (args.outdir / "SUCCESS").write_text(
        f"schema_version\t{SCHEMA_VERSION}\nrows\t{len(rows)}\nstatus\tPASS\n",
        encoding="utf-8",
    )
    print(f"[PASS] Canonical input validated: {len(rows)} rows")
    print(f"[INFO] {report}")


if __name__ == "__main__":
    main()
