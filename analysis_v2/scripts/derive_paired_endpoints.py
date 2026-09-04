#!/usr/bin/env python3
"""Derive baseline-adjusted paired endpoints from canonical v2 input."""

from __future__ import annotations

import argparse
import csv
import hashlib
import subprocess
import sys
from pathlib import Path


KEY_FIELDS = [
    "cohort", "sample_id", "analysis_population", "target_label", "profiler",
    "assembly_arm",
]
OUTPUT_FIELDS = [
    "schema_version", "cohort", "study", "sample_id", "condition",
    "analysis_population", "target_label", "target_taxon", "assembly_arm",
    "profiler", "profile_id", "baseline_profile_id", "spike_fraction_total",
    "spike_fraction_target", "implanted_read_pairs_target", "native_unit",
    "baseline_native_abundance", "observed_native_abundance",
    "baseline_abundance_fraction", "observed_abundance_fraction",
    "baseline_retained_after_dilution", "read_proportional_reference",
    "recovered_spike_signal", "response_ratio", "signed_reference_error",
    "absolute_reference_error", "baseline_detected_native_nonzero",
    "observed_detected_native_nonzero", "source_baseline_profile",
    "source_profile", "source_design",
]


def key(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row[field] for field in KEY_FIELDS)


def render(value: float) -> str:
    return format(value, ".17g")


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
    args.outdir.mkdir(parents=True, exist_ok=True)

    validator = Path(__file__).with_name("validate_canonical_input.py")
    validation_dir = args.outdir / "input_validation"
    subprocess.run([
        sys.executable, str(validator), "--input", str(args.input),
        "--outdir", str(validation_dir),
    ], check=True)

    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    included = [row for row in rows if row["include"] == "1"]
    excluded = [row for row in rows if row["include"] == "0"]
    baselines = {
        key(row): row for row in included
        if float(row["spike_fraction_target"]) == 0
    }

    endpoints: list[dict[str, str]] = []
    for row in included:
        target_fraction = float(row["spike_fraction_target"])
        if target_fraction == 0:
            continue
        baseline = baselines[key(row)]
        total_fraction = float(row["spike_fraction_total"])
        baseline_fraction = float(baseline["abundance_fraction"])
        observed_fraction = float(row["abundance_fraction"])
        retained = (1 - total_fraction) * baseline_fraction
        reference = retained + target_fraction
        recovered = observed_fraction - retained
        ratio = recovered / target_fraction
        signed_error = observed_fraction - reference

        endpoint = {
            field: row[field] for field in (
                "schema_version", "cohort", "study", "sample_id", "condition",
                "analysis_population", "target_label", "target_taxon",
                "assembly_arm", "profiler", "profile_id", "baseline_profile_id",
                "spike_fraction_total", "spike_fraction_target",
                "implanted_read_pairs_target", "native_unit", "source_profile",
                "source_design",
            )
        }
        endpoint.update({
            "baseline_native_abundance": baseline["native_abundance"],
            "observed_native_abundance": row["native_abundance"],
            "baseline_abundance_fraction": render(baseline_fraction),
            "observed_abundance_fraction": render(observed_fraction),
            "baseline_retained_after_dilution": render(retained),
            "read_proportional_reference": render(reference),
            "recovered_spike_signal": render(recovered),
            "response_ratio": render(ratio),
            "signed_reference_error": render(signed_error),
            "absolute_reference_error": render(abs(signed_error)),
            "baseline_detected_native_nonzero": baseline["detected_native_nonzero"],
            "observed_detected_native_nonzero": row["detected_native_nonzero"],
            "source_baseline_profile": baseline["source_profile"],
        })
        endpoints.append(endpoint)

    output = args.outdir / "paired_endpoints.tsv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(endpoints)

    exclusions = args.outdir / "excluded_canonical_rows.tsv"
    with exclusions.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(excluded)

    summary = args.outdir / "endpoint_derivation_summary.tsv"
    summary.write_text(
        "metric\tvalue\n"
        f"canonical_rows\t{len(rows)}\n"
        f"included_baseline_rows\t{len(baselines)}\n"
        f"derived_positive_dose_rows\t{len(endpoints)}\n"
        f"excluded_rows\t{len(excluded)}\n"
        "status\tPASS\n",
        encoding="utf-8",
    )
    checksum = args.outdir / "endpoint_derivation.sha256"
    paths = [args.input, output, exclusions, summary,
             validation_dir / "canonical_input_validation.tsv"]
    checksum.write_text(
        "".join(f"{sha256(path)}  {path.resolve()}\n" for path in paths),
        encoding="utf-8",
    )
    (args.outdir / "SUCCESS").write_text(
        f"endpoint_rows\t{len(endpoints)}\nstatus\tPASS\n", encoding="utf-8"
    )
    print(f"[PASS] Derived {len(endpoints)} paired positive-dose endpoints")
    print(f"[INFO] {output}")


if __name__ == "__main__":
    main()
