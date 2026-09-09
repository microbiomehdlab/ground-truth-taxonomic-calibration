#!/usr/bin/env python3
"""Fail closed unless the frozen downstream analysis policy is complete."""
import argparse
import csv
import hashlib
from pathlib import Path

REQUIRED = {
    "policy_version": "1", "native_detection_rule": "abundance_fraction_gt_0",
    "quantitative_primary_population": "unconditional", "quantitative_nondetection_value": "0",
    "paired_endpoint_pseudocount": "none", "paired_endpoint_transform": "none",
    "negative_baseline_adjusted_response": "retain",
    "response_ratio_domain": "positive_implanted_fraction_only",
    "baseline_zero_handling": "retain", "detected_only_role": "secondary_descriptive",
    "biomarker_log2_pseudocount_fraction": "1e-8",
    "artificial_biomarker_primary_q": "0.05", "artificial_biomarker_sensitivity_q": "0.10",
    "disease_biomarker_primary_q": "0.05", "disease_biomarker_sensitivity_q": "0.10",
    "species_closed_renormalization_role": "sensitivity_only",
    "primary_profiler_scale": "profiler_native_fraction",
}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--policy", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    with args.policy.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or set(rows[0]) != {"field", "value"}:
        raise SystemExit("[ERROR] policy must contain exactly field and value columns")
    fields = [row["field"] for row in rows]
    if len(fields) != len(set(fields)):
        raise SystemExit("[ERROR] duplicate analysis-policy field")
    policy = {row["field"]: row["value"] for row in rows}
    for field, expected in REQUIRED.items():
        if policy.get(field) != expected:
            raise SystemExit(f"[ERROR] frozen policy mismatch: {field}")
    for field in ("detection_multiplicity_family", "artificial_biomarker_multiplicity_family",
                  "disease_biomarker_multiplicity_family", "cross_cohort_multiplicity_family"):
        if not policy.get(field):
            raise SystemExit(f"[ERROR] empty analysis-policy field: {field}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    copied = args.outdir / "analysis_policy.tsv"
    copied.write_bytes(args.policy.read_bytes())
    digest = hashlib.sha256(copied.read_bytes()).hexdigest()
    (args.outdir / "analysis_policy.sha256").write_text(
        f"{digest}  analysis_policy.tsv\n", encoding="utf-8")
    (args.outdir / "SUCCESS").write_text("analysis_policy\tPASS\n", encoding="utf-8")
    print(f"[PASS] Frozen analysis policy validated: {copied}")

if __name__ == "__main__":
    main()
