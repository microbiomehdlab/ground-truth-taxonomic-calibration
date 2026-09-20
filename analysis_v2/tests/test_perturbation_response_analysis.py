#!/usr/bin/env python3
"""Focused fixture for the streamed perturbation-response analysis."""

from __future__ import annotations

import csv
import math
import subprocess
import tempfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis_v2/scripts/analyze_perturbation_response.py"
FIELDS = [
    "observation_id", "cohort", "study", "sample_id", "condition",
    "analysis_population", "assembly_arm", "profiler", "perturbation_id",
    "baseline_id", "dose_level", "dose_index", "nominal_total_fraction",
    "nominal_target_fraction", "recorded_total_fraction",
    "effective_total_fraction", "implanted_target_count", "implanted_targets",
    "implanted_features", "implanted_fraction_by_target", "feature",
    "baseline_abundance_fraction", "observed_abundance_fraction",
    "target_fraction_for_feature", "expected_abundance_fraction",
    "dilution_retained_baseline", "response_signal", "response_delta",
    "signed_bounded_error", "absolute_bounded_error", "quantitative_log2_ratio",
    "baseline_detected", "observed_detected", "detection_transition",
    "is_direct_target", "reference_kind",
]
NOMINAL = {
    "independent": {1: 0.0001, 2: 0.0005},
    "community": {3: 0.001},
}
FEATURE_BY_TARGET = {"T1": "F1", "T2": "F2"}
BASELINE = {"F1": 0.1, "F2": 0.2, "X": 0.1, "H": 0.0, "D": 0.001}


def make_row(cohort: str, sample: str, population: str, dose_index: int,
             target_fractions: dict[str, float], feature: str, signal: float) -> dict[str, str]:
    total = sum(target_fractions.values())
    target_features = [FEATURE_BY_TARGET[target] for target in target_fractions]
    baseline = BASELINE[feature]
    retained = (1 - total) * baseline
    direct = sum(value for target, value in target_fractions.items()
                 if FEATURE_BY_TARGET[target] == feature)
    expected = retained + direct
    observed = retained + signal
    if observed < 0:
        raise AssertionError("fixture generated negative abundance")
    delta = observed - expected
    targets = ";".join(target_fractions)
    fractions = ";".join(f"{target}={value:.17g}"
                         for target, value in target_fractions.items())
    perturbation = f"{cohort}_{sample}_{population}_{dose_index}_{targets}"
    baseline_detected = int(baseline > 0)
    observed_detected = int(observed > 0)
    return {
        "observation_id": perturbation,
        "cohort": cohort, "study": cohort, "sample_id": sample,
        "condition": "Control", "analysis_population": population,
        "assembly_arm": "original", "profiler": "profiler_a",
        "perturbation_id": perturbation, "baseline_id": f"{cohort}_{sample}_base",
        "dose_level": f"dose_{dose_index:02d}", "dose_index": str(dose_index),
        "nominal_total_fraction": f"{NOMINAL[population][dose_index]:.17g}",
        "nominal_target_fraction": f"{next(iter(target_fractions.values())):.17g}",
        "recorded_total_fraction": f"{total:.17g}",
        "effective_total_fraction": f"{total:.17g}",
        "implanted_target_count": str(len(target_fractions)),
        "implanted_targets": targets, "implanted_features": ";".join(target_features),
        "implanted_fraction_by_target": fractions, "feature": feature,
        "baseline_abundance_fraction": f"{baseline:.17g}",
        "observed_abundance_fraction": f"{observed:.17g}",
        "target_fraction_for_feature": f"{direct:.17g}",
        "expected_abundance_fraction": f"{expected:.17g}",
        "dilution_retained_baseline": f"{retained:.17g}",
        "response_signal": f"{signal:.17g}", "response_delta": f"{delta:.17g}",
        "signed_bounded_error": f"{(delta / (observed + expected) if observed + expected else 0):.17g}",
        "absolute_bounded_error": f"{abs(delta / (observed + expected) if observed + expected else 0):.17g}",
        "quantitative_log2_ratio": (
            f"{math.log2(observed / expected):.17g}" if observed > 0 and expected > 0 else ""
        ),
        "baseline_detected": str(baseline_detected),
        "observed_detected": str(observed_detected),
        "detection_transition": (
            "RETAINED" if baseline_detected and observed_detected else
            "LOST" if baseline_detected else "GAINED" if observed_detected else "STABLE_ABSENT"
        ),
        "is_direct_target": str(int(feature in target_features)),
        "reference_kind": "composition_plus_implant",
    }


def fixture_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for cohort in ("alpha", "beta"):
        sample = f"{cohort}_S1"
        for target in ("T1", "T2"):
            for dose_index, fraction in ((1, 0.0001), (2, 0.0005)):
                for feature in BASELINE:
                    if feature == FEATURE_BY_TARGET[target]:
                        signal = fraction
                    elif feature == "X":
                        signal = fraction * (2 if target == "T1" else 3)
                    elif feature == "H":
                        signal = fraction if target == "T1" else 0.0
                    elif feature == "D" and target == "T1":
                        signal = -(1 - fraction) * BASELINE[feature]
                    else:
                        signal = 0.0
                    rows.append(make_row(cohort, sample, "independent", dose_index,
                                         {target: fraction}, feature, signal))

        # Community dose 3 has the same per-target dose as independent dose 1.
        components = {"T1": 0.0001, "T2": 0.0001}
        for feature, signal in (("F1", 0.0001), ("F2", 0.0001), ("X", 0.0005)):
            rows.append(make_row(cohort, sample, "community", 3, components, feature, signal))
    return rows


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    rows = fixture_rows()
    source = root / "responses.tsv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    parquet = root / "responses.parquet"
    pq.write_table(pa.Table.from_pylist(rows), parquet)

    outputs = []
    for index, input_path in enumerate((source, parquet), start=1):
        output = root / f"out{index}"
        subprocess.run([
            "python3", str(SCRIPT), "--responses", str(input_path),
            "--outdir", str(output),
            # The fixture is a read-proportional table, so every existing
            # numeric expectation below must still hold exactly.
            "--reference-scale", "read_proportional",
        ], check=True)
        outputs.append(output)
        assert (output / "SUCCESS").is_file()
        assert (output / "DEVELOPMENT_ONLY.txt").is_file()
        subprocess.run(["sha256sum", "-c", "perturbation_response.sha256"],
                       cwd=output, check=True, stdout=subprocess.DEVNULL)

    operators = read_tsv(outputs[1] / "response_operator.tsv")
    x_t1 = next(row for row in operators if row["cohort"] == "alpha" and
                row["target_label"] == "T1" and row["feature"] == "X")
    assert math.isclose(float(x_t1["operator_slope"]), 2.0, abs_tol=1e-10)
    diagonal = next(row for row in operators if row["cohort"] == "alpha" and
                    row["target_label"] == "T1" and row["feature"] == "F1")
    assert diagonal["feature_role"] == "implanted_target"
    assert math.isclose(float(diagonal["operator_slope"]), 1.0, abs_tol=1e-10)

    superposition = read_tsv(outputs[1] / "superposition_summary.tsv")
    bystander = next(row for row in superposition if row["cohort"] == "alpha" and
                     row["feature_role"] == "bystander")
    assert float(bystander["operator_prediction_mae"]) < 1e-12
    assert float(bystander["composition_only_mae"]) > 0

    certificates = read_tsv(outputs[1] / "reliability_certificates.tsv")
    f1 = next(row for row in certificates if row["holdout_cohort"] == "NONE" and
              row["feature"] == "F1")
    assert f1["eligible_targets"] == "1"  # Its own T1 perturbations were excluded.
    assert f1["overall_reliability"]
    h = next(row for row in certificates if row["holdout_cohort"] == "NONE" and
             row["feature"] == "H")
    assert float(h["unexpected_appearance_rate"]) > 0

    holdout = read_tsv(outputs[1] / "cohort_holdout_summary.tsv")
    assert {row["holdout_cohort"] for row in holdout} == {"alpha", "beta"}
    assert all(int(row["features_shared"]) > 0 for row in holdout)

    malformed = root / "malformed.tsv"
    malformed.write_text("cohort\nalpha\n", encoding="utf-8")
    failed = subprocess.run([
        "python3", str(SCRIPT), "--responses", str(malformed),
        "--outdir", str(root / "bad"), "--reference-scale", "read_proportional",
    ], text=True, capture_output=True)
    assert failed.returncode != 0
    assert "input missing columns" in failed.stderr

    # No silent default: the reference scale must be stated explicitly.
    no_scale = subprocess.run([
        "python3", str(SCRIPT), "--responses", str(source),
        "--outdir", str(root / "noscale"),
    ], text=True, capture_output=True)
    assert no_scale.returncode != 0
    assert "--reference-scale" in no_scale.stderr

    # profiler_scale needs the migrated builder's columns; the legacy fixture
    # lacks them, so it must fail loudly rather than silently mixing scales.
    wrong_scale = subprocess.run([
        "python3", str(SCRIPT), "--responses", str(source),
        "--outdir", str(root / "wrongscale"),
        "--reference-scale", "profiler_scale",
    ], text=True, capture_output=True)
    assert wrong_scale.returncode != 0
    combined = wrong_scale.stdout + wrong_scale.stderr
    assert ("profiler_scale needs column" in combined
            or "requires profiler and" in combined), combined

    # A one-cohort input remains an error by default, because publication-scale
    # analysis requires genuine external holdout validation.
    one_cohort = root / "one_cohort.tsv"
    one_rows = [row for row in rows if row["cohort"] == "alpha"]
    with one_cohort.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(one_rows)
    rejected = subprocess.run([
        "python3", str(SCRIPT), "--responses", str(one_cohort),
        "--outdir", str(root / "single_rejected"),
        "--reference-scale", "read_proportional",
    ], text=True, capture_output=True)
    assert rejected.returncode != 0
    assert "at least two cohorts" in rejected.stderr

    # The explicit development-only mode preserves within-cohort analyses and
    # writes cross-cohort products as auditable header-only tables.
    single_output = root / "single_allowed"
    subprocess.run([
        "python3", str(SCRIPT), "--responses", str(one_cohort),
        "--outdir", str(single_output),
        "--reference-scale", "read_proportional",
        "--cohort-validation", "allow_single_cohort",
    ], check=True)
    assert read_tsv(single_output / "response_operator.tsv")
    assert read_tsv(single_output / "superposition_summary.tsv")
    assert read_tsv(single_output / "reliability_certificates.tsv")
    assert not read_tsv(single_output / "heldout_operator_validation.tsv")
    assert not read_tsv(single_output / "cohort_holdout_summary.tsv")
    assert not read_tsv(single_output / "heldout_validation.tsv")
    assert not read_tsv(single_output / "panel_saturation.tsv")
    validation = {row["metric"]: row["value"] for row in
                  read_tsv(single_output / "perturbation_response_validation.tsv")}
    assert validation["cohort_holdout_status"] == "NOT_APPLICABLE_SINGLE_COHORT"

print("[PASS] perturbation-response operator and reliability fixture")
