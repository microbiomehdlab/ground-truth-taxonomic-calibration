#!/usr/bin/env python3
"""Estimate target-excluded profiler response operators and reliability certificates.

The input is a compact paired-feature table: one row per positive-dose profile and
reported-feature union.  Community profiles must occur only once, with their
component taxa encoded in the two semicolon-delimited mapping fields.  This
script deliberately separates binary detection transitions from quantitative
errors and never uses a feature's own implantation when estimating its
target-agnostic reliability.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import math
import random
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import pyarrow.parquet as pq

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from reference_scale import (  # noqa: E402
    PROFILER_SCALE_ALIASES, add_reference_scale_argument,
    validate_profiler_reference_pairs)


REQUIRED = {
    "cohort", "sample_id", "analysis_population", "profiler",
    "observation_id", "perturbation_id", "baseline_id", "feature",
    "dose_level", "dose_index", "nominal_total_fraction",
    "effective_total_fraction", "implanted_target_count", "implanted_targets",
    "implanted_features", "implanted_fraction_by_target",
    "baseline_abundance_fraction", "observed_abundance_fraction",
    "target_fraction_for_feature", "expected_abundance_fraction",
    "dilution_retained_baseline", "response_signal", "response_delta",
    "signed_bounded_error", "absolute_bounded_error",
    "baseline_detected", "observed_detected", "detection_transition",
    "is_direct_target", "reference_kind",
}

NOMINAL_DOSES = {
    "independent": (0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05),
    "community": (0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1),
}

DETECTION_FIELDS = [
    "cohort", "analysis_population", "profiler", "target_label",
    "feature_role", "dose_level", "dose_fraction_nominal",
    "eligible_contexts", "baseline_absent_contexts",
    "baseline_present_contexts", "remained_absent",
    "unexpected_appearances", "remained_detected", "unexpected_dropouts",
    "appearance_rate", "dropout_rate",
]

QUANTITATIVE_FIELDS = [
    "cohort", "analysis_population", "profiler", "target_label",
    "feature_role", "dose_level", "dose_fraction_nominal",
    "eligible_contexts", "mean_expected_abundance",
    "mean_observed_abundance", "mean_signed_error", "mean_absolute_error",
    "rmse", "mean_log2_response_error", "quantitative_instability_rate",
]

OPERATOR_FIELDS = [
    "holdout_cohort", "training_cohorts", "cohort", "analysis_population",
    "profiler", "target_label", "feature",
    "feature_role", "eligible_contexts", "distinct_samples",
    "distinct_doses", "operator_slope", "operator_intercept",
    "operator_r_squared", "mean_signed_error", "mean_absolute_error",
    "detection_transition_rate",
]

SUPERPOSITION_FIELDS = [
    "cohort", "profiler", "dose_level", "dose_fraction_nominal",
    "feature_role", "eligible_contexts", "distinct_samples",
    "composition_only_mae", "operator_prediction_mae",
    "operator_prediction_rmse", "operator_r_squared", "calibration_slope",
    "improvement_over_composition",
]

HOLDOUT_FIELDS = [
    "holdout_cohort", "training_cohorts", "profiler", "feature_role",
    "eligible_contexts", "distinct_samples", "composition_only_mae",
    "operator_prediction_mae", "operator_prediction_rmse",
    "operator_r_squared", "improvement_over_composition",
]

CERTIFICATE_FIELDS = [
    "holdout_cohort", "training_cohorts", "profiler", "feature",
    "eligible_contexts", "eligible_targets", "eligible_samples",
    "baseline_absent_contexts", "baseline_present_contexts",
    "quantitative_contexts", "unexpected_appearance_rate",
    "unexpected_dropout_rate", "quantitative_instability_rate",
    "mean_absolute_error", "hallucination_reliability",
    "stability_reliability", "quantitative_reliability",
    "overall_reliability", "reliability_ci_lower", "reliability_ci_upper",
    "evaluation_contexts", "evaluation_overall_reliability",
]

COHORT_HOLDOUT_FIELDS = [
    "holdout_cohort", "training_cohorts", "profiler", "features_trained",
    "features_evaluated", "features_shared", "spearman_reliability",
    "mean_absolute_reliability_gap", "mean_training_reliability",
    "mean_evaluation_reliability",
]

TARGET_HOLDOUT_FIELDS = [
    "holdout_cohort", "training_cohorts", "profiler", "validation_type",
    "heldout_target", "feature", "predicted_reliability",
    "observed_reliability", "training_contexts", "evaluation_contexts",
]

PANEL_FIELDS = [
    "holdout_cohort", "training_cohorts", "profiler", "panel_size",
    "panel_replicate", "panel_targets", "features_evaluated",
    "rank_correlation_full_panel", "heldout_reliability_mae",
    "high_risk_taxon_recall",
]

BIOMARKER_FIELDS = [
    "holdout_cohort", "training_cohorts", "analysis_population", "profiler",
    "feature", "source_effect", "source_q_value", "overall_reliability",
    "external_evaluable", "replication_status", "replicated",
]


def render(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return format(value, ".17g")


def divide(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator else None


def mean(total: float, count: int) -> float | None:
    return divide(total, count)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_mapping(raw: str, field_name: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in raw.split(";"):
        item = item.strip()
        if not item:
            continue
        separator = "=" if "=" in item else ":" if ":" in item else None
        if separator is None:
            raise ValueError(f"{field_name} entry lacks '=' or ':': {item!r}")
        key, value = (part.strip() for part in item.split(separator, 1))
        if not key or not value or key in result:
            raise ValueError(f"invalid or duplicate {field_name} entry: {item!r}")
        result[key] = value
    if not result:
        raise ValueError(f"empty {field_name}")
    return result


def number(raw: str, field_name: str, line: int) -> float:
    try:
        value = float(raw)
    except (TypeError, ValueError) as error:
        raise ValueError(f"line {line}: nonnumeric {field_name}: {raw!r}") from error
    if not math.isfinite(value):
        raise ValueError(f"line {line}: nonfinite {field_name}")
    return value


def nominal_dose(population: str, dose_level: str, line: int) -> float:
    grid = NOMINAL_DOSES.get(population)
    if grid is None:
        raise ValueError(f"line {line}: unsupported analysis_population {population!r}")
    if not dose_level.startswith("dose_"):
        raise ValueError(f"line {line}: dose_level must be dose_01, dose_02, ...")
    try:
        rank = int(dose_level[5:])
    except ValueError as error:
        raise ValueError(f"line {line}: invalid dose_level {dose_level!r}") from error
    if rank < 1 or rank > len(grid) or dose_level != f"dose_{rank:02d}":
        raise ValueError(f"line {line}: dose_level outside frozen {population} grid")
    return grid[rank - 1]


@dataclass(frozen=True)
class ResponseRow:
    observation_id: str
    cohort: str
    sample_id: str
    population: str
    profiler: str
    perturbation_id: str
    baseline_id: str
    feature: str
    dose_level: str
    dose_index: int
    nominal_total: float
    actual_total: float
    fractions: dict[str, float]
    target_feature_by_label: dict[str, str]
    implanted_features: frozenset[str]
    baseline: float
    observed: float
    supplied_direct_fraction: float
    supplied_retained_baseline: float
    supplied_expected: float
    supplied_signal: float
    supplied_signed_error: float
    baseline_detected: bool
    observed_detected: bool

    @property
    def targets(self) -> tuple[str, ...]:
        return tuple(sorted(self.fractions))

    @property
    def target_label(self) -> str:
        return self.targets[0] if self.population == "independent" else "COMMUNITY"

    @property
    def is_target(self) -> bool:
        return self.feature in self.implanted_features

    @property
    def role(self) -> str:
        return "implanted_target" if self.is_target else "bystander"

    @property
    def direct_fraction(self) -> float:
        return self.supplied_direct_fraction

    @property
    def retained_baseline(self) -> float:
        """The selected-scale retained baseline exactly as supplied.

        Never reconstructed as ``(1 - effective_total_fraction) * baseline``:
        that is the read-proportional formula only. Under ``profiler_scale``
        MetaPhlAn's retained baseline is ``(1-F)o / D_i``, which the builder
        emits as ``retained_baseline_profiler_scale``.
        """
        return self.supplied_retained_baseline

    @property
    def expected(self) -> float:
        return self.supplied_expected

    @property
    def signal(self) -> float:
        return self.supplied_signal

    @property
    def signed_error(self) -> float:
        return self.supplied_signed_error


def truth(raw: object, field_name: str, line: int) -> bool:
    normalized = str(raw).strip().lower()
    if normalized in {"1", "true", "t", "yes"}:
        return True
    if normalized in {"0", "false", "f", "no"}:
        return False
    raise ValueError(f"line {line}: invalid Boolean {field_name}: {raw!r}")


def apply_reference_scale(fields: list[str], scale: str) -> dict[str, str]:
    """Map canonical field -> source field for the chosen reference scale.

    `read_proportional` keeps the original columns. `profiler_scale` reads the
    `*_profiler_scale` columns instead, which are read-proportional for Bracken
    and genome-equivalent for MetaPhlAn, so Bracken results are unchanged.
    """
    if scale == "read_proportional":
        return {name: name for name in fields}
    missing = sorted(set(PROFILER_SCALE_ALIASES.values()) - set(fields))
    if missing:
        raise SystemExit(
            "[ERROR] --reference-scale profiler_scale needs column(s) "
            + ", ".join(missing)
            + "; rebuild the response input with the migrated builder")
    mapping = {name: name for name in fields}
    mapping.update(PROFILER_SCALE_ALIASES)
    return mapping


def validate_reference_pairs_from_rows(rows, scale: str) -> str:
    """Collect distinct (profiler, reference_type) pairs and validate them.

    Validation is within profiler: a primary table legitimately holds
    read_proportional Bracken rows alongside genome_equivalent MetaPhlAn rows.
    """
    pairs = {(row.get("profiler"), row.get("reference_type")) for row in rows}
    return validate_profiler_reference_pairs(sorted(pairs), scale)


def source_rows(path: Path, scale: str = "read_proportional"
                ) -> tuple[list[str], Iterable[tuple[int, dict[str, object]]]]:
    """Return a schema and a batch-streamed row iterator."""
    if path.suffix.lower() == ".parquet":
        parquet = pq.ParquetFile(path)
        fields = list(parquet.schema_arrow.names)

        remap = apply_reference_scale(fields, scale)

        def parquet_rows() -> Iterable[tuple[int, dict[str, object]]]:
            line = 1
            for batch in parquet.iter_batches(batch_size=65536):
                columns = batch.to_pydict()
                for index in range(batch.num_rows):
                    line += 1
                    yield line, {name: columns[remap[name]][index] for name in fields}

        return fields, parquet_rows()

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])

    # The TSV branch applies the identical remap, so a profiler-scale request
    # against a table that lacks those columns fails instead of silently
    # falling back to the read-proportional values.
    tsv_remap = apply_reference_scale(fields, scale)

    def tsv_rows() -> Iterable[tuple[int, dict[str, object]]]:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for line, row in enumerate(reader, start=2):
                yield line, {name: row[tsv_remap[name]] for name in fields}

    return fields, tsv_rows()


def iter_rows(path: Path, scale: str = "read_proportional") -> Iterable[ResponseRow]:
    # Validate the profiler/reference mapping before any quantitative use.
    probe_fields, probe = source_rows(path, "read_proportional")
    if "profiler" in probe_fields and "reference_type" in probe_fields:
        seen = []
        for _, row in probe:
            seen.append({"profiler": row.get("profiler"),
                         "reference_type": row.get("reference_type")})
        validate_reference_pairs_from_rows(seen, scale)
    elif scale == "profiler_scale":
        raise SystemExit(
            "[ERROR] --reference-scale profiler_scale requires profiler and "
            "reference_type columns; rebuild the response input with the "
            "migrated build_perturbation_response_input.py")
    field_names, rows = source_rows(path, scale)
    fields = set(field_names)
    if not field_names:
        raise ValueError("input is empty")
    if not REQUIRED <= fields:
        raise ValueError("input missing columns: " + ", ".join(sorted(REQUIRED - fields)))
    for line, source in rows:
            source = {key: "" if value is None else str(value) for key, value in source.items()}
            for field_name in REQUIRED - {
                "baseline_abundance_fraction", "observed_abundance_fraction",
                "effective_total_fraction", "nominal_total_fraction",
                "target_fraction_for_feature", "expected_abundance_fraction",
                "dilution_retained_baseline", "response_signal", "response_delta",
                "signed_bounded_error", "absolute_bounded_error",
                "quantitative_log2_ratio",
            }:
                if not source[field_name].strip():
                    raise ValueError(f"line {line}: empty {field_name}")
            population = source["analysis_population"].strip()
            nominal = nominal_dose(population, source["dose_level"].strip(), line)
            supplied_nominal = number(
                source["nominal_total_fraction"], "nominal_total_fraction", line
            )
            if not math.isclose(supplied_nominal, nominal, rel_tol=1e-9, abs_tol=1e-12):
                raise ValueError(f"line {line}: nominal fraction disagrees with dose rank")
            actual_total = number(
                source["effective_total_fraction"], "effective_total_fraction", line
            )
            baseline = number(
                source["baseline_abundance_fraction"], "baseline_abundance_fraction", line
            )
            observed = number(
                source["observed_abundance_fraction"], "observed_abundance_fraction", line
            )
            if not 0 < actual_total < 1:
                raise ValueError(f"line {line}: spike_fraction_total must be between 0 and 1")
            if not 0 <= baseline <= 1 or not 0 <= observed <= 1:
                raise ValueError(f"line {line}: abundance fraction outside [0,1]")
            raw_fractions = parse_mapping(
                source["implanted_fraction_by_target"], "implanted_fraction_by_target"
            )
            fractions = {
                key: number(value, "implanted_fraction_by_target", line)
                for key, value in raw_fractions.items()
            }
            target_labels = tuple(
                item.strip() for item in source["implanted_targets"].split(";") if item.strip()
            )
            feature_list = tuple(
                item.strip() for item in source["implanted_features"].split(";") if item.strip()
            )
            implanted_features = frozenset(feature_list)
            if set(target_labels) != set(fractions):
                raise ValueError(f"line {line}: implanted target list and mapping keys differ")
            if len(feature_list) != len(target_labels):
                raise ValueError(f"line {line}: implanted targets/features have different lengths")
            target_feature_by_label = dict(zip(target_labels, feature_list))
            try:
                dose_index = int(source["dose_index"])
                target_count = int(source["implanted_target_count"])
            except ValueError as error:
                raise ValueError(f"line {line}: invalid integer dose_index/target_count") from error
            if dose_index != int(source["dose_level"][5:]) or target_count != len(target_labels):
                raise ValueError(f"line {line}: inconsistent dose index or target count")
            if any(value <= 0 for value in fractions.values()):
                raise ValueError(f"line {line}: component fractions must be positive")
            if not math.isclose(sum(fractions.values()), actual_total,
                                rel_tol=1e-6, abs_tol=1e-12):
                raise ValueError(f"line {line}: component fractions do not sum to total")
            if population == "independent" and len(fractions) != 1:
                raise ValueError(f"line {line}: independent perturbation must have one target")
            if population == "community" and len(fractions) < 2:
                raise ValueError(f"line {line}: community perturbation needs multiple targets")
            direct = number(source["target_fraction_for_feature"],
                            "target_fraction_for_feature", line)
            expected = number(source["expected_abundance_fraction"],
                              "expected_abundance_fraction", line)
            supplied_retained = number(source["dilution_retained_baseline"],
                                       "dilution_retained_baseline", line)
            signal = number(source["response_signal"], "response_signal", line)
            response_delta = number(source["response_delta"], "response_delta", line)
            signed_error = number(source["signed_bounded_error"], "signed_bounded_error", line)
            absolute_error = number(source["absolute_bounded_error"],
                                    "absolute_bounded_error", line)
            baseline_detected = truth(source["baseline_detected"], "baseline_detected", line)
            observed_detected = truth(source["observed_detected"], "observed_detected", line)
            direct_flag = truth(source["is_direct_target"], "is_direct_target", line)
            transition = source["detection_transition"].strip()
            # Scale-independent identities. The selected retained baseline is
            # authoritative: reconstructing it as (1-F)*baseline is only the
            # read-proportional formula and is wrong for genome-equivalent
            # MetaPhlAn, where the composition is renormalised by D_i.
            # effective_total_fraction remains available as perturbation
            # metadata but must not define the selected retained abundance.
            retained = supplied_retained
            checks = (
                (expected, retained + direct, "expected abundance"),
                (signal, observed - retained, "response signal"),
                (response_delta, observed - expected, "response delta"),
                (signed_error,
                 (observed - expected) / (observed + expected)
                 if observed + expected > 0 else 0.0,
                 "signed bounded error"),
                (absolute_error, abs(signed_error), "absolute bounded error"),
            )
            for supplied, calculated, label in checks:
                if not math.isclose(supplied, calculated, rel_tol=1e-7, abs_tol=1e-12):
                    raise ValueError(f"line {line}: inconsistent {label}")
            if baseline_detected != (baseline > 0) or observed_detected != (observed > 0):
                raise ValueError(f"line {line}: detection flag disagrees with abundance")
            expected_transition = (
                "RETAINED" if baseline_detected and observed_detected else
                "LOST" if baseline_detected else
                "GAINED" if observed_detected else "STABLE_ABSENT"
            )
            if transition != expected_transition:
                raise ValueError(f"line {line}: inconsistent detection transition")
            if direct_flag != (source["feature"].strip() in implanted_features):
                raise ValueError(f"line {line}: direct-target flag disagrees with feature list")
            yield ResponseRow(
                observation_id=source["observation_id"].strip(),
                cohort=source["cohort"].strip(),
                sample_id=source["sample_id"].strip(),
                population=population,
                profiler=source["profiler"].strip(),
                perturbation_id=source["perturbation_id"].strip(),
                baseline_id=source["baseline_id"].strip(),
                feature=source["feature"].strip(),
                dose_level=source["dose_level"].strip(),
                dose_index=dose_index,
                nominal_total=nominal,
                actual_total=actual_total,
                fractions=fractions,
                target_feature_by_label=target_feature_by_label,
                implanted_features=implanted_features,
                baseline=baseline,
                observed=observed,
                supplied_direct_fraction=direct,
                supplied_retained_baseline=supplied_retained,
                supplied_expected=expected,
                supplied_signal=signal,
                supplied_signed_error=signed_error,
                baseline_detected=baseline_detected,
                observed_detected=observed_detected,
            )


@dataclass(frozen=True)
class IndependentDriver:
    """The selected-scale implanted dose of one independent perturbation.

    The driver is read off the perturbation's *direct-target* row, because that
    is the only row whose ``target_fraction_for_feature`` is non-zero. Under
    ``profiler_scale`` the canonical name resolves to
    ``implanted_signal_profiler_scale``, i.e. the read fraction for Bracken and
    ``q_it / D_i`` for MetaPhlAn; under ``read_proportional`` it is the read
    fraction for both profilers.
    """

    observation_id: str
    cohort: str
    sample_id: str
    profiler: str
    target: str
    feature: str
    value: float


def collect_independent_drivers(
    rows: Iterable[ResponseRow],
) -> dict[str, IndependentDriver]:
    """Deterministic pre-pass: one selected-scale driver per independent observation.

    ``row.fractions[target]`` must never be used instead. It is the implanted
    *read* fraction taken from the perturbation metadata and is wrong on the
    genome-equivalent scale, where the operator's regressor is ``q_it / D_i``.
    Non-implanted features of the same perturbation share this driver; their own
    ``direct_fraction`` is structurally zero.
    """
    drivers: dict[str, IndependentDriver] = {}
    observations: dict[str, tuple[str, str, str, str]] = {}
    for row in rows:
        if row.population != "independent":
            continue
        observation = row.observation_id
        target = row.targets[0]
        identity = (row.cohort, row.sample_id, row.profiler, target)
        previous_identity = observations.setdefault(observation, identity)
        if previous_identity != identity:
            raise ValueError(
                f"observation {observation!r} carries conflicting physical identity "
                f"{previous_identity} versus {identity}")
        if row.feature != row.target_feature_by_label[target]:
            continue
        value = row.direct_fraction
        if not math.isfinite(value) or value <= 0:
            raise ValueError(
                f"observation {observation!r} has a non-positive or nonfinite "
                f"selected-scale implanted driver: {value!r}")
        candidate = IndependentDriver(
            observation_id=observation, cohort=row.cohort, sample_id=row.sample_id,
            profiler=row.profiler, target=target, feature=row.feature, value=value)
        previous = drivers.get(observation)
        if previous is None:
            drivers[observation] = candidate
        else:
            raise ValueError(
                f"observation {observation!r} has more than one direct-target "
                f"driver row: {previous.value!r} and {value!r}")
    missing = sorted(set(observations) - set(drivers))
    if missing:
        raise ValueError(
            f"{len(missing)} independent observation(s) lack a direct-target row "
            f"supplying the selected implanted driver, e.g. {missing[0]!r}")
    return drivers


def driver_for(drivers: dict[str, IndependentDriver], row: ResponseRow) -> float:
    """Selected-scale driver for `row`, validated against its physical identity."""
    driver = drivers.get(row.observation_id)
    if driver is None:
        raise ValueError(
            f"no independent driver recorded for observation {row.observation_id!r}")
    target = row.targets[0]
    if (driver.cohort, driver.sample_id, driver.profiler, driver.target) != (
            row.cohort, row.sample_id, row.profiler, target):
        raise ValueError(
            f"driver identity mismatch for observation {row.observation_id!r}: "
            f"{(driver.cohort, driver.sample_id, driver.profiler, driver.target)} "
            f"versus {(row.cohort, row.sample_id, row.profiler, target)}")
    return driver.value


@dataclass
class DetectionAggregate:
    n: int = 0
    absent: int = 0
    present: int = 0
    remained_absent: int = 0
    appearances: int = 0
    remained_detected: int = 0
    dropouts: int = 0

    def add(self, row: ResponseRow) -> None:
        self.n += 1
        if row.baseline > 0:
            self.present += 1
            if row.observed > 0:
                self.remained_detected += 1
            else:
                self.dropouts += 1
        else:
            self.absent += 1
            if row.observed > 0:
                self.appearances += 1
            else:
                self.remained_absent += 1


@dataclass
class QuantitativeAggregate:
    n: int = 0
    expected_sum: float = 0.0
    observed_sum: float = 0.0
    signed_sum: float = 0.0
    absolute_sum: float = 0.0
    squared_sum: float = 0.0
    log2_sum: float = 0.0
    log2_n: int = 0
    unstable: int = 0

    def add(self, row: ResponseRow, threshold: float) -> None:
        self.n += 1
        error = row.signed_error
        self.expected_sum += row.expected
        self.observed_sum += row.observed
        self.signed_sum += error
        self.absolute_sum += abs(error)
        self.squared_sum += error * error
        if row.expected > 0 and row.observed > 0:
            log_error = math.log2(row.observed / row.expected)
            self.log2_sum += log_error
            self.log2_n += 1
            self.unstable += abs(log_error) > threshold


@dataclass
class RegressionAggregate:
    n: int = 0
    sum_x: float = 0.0
    sum_y: float = 0.0
    sum_xx: float = 0.0
    sum_xy: float = 0.0
    sum_yy: float = 0.0
    error_sum: float = 0.0
    absolute_error_sum: float = 0.0
    transitions: int = 0
    sample_bits: int = 0
    dose_bits: int = 0

    def add(self, row: ResponseRow, x: float, sample_bit: int, dose_bit: int) -> None:
        y = row.signal
        self.n += 1
        self.sum_x += x
        self.sum_y += y
        self.sum_xx += x * x
        self.sum_xy += x * y
        self.sum_yy += y * y
        self.error_sum += row.signed_error
        self.absolute_error_sum += abs(row.signed_error)
        self.transitions += (row.baseline > 0) != (row.observed > 0)
        self.sample_bits |= sample_bit
        self.dose_bits |= dose_bit

    def slope_origin(self) -> float | None:
        return divide(self.sum_xy, self.sum_xx)

    def intercept(self) -> float | None:
        denominator = self.n * self.sum_xx - self.sum_x * self.sum_x
        if self.n < 2 or denominator == 0:
            return None
        slope = (self.n * self.sum_xy - self.sum_x * self.sum_y) / denominator
        return (self.sum_y - slope * self.sum_x) / self.n

    def r_squared_origin(self) -> float | None:
        slope = self.slope_origin()
        if slope is None or self.sum_yy == 0:
            return None
        sse = self.sum_yy - 2 * slope * self.sum_xy + slope * slope * self.sum_xx
        return 1.0 - sse / self.sum_yy


@dataclass
class OperatorDesign:
    """Denominator information, including structurally zero feature responses."""

    n: int = 0
    sum_x: float = 0.0
    sum_xx: float = 0.0
    sample_bits: int = 0
    dose_bits: int = 0

    def add(self, x: float, sample_bit: int, dose_bit: int) -> None:
        self.n += 1
        self.sum_x += x
        self.sum_xx += x * x
        self.sample_bits |= sample_bit
        self.dose_bits |= dose_bit


def complete_regression(
    observed: RegressionAggregate | None, design: OperatorDesign,
) -> RegressionAggregate:
    """Add response rows that are implicit zeroes in the sparse input."""
    source = observed or RegressionAggregate()
    return RegressionAggregate(
        n=design.n, sum_x=design.sum_x, sum_y=source.sum_y,
        sum_xx=design.sum_xx, sum_xy=source.sum_xy, sum_yy=source.sum_yy,
        error_sum=source.error_sum,
        absolute_error_sum=source.absolute_error_sum,
        transitions=source.transitions, sample_bits=design.sample_bits,
        dose_bits=design.dose_bits,
    )


@dataclass
class ReliabilityAggregate:
    n: int = 0
    absent: int = 0
    appearances: int = 0
    present: int = 0
    dropouts: int = 0
    quantitative: int = 0
    unstable: int = 0
    absolute_error_sum: float = 0.0
    target_bits: int = 0
    sample_bits: int = 0

    def add(self, row: ResponseRow, threshold: float,
            sample_bit: int, target_bits: int) -> None:
        self.n += 1
        self.target_bits |= target_bits
        self.sample_bits |= sample_bit
        if row.baseline > 0:
            self.present += 1
            self.dropouts += row.observed == 0
        else:
            self.absent += 1
            self.appearances += row.observed > 0
        if row.expected > 0 and row.observed > 0:
            self.quantitative += 1
            log_error = math.log2(row.observed / row.expected)
            self.unstable += abs(log_error) > threshold
            self.absolute_error_sum += abs(row.signed_error)

    def merge(self, other: "ReliabilityAggregate") -> None:
        for name in ("n", "absent", "appearances", "present", "dropouts",
                     "quantitative", "unstable"):
            setattr(self, name, getattr(self, name) + getattr(other, name))
        self.absolute_error_sum += other.absolute_error_sum
        self.target_bits |= other.target_bits
        self.sample_bits |= other.sample_bits

    @staticmethod
    def posterior_failure_rate(events: int, eligible: int) -> float | None:
        return (events + 0.5) / (eligible + 1.0) if eligible else None

    def components(self) -> tuple[float | None, float | None, float | None]:
        risks = (
            self.posterior_failure_rate(self.appearances, self.absent),
            self.posterior_failure_rate(self.dropouts, self.present),
            self.posterior_failure_rate(self.unstable, self.quantitative),
        )
        return tuple(None if risk is None else 1.0 - risk for risk in risks)

    def score(self) -> float | None:
        components = [value for value in self.components() if value is not None]
        return sum(components) / len(components) if components else None

    def credible_interval(self) -> tuple[float | None, float | None]:
        """Normal approximation to the mean of independent Jeffreys posteriors."""
        pairs = (
            (self.absent - self.appearances, self.appearances),
            (self.present - self.dropouts, self.dropouts),
            (self.quantitative - self.unstable, self.unstable),
        )
        moments = []
        for successes, failures in pairs:
            if successes + failures == 0:
                continue
            alpha, beta = successes + 0.5, failures + 0.5
            total = alpha + beta
            moments.append((alpha / total, alpha * beta / (total * total * (total + 1))))
        if not moments:
            return None, None
        estimate = sum(item[0] for item in moments) / len(moments)
        variance = sum(item[1] for item in moments) / (len(moments) ** 2)
        half = 1.959963984540054 * math.sqrt(variance)
        return max(0.0, estimate - half), min(1.0, estimate + half)


def complete_reliability(
    observed: ReliabilityAggregate | None, eligible: int,
    sample_bits: int = 0, target_bits: int = 0,
) -> ReliabilityAggregate:
    """Restore implicit stable absences omitted from the sparse response table."""
    source = observed or ReliabilityAggregate()
    if source.present > eligible:
        raise ValueError("baseline-present count exceeds reliability opportunities")
    result = ReliabilityAggregate(
        n=eligible, absent=eligible - source.present,
        appearances=source.appearances, present=source.present,
        dropouts=source.dropouts, quantitative=source.quantitative,
        unstable=source.unstable, absolute_error_sum=source.absolute_error_sum,
        target_bits=target_bits or source.target_bits,
        sample_bits=sample_bits or source.sample_bits,
    )
    if result.appearances > result.absent:
        raise ValueError("appearance count exceeds baseline-absent opportunities")
    return result


@dataclass
class PredictionAggregate:
    n: int = 0
    composition_abs: float = 0.0
    operator_abs: float = 0.0
    operator_sq: float = 0.0
    observed_sum: float = 0.0
    predicted_sum: float = 0.0
    observed_sq: float = 0.0
    predicted_sq: float = 0.0
    cross: float = 0.0
    sample_bits: int = 0

    def add(self, row: ResponseRow, prediction: float, sample_bit: int) -> None:
        self.n += 1
        self.composition_abs += abs(row.observed - row.expected)
        error = row.observed - prediction
        self.operator_abs += abs(error)
        self.operator_sq += error * error
        self.observed_sum += row.observed
        self.predicted_sum += prediction
        self.observed_sq += row.observed * row.observed
        self.predicted_sq += prediction * prediction
        self.cross += row.observed * prediction
        self.sample_bits |= sample_bit

    def r_squared(self) -> float | None:
        if self.n < 2:
            return None
        sst = self.observed_sq - self.observed_sum * self.observed_sum / self.n
        return None if sst <= 0 else 1.0 - self.operator_sq / sst

    def calibration_slope(self) -> float | None:
        denominator = self.predicted_sq - self.predicted_sum * self.predicted_sum / self.n
        numerator = self.cross - self.predicted_sum * self.observed_sum / self.n
        return divide(numerator, denominator)


SUPERPOSITION_COLUMNS = {
    "read_proportional": {
        "component_signal": "response_signal",
        "community_expected": "expected_abundance_fraction",
        "community_retained": "dilution_retained_baseline",
    },
    "profiler_scale": {
        "component_signal": "response_signal_profiler_scale",
        "community_expected": "expected_abundance_profiler_scale",
        "community_retained": "retained_baseline_profiler_scale",
    },
}


def exact_superposition_rows(path: Path, reference_scale: str) -> list[dict[str, object]]:
    """Match community mixtures to same-sample single-target perturbations.

    Community dose indices 3--7 have per-target fractions corresponding to
    independent indices 1--5.  The ordinal join is intentional: historical
    achieved fractions differ slightly and must never be joined as floats.
    DuckDB performs the many-million-row join out of core.

    Exactly one internally consistent column triple is selected, so a
    profiler-scale run never reads a read-proportional column. The community
    retained baseline is read from its own column and is never derived as
    ``expected - target_fraction_for_feature``: that subtraction silently mixes
    scales and is wrong whenever the two sides use different references.
    """
    try:
        columns = SUPERPOSITION_COLUMNS[reference_scale]
    except KeyError:
        raise ValueError(
            f"unknown reference scale for superposition: {reference_scale!r}") from None
    try:
        import duckdb
    except ImportError as error:
        raise ValueError("DuckDB is required for exact superposition") from error
    connection = duckdb.connect(database=":memory:")
    try:
        connection.execute("PRAGMA threads=2")
        connection.execute("PRAGMA preserve_insertion_order=false")
        if path.suffix.lower() == ".parquet":
            connection.read_parquet(str(path)).create_view("responses")
        else:
            connection.read_csv(str(path), header=True, delimiter="\t",
                                sample_size=-1).create_view("responses")
        available = {row[0] for row in
                     connection.execute("DESCRIBE responses").fetchall()}
        missing = sorted(set(columns.values()) - available)
        if missing:
            raise ValueError(
                f"superposition on the {reference_scale} scale needs column(s) "
                + ", ".join(missing))
        query = """
        WITH independent AS (
          SELECT cohort, sample_id, profiler, feature, CAST(dose_index AS INTEGER) AS dose_index,
                 trim(implanted_targets) AS target_label,
                 CAST({component_signal} AS DOUBLE) AS response_signal
          FROM responses
          WHERE analysis_population = 'independent'
        ),
        independent_observations AS (
          SELECT DISTINCT cohort, sample_id, profiler, CAST(dose_index AS INTEGER) AS dose_index,
                 trim(implanted_targets) AS target_label
          FROM responses
          WHERE analysis_population = 'independent'
        ),
        expanded_community AS (
          SELECT observation_id, cohort, sample_id, profiler, perturbation_id, feature,
                 dose_level, CAST(dose_index AS INTEGER) AS dose_index,
                 CAST(nominal_total_fraction AS DOUBLE) AS nominal_total_fraction,
                 CAST(implanted_target_count AS INTEGER) AS implanted_target_count,
                 CAST(observed_abundance_fraction AS DOUBLE) AS observed,
                 CAST({community_expected} AS DOUBLE) AS expected,
                 CAST({community_retained} AS DOUBLE) AS retained,
                 CASE WHEN CAST(is_direct_target AS BOOLEAN)
                      THEN 'implanted_target' ELSE 'bystander' END AS feature_role,
                 unnest(string_split(implanted_targets, ';')) AS target_label
          FROM responses
          WHERE analysis_population = 'community'
            AND CAST(dose_index AS INTEGER) BETWEEN 3 AND 7
        ),
        components AS (
          SELECT c.*, io.target_label AS matched_target,
                 coalesce(i.response_signal, 0.0) AS component_signal
          FROM expanded_community c
          LEFT JOIN independent_observations io
            ON io.cohort = c.cohort AND io.sample_id = c.sample_id
           AND io.profiler = c.profiler AND io.target_label = trim(c.target_label)
           AND io.dose_index = c.dose_index - 2
          LEFT JOIN independent i
            ON i.cohort = io.cohort AND i.sample_id = io.sample_id
           AND i.profiler = io.profiler AND i.target_label = io.target_label
           AND i.dose_index = io.dose_index AND i.feature = c.feature
        ),
        predictions AS (
          SELECT c.observation_id, c.cohort, c.sample_id, c.profiler, c.perturbation_id,
                 c.feature, c.dose_level, c.dose_index, c.nominal_total_fraction,
                 c.feature_role, max(c.observed) AS observed, max(c.expected) AS expected,
                 max(c.retained) + sum(c.component_signal) AS predicted,
                 count(c.matched_target) AS matched_targets,
                 max(c.implanted_target_count) AS expected_targets
          FROM components c
          GROUP BY c.observation_id, c.cohort, c.sample_id, c.profiler,
                   c.perturbation_id, c.feature, c.dose_level, c.dose_index,
                   c.nominal_total_fraction, c.feature_role
          HAVING count(c.matched_target) = max(c.implanted_target_count)
        )
        SELECT cohort, profiler, dose_level,
               nominal_total_fraction AS dose_fraction_nominal, feature_role,
               count(*) AS eligible_contexts, count(DISTINCT sample_id) AS distinct_samples,
               avg(abs(observed - expected)) AS composition_only_mae,
               avg(abs(observed - predicted)) AS operator_prediction_mae,
               sqrt(avg(pow(observed - predicted, 2))) AS operator_prediction_rmse,
               CASE WHEN var_pop(observed) > 0
                    THEN 1 - sum(pow(observed - predicted, 2)) /
                             (count(*) * var_pop(observed)) ELSE NULL END AS operator_r_squared,
               regr_slope(observed, predicted) AS calibration_slope,
               CASE WHEN avg(abs(observed - expected)) > 0
                    THEN (avg(abs(observed - expected)) - avg(abs(observed - predicted)))
                         / avg(abs(observed - expected))
                    ELSE NULL END AS improvement_over_composition
        FROM predictions
        GROUP BY cohort, profiler, dose_level, nominal_total_fraction, feature_role
        ORDER BY cohort, profiler, dose_level, feature_role
        """.format(**columns)
        result = connection.execute(query)
        fields = [item[0] for item in result.description]
        return [
            {
                field: (render(value) if isinstance(value, float) else value)
                for field, value in zip(fields, row)
            }
            for row in result.fetchall()
        ]
    finally:
        connection.close()


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0] * len(values)
    position = 0
    while position < len(order):
        end = position + 1
        while end < len(order) and values[order[end]] == values[order[position]]:
            end += 1
        average_rank = (position + 1 + end) / 2.0
        for index in order[position:end]:
            result[index] = average_rank
        position = end
    return result


def correlation(left: list[float], right: list[float]) -> float | None:
    if len(left) < 2 or len(left) != len(right):
        return None
    lm, rm = sum(left) / len(left), sum(right) / len(right)
    cross = sum((x - lm) * (y - rm) for x, y in zip(left, right))
    lv = sum((x - lm) ** 2 for x in left)
    rv = sum((y - rm) ** 2 for y in right)
    return divide(cross, math.sqrt(lv * rv)) if lv and rv else None


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    temporary.replace(path)


def combine_reliability(
    source: dict[tuple[str, str, str], ReliabilityAggregate],
    cohorts: Iterable[str], profiler: str, feature: str,
) -> ReliabilityAggregate:
    result = ReliabilityAggregate()
    for cohort in cohorts:
        aggregate = source.get((cohort, profiler, feature))
        if aggregate:
            result.merge(aggregate)
    return result


def reliability_score_for_targets(
    source: dict[tuple[str, str, str, str], ReliabilityAggregate],
    cohorts: Iterable[str], profiler: str, feature: str, targets: Iterable[str],
) -> tuple[float | None, int]:
    result = ReliabilityAggregate()
    for cohort in cohorts:
        for target in targets:
            aggregate = source.get((cohort, profiler, feature, target))
            if aggregate:
                result.merge(aggregate)
    return result.score(), result.n


def build_biomarker_replication(
    disease_path: Path, certificate_path: Path, cohorts: set[str],
    feature_aliases: Path | None = None,
) -> list[dict[str, object]]:
    """Link baseline biomarkers to target-excluded, leave-cohort-out reliability."""
    try:
        import duckdb
    except ImportError as error:
        raise ValueError("DuckDB is required for biomarker replication") from error
    connection = duckdb.connect(database=":memory:")
    try:
        connection.execute("PRAGMA threads=2")
        connection.read_csv(str(disease_path), header=True, delimiter="\t",
                            sample_size=-1).create_view("disease")
        connection.read_csv(str(certificate_path), header=True, delimiter="\t",
                            sample_size=-1).create_view("certificates")
        connection.execute("""
            CREATE TABLE feature_aliases (
                profiler VARCHAR, source_feature VARCHAR,
                canonical_feature VARCHAR, rationale VARCHAR
            )
        """)
        if feature_aliases:
            # The response runner passes the spike-panel CSV (canonical/alias/tool),
            # while other callers may pass the native TSV mapping schema.
            with feature_aliases.open(newline="", encoding="utf-8") as handle:
                first = handle.readline()
                handle.seek(0)
                delimiter = "," if "," in first and "\t" not in first else "\t"
                reader = csv.DictReader(handle, delimiter=delimiter)
                fields = set(reader.fieldnames or ())
                if {"canonical", "alias", "tool"} <= fields:
                    aliases = [(r["tool"].strip(), r["alias"].strip(),
                                r["canonical"].strip(), "spike-panel alias") for r in reader]
                elif {"profiler", "source_feature", "canonical_feature"} <= fields:
                    aliases = [(r["profiler"].strip(), r["source_feature"].strip(),
                                r["canonical_feature"].strip(), r.get("rationale", ""))
                               for r in reader]
                else:
                    raise ValueError("unrecognized feature alias schema")
            if any(not all(item[:3]) for item in aliases):
                raise ValueError("feature alias mapping contains blank identifiers")
            seen_aliases = {}
            for profiler, source, canonical, _ in aliases:
                key = profiler, source
                if key in seen_aliases and seen_aliases[key] != canonical:
                    raise ValueError(f"conflicting feature alias: {key}")
                seen_aliases[key] = canonical
            if aliases:
                unique_aliases = {(p, s): (p, s, c, r) for p, s, c, r in aliases}
                connection.executemany("INSERT INTO feature_aliases VALUES (?, ?, ?, ?)",
                                       list(unique_aliases.values()))
        cohort_sql = ",".join("'" + item.replace("'", "''") + "'" for item in sorted(cohorts))
        query = f"""
        WITH baseline AS (
          SELECT d.cohort, d.analysis_population, d.profiler,
                 coalesce(a.canonical_feature, d.feature) AS feature,
                 median(TRY_CAST(effect AS DOUBLE)) AS effect,
                 median(TRY_CAST(q_value AS DOUBLE)) AS q_value
          FROM disease d
          LEFT JOIN feature_aliases a
            ON a.profiler=d.profiler AND a.source_feature=d.feature
          WHERE CAST(d.include AS INTEGER) = 1 AND d.dose_level = 'baseline'
            AND d.contrast = 'CRC_vs_Control' AND d.cohort IN ({cohort_sql})
          GROUP BY d.cohort, d.analysis_population, d.profiler,
                   coalesce(a.canonical_feature, d.feature)
        ), sources AS (
          SELECT * FROM baseline WHERE q_value <= 0.05
        ), external AS (
          SELECT s.cohort AS holdout_cohort, s.analysis_population, s.profiler,
                 s.feature, s.effect AS source_effect, s.q_value AS source_q_value,
                 count(b.cohort) AS external_evaluable,
                 sum(CASE WHEN b.q_value <= 0.05 AND sign(b.effect) = sign(s.effect)
                          THEN 1 ELSE 0 END) AS replicated_count
          FROM sources s
          LEFT JOIN baseline b ON b.cohort <> s.cohort
            AND b.analysis_population = s.analysis_population
            AND b.profiler = s.profiler AND b.feature = s.feature
          GROUP BY s.cohort, s.analysis_population, s.profiler, s.feature,
                   s.effect, s.q_value
        )
        SELECT e.holdout_cohort, c.training_cohorts, e.analysis_population,
               e.profiler, e.feature, e.source_effect, e.source_q_value,
               TRY_CAST(c.overall_reliability AS DOUBLE) AS overall_reliability,
               e.external_evaluable,
               CASE WHEN e.external_evaluable = 0 THEN 'NOT_EVALUABLE'
                    WHEN e.replicated_count > 0 THEN 'DIRECTIONALLY_REPLICATED'
                    ELSE 'NOT_SIGNIFICANT_ELSEWHERE' END AS replication_status,
               CASE WHEN e.replicated_count > 0 THEN 1 ELSE 0 END AS replicated
        FROM external e
        JOIN certificates c ON c.holdout_cohort = e.holdout_cohort
          AND c.profiler = e.profiler AND c.feature = e.feature
        WHERE TRY_CAST(c.overall_reliability AS DOUBLE) IS NOT NULL
        ORDER BY e.holdout_cohort, e.analysis_population, e.profiler, e.feature
        """
        result = connection.execute(query)
        fields = [item[0] for item in result.description]
        return [
            {field: render(value) if isinstance(value, float) else value
             for field, value in zip(fields, record)}
            for record in result.fetchall()
        ]
    finally:
        connection.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", "--responses", dest="input", type=Path, required=True,
                        help="Compact paired-feature response TSV.")
    parser.add_argument("--outdir", type=Path, required=True)
    add_reference_scale_argument(parser)
    parser.add_argument("--log2-instability-threshold", type=float, default=1.0,
                        help="Absolute conditional log2 error defining instability.")
    parser.add_argument("--panel-replicates", type=int, default=50)
    parser.add_argument(
        "--cohort-validation",
        choices=("require_holdout", "allow_single_cohort"),
        default="require_holdout",
        help=("Require cross-cohort holdout validation (default), or explicitly "
              "allow a single-cohort development audit. In the latter case, "
              "cross-cohort output tables are written header-only."),
    )
    parser.add_argument("--disease-results", type=Path,
                        help="Optional primary disease DA table for replication analysis.")
    parser.add_argument("--feature-aliases", type=Path,
                        help="Alias table already applied while building response input.")
    args = parser.parse_args()
    if not args.input.is_file() or not args.input.stat().st_size:
        raise SystemExit("[ERROR] input is missing or empty")
    if not math.isfinite(args.log2_instability_threshold) or args.log2_instability_threshold <= 0:
        raise SystemExit("[ERROR] log2 instability threshold must be positive")

    try:
        args.outdir.mkdir(parents=True, exist_ok=True)
        for stale in (args.outdir / "SUCCESS", args.outdir / "perturbation_response.sha256"):
            if stale.exists():
                stale.unlink()

        detection: dict[tuple[str, ...], DetectionAggregate] = defaultdict(DetectionAggregate)
        quantitative: dict[tuple[str, ...], QuantitativeAggregate] = defaultdict(QuantitativeAggregate)
        operators: dict[tuple[str, str, str, str], RegressionAggregate] = defaultdict(RegressionAggregate)
        operator_designs: dict[tuple[str, str, str], OperatorDesign] = defaultdict(OperatorDesign)
        reliability_observed: dict[tuple[str, str, str], ReliabilityAggregate] = defaultdict(ReliabilityAggregate)
        reliability_target_observed: dict[tuple[str, str, str, str], ReliabilityAggregate] = defaultdict(ReliabilityAggregate)
        observation_counts: dict[tuple[str, ...], int] = defaultdict(int)
        observation_target_counts: dict[tuple[str, ...], int] = {}
        independent_counts: dict[tuple[str, str], int] = defaultdict(int)
        independent_sample_bits: dict[tuple[str, str], int] = defaultdict(int)
        independent_target_bits: dict[tuple[str, str], int] = defaultdict(int)
        direct_opportunities: dict[tuple[str, str, str], int] = defaultdict(int)
        features_by_profiler: dict[str, set[str]] = defaultdict(set)
        target_features: dict[tuple[str, str], str] = {}
        seen_observations: set[str] = set()
        cohorts: set[str] = set()
        profilers: set[str] = set()
        features: set[str] = set()
        input_rows = 0
        target_rows = 0
        bystander_rows = 0
        target_pairs: set[tuple[str, str, str, str]] = set()
        sample_indices: dict[tuple[str, str], int] = {}
        target_indices: dict[str, int] = {}

        # Deterministic pre-pass: the selected-scale implanted driver of every
        # independent perturbation, read from its direct-target row. Every
        # feature of that perturbation -- including non-implanted ones -- is
        # regressed on this one value.
        drivers = collect_independent_drivers(
            iter_rows(args.input, args.reference_scale))

        for row in iter_rows(args.input, args.reference_scale):
            input_rows += 1
            cohorts.add(row.cohort)
            profilers.add(row.profiler)
            features.add(row.feature)
            features_by_profiler[row.profiler].add(row.feature)
            sample_key = (row.cohort, row.sample_id)
            if sample_key not in sample_indices:
                sample_indices[sample_key] = len(sample_indices)
            sample_bit = 1 << sample_indices[sample_key]
            for target in row.targets:
                if target not in target_indices:
                    target_indices[target] = len(target_indices)
            row_target_bits = sum(1 << target_indices[target] for target in row.targets)
            for target, target_feature in row.target_feature_by_label.items():
                alias_key = (row.profiler, target)
                previous = target_features.setdefault(alias_key, target_feature)
                if previous != target_feature:
                    raise ValueError(f"inconsistent feature alias for {alias_key}")
            if row.observation_id not in seen_observations:
                seen_observations.add(row.observation_id)
                count_key = (
                    row.cohort, row.population, row.profiler, row.target_label,
                    row.dose_level, render(row.nominal_total),
                )
                observation_counts[count_key] += 1
                previous_target_count = observation_target_counts.setdefault(
                    count_key, len(row.targets))
                if previous_target_count != len(row.targets):
                    raise ValueError(f"inconsistent target count within {count_key}")
                if row.population == "independent":
                    target = row.targets[0]
                    design_key = (row.cohort, row.profiler, target)
                    operator_designs[design_key].add(
                        driver_for(drivers, row), sample_bit, 1 << (row.dose_index - 1)
                    )
                    group_key = (row.cohort, row.profiler)
                    independent_counts[group_key] += 1
                    independent_sample_bits[group_key] |= sample_bit
                    independent_target_bits[group_key] |= row_target_bits
                    direct_opportunities[(row.cohort, row.profiler,
                                          row.target_feature_by_label[target])] += 1
            target_rows += row.is_target
            bystander_rows += not row.is_target
            summary_key = (
                row.cohort, row.population, row.profiler, row.target_label,
                row.role, row.dose_level, render(row.nominal_total),
            )
            detection[summary_key].add(row)
            quantitative[summary_key].add(row, args.log2_instability_threshold)
            if row.population == "independent":
                target = row.targets[0]
                operator_key = (row.cohort, row.profiler, target, row.feature)
                operators[operator_key].add(
                    row, driver_for(drivers, row), sample_bit, 1 << (row.dose_index - 1)
                )
                if row.is_target:
                    target_pairs.add(operator_key)
            if row.population == "independent" and not row.is_target:
                reliability_observed[(row.cohort, row.profiler, row.feature)].add(
                    row, args.log2_instability_threshold, sample_bit, row_target_bits
                )
                reliability_target_observed[(row.cohort, row.profiler, row.feature,
                                             row.targets[0])].add(
                    row, args.log2_instability_threshold, sample_bit, row_target_bits
                )
        if input_rows == 0:
            raise ValueError("input contains no data rows")
        single_cohort = len(cohorts) == 1
        if single_cohort and args.cohort_validation == "require_holdout":
            raise ValueError("at least two cohorts are required for holdout validation")

        reliability: dict[tuple[str, str, str], ReliabilityAggregate] = {}
        reliability_by_target: dict[tuple[str, str, str, str], ReliabilityAggregate] = {}
        for cohort in cohorts:
            for profiler in profilers:
                group_key = (cohort, profiler)
                for feature in features_by_profiler[profiler]:
                    eligible = independent_counts[group_key] - direct_opportunities.get(
                        (cohort, profiler, feature), 0)
                    if eligible > 0:
                        reliability[(cohort, profiler, feature)] = complete_reliability(
                            reliability_observed.get((cohort, profiler, feature)), eligible,
                            independent_sample_bits[group_key], independent_target_bits[group_key],
                        )
                for target in sorted({key[2] for key in operator_designs
                                      if key[:2] == group_key}):
                    design = operator_designs[(cohort, profiler, target)]
                    direct_feature = target_features[(profiler, target)]
                    target_bit = 1 << target_indices[target]
                    for feature in features_by_profiler[profiler] - {direct_feature}:
                        reliability_by_target[(cohort, profiler, feature, target)] = complete_reliability(
                            reliability_target_observed.get((cohort, profiler, feature, target)),
                            design.n, design.sample_bits, target_bit,
                        )
        # Reconstitute cohort-level certificates from the completed target-level
        # cells so an implanted feature's own target is genuinely excluded.
        reliability = {}
        for (cohort, profiler, feature, _target), aggregate in reliability_by_target.items():
            destination = reliability.setdefault(
                (cohort, profiler, feature), ReliabilityAggregate())
            destination.merge(aggregate)

        detection_rows = []
        for key in sorted(detection):
            aggregate = detection[key]
            cohort, population, profiler, target_label, role, dose_level, nominal = key
            physical_n = observation_counts[(cohort, population, profiler, target_label,
                                             dose_level, nominal)]
            direct_count = observation_target_counts[(cohort, population, profiler,
                                                       target_label, dose_level, nominal)]
            eligible = (physical_n * direct_count if role == "implanted_target" else
                        physical_n * (len(features_by_profiler[profiler]) - direct_count))
            implicit_absent = eligible - aggregate.n
            if implicit_absent < 0:
                raise ValueError(
                    f"explicit detection rows exceed physical opportunities for {key}: "
                    f"rows={aggregate.n}, eligible={eligible}, observations={physical_n}"
                )
            detection_rows.append(dict(zip(DETECTION_FIELDS[:7], key)) | {
                "eligible_contexts": eligible,
                "baseline_absent_contexts": aggregate.absent + implicit_absent,
                "baseline_present_contexts": aggregate.present,
                "remained_absent": aggregate.remained_absent + implicit_absent,
                "unexpected_appearances": aggregate.appearances,
                "remained_detected": aggregate.remained_detected,
                "unexpected_dropouts": aggregate.dropouts,
                "appearance_rate": render(divide(aggregate.appearances,
                                                  aggregate.absent + implicit_absent)),
                "dropout_rate": render(divide(aggregate.dropouts, aggregate.present)),
            })

        quantitative_rows = []
        for key in sorted(quantitative):
            aggregate = quantitative[key]
            quantitative_rows.append(dict(zip(QUANTITATIVE_FIELDS[:7], key)) | {
                "eligible_contexts": aggregate.n,
                "mean_expected_abundance": render(mean(aggregate.expected_sum, aggregate.n)),
                "mean_observed_abundance": render(mean(aggregate.observed_sum, aggregate.n)),
                "mean_signed_error": render(mean(aggregate.signed_sum, aggregate.n)),
                "mean_absolute_error": render(mean(aggregate.absolute_sum, aggregate.n)),
                "rmse": render(math.sqrt(aggregate.squared_sum / aggregate.n)),
                "mean_log2_response_error": render(mean(aggregate.log2_sum, aggregate.log2_n)),
                "quantitative_instability_rate": render(divide(aggregate.unstable, aggregate.log2_n)),
            })

        operator_rows = []
        slopes: dict[tuple[str, str, str, str], float] = {}
        complete_operators: dict[tuple[str, str, str, str], RegressionAggregate] = {}
        for design_key in sorted(operator_designs):
            cohort, profiler, target = design_key
            for feature in sorted(features_by_profiler[profiler]):
                key = (cohort, profiler, target, feature)
                aggregate = complete_regression(operators.get(key), operator_designs[design_key])
                complete_operators[key] = aggregate
                slope = aggregate.slope_origin()
                if slope is not None:
                    slopes[key] = slope
                operator_rows.append({
                    "holdout_cohort": "NONE", "training_cohorts": cohort,
                    "cohort": cohort, "analysis_population": "independent",
                    "profiler": profiler, "target_label": target, "feature": feature,
                    "feature_role": "implanted_target" if feature == target_features[(profiler, target)] else "bystander",
                    "eligible_contexts": aggregate.n,
                    "distinct_samples": aggregate.sample_bits.bit_count(),
                    "distinct_doses": aggregate.dose_bits.bit_count(),
                    "operator_slope": render(slope),
                    "operator_intercept": render(aggregate.intercept()),
                    "operator_r_squared": render(aggregate.r_squared_origin()),
                    "mean_signed_error": render(mean(aggregate.error_sum, aggregate.n)),
                    "mean_absolute_error": render(mean(aggregate.absolute_error_sum, aggregate.n)),
                    "detection_transition_rate": render(divide(aggregate.transitions, aggregate.n)),
                })

        # Exact, same-sample independent-to-community superposition is an
        # out-of-core ordinal join rather than an in-memory feature lookup.
        superposition_rows = exact_superposition_rows(args.input, args.reference_scale)

        # Cross-cohort operator prediction on unseen independent perturbations.
        cross_slopes: dict[tuple[str, str, str, str], float] = {}
        for holdout in sorted(cohorts):
            for profiler in sorted(profilers):
                relevant = {(target, feature) for cohort, p, target, feature in complete_operators
                            if cohort != holdout and p == profiler}
                for target, feature in relevant:
                    combined = RegressionAggregate()
                    for cohort in sorted(cohorts - {holdout}):
                        source = complete_operators.get((cohort, profiler, target, feature))
                        if not source:
                            continue
                        for name in ("n", "sum_x", "sum_y", "sum_xx", "sum_xy", "sum_yy",
                                     "error_sum", "absolute_error_sum", "transitions"):
                            setattr(combined, name, getattr(combined, name) + getattr(source, name))
                    value = combined.slope_origin()
                    if value is not None:
                        cross_slopes[(holdout, profiler, target, feature)] = value

        holdout_predictions: dict[tuple[str, str, str], PredictionAggregate] = defaultdict(PredictionAggregate)
        for row in iter_rows(args.input, args.reference_scale):
            if row.population != "independent":
                continue
            target = row.targets[0]
            slope = cross_slopes.get((row.cohort, row.profiler, target, row.feature))
            if slope is None:
                continue
            # Selected retained baseline + selected driver * fitted slope.
            # Neither term may be reconstructed from read-proportional metadata.
            prediction = row.retained_baseline + driver_for(drivers, row) * slope
            sample_bit = 1 << sample_indices[(row.cohort, row.sample_id)]
            holdout_predictions[(row.cohort, row.profiler, row.role)].add(
                row, prediction, sample_bit
            )
        holdout_rows = []
        for key in sorted(holdout_predictions):
            holdout, profiler, role = key
            aggregate = holdout_predictions[key]
            composition_mae = aggregate.composition_abs / aggregate.n
            operator_mae = aggregate.operator_abs / aggregate.n
            holdout_rows.append({
                "holdout_cohort": holdout,
                "training_cohorts": ",".join(sorted(cohorts - {holdout})),
                "profiler": profiler,
                "feature_role": role,
                "eligible_contexts": aggregate.n,
                "distinct_samples": aggregate.sample_bits.bit_count(),
                "composition_only_mae": render(composition_mae),
                "operator_prediction_mae": render(operator_mae),
                "operator_prediction_rmse": render(math.sqrt(aggregate.operator_sq / aggregate.n)),
                "operator_r_squared": render(aggregate.r_squared()),
                "improvement_over_composition": render(composition_mae - operator_mae),
            })

        certificate_rows: list[dict[str, object]] = []
        certificate_scores: dict[tuple[str, str, str], tuple[float | None, float | None]] = {}
        holdout_labels: list[str | None] = [None]
        if not single_cohort:
            holdout_labels.extend(sorted(cohorts))
        for holdout in holdout_labels:
            training = sorted(cohorts if holdout is None else cohorts - {holdout})
            holdout_label = "NONE" if holdout is None else holdout
            for profiler in sorted(profilers):
                profiler_features = sorted({feature for cohort, p, feature in reliability
                                            if p == profiler and cohort in training})
                for feature in profiler_features:
                    trained = combine_reliability(reliability, training, profiler, feature)
                    evaluated = (ReliabilityAggregate() if holdout is None else
                                 combine_reliability(reliability, [holdout], profiler, feature))
                    components = trained.components()
                    score = trained.score()
                    evaluation_score = evaluated.score()
                    lower, upper = trained.credible_interval()
                    certificate_scores[(holdout_label, profiler, feature)] = (score, evaluation_score)
                    certificate_rows.append({
                        "holdout_cohort": holdout_label,
                        "training_cohorts": ",".join(training),
                        "profiler": profiler,
                        "feature": feature,
                        "eligible_contexts": trained.n,
                        "eligible_targets": trained.target_bits.bit_count(),
                        "eligible_samples": trained.sample_bits.bit_count(),
                        "baseline_absent_contexts": trained.absent,
                        "baseline_present_contexts": trained.present,
                        "quantitative_contexts": trained.quantitative,
                        "unexpected_appearance_rate": render(divide(trained.appearances, trained.absent)),
                        "unexpected_dropout_rate": render(divide(trained.dropouts, trained.present)),
                        "quantitative_instability_rate": render(divide(trained.unstable, trained.quantitative)),
                        "mean_absolute_error": render(mean(trained.absolute_error_sum, trained.quantitative)),
                        "hallucination_reliability": render(components[0]),
                        "stability_reliability": render(components[1]),
                        "quantitative_reliability": render(components[2]),
                        "overall_reliability": render(score),
                        "reliability_ci_lower": render(lower),
                        "reliability_ci_upper": render(upper),
                        "evaluation_contexts": evaluated.n if holdout is not None else "",
                        "evaluation_overall_reliability": render(evaluation_score),
                    })

        cohort_holdout_rows = []
        for holdout in ([] if single_cohort else sorted(cohorts)):
            for profiler in sorted(profilers):
                trained_features = {
                    feature for h, p, feature in certificate_scores
                    if h == holdout and p == profiler
                }
                evaluated_features = {
                    feature for cohort, p, feature in reliability
                    if cohort == holdout and p == profiler
                }
                shared = sorted(trained_features & evaluated_features)
                pairs = [certificate_scores[(holdout, profiler, feature)] for feature in shared]
                pairs = [(left, right) for left, right in pairs if left is not None and right is not None]
                left = [pair[0] for pair in pairs]
                right = [pair[1] for pair in pairs]
                cohort_holdout_rows.append({
                    "holdout_cohort": holdout,
                    "training_cohorts": ",".join(sorted(cohorts - {holdout})),
                    "profiler": profiler,
                    "features_trained": len(trained_features),
                    "features_evaluated": len(evaluated_features),
                    "features_shared": len(pairs),
                    "spearman_reliability": render(correlation(ranks(left), ranks(right))),
                    "mean_absolute_reliability_gap": render(
                        mean(sum(abs(x - y) for x, y in pairs), len(pairs))
                    ),
                    "mean_training_reliability": render(mean(sum(left), len(left))),
                    "mean_evaluation_reliability": render(mean(sum(right), len(right))),
                })

        # Leave one implanted taxon out: predict its perturbation behaviour
        # using only the other targets, never the held-out target itself.
        target_holdout_rows: list[dict[str, object]] = []
        all_targets = sorted(target_indices)
        for holdout in ([] if single_cohort else sorted(cohorts)):
            training_cohorts = sorted(cohorts - {holdout})
            for profiler in sorted(profilers):
                for heldout_target in all_targets:
                    training_targets = [item for item in all_targets if item != heldout_target]
                    direct_feature = target_features.get((profiler, heldout_target))
                    for feature in sorted(features_by_profiler[profiler]):
                        if feature == direct_feature:
                            continue
                        predicted, training_n = reliability_score_for_targets(
                            reliability_by_target, training_cohorts, profiler, feature,
                            training_targets,
                        )
                        observed, evaluation_n = reliability_score_for_targets(
                            reliability_by_target, [holdout], profiler, feature,
                            [heldout_target],
                        )
                        if predicted is None or observed is None:
                            continue
                        target_holdout_rows.append({
                            "holdout_cohort": holdout,
                            "training_cohorts": ",".join(training_cohorts),
                            "profiler": profiler,
                            "validation_type": "leave_one_target_and_cohort_out",
                            "heldout_target": heldout_target, "feature": feature,
                            "predicted_reliability": render(predicted),
                            "observed_reliability": render(observed),
                            "training_contexts": training_n,
                            "evaluation_contexts": evaluation_n,
                        })

        # How many independently implanted sentinel taxa are needed before the
        # feature-risk ranking saturates? Subsets are deterministic and evaluated
        # in an untouched cohort.
        panel_rows: list[dict[str, object]] = []
        randomizer = random.Random(20260915)
        panel_sizes = [size for size in (1, 2, 3, 5, 7, 10) if size <= len(all_targets)]
        for holdout in ([] if single_cohort else sorted(cohorts)):
            training_cohorts = sorted(cohorts - {holdout})
            for profiler in sorted(profilers):
                feature_pool = sorted(features_by_profiler[profiler])
                full_scores = {
                    feature: reliability_score_for_targets(
                        reliability_by_target, training_cohorts, profiler, feature,
                        all_targets,
                    )[0] for feature in feature_pool
                }
                evaluation_scores = {
                    feature: reliability_score_for_targets(
                        reliability_by_target, [holdout], profiler, feature, all_targets,
                    )[0] for feature in feature_pool
                }
                valid_full = sorted(feature for feature in feature_pool
                                    if full_scores[feature] is not None and
                                    evaluation_scores[feature] is not None)
                if not valid_full:
                    continue
                ordered_risk = sorted(valid_full, key=lambda item: evaluation_scores[item])
                high_risk = set(ordered_risk[:max(1, math.ceil(len(ordered_risk) * 0.2))])
                for panel_size in panel_sizes:
                    combinations = list(itertools.combinations(all_targets, panel_size))
                    if len(combinations) > args.panel_replicates:
                        combinations = randomizer.sample(combinations, args.panel_replicates)
                    for replicate, subset in enumerate(sorted(combinations), start=1):
                        subset_scores = {
                            feature: reliability_score_for_targets(
                                reliability_by_target, training_cohorts, profiler,
                                feature, subset,
                            )[0] for feature in valid_full
                        }
                        usable = [feature for feature in valid_full
                                  if subset_scores[feature] is not None]
                        subset_risk = set(sorted(
                            usable, key=lambda item: subset_scores[item]
                        )[:max(1, math.ceil(len(usable) * 0.2))])
                        panel_rows.append({
                            "holdout_cohort": holdout,
                            "training_cohorts": ",".join(training_cohorts),
                            "profiler": profiler, "panel_size": panel_size,
                            "panel_replicate": replicate,
                            "panel_targets": ";".join(subset),
                            "features_evaluated": len(usable),
                            "rank_correlation_full_panel": render(correlation(
                                ranks([subset_scores[item] for item in usable]),
                                ranks([full_scores[item] for item in usable]),
                            )),
                            "heldout_reliability_mae": render(mean(sum(
                                abs(subset_scores[item] - evaluation_scores[item])
                                for item in usable), len(usable))),
                            "high_risk_taxon_recall": render(divide(
                                len(high_risk & subset_risk), len(high_risk))),
                        })

        outputs = [
            ("detection_response_summary.tsv", DETECTION_FIELDS, detection_rows),
            ("quantitative_response_summary.tsv", QUANTITATIVE_FIELDS, quantitative_rows),
            ("response_operator.tsv", OPERATOR_FIELDS, operator_rows),
            ("superposition_summary.tsv", SUPERPOSITION_FIELDS, superposition_rows),
            ("heldout_operator_validation.tsv", HOLDOUT_FIELDS, holdout_rows),
            ("reliability_certificates.tsv", CERTIFICATE_FIELDS, certificate_rows),
            ("cohort_holdout_summary.tsv", COHORT_HOLDOUT_FIELDS, cohort_holdout_rows),
            ("heldout_validation.tsv", TARGET_HOLDOUT_FIELDS, target_holdout_rows),
            ("panel_saturation.tsv", PANEL_FIELDS, panel_rows),
        ]
        output_paths: list[Path] = []
        for filename, fields, rows_to_write in outputs:
            path = args.outdir / filename
            write_tsv(path, fields, rows_to_write)
            output_paths.append(path)

        if args.disease_results:
            if not args.disease_results.is_file() or not args.disease_results.stat().st_size:
                raise ValueError("disease results are missing or empty")
            biomarker_path = args.outdir / "biomarker_replication.tsv"
            biomarker_rows = build_biomarker_replication(
                args.disease_results, args.outdir / "reliability_certificates.tsv", cohorts,
                args.feature_aliases,
            )
            write_tsv(biomarker_path, BIOMARKER_FIELDS, biomarker_rows)
            output_paths.append(biomarker_path)

        validation = args.outdir / "perturbation_response_validation.tsv"
        validation_rows = [
            {"metric": "input_rows", "value": input_rows},
            {"metric": "input_sha256", "value": sha256(args.input)},
            {"metric": "cohorts", "value": len(cohorts)},
            {"metric": "cohort_holdout_status", "value": (
                "NOT_APPLICABLE_SINGLE_COHORT" if single_cohort else "PASS"
            )},
            {"metric": "profilers", "value": len(profilers)},
            {"metric": "features", "value": len(features)},
            {"metric": "implanted_target_rows", "value": target_rows},
            {"metric": "target_excluded_bystander_rows", "value": bystander_rows},
            {"metric": "operator_rows", "value": len(operator_rows)},
            {"metric": "superposition_summary_rows", "value": len(superposition_rows)},
            {"metric": "reliability_certificate_rows", "value": len(certificate_rows)},
            {"metric": "status", "value": "PASS"},
        ]
        write_tsv(validation, ["metric", "value"], validation_rows)
        output_paths.append(validation)
        settings = args.outdir / "perturbation_response_settings.tsv"
        write_tsv(settings, ["setting", "value"], [
            {"setting": "schema", "value": "paired-feature-response-v1"},
            {"setting": "nominal_dose_mapping", "value": "population_specific_dose_rank"},
            {"setting": "operator_fit", "value": "through_origin"},
            {"setting": "reliability_target_exclusion", "value": "all_implanted_features"},
            {"setting": "binary_rate_prior", "value": "Jeffreys_beta_0.5_0.5"},
            {"setting": "log2_instability_threshold", "value": render(args.log2_instability_threshold)},
            {"setting": "cohort_validation", "value": args.cohort_validation},
            {"setting": "status", "value": "DEVELOPMENT_ONLY"},
        ])
        output_paths.append(settings)
        development = args.outdir / "DEVELOPMENT_ONLY.txt"
        development.write_text(
            "status=DEVELOPMENT_ONLY\n"
            "use_for_manuscript=NO\n"
            "reason=exploratory target-agnostic response methodology requires prespecified final validation\n",
            encoding="utf-8",
        )
        output_paths.append(development)
        seal = args.outdir / "perturbation_response.sha256"
        seal.write_text(
            "".join(f"{sha256(path)}  {path.name}\n" for path in sorted(output_paths)),
            encoding="utf-8",
        )
        success = args.outdir / "SUCCESS"
        success.write_text(
            f"status\tPASS\ninput_rows\t{input_rows}\noperator_rows\t{len(operator_rows)}\n"
            f"certificate_rows\t{len(certificate_rows)}\n",
            encoding="utf-8",
        )
        print(f"[PASS] Perturbation response analysis: {args.outdir}")
        print(f"[INFO] Input rows: {input_rows}; operator rows: {len(operator_rows)}; "
              f"certificates: {len(certificate_rows)}")
    except (ValueError, OSError, KeyError) as error:
        raise SystemExit(f"[ERROR] {error}") from error


if __name__ == "__main__":
    main()
