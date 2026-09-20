#!/usr/bin/env python3
"""Audit whether MetaPhlAn's genome-size-dependent recovery bias is removed on
the genome-equivalent profiler scale.

The read-proportional reference treats the implanted signal as a read fraction,
so a marker-length-normalised profiler should under-report large genomes in
proportion to their length:

    log2(observed / expected_read) ~= constant - log2(G_t)

i.e. a slope near **-1** against log2 of the implanted target's genome size.
Replacing the implanted read fraction with the genome-equivalent implanted
signal q_it/D_i should remove that dependence, giving a slope near **0**.

This script only measures that change. It never assumes, forces or rewards
either result, and it never estimates a genome size from a recovery outcome:
lengths come exclusively from the FASTA-measured `target_genome_sizes.tsv`.

Statistical unit
----------------
The implanted taxon, not the observation row. Each target contributes one
median log2 recovery ratio, and the regression is fitted across the ten target
points. Observation rows are repeated measures of the same ten taxa; regressing
them directly would be pseudoreplication. The bootstrap therefore resamples
**targets**, and is descriptive uncertainty over the ten implanted taxa only.

DEVELOPMENT_ONLY: Yachida, MetaPhlAn 4.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import random
from collections import defaultdict
from pathlib import Path

# The frozen implanted panel. Any deviation is a hard failure.
EXPECTED_LABELS = ("Bfrag", "Csym", "Dpne", "Fnuc", "Hhat",
                   "Pmic", "Pana", "Psto", "Porp", "Pint")
EXPECTED_TARGET_COUNT = len(EXPECTED_LABELS)

# Plausibility window already enforced by build_target_genome_sizes.py.
MIN_GENOME_BP = 100_000
MAX_GENOME_BP = 50_000_000

METAPHLAN = "metaphlan4"
BRACKEN = "kraken2_bracken"
REQUIRED_PROFILERS = (BRACKEN, METAPHLAN)

# This development audit is Yachida-only and proves it rather than asserting it.
# The spelling is the canonical lower-case cohort identifier used throughout the
# pipeline (build_perturbation_response_input.py tests `cohort = 'yachida'`).
DEVELOPMENT_COHORT = "yachida"

# Selected-reference contract, per arm. `reference_type` is the ORIGINAL
# row-level provenance of the source profile and is NOT the estimand: a
# MetaPhlAn row keeps `genome_equivalent` even in the read-proportional arm.
# The estimand is fixed by (reference_scale, selected_reference_type).
ARM_SELECTION = {
    "read_proportional_sensitivity": ("read_proportional", "read_proportional"),
    "genome_equivalent_primary": ("profiler_scale", "profiler_scale_primary"),
}
# The selected estimand each arm must resolve to for MetaPhlAn.
ARM_SELECTED_ESTIMAND = {
    "read_proportional_sensitivity": "read_proportional",
    "genome_equivalent_primary": "genome_equivalent",
}
METAPHLAN_ROW_REFERENCE = "genome_equivalent"

# Supplied absolute relative error must reproduce abs(ratio - 1). The window is
# sized for 17-significant-digit TSV round-tripping, not for rescuing a wrong
# value: the supplied number is validated, never recomputed and substituted.
ERROR_TOLERANCE_POLICY = (
    "abs(actual - expected) <= max(1e-12, 1e-9 * max(abs(actual), "
    "abs(expected), 1)); expected = abs(observed_over_expected - 1)")
ERROR_ABSOLUTE_TOLERANCE = 1e-12
ERROR_RELATIVE_TOLERANCE = 1e-9

BOOTSTRAP_SEED = 20260920
DEFAULT_BOOTSTRAP_REPLICATES = 10000
MIN_VALID_BOOTSTRAP_FRACTION = 0.95
CI_LOWER_Q = 0.025
CI_UPPER_Q = 0.975

SCOPES = ("independent", "community", "pooled")
POPULATION_SCOPES = ("independent", "community")

# arm -> (ratio column, absolute-relative-error column, provenance columns,
#         expected slope against log2 genome size). The provenance quadruple is
# (row reference type, selected reference type, reference scale, selected
# estimand); the ambiguous legacy `*_reference_type` alias is deliberately not
# used. The comparator's `*_selected_estimand` is redundant with the pair
# (scale, selected type), but a contradiction between them means the upstream
# provenance is untrustworthy, so it fails closed instead of being ignored.
ARMS = {
    "read_proportional_sensitivity": (
        "sensitivity_observed_over_expected",
        "sensitivity_absolute_relative_error",
        ("sensitivity_row_reference_type",
         "sensitivity_selected_reference_type",
         "sensitivity_reference_scale",
         "sensitivity_selected_estimand"),
        -1.0,
    ),
    "genome_equivalent_primary": (
        "primary_observed_over_expected",
        "primary_absolute_relative_error",
        ("primary_row_reference_type",
         "primary_selected_reference_type",
         "primary_reference_scale",
         "primary_selected_estimand"),
        0.0,
    ),
}
ARM_ORDER = ("read_proportional_sensitivity", "genome_equivalent_primary")

# Physical observation identity, as written by compare_target_recovery_references.py.
OBSERVATION_KEY = (
    "cohort", "study", "sample_id", "condition", "analysis_population",
    "profiler", "baseline_id", "target_feature", "dose_index",
    "nominal_target_fraction",
)

TARGET_LEVEL_FIELDS = [
    "scope", "arm", "target_label", "target_feature", "genome_size_bp",
    "log2_genome_size_bp", "eligible_observations", "excluded_observations",
    "exclusion_reasons", "median_observed_over_expected",
    "median_log2_observed_over_expected", "median_absolute_relative_error",
    # `reference_type` is the SELECTED estimand. The original row provenance and
    # the selection labels are kept beside it so nothing is lost.
    "reference_type", "row_reference_type", "selected_reference_type",
    "reference_scale",
]
REGRESSION_FIELDS = [
    "scope", "arm", "n_targets", "slope", "intercept", "r_squared",
    "pearson_correlation", "expected_slope", "absolute_distance_from_expected",
    "bootstrap_replicates_requested", "bootstrap_replicates_valid",
    "bootstrap_replicates_discarded", "bootstrap_ci_lower", "bootstrap_ci_upper",
    "bootstrap_seed",
]
PAIRED_CHANGE_FIELDS = [
    "scope", "target_label", "target_feature", "genome_size_bp",
    "sensitivity_median_log2_ratio", "primary_median_log2_ratio",
    "primary_minus_sensitivity_log2_ratio",
    "sensitivity_median_absolute_relative_error",
    "primary_median_absolute_relative_error",
    "primary_minus_sensitivity_absolute_relative_error",
]


class AuditError(Exception):
    """A fail-closed audit gate."""


# --------------------------------------------------------------------- utility
def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def render(value: float | None) -> str:
    """Undefined quantities are written as an empty cell, never as 0 or NaN."""
    if value is None or not math.isfinite(value):
        return ""
    return format(value, ".17g")


def median(values: list[float]) -> float:
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2.0


def percentile(ordered: list[float], quantile: float) -> float:
    """Linear-interpolation percentile of an already sorted list.

    Implemented here so the audit adds no SciPy/NumPy dependency and stays
    bit-identical across runs and platforms.
    """
    if not ordered:
        raise AuditError("percentile of an empty sample")
    if len(ordered) == 1:
        return ordered[0]
    position = quantile * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[int(position)]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def require_file(path: Path, label: str) -> None:
    if not path.is_file() or not path.stat().st_size:
        raise AuditError(f"missing or empty {label}: {path}")


def read_tsv(path: Path, label: str) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = [row for row in reader]
    if not fields:
        raise AuditError(f"{label} has no header: {path}")
    if not rows:
        raise AuditError(f"{label} has no data rows: {path}")
    return fields, rows


def read_csv_rows(path: Path, label: str) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        rows = [row for row in reader]
    if not fields:
        raise AuditError(f"{label} has no header: {path}")
    if not rows:
        raise AuditError(f"{label} has no data rows: {path}")
    return fields, rows


# ------------------------------------------------------------ exact target map
def load_genome_sizes(path: Path, panel: dict[str, str]) -> dict[str, float]:
    """target_label -> FASTA-measured genome size, validated against the panel."""
    fields, rows = read_tsv(path, "target genome-size table")
    for required in ("target_label", "target_feature", "genome_size_bp"):
        if required not in fields:
            raise AuditError(
                f"target genome-size table lacks column {required!r}")
    sizes: dict[str, float] = {}
    for line, row in enumerate(rows, start=2):
        label = (row.get("target_label") or "").strip()
        if not label:
            raise AuditError(f"target genome-size table line {line}: blank target_label")
        raw = (row.get("genome_size_bp") or "").strip()
        try:
            size = float(raw)
        except (TypeError, ValueError) as error:
            raise AuditError(
                f"target genome-size table line {line}: nonnumeric "
                f"genome_size_bp {raw!r} for {label}") from error
        if not math.isfinite(size):
            raise AuditError(
                f"target genome-size table line {line}: nonfinite genome_size_bp "
                f"for {label}")
        if size <= 0:
            raise AuditError(
                f"target genome-size table line {line}: non-positive "
                f"genome_size_bp for {label}")
        if not MIN_GENOME_BP <= size <= MAX_GENOME_BP:
            raise AuditError(
                f"{label}: genome size {render(size)} bp is outside the "
                f"plausibility bounds [{MIN_GENOME_BP}, {MAX_GENOME_BP}] "
                "enforced by build_target_genome_sizes.py")
        if label in sizes and sizes[label] != size:
            raise AuditError(
                f"{label}: inconsistent duplicate genome sizes "
                f"{render(sizes[label])} and {render(size)}")
        sizes[label] = size
        # The genome-size row must belong to the same frozen panel entry,
        # matched by exact string equality on the taxon name.
        feature = (row.get("target_feature") or "").strip()
        if label in panel and feature != panel[label]:
            raise AuditError(
                f"{label}: genome-size table target_feature {feature!r} does not "
                f"exactly equal the spike-panel taxon_name {panel[label]!r}")
    return sizes


def load_spike_panel(path: Path) -> dict[str, str]:
    """target_label -> canonical taxon name, from the frozen spike panel."""
    fields, rows = read_tsv(path, "spike panel")
    for required in ("label", "taxon_name"):
        if required not in fields:
            raise AuditError(f"spike panel lacks column {required!r}")
    panel: dict[str, str] = {}
    for line, row in enumerate(rows, start=2):
        label = (row.get("label") or "").strip()
        taxon = (row.get("taxon_name") or "").strip()
        if not label or not taxon:
            raise AuditError(f"spike panel line {line}: blank label or taxon_name")
        if label in panel:
            raise AuditError(f"spike panel: duplicate label {label!r}")
        panel[label] = taxon
    return panel


def load_metaphlan_aliases(path: Path, panel: dict[str, str]) -> dict[str, str]:
    """target_label -> the exact MetaPhlAn feature string for that taxon.

    `examples/spike_taxon_aliases.csv` is the same frozen table that produced the
    `target_feature` values in the upstream comparison, so this is an exact
    equality join on (canonical taxon name, tool). Nothing is normalised,
    lower-cased, prefix-matched or guessed.
    """
    fields, rows = read_csv_rows(path, "feature alias table")
    for required in ("canonical", "alias", "tool"):
        if required not in fields:
            raise AuditError(f"feature alias table lacks column {required!r}")
    alias_by_canonical: dict[str, str] = {}
    for line, row in enumerate(rows, start=2):
        if (row.get("tool") or "").strip() != METAPHLAN:
            continue
        canonical = (row.get("canonical") or "").strip()
        alias = (row.get("alias") or "").strip()
        if not canonical or not alias:
            raise AuditError(
                f"feature alias table line {line}: blank canonical or alias")
        previous = alias_by_canonical.setdefault(canonical, alias)
        if previous != alias:
            raise AuditError(
                f"feature alias table: {canonical!r} has conflicting "
                f"{METAPHLAN} aliases {previous!r} and {alias!r}")
    aliases: dict[str, str] = {}
    for label, taxon in panel.items():
        if taxon not in alias_by_canonical:
            raise AuditError(
                f"{label}: no exact {METAPHLAN} alias for the spike-panel taxon "
                f"{taxon!r} in the feature alias table")
        aliases[label] = alias_by_canonical[taxon]
    return aliases


def build_feature_index(panel: dict[str, str],
                        aliases: dict[str, str]) -> dict[str, str]:
    """Exact feature string -> target_label.

    Both the canonical spike-panel taxon name and the frozen MetaPhlAn alias are
    accepted, because the upstream comparison carries the profiler alias while
    the genome-size table carries the canonical name. The index must be
    injective: if one feature string could denote two implanted taxa the audit
    fails rather than choosing one.
    """
    index: dict[str, str] = {}
    for source in (panel, aliases):
        for label, feature in source.items():
            owner = index.setdefault(feature, label)
            if owner != label:
                raise AuditError(
                    f"feature {feature!r} maps to more than one target label: "
                    f"{owner} and {label}")
    return index


# ------------------------------------------------------------------- gates
def check_bracken_validation(path: Path) -> str:
    """The upstream Bracken-identity proof is a required input gate."""
    _, rows = read_tsv(path, "comparison validation table")
    values = {}
    for row in rows:
        metric = (row.get("metric") or "").strip()
        if not metric:
            continue
        if metric in values:
            raise AuditError(
                f"comparison validation table has duplicate metric {metric!r}")
        values[metric] = (row.get("value") or "").strip()
    if "bracken_identical" not in values:
        raise AuditError(
            "comparison validation table has no bracken_identical metric; the "
            "Bracken-identity gate cannot be inferred from a filename")
    if values["bracken_identical"] != "PASS":
        raise AuditError(
            "Bracken identity gate is "
            f"{values['bracken_identical']!r}, not PASS")
    return values["bracken_identical"]


def arm_provenance(row: dict[str, str], arm: str,
                   columns: tuple[str, str, str, str],
                   line: int) -> tuple[str, str, str, str]:
    """Validate one arm's explicit reference provenance for a MetaPhlAn row.

    Returns (row reference type, selected reference type, reference scale,
    selected estimand). The original row type is NOT the estimand: a MetaPhlAn
    row keeps ``genome_equivalent`` in the read-proportional arm because it
    records the source profile, so the estimand is taken from
    (reference_scale, selected_reference_type) and cross-checked against the
    arm's contract.
    """
    row_column, selected_column, scale_column, estimand_column = columns
    row_type = row[row_column].strip()
    selected = row[selected_column].strip()
    scale = row[scale_column].strip()
    if row_type != METAPHLAN_ROW_REFERENCE:
        raise AuditError(
            f"reference comparison line {line}: {row_column} is {row_type!r}; "
            f"MetaPhlAn source rows must be {METAPHLAN_ROW_REFERENCE!r}")
    required = ARM_SELECTION[arm]
    if (scale, selected) != required:
        raise AuditError(
            f"reference comparison line {line}: the {arm} arm carries "
            f"reference_scale={scale!r} and selected_reference_type={selected!r}; "
            f"{required[0]!r} and {required[1]!r} are required")
    estimand = ARM_SELECTED_ESTIMAND[arm]
    # The comparator also states the estimand outright. Redundant, but a
    # disagreement means one of the two provenance paths is wrong, and guessing
    # which would defeat the point of carrying provenance at all.
    supplied = row[estimand_column].strip()
    if supplied != estimand:
        raise AuditError(
            f"reference comparison line {line}: {estimand_column} is "
            f"{supplied!r} but reference_scale={scale!r} with "
            f"selected_reference_type={selected!r} derives {estimand!r}; the "
            "comparator's reference provenance contradicts itself")
    return row_type, selected, scale, estimand


def ratio_problem(raw: str) -> tuple[float | None, str | None]:
    """Classify a recovery ratio that is present but possibly unusable.

    A *missing* cell is a structural defect and is rejected by the caller before
    this point; only a *present* value can be excluded. No pseudocount is ever
    added to rescue a zero, negative or nonfinite ratio.
    """
    text = (raw or "").strip()
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None, "nonnumeric"
    if math.isnan(value) or math.isinf(value):
        return None, "nonfinite"
    if value <= 0:
        return None, "nonpositive"
    return value, None


# ------------------------------------------------------------------ regression
def fit(points: list[tuple[float, float]]) -> dict[str, float | None]:
    """Ordinary least squares of y on x across the target-level points."""
    n = len(points)
    if n < 2:
        raise AuditError(f"cannot fit a slope through {n} target point(s)")
    mean_x = sum(x for x, _ in points) / n
    mean_y = sum(y for _, y in points) / n
    sxx = sum((x - mean_x) ** 2 for x, _ in points)
    syy = sum((y - mean_y) ** 2 for _, y in points)
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in points)
    if sxx == 0:
        raise AuditError(
            "all targets share one genome size; the slope is undefined")
    slope = sxy / sxx
    intercept = mean_y - slope * mean_x
    # A zero-variance response is legitimate (a perfectly corrected arm can be
    # exactly flat). The slope is still 0; R^2 and the correlation are undefined
    # and are reported as such rather than as 0, 1 or NaN.
    r_squared = None
    correlation = None
    if syy > 0:
        residual = sum((y - (intercept + slope * x)) ** 2 for x, y in points)
        r_squared = 1.0 - residual / syy
        correlation = sxy / math.sqrt(sxx * syy)
    return {"slope": slope, "intercept": intercept, "r_squared": r_squared,
            "pearson_correlation": correlation}


def bootstrap_slope_ci(points: list[tuple[float, float]], replicates: int,
                       seed: int) -> dict[str, object]:
    """Target-cluster percentile bootstrap: resample the ten taxa, not the rows.

    A fresh generator seeded with `seed` is used for every regression, so the
    same target resamples are applied to both arms. That makes the two arms'
    intervals directly comparable and the whole audit byte-identical on rerun.
    """
    generator = random.Random(seed)
    n = len(points)
    slopes: list[float] = []
    discarded = 0
    for _ in range(replicates):
        sample = [points[generator.randrange(n)] for _ in range(n)]
        if len({x for x, _ in sample}) < 2:
            discarded += 1
            continue
        slopes.append(fit(sample)["slope"])
    valid = len(slopes)
    if valid < math.ceil(MIN_VALID_BOOTSTRAP_FRACTION * replicates):
        raise AuditError(
            f"only {valid} of {replicates} bootstrap replicates were valid; at "
            f"least {MIN_VALID_BOOTSTRAP_FRACTION:.0%} are required")
    slopes.sort()
    return {
        "bootstrap_replicates_requested": replicates,
        "bootstrap_replicates_valid": valid,
        "bootstrap_replicates_discarded": discarded,
        "bootstrap_ci_lower": percentile(slopes, CI_LOWER_Q),
        "bootstrap_ci_upper": percentile(slopes, CI_UPPER_Q),
        "bootstrap_seed": seed,
    }


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


# ----------------------------------------------------------------------- main
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--comparison", type=Path, required=True,
        help="reference_comparison/target_recovery_reference_comparison.tsv")
    parser.add_argument(
        "--comparison-validation", type=Path, required=True,
        help="reference_comparison/reference_comparison_validation.tsv; its "
             "bracken_identical metric must read PASS.")
    parser.add_argument("--target-genome-sizes", type=Path, required=True)
    parser.add_argument("--spike-panel", type=Path, required=True)
    parser.add_argument(
        "--feature-aliases", type=Path, required=True,
        help="Frozen canonical/alias/tool table (examples/spike_taxon_aliases.csv) "
             "used for the exact taxon-name-to-profiler-feature join.")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--bootstrap-replicates", type=int,
                        default=DEFAULT_BOOTSTRAP_REPLICATES)
    return parser.parse_args()


def run(args: argparse.Namespace) -> None:
    inputs = {
        "comparison": args.comparison,
        "comparison_validation": args.comparison_validation,
        "target_genome_sizes": args.target_genome_sizes,
        "spike_panel": args.spike_panel,
        "feature_aliases": args.feature_aliases,
    }
    for label, path in inputs.items():
        require_file(path, label.replace("_", " "))
    if args.bootstrap_replicates < 2:
        raise AuditError("--bootstrap-replicates must be at least 2")
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise AuditError(f"OUTDIR must be new or empty: {args.outdir}")

    bracken_gate = check_bracken_validation(args.comparison_validation)

    panel = load_spike_panel(args.spike_panel)
    unexpected_panel = sorted(set(panel) - set(EXPECTED_LABELS))
    missing_panel = sorted(set(EXPECTED_LABELS) - set(panel))
    if unexpected_panel or missing_panel:
        raise AuditError(
            "spike panel labels do not match the frozen set; "
            f"unexpected={unexpected_panel or 'none'} "
            f"missing={missing_panel or 'none'}")
    sizes = load_genome_sizes(args.target_genome_sizes, panel)
    unexpected_size = sorted(set(sizes) - set(EXPECTED_LABELS))
    missing_size = sorted(set(EXPECTED_LABELS) - set(sizes))
    if unexpected_size or missing_size:
        raise AuditError(
            "target genome-size labels do not match the frozen set; "
            f"unexpected={unexpected_size or 'none'} "
            f"missing={missing_size or 'none'}")
    aliases = load_metaphlan_aliases(args.feature_aliases, panel)
    feature_index = build_feature_index(panel, aliases)

    fields, rows = read_tsv(args.comparison, "reference comparison")
    required_columns = set(OBSERVATION_KEY)
    for ratio_column, error_column, type_columns, _ in ARMS.values():
        required_columns |= {ratio_column, error_column, *type_columns}
    missing_columns = sorted(required_columns - set(fields))
    if missing_columns:
        raise AuditError(
            "reference comparison lacks column(s): " + ", ".join(missing_columns))

    seen_keys: set[tuple[str, ...]] = set()
    profilers_present: set[str] = set()
    bracken_rows = 0
    metaphlan_rows = 0
    # (scope, target_label) -> arm -> list of (ratio, absolute_relative_error)
    observations: dict[tuple[str, str], dict[str, list[tuple[float, float]]]] = \
        defaultdict(lambda: {arm: [] for arm in ARMS})
    excluded: dict[tuple[str, str], int] = defaultdict(int)
    exclusion_reasons: dict[tuple[str, str], dict[str, int]] = \
        defaultdict(lambda: defaultdict(int))
    observed_features: dict[str, set[str]] = defaultdict(set)
    # (arm, target_label) -> set of (row type, selected type, scale, estimand)
    provenance: dict[tuple[str, str], set[tuple[str, str, str, str]]] = \
        defaultdict(set)
    cohorts_present: set[str] = set()

    for line, row in enumerate(rows, start=2):
        # Defect 3: every component of the physical identity must be present and
        # nonblank BEFORE any duplicate check, so two rows can never collide
        # merely because a field was left empty. Blanks are never replaced by a
        # sentinel.
        for name in OBSERVATION_KEY:
            if not (row.get(name) or "").strip():
                raise AuditError(
                    f"reference comparison line {line}: blank {name}; every "
                    "physical observation identity field is required")
        key = tuple(row[name].strip() for name in OBSERVATION_KEY)
        if key in seen_keys:
            raise AuditError(
                f"reference comparison line {line}: duplicate physical "
                f"observation key {key}")
        seen_keys.add(key)
        profiler = row["profiler"].strip()
        profilers_present.add(profiler)
        # Defect 2: the cohort of every row is inspected, including Bracken's.
        # Non-Yachida rows are never silently filtered out.
        cohorts_present.add(row["cohort"].strip())
        if profiler != METAPHLAN:
            # Bracken is required as an input gate but is never regressed: its
            # two arms are numerically identical by construction.
            if profiler == BRACKEN:
                bracken_rows += 1
            continue
        metaphlan_rows += 1

        population = (row.get("analysis_population") or "").strip()
        if population not in POPULATION_SCOPES:
            raise AuditError(
                f"reference comparison line {line}: unsupported "
                f"analysis_population {population!r}")
        feature = (row.get("target_feature") or "").strip()
        if feature not in feature_index:
            raise AuditError(
                f"reference comparison line {line}: MetaPhlAn target_feature "
                f"{feature!r} does not match exactly one implanted target; the "
                "audit never falls back to fuzzy, prefix or substring matching")
        label = feature_index[feature]
        observed_features[label].add(feature)

        # Paired eligibility. An observation is used only when BOTH arms yield a
        # usable ratio, because the arms are compared as paired measurements of
        # the same physical observation. No pseudocount is ever added.
        # A blank cell in either arm is a missing paired value, not an
        # excludable observation: the two arms are then not the same physical
        # measurement and the comparison is structurally incomplete.
        for arm in ARM_ORDER:
            ratio_column, error_column, type_columns, _ = ARMS[arm]
            for column in (ratio_column, error_column, *type_columns):
                if not (row.get(column) or "").strip():
                    raise AuditError(
                        f"reference comparison line {line}: missing {column}; "
                        "the paired arms are incomplete")

        parsed: dict[str, float] = {}
        reasons: list[str] = []
        for arm in ARM_ORDER:
            ratio_column, _, type_columns, _ = ARMS[arm]
            value, problem = ratio_problem(row.get(ratio_column, ""))
            if problem is not None:
                reasons.append(f"{arm}_{problem}")
            else:
                parsed[arm] = value
            provenance[(arm, label)].add(
                arm_provenance(row, arm, type_columns, line))
        if reasons:
            for scope in (population, "pooled"):
                excluded[(scope, label)] += 1
                for reason in reasons:
                    exclusion_reasons[(scope, label)][reason] += 1
            continue

        # Defect 4: the supplied absolute relative error is validated, never
        # recomputed and silently substituted. Rows already paired-excluded
        # above never reach this check.
        errors: dict[str, float] = {}
        for arm in ARM_ORDER:
            _, error_column, _, _ = ARMS[arm]
            raw = row[error_column].strip()
            try:
                error = float(raw)
            except (TypeError, ValueError) as problem:
                raise AuditError(
                    f"reference comparison line {line}: nonnumeric "
                    f"{error_column} {raw!r}") from problem
            if not math.isfinite(error):
                raise AuditError(
                    f"reference comparison line {line}: nonfinite {error_column}")
            if error < 0:
                raise AuditError(
                    f"reference comparison line {line}: negative {error_column} "
                    f"{render(error)}; an absolute relative error cannot be "
                    "negative")
            expected_error = abs(parsed[arm] - 1.0)
            tolerance = max(ERROR_ABSOLUTE_TOLERANCE,
                            ERROR_RELATIVE_TOLERANCE
                            * max(abs(error), abs(expected_error), 1.0))
            if abs(error - expected_error) > tolerance:
                raise AuditError(
                    f"reference comparison line {line}: {error_column} "
                    f"{render(error)} is inconsistent with "
                    f"abs(observed_over_expected - 1) = {render(expected_error)} "
                    f"beyond the {render(tolerance)} round-trip tolerance")
            errors[arm] = error

        for scope in (population, "pooled"):
            for arm in ARM_ORDER:
                observations[(scope, label)][arm].append((parsed[arm], errors[arm]))

    missing_profilers = sorted(set(REQUIRED_PROFILERS) - profilers_present)
    if missing_profilers:
        raise AuditError(
            "reference comparison must contain both profilers; missing "
            + ", ".join(missing_profilers))
    unknown_profilers = sorted(profilers_present - set(REQUIRED_PROFILERS))
    if unknown_profilers:
        raise AuditError(
            "reference comparison has unknown profiler(s): "
            + ", ".join(unknown_profilers))
    if not metaphlan_rows:
        raise AuditError("reference comparison has no MetaPhlAn observations")

    if cohorts_present != {DEVELOPMENT_COHORT}:
        raise AuditError(
            "this development audit is restricted to the "
            f"{DEVELOPMENT_COHORT!r} cohort; the comparison carries "
            + ", ".join(repr(name) for name in sorted(cohorts_present))
            + ". Non-matching rows are never silently filtered out")

    observed_labels = set(observed_features)
    unexpected = sorted(observed_labels - set(EXPECTED_LABELS))
    missing = sorted(set(EXPECTED_LABELS) - observed_labels)
    if unexpected or missing:
        raise AuditError(
            "MetaPhlAn target labels do not match the frozen panel; "
            f"unexpected={unexpected or 'none'} missing={missing or 'none'}")
    for label in sorted(observed_features):
        if len(observed_features[label]) != 1:
            raise AuditError(
                f"target label {label} maps to multiple MetaPhlAn features: "
                + ", ".join(sorted(observed_features[label])))

    # ---------------------------------------------------------- aggregation
    target_rows: list[dict[str, object]] = []
    medians: dict[tuple[str, str, str], float] = {}
    error_medians: dict[tuple[str, str, str], float] = {}
    for scope in SCOPES:
        for label in EXPECTED_LABELS:
            eligible = observations.get((scope, label))
            usable = eligible[ARM_ORDER[0]] if eligible else []
            if not usable:
                raise AuditError(
                    f"target {label} has no eligible MetaPhlAn observation in "
                    f"scope {scope}; the target is entirely unusable")
            feature = sorted(observed_features[label])[0]
            size = sizes[label]
            for arm in ARM_ORDER:
                values = eligible[arm]
                ratios = [ratio for ratio, _ in values]
                errors = [error for _, error in values]
                median_ratio = median(ratios)
                median_log2 = median([math.log2(ratio) for ratio in ratios])
                median_error = median(errors)
                medians[(scope, arm, label)] = median_log2
                error_medians[(scope, arm, label)] = median_error
                records = sorted(provenance[(arm, label)])
                if len(records) != 1:
                    raise AuditError(
                        f"{arm} arm has inconsistent selected-reference "
                        f"provenance for {label}: "
                        + "; ".join("/".join(item) for item in records))
                row_type, selected_type, scale, estimand = records[0]
                reasons = exclusion_reasons.get((scope, label), {})
                target_rows.append({
                    "scope": scope, "arm": arm, "target_label": label,
                    "target_feature": feature,
                    "genome_size_bp": render(size),
                    "log2_genome_size_bp": render(math.log2(size)),
                    "eligible_observations": len(values),
                    "excluded_observations": excluded.get((scope, label), 0),
                    "exclusion_reasons": ";".join(
                        f"{name}={count}" for name, count in sorted(reasons.items())),
                    "median_observed_over_expected": render(median_ratio),
                    "median_log2_observed_over_expected": render(median_log2),
                    "median_absolute_relative_error": render(median_error),
                    # The SELECTED estimand, not the source row's provenance.
                    "reference_type": estimand,
                    "row_reference_type": row_type,
                    "selected_reference_type": selected_type,
                    "reference_scale": scale,
                })

    # ----------------------------------------------------------- regression
    regression_rows: list[dict[str, object]] = []
    for scope in SCOPES:
        for arm in ARM_ORDER:
            expected_slope = ARMS[arm][3]
            points = [(math.log2(sizes[label]), medians[(scope, arm, label)])
                      for label in EXPECTED_LABELS]
            if len(points) < EXPECTED_TARGET_COUNT:
                raise AuditError(
                    f"scope {scope} arm {arm} has only {len(points)} target "
                    f"points; {EXPECTED_TARGET_COUNT} are required")
            estimate = fit(points)
            interval = bootstrap_slope_ci(points, args.bootstrap_replicates,
                                          BOOTSTRAP_SEED)
            regression_rows.append({
                "scope": scope, "arm": arm, "n_targets": len(points),
                "slope": render(estimate["slope"]),
                "intercept": render(estimate["intercept"]),
                "r_squared": render(estimate["r_squared"]),
                "pearson_correlation": render(estimate["pearson_correlation"]),
                "expected_slope": render(expected_slope),
                "absolute_distance_from_expected":
                    render(abs(estimate["slope"] - expected_slope)),
                "bootstrap_replicates_requested":
                    interval["bootstrap_replicates_requested"],
                "bootstrap_replicates_valid": interval["bootstrap_replicates_valid"],
                "bootstrap_replicates_discarded":
                    interval["bootstrap_replicates_discarded"],
                "bootstrap_ci_lower": render(interval["bootstrap_ci_lower"]),
                "bootstrap_ci_upper": render(interval["bootstrap_ci_upper"]),
                "bootstrap_seed": interval["bootstrap_seed"],
            })

    # -------------------------------------------------------- paired change
    paired_rows: list[dict[str, object]] = []
    for scope in SCOPES:
        for label in EXPECTED_LABELS:
            sensitivity = medians[(scope, "read_proportional_sensitivity", label)]
            primary = medians[(scope, "genome_equivalent_primary", label)]
            sensitivity_error = error_medians[
                (scope, "read_proportional_sensitivity", label)]
            primary_error = error_medians[(scope, "genome_equivalent_primary", label)]
            paired_rows.append({
                "scope": scope, "target_label": label,
                "target_feature": sorted(observed_features[label])[0],
                "genome_size_bp": render(sizes[label]),
                "sensitivity_median_log2_ratio": render(sensitivity),
                "primary_median_log2_ratio": render(primary),
                "primary_minus_sensitivity_log2_ratio": render(primary - sensitivity),
                "sensitivity_median_absolute_relative_error": render(sensitivity_error),
                "primary_median_absolute_relative_error": render(primary_error),
                "primary_minus_sensitivity_absolute_relative_error":
                    render(primary_error - sensitivity_error),
            })

    # -------------------------------------------------------------- outputs
    args.outdir.mkdir(parents=True, exist_ok=True)
    target_path = args.outdir / "metaphlan_genome_size_target_level.tsv"
    regression_path = args.outdir / "metaphlan_genome_size_regression.tsv"
    paired_path = args.outdir / "metaphlan_genome_size_paired_change.tsv"
    validation_path = args.outdir / "metaphlan_genome_size_audit_validation.tsv"
    write_tsv(target_path, TARGET_LEVEL_FIELDS, target_rows)
    write_tsv(regression_path, REGRESSION_FIELDS, regression_rows)
    write_tsv(paired_path, PAIRED_CHANGE_FIELDS, paired_rows)

    total_excluded = sum(excluded[(scope, label)] for scope in POPULATION_SCOPES
                         for label in EXPECTED_LABELS)
    reason_totals: dict[str, int] = defaultdict(int)
    for scope in POPULATION_SCOPES:
        for label in EXPECTED_LABELS:
            for reason, count in exclusion_reasons.get((scope, label), {}).items():
                reason_totals[reason] += count
    validation_rows = [
        {"metric": "comparison_path", "value": str(args.comparison.resolve())},
        {"metric": "comparison_sha256", "value": sha256(args.comparison)},
        {"metric": "comparison_validation_path",
         "value": str(args.comparison_validation.resolve())},
        {"metric": "comparison_validation_sha256",
         "value": sha256(args.comparison_validation)},
        {"metric": "target_genome_sizes_path",
         "value": str(args.target_genome_sizes.resolve())},
        {"metric": "target_genome_sizes_sha256",
         "value": sha256(args.target_genome_sizes)},
        {"metric": "spike_panel_path", "value": str(args.spike_panel.resolve())},
        {"metric": "spike_panel_sha256", "value": sha256(args.spike_panel)},
        {"metric": "feature_aliases_path",
         "value": str(args.feature_aliases.resolve())},
        {"metric": "feature_aliases_sha256", "value": sha256(args.feature_aliases)},
        {"metric": "target_count", "value": len(observed_labels)},
        {"metric": "expected_target_count", "value": EXPECTED_TARGET_COUNT},
        {"metric": "metaphlan_observations", "value": metaphlan_rows},
        {"metric": "bracken_observations", "value": bracken_rows},
        {"metric": "bracken_rows_in_regression", "value": 0},
        {"metric": "bracken_identical", "value": bracken_gate},
        {"metric": "observed_cohorts",
         "value": ",".join(sorted(cohorts_present))},
        {"metric": "validated_cohort_count", "value": len(cohorts_present)},
        {"metric": "required_cohort", "value": DEVELOPMENT_COHORT},
        {"metric": "observation_key_fields", "value": ",".join(OBSERVATION_KEY)},
        {"metric": "observation_key_blank_policy",
         "value": "every component required and nonblank; no sentinel"},
        {"metric": "selected_reference_sensitivity",
         "value": "/".join(ARM_SELECTION["read_proportional_sensitivity"])
                  + " -> " + ARM_SELECTED_ESTIMAND["read_proportional_sensitivity"]},
        {"metric": "selected_reference_primary",
         "value": "/".join(ARM_SELECTION["genome_equivalent_primary"])
                  + " -> " + ARM_SELECTED_ESTIMAND["genome_equivalent_primary"]},
        {"metric": "row_reference_type_metaphlan",
         "value": METAPHLAN_ROW_REFERENCE},
        {"metric": "absolute_relative_error_policy",
         "value": "validated against abs(observed_over_expected - 1); never "
                  "recomputed and substituted"},
        {"metric": "absolute_relative_error_tolerance",
         "value": ERROR_TOLERANCE_POLICY},
        {"metric": "paired_excluded_observations", "value": total_excluded},
        {"metric": "exclusion_reasons",
         "value": ";".join(f"{name}={count}"
                           for name, count in sorted(reason_totals.items())) or "none"},
        {"metric": "pseudocount_applied", "value": "NONE"},
        {"metric": "eligibility_rule",
         "value": "paired: an observation is used only when both arms give a "
                  "positive finite ratio"},
        {"metric": "statistical_unit", "value": "implanted_target"},
        {"metric": "bootstrap_unit", "value": "implanted_target"},
        {"metric": "bootstrap_seed", "value": BOOTSTRAP_SEED},
        {"metric": "bootstrap_replicates_requested",
         "value": args.bootstrap_replicates},
        {"metric": "bootstrap_interpretation",
         "value": "descriptive uncertainty over the ten implanted taxa only"},
        {"metric": "expected_slope_read_proportional_sensitivity", "value": "-1"},
        {"metric": "expected_slope_genome_equivalent_primary", "value": "0"},
        {"metric": "genome_size_source", "value": "measured spike FASTAs"},
        {"metric": "fitted_genome_size_constant_used", "value": "NONE"},
        {"metric": "cohort_scope", "value": "yachida"},
        {"metric": "audit_status", "value": "PASS"},
        {"metric": "development_status", "value": "DEVELOPMENT_ONLY"},
    ]
    write_tsv(validation_path, ["metric", "value"], validation_rows)

    development = args.outdir / "DEVELOPMENT_ONLY.txt"
    development.write_text(
        "status=DEVELOPMENT_ONLY\n"
        "cohort=yachida\n"
        "profiler=metaphlan4\n"
        "use_for_manuscript=NO\n"
        "reason=Yachida-only genome-size residual audit; requires replication in "
        "Feng and Zeller and a final three-cohort execution before any "
        "manuscript claim\n"
        "requires=feng_replication,zeller_replication,three_cohort_validation\n",
        encoding="utf-8")

    outputs = [target_path, regression_path, paired_path, validation_path,
               development]
    checksum_path = args.outdir / "metaphlan_genome_size_audit.sha256"
    with checksum_path.open("w", encoding="utf-8") as handle:
        for path in list(inputs.values()) + outputs:
            handle.write(f"{sha256(path)}  {path.resolve()}\n")

    (args.outdir / "SUCCESS").write_text(
        f"status\tPASS\ntargets\t{len(observed_labels)}\n"
        f"metaphlan_observations\t{metaphlan_rows}\n"
        f"bootstrap_seed\t{BOOTSTRAP_SEED}\n", encoding="utf-8")
    print(f"[PASS] MetaPhlAn genome-size residual audit: {args.outdir}")
    for row in regression_rows:
        print(f"[INFO] {row['scope']:<12} {row['arm']:<30} "
              f"slope={row['slope']} expected={row['expected_slope']} "
              f"ci=[{row['bootstrap_ci_lower']}, {row['bootstrap_ci_upper']}]")


def main() -> None:
    args = parse_args()
    try:
        run(args)
    except (AuditError, OSError, KeyError) as error:
        raise SystemExit(f"[ERROR] {error}") from error


if __name__ == "__main__":
    main()
