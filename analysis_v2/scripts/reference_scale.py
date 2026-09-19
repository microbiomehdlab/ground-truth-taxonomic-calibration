#!/usr/bin/env python3
"""Profiler-aware selection and validation of the quantitative reference scale.

`build_perturbation_response_input.py` emits both references side by side:

* read-proportional -- `expected_abundance_fraction` and friends, unchanged;
* profiler scale    -- `*_profiler_scale`, read-proportional for Bracken and
                       genome-equivalent for MetaPhlAn.

A **primary** table therefore normally and correctly holds two row-level
reference types at once:

    kraken2_bracken -> read_proportional
    metaphlan4      -> genome_equivalent

That is not a mixed-reference error. Validation is performed **within profiler**,
never globally across profilers.
"""

from __future__ import annotations

CHOICES = ("profiler_scale", "read_proportional")

# Primary profiler -> required row-level reference_type under profiler_scale.
PRIMARY_REFERENCE_BY_PROFILER = {
    "kraken2_bracken": "read_proportional",
    "metaphlan4": "genome_equivalent",
}

# Canonical name -> profiler-scale source column.
PROFILER_SCALE_ALIASES = {
    "expected_abundance_fraction": "expected_abundance_profiler_scale",
    "dilution_retained_baseline": "retained_baseline_profiler_scale",
    "response_signal": "response_signal_profiler_scale",
    "response_delta": "response_delta_profiler_scale",
    "quantitative_log2_ratio": "quantitative_log2_ratio_profiler_scale",
    # The implanted signal must move with the response signal, otherwise a
    # recovery ratio would divide a profiler-scale numerator by a read fraction.
    "target_fraction_for_feature": "implanted_signal_profiler_scale",
    # Bounded errors are derived from the expected value, so they must move with
    # it; leaving them read-proportional made the analyzer reject genuine
    # profiler-scale builder output.
    "signed_bounded_error": "signed_bounded_error_profiler_scale",
    "absolute_bounded_error": "absolute_bounded_error_profiler_scale",
}


def add_reference_scale_argument(parser) -> None:
    parser.add_argument(
        "--reference-scale", choices=list(CHOICES), default=None, required=True,
        help="profiler_scale (primary: read-proportional for Bracken, "
             "genome-equivalent for MetaPhlAn) or read_proportional "
             "(labelled sensitivity). No default.")


def select_expression(columns: set[str], scale: str) -> str:
    """SELECT list aliasing the chosen reference onto the canonical names."""
    if scale not in CHOICES:
        raise SystemExit(f"[FAIL] unknown --reference-scale {scale!r}")
    if scale == "read_proportional":
        return "*"
    missing = sorted(set(PROFILER_SCALE_ALIASES.values()) - columns)
    if missing:
        raise SystemExit(
            "[FAIL] --reference-scale profiler_scale needs column(s) "
            + ", ".join(missing)
            + "; rebuild the response input with the migrated "
              "build_perturbation_response_input.py")
    replaced = set(PROFILER_SCALE_ALIASES)
    kept = [c for c in sorted(columns) if c not in replaced]
    aliased = [f"{src} AS {dst}" for dst, src in sorted(PROFILER_SCALE_ALIASES.items())]
    return ", ".join(kept + aliased)


def validate_profiler_reference_pairs(pairs, scale: str) -> str:
    """Validate (profiler, reference_type) pairs within profiler.

    `pairs` is an iterable of 2-tuples. Returns the selected reference label to
    stamp on analysis outputs. Row-level provenance is never rewritten.
    """
    if scale not in CHOICES:
        raise SystemExit(f"[FAIL] unknown --reference-scale {scale!r}")

    observed: dict[str, set[str]] = {}
    for profiler, reference in pairs:
        profiler = (profiler or "").strip()
        reference = (reference or "").strip()
        if not profiler:
            raise SystemExit("[FAIL] response table has a blank profiler")
        if not reference:
            raise SystemExit(
                f"[FAIL] response table has a blank reference_type for {profiler}")
        observed.setdefault(profiler, set()).add(reference)

    if not observed:
        raise SystemExit(
            "[FAIL] response table has no profiler/reference_type values; rebuild "
            "it with the migrated build_perturbation_response_input.py")

    if scale == "read_proportional":
        # The preserved read-reference fields are always available, including in
        # a primary table whose MetaPhlAn rows are labelled genome_equivalent.
        # Row-level provenance stands; only the selected scale is reported.
        return "read_proportional"

    for profiler, references in sorted(observed.items()):
        if profiler not in PRIMARY_REFERENCE_BY_PROFILER:
            raise SystemExit(
                f"[FAIL] unknown profiler {profiler!r}; expected one of "
                + ", ".join(sorted(PRIMARY_REFERENCE_BY_PROFILER)))
        if len(references) > 1:
            raise SystemExit(
                f"[FAIL] profiler {profiler} carries multiple reference types: "
                + ", ".join(sorted(references)))
        expected = PRIMARY_REFERENCE_BY_PROFILER[profiler]
        actual = next(iter(references))
        if actual != expected:
            raise SystemExit(
                f"[FAIL] profiler {profiler} has reference_type {actual!r}; the "
                f"primary scale requires {expected!r}")
    return "profiler_scale_primary"


def require_valid_reference(con, view: str, scale: str) -> str:
    """DuckDB wrapper around `validate_profiler_reference_pairs`."""
    columns = {row[0] for row in con.execute(f"DESCRIBE {view}").fetchall()}
    for required in ("profiler", "reference_type"):
        if required not in columns:
            raise SystemExit(
                f"[FAIL] response table lacks a {required} column; rebuild it "
                "with the migrated build_perturbation_response_input.py")
    pairs = con.execute(
        f"SELECT DISTINCT profiler, reference_type FROM {view} "
        "ORDER BY profiler, reference_type").fetchall()
    return validate_profiler_reference_pairs(pairs, scale)
