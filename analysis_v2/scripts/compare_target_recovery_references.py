#!/usr/bin/env python3
"""Compare profiler-scale target recovery with read-reference sensitivity.

The comparison is paired at the physical observation level. Bracken must be
numerically identical between arms; MetaPhlAn differences are summarized as
changes caused by using the genome-equivalent rather than read-proportional
implanted signal.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from collections import defaultdict
from pathlib import Path


KEY = (
    "cohort", "study", "sample_id", "condition", "analysis_population",
    "profiler", "baseline_id", "target_feature", "dose_index",
    "nominal_target_fraction",
)
NUMERIC = (
    "target_fraction_for_feature", "response_signal_selected",
    "implanted_signal_selected", "observed_over_expected",
    "absolute_relative_error",
)
CLASS_RANK = {"Good": 0, "Average": 1, "Poor / missed": 2}

# Provenance carried by summarize_target_recovery.py on every row:
#   reference_type          -- ORIGINAL row-level scale of the source profile
#   reference_scale         -- the --reference-scale the arm was summarised with
#   selected_reference_type -- the label stamped for that selection
# These are three different things. The original row type of a MetaPhlAn row
# stays 'genome_equivalent' even in the read-proportional arm, because it
# records where the row came from, not which estimand was selected. Carrying
# only `reference_type` therefore mislabels the sensitivity arm.
PROVENANCE = ("reference_type", "reference_scale", "selected_reference_type")

# (reference_scale, selected_reference_type) -> selected estimand, where None
# means "the row's own reference type". Under profiler_scale the selected
# estimand IS the row's native reference: read-proportional for Bracken,
# genome-equivalent for MetaPhlAn. Under read_proportional it is
# read-proportional for both profilers.
SELECTION_CONTRACT = {
    ("profiler_scale", "profiler_scale_primary"): None,
    ("read_proportional", "read_proportional"): "read_proportional",
}
# arm -> the (reference_scale, selected_reference_type) pair it must carry.
ARM_SELECTION = {
    "primary": ("profiler_scale", "profiler_scale_primary"),
    "sensitivity": ("read_proportional", "read_proportional"),
}
# profiler -> required ORIGINAL row reference type, and the selected estimand
# each arm must resolve to.
PROFILER_ROW_REFERENCE = {
    "kraken2_bracken": "read_proportional",
    "metaphlan4": "genome_equivalent",
}
ARM_SELECTED_ESTIMAND = {
    ("kraken2_bracken", "primary"): "read_proportional",
    ("kraken2_bracken", "sensitivity"): "read_proportional",
    ("metaphlan4", "primary"): "genome_equivalent",
    ("metaphlan4", "sensitivity"): "read_proportional",
}
OBSERVATION_FIELDS = list(KEY) + [
    # Explicit, unambiguous provenance for each arm.
    "primary_row_reference_type", "sensitivity_row_reference_type",
    "primary_selected_reference_type", "sensitivity_selected_reference_type",
    "primary_reference_scale", "sensitivity_reference_scale",
    "primary_selected_estimand", "sensitivity_selected_estimand",
    # Backward-compatible aliases of the ORIGINAL row type. Ambiguous by name;
    # new consumers must use the explicit fields above.
    "primary_reference_type", "sensitivity_reference_type",
    "primary_implanted_signal", "sensitivity_implanted_signal",
    "primary_response_signal", "sensitivity_response_signal",
    "primary_observed_over_expected", "sensitivity_observed_over_expected",
    "primary_absolute_relative_error", "sensitivity_absolute_relative_error",
    "absolute_relative_error_change_primary_minus_sensitivity",
    "primary_recovery_class", "sensitivity_recovery_class",
    "class_transition", "primary_class_improved", "primary_class_worsened",
]
SUMMARY_FIELDS = [
    "profiler", "cohort", "analysis_population", "condition",
    "target_feature", "nominal_target_fraction", "observations",
    "primary_good_fraction", "sensitivity_good_fraction",
    "primary_average_fraction", "sensitivity_average_fraction",
    "primary_poor_fraction", "sensitivity_poor_fraction",
    "class_changed_fraction", "primary_improved_fraction",
    "primary_worsened_fraction", "median_primary_observed_over_expected",
    "median_sensitivity_observed_over_expected",
    "median_primary_absolute_relative_error",
    "median_sensitivity_absolute_relative_error",
    "median_absolute_relative_error_change_primary_minus_sensitivity",
]
OVERALL_FIELDS = [
    "profiler", "analysis_population", "observations",
    "primary_good_fraction", "sensitivity_good_fraction",
    "primary_average_fraction", "sensitivity_average_fraction",
    "primary_poor_fraction", "sensitivity_poor_fraction",
    "class_changed_fraction", "primary_improved_fraction",
    "primary_worsened_fraction", "median_primary_observed_over_expected",
    "median_sensitivity_observed_over_expected",
    "median_primary_absolute_relative_error",
    "median_sensitivity_absolute_relative_error",
    "median_absolute_relative_error_change_primary_minus_sensitivity",
]


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_rows(path: Path) -> tuple[list[str], dict[tuple[str, ...], dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        required = set(KEY) | set(NUMERIC) | {"recovery_class"} | set(PROVENANCE)
        missing = sorted(required - set(fields))
        if missing:
            raise ValueError(f"{path} lacks columns: {', '.join(missing)}")
        result: dict[tuple[str, ...], dict[str, str]] = {}
        for line, row in enumerate(reader, start=2):
            key = tuple(row[field] for field in KEY)
            if key in result:
                raise ValueError(f"{path}:{line}: duplicate observation key {key}")
            if row["recovery_class"] not in CLASS_RANK:
                raise ValueError(f"{path}:{line}: unknown recovery class")
            result[key] = row
    if not result:
        raise ValueError(f"{path} contains no observations")
    return fields, result


def selected_estimand(row: dict[str, str], arm: str, profiler: str,
                      key: tuple[str, ...]) -> str:
    """Resolve and validate one arm's selected quantitative reference.

    Fails closed rather than letting an arm inherit a label that describes the
    source row instead of the estimand the arm actually measures.
    """
    values = {name: (row.get(name) or "").strip() for name in PROVENANCE}
    for name, value in values.items():
        if not value:
            raise ValueError(f"{arm} arm has a blank {name} for {key}")
    pair = (values["reference_scale"], values["selected_reference_type"])
    if pair not in SELECTION_CONTRACT:
        raise ValueError(
            f"{arm} arm has an unknown reference selection "
            f"{pair} for {key}")
    if pair != ARM_SELECTION[arm]:
        raise ValueError(
            f"{arm} arm must be summarised as {ARM_SELECTION[arm]}, not {pair}, "
            f"for {key}")
    expected_row = PROFILER_ROW_REFERENCE.get(profiler)
    if expected_row is None:
        raise ValueError(f"unknown profiler {profiler!r} for {key}")
    if values["reference_type"] != expected_row:
        raise ValueError(
            f"{arm} arm has row reference_type {values['reference_type']!r} for "
            f"{profiler}; the source rows must be {expected_row!r} for {key}")
    resolved = SELECTION_CONTRACT[pair] or values["reference_type"]
    required = ARM_SELECTED_ESTIMAND[(profiler, arm)]
    if resolved != required:
        raise ValueError(
            f"{arm} arm resolves to the {resolved!r} estimand for {profiler}; "
            f"{required!r} is required for {key}")
    return resolved


def number(row: dict[str, str], field: str, key: tuple[str, ...]) -> float:
    try:
        value = float(row[field])
    except (TypeError, ValueError) as error:
        raise ValueError(f"nonnumeric {field} for {key}") from error
    if not math.isfinite(value):
        raise ValueError(f"nonfinite {field} for {key}")
    return value


def render(value: float) -> str:
    return format(value, ".17g")


def median(values: list[float]) -> float:
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def aggregate(rows: list[dict[str, object]]) -> dict[str, object]:
    n = len(rows)
    fraction = lambda count: render(count / n)
    return {
        "observations": n,
        "primary_good_fraction": fraction(sum(r["primary_recovery_class"] == "Good" for r in rows)),
        "sensitivity_good_fraction": fraction(sum(r["sensitivity_recovery_class"] == "Good" for r in rows)),
        "primary_average_fraction": fraction(sum(r["primary_recovery_class"] == "Average" for r in rows)),
        "sensitivity_average_fraction": fraction(sum(r["sensitivity_recovery_class"] == "Average" for r in rows)),
        "primary_poor_fraction": fraction(sum(r["primary_recovery_class"] == "Poor / missed" for r in rows)),
        "sensitivity_poor_fraction": fraction(sum(r["sensitivity_recovery_class"] == "Poor / missed" for r in rows)),
        "class_changed_fraction": fraction(sum(r["primary_recovery_class"] != r["sensitivity_recovery_class"] for r in rows)),
        "primary_improved_fraction": fraction(sum(int(r["primary_class_improved"]) for r in rows)),
        "primary_worsened_fraction": fraction(sum(int(r["primary_class_worsened"]) for r in rows)),
        "median_primary_observed_over_expected": render(median([float(r["primary_observed_over_expected"]) for r in rows])),
        "median_sensitivity_observed_over_expected": render(median([float(r["sensitivity_observed_over_expected"]) for r in rows])),
        "median_primary_absolute_relative_error": render(median([float(r["primary_absolute_relative_error"]) for r in rows])),
        "median_sensitivity_absolute_relative_error": render(median([float(r["sensitivity_absolute_relative_error"]) for r in rows])),
        "median_absolute_relative_error_change_primary_minus_sensitivity": render(median([float(r["absolute_relative_error_change_primary_minus_sensitivity"]) for r in rows])),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--sensitivity", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    for path in (args.primary, args.sensitivity):
        if not path.is_file() or not path.stat().st_size:
            raise SystemExit(f"[ERROR] missing or empty input: {path}")
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise SystemExit("[ERROR] OUTDIR must be new or empty")
    args.outdir.mkdir(parents=True, exist_ok=True)

    try:
        _, primary = read_rows(args.primary)
        _, sensitivity = read_rows(args.sensitivity)
        if primary.keys() != sensitivity.keys():
            only_primary = len(primary.keys() - sensitivity.keys())
            only_sensitivity = len(sensitivity.keys() - primary.keys())
            raise ValueError("observation keys differ between arms: "
                             f"primary_only={only_primary}, "
                             f"sensitivity_only={only_sensitivity}")

        observation_rows: list[dict[str, object]] = []
        grouped: dict[tuple[str, ...], list[dict[str, object]]] = defaultdict(list)
        bracken_rows = 0
        metaphlan_rows = 0
        for key in sorted(primary):
            left, right = primary[key], sensitivity[key]
            profiler = left["profiler"]
            if right["profiler"] != profiler:
                raise ValueError(f"profiler mismatch for {key}")
            if profiler == "kraken2_bracken":
                for field in NUMERIC:
                    if number(left, field, key) != number(right, field, key):
                        raise ValueError(f"Bracken differs in {field} for {key}")
                if left["recovery_class"] != right["recovery_class"]:
                    raise ValueError(f"Bracken recovery class differs for {key}")
                bracken_rows += 1
            elif profiler == "metaphlan4":
                metaphlan_rows += 1
            else:
                raise ValueError(f"unknown profiler {profiler!r} for {key}")

            primary_estimand = selected_estimand(left, "primary", profiler, key)
            sensitivity_estimand = selected_estimand(
                right, "sensitivity", profiler, key)

            p_class = left["recovery_class"]
            s_class = right["recovery_class"]
            p_error = number(left, "absolute_relative_error", key)
            s_error = number(right, "absolute_relative_error", key)
            row: dict[str, object] = dict(zip(KEY, key))
            row.update({
                "primary_row_reference_type": left["reference_type"],
                "sensitivity_row_reference_type": right["reference_type"],
                "primary_selected_reference_type": left["selected_reference_type"],
                "sensitivity_selected_reference_type":
                    right["selected_reference_type"],
                "primary_reference_scale": left["reference_scale"],
                "sensitivity_reference_scale": right["reference_scale"],
                "primary_selected_estimand": primary_estimand,
                "sensitivity_selected_estimand": sensitivity_estimand,
                "primary_reference_type": left["reference_type"],
                "sensitivity_reference_type": right["reference_type"],
                "primary_implanted_signal": left["implanted_signal_selected"],
                "sensitivity_implanted_signal": right["implanted_signal_selected"],
                "primary_response_signal": left["response_signal_selected"],
                "sensitivity_response_signal": right["response_signal_selected"],
                "primary_observed_over_expected": left["observed_over_expected"],
                "sensitivity_observed_over_expected": right["observed_over_expected"],
                "primary_absolute_relative_error": left["absolute_relative_error"],
                "sensitivity_absolute_relative_error": right["absolute_relative_error"],
                "absolute_relative_error_change_primary_minus_sensitivity": render(p_error-s_error),
                "primary_recovery_class": p_class,
                "sensitivity_recovery_class": s_class,
                "class_transition": f"{s_class} -> {p_class}",
                "primary_class_improved": int(CLASS_RANK[p_class] < CLASS_RANK[s_class]),
                "primary_class_worsened": int(CLASS_RANK[p_class] > CLASS_RANK[s_class]),
            })
            observation_rows.append(row)
            group_key = (profiler, left["cohort"], left["analysis_population"],
                         left["condition"], left["target_feature"],
                         left["nominal_target_fraction"])
            grouped[group_key].append(row)

        if not bracken_rows or not metaphlan_rows:
            raise ValueError("both Bracken and MetaPhlAn observations are required")

        summary_rows: list[dict[str, object]] = []
        for group_key in sorted(grouped):
            rows = grouped[group_key]
            summary = dict(zip(SUMMARY_FIELDS[:6], group_key))
            summary.update(aggregate(rows))
            summary_rows.append(summary)

        overall_groups: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
        for row in observation_rows:
            overall_groups[(str(row["profiler"]),
                            str(row["analysis_population"]))].append(row)
        overall_rows: list[dict[str, object]] = []
        for group_key in sorted(overall_groups):
            row = dict(zip(OVERALL_FIELDS[:2], group_key))
            row.update(aggregate(overall_groups[group_key]))
            overall_rows.append(row)

        observations_path = args.outdir / "target_recovery_reference_comparison.tsv"
        summary_path = args.outdir / "target_recovery_reference_summary.tsv"
        overall_path = args.outdir / "target_recovery_reference_overall.tsv"
        write_tsv(observations_path, OBSERVATION_FIELDS, observation_rows)
        write_tsv(summary_path, SUMMARY_FIELDS, summary_rows)
        write_tsv(overall_path, OVERALL_FIELDS, overall_rows)
        validation_path = args.outdir / "reference_comparison_validation.tsv"
        write_tsv(validation_path, ["metric", "value"], [
            {"metric": "paired_observations", "value": len(observation_rows)},
            {"metric": "bracken_observations", "value": bracken_rows},
            {"metric": "metaphlan_observations", "value": metaphlan_rows},
            {"metric": "bracken_identical", "value": "PASS"},
            {"metric": "primary_selection",
             "value": "/".join(ARM_SELECTION["primary"])},
            {"metric": "sensitivity_selection",
             "value": "/".join(ARM_SELECTION["sensitivity"])},
            {"metric": "selected_estimand_contract",
             "value": ";".join(f"{profiler}/{arm}={estimand}" for
                               (profiler, arm), estimand
                               in sorted(ARM_SELECTED_ESTIMAND.items()))},
            {"metric": "comparison_status", "value": "DEVELOPMENT_ONLY"},
        ])
        checksum_path = args.outdir / "reference_comparison.sha256"
        with checksum_path.open("w", encoding="utf-8") as handle:
            for path in (args.primary, args.sensitivity, observations_path,
                         summary_path, overall_path, validation_path):
                handle.write(f"{digest(path)}  {path.resolve()}\n")
        (args.outdir / "SUCCESS").write_text(
            f"status\tPASS\npaired_observations\t{len(observation_rows)}\n",
            encoding="utf-8",
        )
        print(f"[PASS] Reference-scale comparison: {len(observation_rows)} paired observations")
        print(f"[PASS] Bracken identical: {bracken_rows} observations")
    except (OSError, ValueError, KeyError) as error:
        raise SystemExit(f"[ERROR] {error}") from error


if __name__ == "__main__":
    main()
