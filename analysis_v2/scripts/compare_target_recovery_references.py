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
OBSERVATION_FIELDS = list(KEY) + [
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
        required = set(KEY) | set(NUMERIC) | {"recovery_class", "reference_type"}
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

            p_class = left["recovery_class"]
            s_class = right["recovery_class"]
            p_error = number(left, "absolute_relative_error", key)
            s_error = number(right, "absolute_relative_error", key)
            row: dict[str, object] = dict(zip(KEY, key))
            row.update({
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
