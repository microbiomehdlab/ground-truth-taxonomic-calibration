#!/usr/bin/env python3
"""Summarize direct-target recovery and paired biomarker-discovery thresholds.

This intentionally scans only direct implanted-target rows from the large
response Parquet.  Recovery is evaluated against the exact target fraction:
good <=10% relative error, average <=50%, and poor/missed otherwise.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

try:
    import duckdb
except ImportError as error:
    raise SystemExit("[ERROR] duckdb is required") from error

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reference_scale import (  # noqa: E402
    add_reference_scale_argument, require_valid_reference, select_expression)


def sql_path(path: Path) -> str:
    return str(path.resolve()).replace("'", "''")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--responses", type=Path, required=True)
    parser.add_argument("--biomarker-metrics", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    add_reference_scale_argument(parser)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    if not args.responses.is_file() or not args.responses.stat().st_size:
        raise SystemExit(f"[ERROR] Missing responses: {args.responses}")
    if args.biomarker_metrics and (
        not args.biomarker_metrics.is_file() or not args.biomarker_metrics.stat().st_size
    ):
        raise SystemExit(f"[ERROR] Missing biomarker metrics: {args.biomarker_metrics}")
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise SystemExit("[ERROR] OUTDIR must be new or empty")
    args.outdir.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect(database=":memory:")
    con.execute(f"PRAGMA threads={max(1, args.threads)}")
    con.execute(
        f"CREATE VIEW responses_raw AS SELECT * FROM read_parquet('{sql_path(args.responses)}')"
    )
    available = {row[0] for row in con.execute("DESCRIBE responses_raw").fetchall()}
    reference_label = require_valid_reference(con, "responses_raw",
                                              args.reference_scale)
    # Alias the chosen reference onto the canonical column names so the
    # recovery SQL below is unchanged and Bracken results stay identical.
    con.execute("CREATE VIEW responses AS SELECT "
                + select_expression(available, args.reference_scale)
                + " FROM responses_raw")
    required = {
        "cohort", "study", "sample_id", "condition", "analysis_population",
        "profiler", "baseline_id", "dose_index", "nominal_target_fraction",
        "target_fraction_for_feature", "feature", "baseline_abundance_fraction",
        "observed_abundance_fraction", "response_signal", "observed_detected",
        "is_direct_target",
    }
    columns = {row[0] for row in con.execute("DESCRIBE responses").fetchall()}
    missing = sorted(required - columns)
    if missing:
        raise SystemExit("[ERROR] Response table lacks: " + ", ".join(missing))

    # `responses` already aliases the selected scale onto the canonical names,
    # so response_signal and target_fraction_for_feature are both on the
    # selected scale. Naming them explicitly keeps that impossible to misread.
    con.execute("""
        CREATE VIEW target_rows AS
        SELECT *,
          response_signal AS response_signal_selected,
          target_fraction_for_feature AS implanted_signal_selected,
          CASE
            WHEN observed_detected=0 THEN 'Poor / missed'
            WHEN abs(response_signal-target_fraction_for_feature)
                   / target_fraction_for_feature <= 0.10 THEN 'Good'
            WHEN abs(response_signal-target_fraction_for_feature)
                   / target_fraction_for_feature <= 0.50 THEN 'Average'
            ELSE 'Poor / missed'
          END AS recovery_class,
          response_signal/target_fraction_for_feature AS observed_over_expected,
          abs(response_signal-target_fraction_for_feature)
            / target_fraction_for_feature AS absolute_relative_error
        FROM responses
        WHERE is_direct_target=1 AND target_fraction_for_feature>0
    """)
    # A direct implanted target must carry a usable implanted signal on the
    # selected scale; a missing or non-positive one would silently drop rows.
    unusable = con.execute("""
        SELECT count(*) FROM responses
        WHERE is_direct_target=1
          AND (target_fraction_for_feature IS NULL
               OR NOT isfinite(target_fraction_for_feature)
               OR target_fraction_for_feature <= 0)
    """).fetchone()[0]
    if unusable:
        raise SystemExit(
            f"[ERROR] {unusable} direct implanted-target row(s) have a missing, "
            "non-finite or non-positive implanted signal on the selected scale "
            f"({args.reference_scale})")

    outputs: list[Path] = []

    def copy(query: str, name: str) -> Path:
        path = args.outdir / name
        # Every quantitative output records which reference produced it.
        # Row-level reference_type is preserved by the source query where it
        # exists; only the selected scale is added.
        selected_type = ("read_proportional" if args.reference_scale == "read_proportional"
                         else "profiler_scale_primary")
        stamped = (f"SELECT *, '{args.reference_scale}' AS reference_scale, "
                   f"'{selected_type}' AS selected_reference_type FROM ({query})")
        con.execute(
            f"COPY ({stamped}) TO '{sql_path(path)}' (HEADER, DELIMITER '\\t')"
        )
        outputs.append(path)
        return path

    copy("""
        SELECT cohort, study, sample_id, condition, analysis_population, profiler,
               baseline_id, feature AS target_feature, dose_index,
               nominal_target_fraction, target_fraction_for_feature,
               baseline_abundance_fraction, observed_abundance_fraction,
               response_signal, response_signal_selected, implanted_signal_selected,
               observed_over_expected, absolute_relative_error,
               observed_detected, recovery_class,
               -- Original row-level provenance, never rewritten.
               reference_type
        FROM target_rows ORDER BY cohort, analysis_population, profiler,
                                  feature, sample_id, dose_index
    """, "target_recovery_observations.tsv")

    copy("""
        SELECT cohort, analysis_population, condition, profiler, feature AS target_feature,
               nominal_target_fraction, count(*) AS samples,
               avg(CASE WHEN recovery_class='Good' THEN 1 ELSE 0 END) AS good_fraction,
               avg(CASE WHEN recovery_class='Average' THEN 1 ELSE 0 END) AS average_fraction,
               avg(CASE WHEN recovery_class='Poor / missed' THEN 1 ELSE 0 END) AS poor_fraction,
               median(observed_over_expected) AS median_observed_over_expected,
               median(absolute_relative_error) AS median_absolute_relative_error
        FROM target_rows
        GROUP BY ALL ORDER BY cohort, analysis_population, profiler, feature,
                                  nominal_target_fraction, condition
    """, "target_recovery_summary.tsv")

    copy("""
        WITH unique_baselines AS (
          SELECT DISTINCT cohort, analysis_population, condition, profiler, sample_id, baseline_id,
                 baseline_abundance_fraction,
                 CAST(baseline_abundance_fraction>0 AS INTEGER) AS baseline_detected
          FROM target_rows
          WHERE analysis_population='community'
            AND feature='Fusobacterium nucleatum'
        )
        SELECT cohort, analysis_population, condition, profiler, count(*) AS samples,
               sum(baseline_detected) AS detected_samples,
               avg(baseline_detected) AS baseline_prevalence,
               median(baseline_abundance_fraction) AS median_baseline_abundance
        FROM unique_baselines GROUP BY ALL ORDER BY cohort, analysis_population, profiler, condition
    """, "fnuc_baseline_prevalence.tsv")

    copy("""
        WITH counts AS (
        SELECT cohort, analysis_population, condition, profiler, nominal_target_fraction,
                 recovery_class, count(*) AS samples
          FROM target_rows
          WHERE analysis_population='community'
            AND feature='Fusobacterium nucleatum'
          GROUP BY ALL
        )
        SELECT *, samples / sum(samples) OVER (
          PARTITION BY cohort, analysis_population, condition, profiler, nominal_target_fraction
        ) AS fraction
        FROM counts ORDER BY cohort, analysis_population, profiler, condition,
                             nominal_target_fraction, recovery_class
    """, "fnuc_recovery_classes.tsv")

    biomarker_rows = 0
    if args.biomarker_metrics:
        metrics = sql_path(args.biomarker_metrics)
        con.execute(
            f"CREATE VIEW biomarker_raw AS SELECT * FROM read_csv('{metrics}', "
            "delim='\\t', header=true, all_varchar=true)"
        )
        biomarker_columns = {row[0] for row in con.execute("DESCRIBE biomarker_raw").fetchall()}
        needed = {
            "cohort", "analysis_population", "target_label", "profiler", "contrast",
            "spike_fraction_target", "q_threshold", "target_called", "target_effect",
            "target_q_value",
        }
        absent = sorted(needed - biomarker_columns)
        if absent:
            raise SystemExit("[ERROR] Biomarker metrics lack: " + ", ".join(absent))
        path = copy("""
            WITH selected AS (
              SELECT cohort, analysis_population, target_label, profiler, contrast,
                     CAST(spike_fraction_target AS DOUBLE) AS spike_fraction_target,
                     CAST(q_threshold AS DOUBLE) AS q_threshold,
                     CAST(target_called AS INTEGER) AS target_called,
                     CAST(target_effect AS DOUBLE) AS target_effect,
                     CAST(target_q_value AS DOUBLE) AS target_q_value
              FROM biomarker_raw
              WHERE contrast='spiked_vs_matched_baseline__pooled_conditions'
                AND abs(CAST(q_threshold AS DOUBLE)-0.05)<1e-8
            ), contexts AS (
              SELECT cohort, analysis_population, target_label, profiler,
                     min(CASE WHEN target_called=1 THEN spike_fraction_target END)
                       AS minimum_called_fraction,
                     max(target_called) AS ever_called,
                     count(*) AS tested_doses
              FROM selected GROUP BY ALL
            )
            SELECT * FROM contexts ORDER BY cohort, analysis_population, profiler, target_label
        """, "biomarker_detection_thresholds.tsv")
        biomarker_rows = con.execute(
            f"SELECT count(*) FROM read_csv('{sql_path(path)}', delim='\\t', header=true)"
        ).fetchone()[0]

    summary = args.outdir / "target_recovery_diagnostics.tsv"
    direct_rows = con.execute("SELECT count(*) FROM target_rows").fetchone()[0]
    independent_rows = con.execute(
        "SELECT count(*) FROM target_rows WHERE analysis_population='independent'"
    ).fetchone()[0]
    with summary.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"direct_target_rows\t{direct_rows}\n")
        handle.write(f"independent_target_rows\t{independent_rows}\n")
        handle.write(f"biomarker_threshold_contexts\t{biomarker_rows}\n")
        handle.write("recovery_good_relative_error\t0.10\n")
        handle.write("recovery_average_relative_error\t0.50\n")
        handle.write("biomarker_q_threshold\t0.05\n")
        handle.write("status\tDEVELOPMENT_ONLY\n")
    outputs.append(summary)

    checksum = args.outdir / "target_recovery.sha256"
    inputs = [args.responses] + ([args.biomarker_metrics] if args.biomarker_metrics else [])
    with checksum.open("w", encoding="utf-8") as handle:
        for path in inputs + outputs:
            handle.write(f"{digest(path)}  {path.resolve()}\n")
    (args.outdir / "SUCCESS").write_text(
        f"status\tPASS\ndirect_target_rows\t{direct_rows}\n", encoding="utf-8"
    )
    con.close()
    print(f"[PASS] Target-recovery summaries: {direct_rows} direct target rows")


if __name__ == "__main__":
    main()
