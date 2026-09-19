#!/usr/bin/env python3
"""Build a physical-profile, feature-level perturbation response table.

The canonical input is target-expanded: one community profile occurs once for
each implanted member.  This builder deliberately collapses those rows to one
physical observation, reconstructs its implantation vector, and only then
joins the sparse native-profiler abundance table.  The output is Parquet
because the all-feature paired table is too large for a practical TSV.

The read-proportional reference for reported feature j is

    expected_j = (1 - F) * baseline_j + implanted_fraction_j

where F is the sum of the target fractions in the reconstructed input vector.
For legacy Feng/Zeller profiles those fractions are nominal.  Yachida retains
the exact achieved integer-allocation fractions.  MetaPhlAn output remains a
marker-based compositional measurement; it is never pooled with Bracken.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
from pathlib import Path

try:
    import duckdb
except ImportError as error:  # pragma: no cover - exercised by deployment preflight
    raise SystemExit("[ERROR] duckdb is required to build the response table") from error


def sql_path(path: Path) -> str:
    return str(path.resolve()).replace("'", "''")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


# Prespecified relative tolerance between the recorded total implanted fraction
# and the sum of its members (METHODS_DECISION_LOG, 2026-09-14).
ACHIEVED_FRACTION_TOLERANCE = 0.05


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-manifest", required=True, type=Path)
    parser.add_argument("--endpoints", required=True, type=Path)
    parser.add_argument("--abundance", required=True, type=Path)
    parser.add_argument("--feature-aliases", type=Path,
                        help="Frozen TSV mapping profiler/source_feature to canonical_feature.")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--memory-limit", default="16GB")
    parser.add_argument("--expected-community-targets", type=int, default=10)
    parser.add_argument(
        "--metaphlan-reference", choices=["genome_equivalent", "read_proportional"],
        help="Quantitative reference scale for MetaPhlAn rows. Required whenever "
             "MetaPhlAn rows are present; there is no default.")
    parser.add_argument(
        "--target-genome-sizes", type=Path,
        help="TSV of target_label, genome_size_bp from the exact implanted FASTA.")
    parser.add_argument(
        "--effective-genome-sizes", type=Path,
        help="TSV of cohort, sample_id, effective_genome_size_bp (audited G_eff).")
    parser.add_argument(
        "--temp-directory", type=Path,
        help="DuckDB spill directory (default: OUTDIR/.work/tmp).",
    )
    parser.add_argument(
        "--keep-work-files", action="store_true",
        help="Retain the temporary DuckDB database and normalized abundance Parquet.",
    )
    return parser.parse_args()


def require_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty {label}: {path}")


def scalar(con: duckdb.DuckDBPyConnection, query: str):
    return con.execute(query).fetchone()[0]


def main() -> None:
    args = parse_args()
    try:
        for path, label in (
            (args.profile_manifest, "profile manifest"),
            (args.endpoints, "paired endpoints"),
            (args.abundance, "abundance table"),
        ):
            require_file(path, label)
        if args.feature_aliases:
            require_file(args.feature_aliases, "feature aliases")
        if args.threads < 1 or args.expected_community_targets < 2:
            raise ValueError("threads and expected community targets must be positive")
        if args.outdir.exists() and any(args.outdir.iterdir()):
            raise ValueError("OUTDIR must be new or empty")

        args.outdir.mkdir(parents=True, exist_ok=True)
        work = args.outdir / ".work"
        temp = args.temp_directory or work / "tmp"
        work.mkdir()
        temp.mkdir(parents=True, exist_ok=True)
        db_path = work / "response_input.duckdb"
        abundance_parquet = work / "abundance.parquet"

        con = duckdb.connect(str(db_path))
        con.execute(f"PRAGMA threads={args.threads}")
        con.execute(f"PRAGMA memory_limit='{args.memory_limit}'")
        con.execute(f"PRAGMA temp_directory='{sql_path(temp)}'")
        con.execute("PRAGMA preserve_insertion_order=false")

        manifest = sql_path(args.profile_manifest)
        endpoints = sql_path(args.endpoints)
        abundance = sql_path(args.abundance)
        con.execute(
            f"CREATE VIEW manifest_raw AS SELECT * FROM read_csv("
            f"'{manifest}', delim='\\t', header=true, all_varchar=true)"
        )
        con.execute(
            f"CREATE VIEW endpoints_raw AS SELECT * FROM read_csv("
            f"'{endpoints}', delim='\\t', header=true, all_varchar=true)"
        )

        if args.feature_aliases:
            aliases = sql_path(args.feature_aliases)
            con.execute(
                f"CREATE TABLE feature_aliases AS SELECT profiler, source_feature, "
                f"canonical_feature, rationale FROM read_csv('{aliases}', delim='\\t', "
                f"header=true, all_varchar=true)"
            )
            invalid_aliases = scalar(con, """
                SELECT count(*) FROM (
                    SELECT profiler, source_feature, count(*) AS n,
                           min(canonical_feature) AS canonical_feature
                    FROM feature_aliases GROUP BY profiler, source_feature
                ) WHERE n != 1 OR profiler='' OR source_feature=''
                   OR canonical_feature='' OR source_feature=canonical_feature
            """)
            if invalid_aliases:
                raise ValueError("feature alias table contains duplicate or invalid mappings")
        else:
            con.execute("""
                CREATE TABLE feature_aliases (
                    profiler VARCHAR, source_feature VARCHAR,
                    canonical_feature VARCHAR, rationale VARCHAR
                )
            """)

        # Full canonical keys make this a 1:1 analytical-context join.  Names
        # `recorded_total_fraction` and `target_fraction` are intentionally
        # distinct: SQL identifiers are case-insensitive, so F/f aliases are unsafe.
        con.execute("""
            CREATE TABLE target_contexts AS
            SELECT
                m.cohort, m.study, m.analysis_population, m.sample_id,
                m.condition, m.assembly_arm, m.profiler, m.profile_id,
                m.baseline_profile_id, m.source_profile, e.source_baseline_profile,
                m.dose_level, m.target_label,
                coalesce(a.canonical_feature, m.target_feature) AS target_feature,
                CAST(e.spike_fraction_total AS DOUBLE) AS recorded_total_fraction,
                CAST(e.spike_fraction_target AS DOUBLE) AS target_fraction
            FROM manifest_raw AS m
            JOIN endpoints_raw AS e USING (
                cohort, study, analysis_population, sample_id, condition,
                target_label, assembly_arm, profiler, profile_id,
                baseline_profile_id, source_profile
            )
            LEFT JOIN feature_aliases AS a
              ON a.profiler=m.profiler AND a.source_feature=m.target_feature
            WHERE CAST(m.spike_fraction_target AS DOUBLE) > 0
              AND CAST(e.spike_fraction_target AS DOUBLE) > 0
        """)
        # --- profiler-scale reference gate ---------------------------------
        # Bracken reports a read fraction; MetaPhlAn reports a
        # genome-equivalent-like composition. Defaulting is what let the
        # estimand mismatch in, so MetaPhlAn rows require an explicit choice.
        has_metaphlan = scalar(
            con, "SELECT count(*) FROM target_contexts WHERE profiler='metaphlan4'")
        use_ge = False
        if has_metaphlan:
            if args.metaphlan_reference is None:
                raise ValueError(
                    "input contains MetaPhlAn rows; --metaphlan-reference "
                    "{genome_equivalent|read_proportional} is required")
            use_ge = args.metaphlan_reference == "genome_equivalent"
        if use_ge:
            if not (args.target_genome_sizes and args.effective_genome_sizes):
                raise ValueError(
                    "--metaphlan-reference genome_equivalent requires both "
                    "--target-genome-sizes and --effective-genome-sizes")
            require_file(args.target_genome_sizes, "target genome size table")
            require_file(args.effective_genome_sizes, "effective genome size table")
            con.execute(f"""
                CREATE TABLE target_genome_sizes AS
                SELECT target_label,
                       CAST(genome_size_bp AS DOUBLE) AS target_genome_size_bp
                FROM read_csv('{sql_path(args.target_genome_sizes)}',
                              delim='\t', header=true, all_varchar=true)
            """)
            con.execute(f"""
                CREATE TABLE effective_genome_sizes AS
                SELECT cohort, sample_id,
                       CAST(effective_genome_size_bp AS DOUBLE)
                           AS effective_community_genome_size_bp
                FROM read_csv('{sql_path(args.effective_genome_sizes)}',
                              delim='\t', header=true, all_varchar=true)
            """)
            for table, column in (("target_genome_sizes", "target_genome_size_bp"),
                                  ("effective_genome_sizes",
                                   "effective_community_genome_size_bp")):
                bad = scalar(con, f"SELECT count(*) FROM {table} "
                                  f"WHERE {column} IS NULL OR {column} <= 0")
                if bad:
                    raise ValueError(f"{bad} non-positive value(s) in {table}")
            # G_eff is sample-wide, so a duplicate (cohort, sample_id) with a
            # different value is unresolvable here.
            if scalar(con, """
                SELECT count(*) FROM (
                  SELECT cohort, sample_id
                  FROM effective_genome_sizes GROUP BY cohort, sample_id
                  HAVING count(DISTINCT effective_community_genome_size_bp) > 1)
            """):
                raise ValueError("conflicting effective_genome_size_bp for a sample")
            missing_target = scalar(con, """
                SELECT count(DISTINCT t.target_label) FROM target_contexts t
                LEFT JOIN target_genome_sizes g USING (target_label)
                WHERE t.profiler='metaphlan4' AND g.target_genome_size_bp IS NULL
            """)
            if missing_target:
                raise ValueError(
                    f"{missing_target} implanted target(s) lack a genome size")
            missing_geff = scalar(con, """
                SELECT count(*) FROM (
                  SELECT DISTINCT t.cohort, t.sample_id FROM target_contexts t
                  LEFT JOIN effective_genome_sizes e
                    ON e.cohort=t.cohort AND e.sample_id=t.sample_id
                  WHERE t.profiler='metaphlan4'
                    AND e.effective_community_genome_size_bp IS NULL)
            """)
            if missing_geff:
                raise ValueError(
                    f"{missing_geff} MetaPhlAn sample(s) lack an effective genome size")

        joined = scalar(con, "SELECT count(*) FROM target_contexts")
        endpoint_positive = scalar(
            con, "SELECT count(*) FROM endpoints_raw "
                 "WHERE CAST(spike_fraction_target AS DOUBLE) > 0"
        )
        if joined != endpoint_positive:
            raise ValueError(
                f"manifest/endpoints positive join is not 1:1: {joined} versus "
                f"{endpoint_positive}"
            )
        duplicate_contexts = scalar(con, """
            SELECT count(*) FROM (
                SELECT cohort, sample_id, analysis_population, profiler,
                       source_profile, target_label, count(*) AS n
                FROM target_contexts GROUP BY ALL HAVING n != 1
            )
        """)
        if duplicate_contexts:
            raise ValueError("duplicate target contexts after manifest/endpoints join")

        con.execute("""
            CREATE TABLE observation_groups AS
            SELECT
                cohort, study, sample_id, condition, analysis_population,
                assembly_arm, profiler, source_profile AS perturbation_id,
                source_baseline_profile AS baseline_id, dose_level,
                CAST(replace(dose_level, 'dose_', '') AS INTEGER) AS dose_index,
                max(recorded_total_fraction) AS recorded_total_fraction,
                sum(target_fraction) AS effective_total_fraction,
                count(*) AS implanted_target_count,
                string_agg(target_label, ';' ORDER BY target_label) AS implanted_targets,
                string_agg(target_feature, ';' ORDER BY target_label) AS implanted_features,
                string_agg(
                    target_label || '=' || printf('%.17g', target_fraction),
                    ';' ORDER BY target_label
                ) AS implanted_fraction_by_target
            FROM target_contexts
            GROUP BY cohort, study, sample_id, condition, analysis_population,
                     assembly_arm, profiler, source_profile,
                     source_baseline_profile, dose_level
        """)
        con.execute("""
            CREATE TABLE observations AS
            SELECT
                row_number() OVER (
                    ORDER BY cohort, analysis_population, profiler, sample_id,
                             dose_index, perturbation_id
                ) AS observation_id,
                *,
                CASE
                    WHEN analysis_population = 'independent' AND dose_index = 1 THEN 0.0001
                    WHEN analysis_population = 'independent' AND dose_index = 2 THEN 0.0005
                    WHEN analysis_population = 'independent' AND dose_index = 3 THEN 0.001
                    WHEN analysis_population = 'independent' AND dose_index = 4 THEN 0.005
                    WHEN analysis_population = 'independent' AND dose_index = 5 THEN 0.01
                    WHEN analysis_population = 'independent' AND dose_index = 6 THEN 0.05
                    WHEN analysis_population = 'community' AND dose_index = 1 THEN 0.0001
                    WHEN analysis_population = 'community' AND dose_index = 2 THEN 0.0005
                    WHEN analysis_population = 'community' AND dose_index = 3 THEN 0.001
                    WHEN analysis_population = 'community' AND dose_index = 4 THEN 0.005
                    WHEN analysis_population = 'community' AND dose_index = 5 THEN 0.01
                    WHEN analysis_population = 'community' AND dose_index = 6 THEN 0.05
                    WHEN analysis_population = 'community' AND dose_index = 7 THEN 0.1
                END AS nominal_total_fraction,
                CASE
                    WHEN analysis_population = 'independent' AND dose_index = 1 THEN 0.0001
                    WHEN analysis_population = 'independent' AND dose_index = 2 THEN 0.0005
                    WHEN analysis_population = 'independent' AND dose_index = 3 THEN 0.001
                    WHEN analysis_population = 'independent' AND dose_index = 4 THEN 0.005
                    WHEN analysis_population = 'independent' AND dose_index = 5 THEN 0.01
                    WHEN analysis_population = 'independent' AND dose_index = 6 THEN 0.05
                    WHEN analysis_population = 'community' AND dose_index = 1 THEN 0.00001
                    WHEN analysis_population = 'community' AND dose_index = 2 THEN 0.00005
                    WHEN analysis_population = 'community' AND dose_index = 3 THEN 0.0001
                    WHEN analysis_population = 'community' AND dose_index = 4 THEN 0.0005
                    WHEN analysis_population = 'community' AND dose_index = 5 THEN 0.001
                    WHEN analysis_population = 'community' AND dose_index = 6 THEN 0.005
                    WHEN analysis_population = 'community' AND dose_index = 7 THEN 0.01
                END AS nominal_target_fraction,
                CASE WHEN cohort = 'yachida'
                     THEN 'EXACT_ACHIEVED_ALLOCATION'
                     ELSE 'LEGACY_NOMINAL_FRACTION' END AS reference_kind
            FROM observation_groups
        """)
        malformed = scalar(con, """
            SELECT count(*) FROM observations
            WHERE (analysis_population = 'community'
                   AND implanted_target_count != {expected_community_targets})
               OR (analysis_population = 'independent' AND implanted_target_count != 1)
               OR nominal_total_fraction IS NULL OR nominal_target_fraction IS NULL
               OR effective_total_fraction <= 0 OR effective_total_fraction >= 1
        """.format(expected_community_targets=args.expected_community_targets))
        if malformed:
            raise ValueError(f"{malformed} malformed physical perturbation observations")

        # Q_i is only correct if the implanted members reconstruct the recorded
        # total. Tolerance is the prespecified 5% relative achieved-dose window
        # (METHODS_DECISION_LOG, 2026-09-14): integer read allocation produces
        # several-percent deviations at the smallest community doses.
        mismatched = scalar(con, f"""
            SELECT count(*) FROM observations
            WHERE recorded_total_fraction > 0
              AND abs(recorded_total_fraction-effective_total_fraction)
                  / recorded_total_fraction > {ACHIEVED_FRACTION_TOLERANCE}
        """)
        if mismatched:
            worst = scalar(con, """
                SELECT max(abs(recorded_total_fraction-effective_total_fraction)
                           / nullif(recorded_total_fraction, 0)) FROM observations
            """)
            raise ValueError(
                f"{mismatched} observation(s) whose implanted member fractions do "
                f"not reconstruct the recorded total within "
                f"{ACHIEVED_FRACTION_TOLERANCE:.0%} (worst {worst:.3%}); "
                "implanted membership is incomplete or mis-joined")

        # Q_i = sum_k(f_ik * G_eff,i / G_k) over EVERY implanted member of the
        # perturbation, computed once per physical perturbed profile.
        if use_ge:
            con.execute("""
                CREATE TABLE observation_genome_equivalents AS
                SELECT o.observation_id,
                       e.effective_community_genome_size_bp,
                       sum(t.target_fraction
                           * e.effective_community_genome_size_bp
                           / g.target_genome_size_bp)
                           AS total_implanted_genome_equivalent_fraction
                FROM observations o
                JOIN target_contexts t
                  ON o.cohort=t.cohort AND o.analysis_population=t.analysis_population
                 AND o.sample_id=t.sample_id AND o.profiler=t.profiler
                 AND o.perturbation_id=t.source_profile
                JOIN target_genome_sizes g ON g.target_label=t.target_label
                JOIN effective_genome_sizes e
                  ON e.cohort=o.cohort AND e.sample_id=o.sample_id
                WHERE o.profiler='metaphlan4'
                GROUP BY o.observation_id, e.effective_community_genome_size_bp
            """)
            incomplete = scalar(con, """
                SELECT count(*) FROM observations o
                LEFT JOIN observation_genome_equivalents q USING (observation_id)
                WHERE o.profiler='metaphlan4'
                  AND q.total_implanted_genome_equivalent_fraction IS NULL
            """)
            if incomplete:
                raise ValueError(
                    f"{incomplete} MetaPhlAn perturbation(s) have incomplete "
                    "implanted membership or unresolved genome sizes")
            # The renormalisation denominator must stay positive.
            if scalar(con, """
                SELECT count(*) FROM observations o
                JOIN observation_genome_equivalents q USING (observation_id)
                WHERE (1.0-o.effective_total_fraction)
                      + q.total_implanted_genome_equivalent_fraction <= 0
            """):
                raise ValueError("non-positive genome-equivalent renormalisation denominator")
        else:
            con.execute("""
                CREATE TABLE observation_genome_equivalents (
                    observation_id BIGINT,
                    effective_community_genome_size_bp DOUBLE,
                    total_implanted_genome_equivalent_fraction DOUBLE)
            """)
            con.execute("""
                CREATE TABLE target_genome_sizes (
                    target_label VARCHAR, target_genome_size_bp DOUBLE)
            """)
        # Present in both modes so the feature-level join needs no branch.
        con.execute("CREATE VIEW target_genome_sizes_or_empty AS "
                    "SELECT * FROM target_genome_sizes")

        con.execute("""
            CREATE TABLE observation_targets AS
            SELECT o.observation_id, t.target_label, t.target_feature,
                   t.target_fraction
            FROM observations AS o
            JOIN target_contexts AS t
              ON o.cohort = t.cohort
             AND o.analysis_population = t.analysis_population
             AND o.sample_id = t.sample_id
             AND o.profiler = t.profiler
             AND o.perturbation_id = t.source_profile
            ORDER BY o.observation_id, t.target_label
        """)

        # Parse the 6+ GB TSV once.  Subsequent joins use compressed Parquet and
        # can spill safely without retaining all abundance rows in Python memory.
        con.execute(f"""
            COPY (
                SELECT raw.profiler, raw.source_profile,
                       coalesce(alias.canonical_feature, raw.feature) AS feature,
                       sum(CAST(raw.abundance_fraction AS DOUBLE)) AS abundance_fraction
                FROM read_csv('{abundance}', delim='\\t', header=true,
                              all_varchar=true) AS raw
                LEFT JOIN feature_aliases AS alias
                  ON alias.profiler=raw.profiler AND alias.source_feature=raw.feature
                GROUP BY raw.profiler, raw.source_profile,
                         coalesce(alias.canonical_feature, raw.feature)
            ) TO '{sql_path(abundance_parquet)}'
              (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 131072)
        """)
        con.execute(
            f"CREATE VIEW abundance AS SELECT * FROM read_parquet("
            f"'{sql_path(abundance_parquet)}')"
        )
        missing_perturbed = scalar(con, """
            SELECT count(*) FROM observations o
            WHERE NOT EXISTS (
                SELECT 1 FROM abundance a
                WHERE a.profiler=o.profiler AND a.source_profile=o.perturbation_id
            )
        """)
        missing_baseline = scalar(con, """
            SELECT count(*) FROM observations o
            WHERE NOT EXISTS (
                SELECT 1 FROM abundance a
                WHERE a.profiler=o.profiler AND a.source_profile=o.baseline_id
            )
        """)
        if missing_perturbed or missing_baseline:
            raise ValueError(
                "abundance table lacks physical profiles: "
                f"perturbed={missing_perturbed}, baseline={missing_baseline}"
            )

        response_path = args.outdir / "paired_feature_responses.parquet"
        # The third UNION arm guarantees that an undetected implanted feature
        # remains represented, so a direct-target detection failure is auditable.
        con.execute(f"""
            COPY (
                WITH sparse_union AS (
                    SELECT o.observation_id, a.feature,
                           NULL::DOUBLE AS baseline_abundance_fraction,
                           a.abundance_fraction AS observed_abundance_fraction
                    FROM observations o JOIN abundance a
                      ON o.profiler=a.profiler
                     AND o.perturbation_id=a.source_profile
                    UNION ALL
                    SELECT o.observation_id, a.feature,
                           a.abundance_fraction AS baseline_abundance_fraction,
                           NULL::DOUBLE AS observed_abundance_fraction
                    FROM observations o JOIN abundance a
                      ON o.profiler=a.profiler
                     AND o.baseline_id=a.source_profile
                    UNION ALL
                    SELECT observation_id, target_feature,
                           NULL::DOUBLE, NULL::DOUBLE
                    FROM observation_targets
                ), paired AS (
                    SELECT observation_id, feature,
                           coalesce(max(baseline_abundance_fraction), 0.0)
                               AS baseline_abundance_fraction,
                           coalesce(max(observed_abundance_fraction), 0.0)
                               AS observed_abundance_fraction
                    FROM sparse_union GROUP BY observation_id, feature
                ), target_by_feature AS (
                    SELECT ot.observation_id, ot.target_feature AS feature,
                           sum(ot.target_fraction) AS target_fraction_for_feature,
                           sum(ot.target_fraction
                               * coalesce(q.effective_community_genome_size_bp, 0.0)
                               / coalesce(g.target_genome_size_bp, 1.0))
                               AS target_ge_fraction_for_feature,
                           max(g.target_genome_size_bp) AS target_genome_size_bp
                    FROM observation_targets ot
                    LEFT JOIN observation_genome_equivalents q
                           ON q.observation_id = ot.observation_id
                    LEFT JOIN target_genome_sizes_or_empty g
                           ON g.target_label = ot.target_label
                    GROUP BY ot.observation_id, ot.target_feature
                ), base AS (
                    SELECT o.*, p.feature, p.baseline_abundance_fraction,
                           p.observed_abundance_fraction,
                           coalesce(t.target_fraction_for_feature, 0.0)
                               AS target_fraction_for_feature,
                           t.target_genome_size_bp,
                           q.effective_community_genome_size_bp,
                           q.total_implanted_genome_equivalent_fraction,
                           coalesce(t.target_ge_fraction_for_feature, 0.0)
                               AS implanted_genome_equivalent_fraction,
                           -- Read-proportional reference. Kept under its
                           -- original names so existing consumers retain their
                           -- exact meaning.
                           (1.0-o.effective_total_fraction)
                               * p.baseline_abundance_fraction
                               AS dilution_retained_baseline,
                           (1.0-o.effective_total_fraction)
                               * p.baseline_abundance_fraction
                             + coalesce(t.target_fraction_for_feature, 0.0)
                               AS expected_abundance_fraction,
                           -- Profiler scale. Bracken reuses the read-proportional
                           -- value; MetaPhlAn renormalises the composition, which
                           -- is what makes the non-implanted expectation
                           -- e_MP,ij = (1-F)o / ((1-F)+Q) rather than (1-F)o.
                           CASE WHEN q.total_implanted_genome_equivalent_fraction IS NULL
                                THEN (1.0-o.effective_total_fraction)
                                     * p.baseline_abundance_fraction
                                ELSE (1.0-o.effective_total_fraction)
                                     * p.baseline_abundance_fraction
                                     / ((1.0-o.effective_total_fraction)
                                        + q.total_implanted_genome_equivalent_fraction)
                           END AS retained_baseline_profiler_scale,
                           CASE WHEN q.total_implanted_genome_equivalent_fraction IS NULL
                                THEN (1.0-o.effective_total_fraction)
                                     * p.baseline_abundance_fraction
                                   + coalesce(t.target_fraction_for_feature, 0.0)
                                ELSE ((1.0-o.effective_total_fraction)
                                      * p.baseline_abundance_fraction
                                      + coalesce(t.target_ge_fraction_for_feature, 0.0))
                                     / ((1.0-o.effective_total_fraction)
                                        + q.total_implanted_genome_equivalent_fraction)
                           END AS expected_abundance_profiler_scale,
                           -- Implanted signal on the SAME scale as the
                           -- expected value, so that
                           --   expected_ps = retained_ps + implanted_signal_ps
                           -- holds exactly for every row. Bracken: the read
                           -- fraction. MetaPhlAn: q_it / D_i. Non-implanted: 0.
                           CASE WHEN q.total_implanted_genome_equivalent_fraction IS NULL
                                THEN coalesce(t.target_fraction_for_feature, 0.0)
                                ELSE coalesce(t.target_ge_fraction_for_feature, 0.0)
                                     / ((1.0-o.effective_total_fraction)
                                        + q.total_implanted_genome_equivalent_fraction)
                           END AS implanted_signal_profiler_scale,
                           CASE WHEN q.total_implanted_genome_equivalent_fraction IS NULL
                                THEN 'read_proportional' ELSE 'genome_equivalent'
                           END AS reference_type
                    FROM paired p JOIN observations o USING (observation_id)
                    LEFT JOIN target_by_feature t USING (observation_id, feature)
                    LEFT JOIN observation_genome_equivalents q USING (observation_id)
                ), scored AS (
                    SELECT *,
                           expected_abundance_fraction AS read_proportional_reference,
                           target_fraction_for_feature
                               AS implanted_signal_read_proportional,
                           observed_abundance_fraction
                               - retained_baseline_profiler_scale
                               AS response_signal_profiler_scale,
                           observed_abundance_fraction
                               - expected_abundance_profiler_scale
                               AS response_delta_profiler_scale,
                           CASE WHEN observed_abundance_fraction > 0
                                          AND expected_abundance_profiler_scale > 0
                                THEN log2(observed_abundance_fraction
                                          / expected_abundance_profiler_scale)
                                ELSE NULL END
                               AS quantitative_log2_ratio_profiler_scale,
                           CASE WHEN observed_abundance_fraction
                                          + expected_abundance_profiler_scale > 0
                                THEN (observed_abundance_fraction
                                      - expected_abundance_profiler_scale)
                                     / (observed_abundance_fraction
                                        + expected_abundance_profiler_scale)
                                ELSE 0.0 END
                               AS signed_bounded_error_profiler_scale,
                           CASE WHEN observed_abundance_fraction
                                          + expected_abundance_profiler_scale > 0
                                THEN abs(observed_abundance_fraction
                                         - expected_abundance_profiler_scale)
                                     / (observed_abundance_fraction
                                        + expected_abundance_profiler_scale)
                                ELSE 0.0 END
                               AS absolute_bounded_error_profiler_scale,
                           observed_abundance_fraction-dilution_retained_baseline
                               AS response_signal,
                           observed_abundance_fraction-expected_abundance_fraction
                               AS response_delta,
                           CASE WHEN observed_abundance_fraction
                                          + expected_abundance_fraction > 0
                                THEN (observed_abundance_fraction
                                          - expected_abundance_fraction)
                                     / (observed_abundance_fraction
                                          + expected_abundance_fraction)
                                ELSE 0.0 END AS signed_bounded_error,
                           CASE WHEN observed_abundance_fraction
                                          + expected_abundance_fraction > 0
                                THEN abs(observed_abundance_fraction
                                          - expected_abundance_fraction)
                                     / (observed_abundance_fraction
                                          + expected_abundance_fraction)
                                ELSE 0.0 END AS absolute_bounded_error,
                           CASE WHEN observed_abundance_fraction > 0
                                          AND expected_abundance_fraction > 0
                                THEN log2(observed_abundance_fraction
                                          / expected_abundance_fraction)
                                ELSE NULL END AS quantitative_log2_ratio,
                           CAST(baseline_abundance_fraction > 0 AS INTEGER)
                               AS baseline_detected,
                           CAST(observed_abundance_fraction > 0 AS INTEGER)
                               AS observed_detected,
                           CAST(target_fraction_for_feature > 0 AS INTEGER)
                               AS is_direct_target
                    FROM base
                )
                SELECT
                    observation_id, cohort, study, sample_id, condition,
                    analysis_population, assembly_arm, profiler,
                    perturbation_id, baseline_id, dose_level, dose_index,
                    nominal_total_fraction, nominal_target_fraction,
                    recorded_total_fraction, effective_total_fraction,
                    implanted_target_count, implanted_targets, implanted_features,
                    implanted_fraction_by_target, feature,
                    baseline_abundance_fraction, observed_abundance_fraction,
                    target_fraction_for_feature, dilution_retained_baseline,
                    expected_abundance_fraction, response_signal, response_delta,
                    signed_bounded_error, absolute_bounded_error,
                    quantitative_log2_ratio, baseline_detected, observed_detected,
                    CASE
                        WHEN baseline_detected=1 AND observed_detected=1 THEN 'RETAINED'
                        WHEN baseline_detected=1 AND observed_detected=0 THEN 'LOST'
                        WHEN baseline_detected=0 AND observed_detected=1 THEN 'GAINED'
                        ELSE 'STABLE_ABSENT'
                    END AS detection_transition,
                    is_direct_target, reference_kind,
                    -- Explicit reference-scale columns. expected_abundance_fraction
                    -- above remains the READ-PROPORTIONAL reference for every
                    -- profiler; read_proportional_reference is its unambiguous
                    -- alias. The *_profiler_scale columns carry the primary
                    -- publication scale: read-proportional for Bracken,
                    -- genome-equivalent for MetaPhlAn.
                    reference_type, read_proportional_reference,
                    implanted_signal_read_proportional,
                    implanted_signal_profiler_scale,
                    retained_baseline_profiler_scale,
                    expected_abundance_profiler_scale,
                    response_signal_profiler_scale,
                    response_delta_profiler_scale,
                    quantitative_log2_ratio_profiler_scale,
                    signed_bounded_error_profiler_scale,
                    absolute_bounded_error_profiler_scale,
                    effective_community_genome_size_bp,
                    target_genome_size_bp,
                    implanted_genome_equivalent_fraction,
                    total_implanted_genome_equivalent_fraction
                FROM scored
            ) TO '{sql_path(response_path)}'
              (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 131072)
        """)
        con.execute(
            f"CREATE VIEW responses AS SELECT * FROM read_parquet("
            f"'{sql_path(response_path)}')"
        )

        invalid_responses = scalar(con, """
            SELECT count(*) FROM responses
            WHERE baseline_abundance_fraction < 0
               OR observed_abundance_fraction < 0
               OR expected_abundance_fraction < 0
               OR absolute_bounded_error < 0 OR absolute_bounded_error > 1
               OR signed_bounded_error < -1 OR signed_bounded_error > 1
               OR detection_transition NOT IN ('RETAINED','LOST','GAINED','STABLE_ABSENT')
        """)
        if invalid_responses:
            raise ValueError(f"{invalid_responses} invalid feature-response rows")

        observations_path = args.outdir / "observations.tsv"
        targets_path = args.outdir / "observation_targets.tsv"
        con.execute(
            f"COPY (SELECT * FROM observations ORDER BY observation_id) TO "
            f"'{sql_path(observations_path)}' (HEADER, DELIMITER '\\t')"
        )
        con.execute(
            f"COPY (SELECT * FROM observation_targets ORDER BY observation_id, target_label) TO "
            f"'{sql_path(targets_path)}' (HEADER, DELIMITER '\\t')"
        )

        # The profiler-scale expectation must decompose exactly into retained
        # baseline plus implanted signal; otherwise a recovery ratio built from
        # these columns would be internally inconsistent.
        identity_violations = scalar(con, """
            SELECT count(*) FROM responses
            WHERE abs(expected_abundance_profiler_scale
                      - (retained_baseline_profiler_scale
                         + implanted_signal_profiler_scale))
                  > 1e-12 + 1e-9 * abs(expected_abundance_profiler_scale)
        """)
        if identity_violations:
            raise ValueError(
                f"{identity_violations} row(s) violate "
                "expected_abundance_profiler_scale = "
                "retained_baseline_profiler_scale + implanted_signal_profiler_scale")

        # Bounded errors are ratios of like quantities and must stay in range.
        bounded_violations = scalar(con, """
            SELECT count(*) FROM responses
            WHERE NOT isfinite(signed_bounded_error_profiler_scale)
               OR NOT isfinite(absolute_bounded_error_profiler_scale)
               OR signed_bounded_error_profiler_scale < -1
               OR signed_bounded_error_profiler_scale > 1
               OR absolute_bounded_error_profiler_scale < 0
               OR absolute_bounded_error_profiler_scale > 1
               OR abs(absolute_bounded_error_profiler_scale
                      - abs(signed_bounded_error_profiler_scale)) > 1e-12
        """)
        if bounded_violations:
            raise ValueError(
                f"{bounded_violations} row(s) have out-of-range or inconsistent "
                "profiler-scale bounded errors")

        n_observations = scalar(con, "SELECT count(*) FROM observations")
        n_responses = scalar(con, "SELECT count(*) FROM responses")
        n_direct = scalar(con, "SELECT count(*) FROM responses WHERE is_direct_target=1")
        max_fraction_difference = scalar(con, """
            SELECT max(abs(recorded_total_fraction-effective_total_fraction))
            FROM observations
        """)
        cohort_rows = con.execute("""
            SELECT cohort, analysis_population, profiler,
                   count(*) AS physical_observations,
                   count(DISTINCT sample_id) AS biological_samples
            FROM observations GROUP BY ALL ORDER BY ALL
        """).fetchall()

        summary_path = args.outdir / "response_input_summary.tsv"
        with summary_path.open("w", encoding="utf-8") as handle:
            handle.write("metric\tvalue\n")
            handle.write(f"target_expanded_positive_contexts\t{joined}\n")
            handle.write(f"physical_perturbation_observations\t{n_observations}\n")
            handle.write(f"paired_feature_rows\t{n_responses}\n")
            handle.write(f"direct_target_rows\t{n_direct}\n")
            handle.write(f"bystander_rows\t{n_responses-n_direct}\n")
            handle.write(
                "max_recorded_vs_allocated_total_fraction_difference\t"
                f"{max_fraction_difference:.17g}\n"
            )
            handle.write("community_rows_collapsed_to_physical_profiles\tYES\n")
            handle.write(
                f"feature_equivalence_aliases\t"
                f"{scalar(con, 'SELECT count(*) FROM feature_aliases')}\n"
            )
            handle.write("aliased_abundances_summed_before_pairing\tYES\n")
            handle.write("missing_sparse_features_interpreted_within_profiler_universe\tZERO\n")
            handle.write("empirical_sham_null_available\tNO\n")
            handle.write("status\tDEVELOPMENT_ONLY\n")
        cohort_summary = args.outdir / "observation_coverage.tsv"
        with cohort_summary.open("w", encoding="utf-8") as handle:
            handle.write(
                "cohort\tanalysis_population\tprofiler\tphysical_observations\t"
                "biological_samples\n"
            )
            for row in cohort_rows:
                handle.write("\t".join(map(str, row)) + "\n")

        note = args.outdir / "DEVELOPMENT_ONLY.txt"
        note.write_text(
            "status=DEVELOPMENT_ONLY\n"
            "legacy_feng_zeller_fractions=nominal\n"
            "yachida_fractions=exact_achieved_allocations\n"
            "empirical_sham_null=UNAVAILABLE\n"
            "automatic_feature_removal=NO\n"
            "false_positive_probability=NOT_ESTIMATED\n",
            encoding="utf-8",
        )
        manifest_path = args.outdir / "response_input.sha256"
        inputs = [args.profile_manifest, args.endpoints, args.abundance]
        if args.feature_aliases:
            inputs.append(args.feature_aliases)
        outputs = [response_path, observations_path, targets_path, summary_path,
                   cohort_summary, note]
        with manifest_path.open("w", encoding="utf-8") as handle:
            for path in inputs + outputs:
                handle.write(f"{sha256(path)}  {path.resolve()}\n")
        (args.outdir / "SUCCESS").write_text(
            f"status=PASS\nphysical_observations={n_observations}\n"
            f"paired_feature_rows={n_responses}\n",
            encoding="utf-8",
        )
        con.close()
        if not args.keep_work_files:
            shutil.rmtree(work)
        print(
            f"[PASS] Physical perturbation response input: {n_observations} profiles, "
            f"{n_responses} feature rows"
        )
        print(f"[INFO] {response_path}")
    except (OSError, ValueError, duckdb.Error) as error:
        raise SystemExit(f"[ERROR] {error}") from error


if __name__ == "__main__":
    main()
