#!/usr/bin/env python3
"""Descriptive bridge between unspiked CRC replication and local spike reliability.

One observation is a source-cohort significant CRC call evaluated in a *different*
cohort. Reliability is measured from target-excluded spike responses in that
destination cohort. This is an association audit, not a causal or predictive model.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
from pathlib import Path

import duckdb

COHORTS = ("feng", "yachida", "zeller")
PROFILERS = ("kraken2_bracken", "metaphlan4")
FIELDS = ("source_cohort", "destination_cohort", "profiler", "feature", "source_effect",
          "source_q", "source_n", "destination_effect", "destination_q", "destination_n",
          "replication_status", "replicated", "destination_reliability",
          "destination_reliability_contexts", "training_reliability", "training_contexts")


def quote(path: Path) -> str:
    return str(path.resolve()).replace("'", "''")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def alias_rows(path: Path) -> list[tuple[str, str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not {"canonical", "alias", "tool"} <= set(reader.fieldnames or ()):
            raise ValueError("alias CSV needs canonical, alias, tool columns")
        rows = []
        seen = {}
        for row in reader:
            profiler, source, canonical = (row[key].strip() for key in ("tool", "alias", "canonical"))
            if not all((profiler, source, canonical)):
                raise ValueError("blank alias mapping")
            key = profiler, source
            if key in seen and seen[key] != canonical:
                raise ValueError("conflicting alias mapping: " + repr(key))
            seen[key] = canonical
        return [(p, s, c) for (p, s), c in sorted(seen.items())]


def make_rows(calls: Path, certificates: Path, aliases: Path) -> list[dict]:
    con = duckdb.connect(database=":memory:")
    try:
        con.execute("PRAGMA threads=2")
        con.execute("CREATE TABLE aliases(profiler VARCHAR, source_feature VARCHAR, canonical_feature VARCHAR)")
        mappings = alias_rows(aliases)
        if mappings:
            con.executemany("INSERT INTO aliases VALUES (?, ?, ?)", mappings)
        call_path, cert_path = quote(calls), quote(certificates)
        con.execute(f"CREATE VIEW disease AS SELECT * FROM read_csv('{call_path}', delim='\\t', header=true, all_varchar=true)")
        con.execute(f"CREATE VIEW cert AS SELECT * FROM read_csv('{cert_path}', delim='\\t', header=true, all_varchar=true)")
        needed_calls = {"cohort", "analysis_population", "assembly_arm", "profiler", "contrast",
                        "include", "model_spec", "spike_fraction_target", "feature", "effect",
                        "q_value", "n_samples"}
        needed_cert = {"holdout_cohort", "profiler", "feature", "overall_reliability",
                       "eligible_contexts", "evaluation_overall_reliability", "evaluation_contexts"}
        for view, needed in (("disease", needed_calls), ("cert", needed_cert)):
            found = {r[1] for r in con.execute(f"PRAGMA table_info('{view}')").fetchall()}
            if not needed <= found:
                raise ValueError(f"{view} lacks columns: {sorted(needed - found)}")
        sql = """
        WITH baseline_raw AS (
          SELECT d.cohort, d.profiler,
                 coalesce(a.canonical_feature, d.feature) AS feature,
                 try_cast(d.effect AS DOUBLE) AS effect,
                 try_cast(d.q_value AS DOUBLE) AS q,
                 try_cast(d.n_samples AS INTEGER) AS n
          FROM disease d LEFT JOIN aliases a
            ON a.profiler=d.profiler AND a.source_feature=d.feature
          WHERE d.analysis_population='community' AND d.assembly_arm='original'
            AND d.contrast='CRC_vs_Control' AND d.include='1'
            AND d.model_spec='primary_age_sex'
            AND try_cast(d.spike_fraction_target AS DOUBLE)=0
            AND d.cohort IN ('feng','yachida','zeller')
            AND d.profiler IN ('kraken2_bracken','metaphlan4')
        ), baseline AS (
          SELECT cohort, profiler, feature, min(effect) AS effect, min(q) AS q,
                 min(n) AS n, count(DISTINCT (effect,q,n)) AS distinct_fits
          FROM baseline_raw GROUP BY 1,2,3
        ), certificate AS (
          SELECT c.holdout_cohort, c.profiler,
                 coalesce(a.canonical_feature, c.feature) AS feature,
                 try_cast(evaluation_overall_reliability AS DOUBLE) AS destination_reliability,
                 try_cast(evaluation_contexts AS INTEGER) AS destination_contexts,
                 try_cast(overall_reliability AS DOUBLE) AS training_reliability,
                 try_cast(eligible_contexts AS INTEGER) AS training_contexts
          FROM cert c LEFT JOIN aliases a
            ON a.profiler=c.profiler AND a.source_feature=c.feature
          WHERE c.holdout_cohort IN ('feng','yachida','zeller')
        )
        SELECT s.cohort AS source_cohort, dest.cohort AS destination_cohort,
               s.profiler, s.feature, s.effect AS source_effect, s.q AS source_q,
               s.n AS source_n, d.effect AS destination_effect,
               d.q AS destination_q, d.n AS destination_n,
               CASE WHEN d.feature IS NULL THEN 'NOT_EVALUABLE'
                    WHEN sign(s.effect) != sign(d.effect) THEN 'DIRECTION_REVERSED'
                    WHEN d.q <= 0.05 THEN 'REPLICATED'
                    ELSE 'SAME_DIRECTION_NOT_SIGNIFICANT' END AS replication_status,
               CASE WHEN d.feature IS NOT NULL AND sign(s.effect)=sign(d.effect)
                              AND d.q <= 0.05 THEN 1 ELSE 0 END AS replicated,
               c.destination_reliability,
               c.destination_contexts AS destination_reliability_contexts,
               c.training_reliability, c.training_contexts,
               s.distinct_fits, coalesce(d.distinct_fits,1) AS destination_distinct_fits
        FROM baseline s
        CROSS JOIN (SELECT unnest(['feng','yachida','zeller']) AS cohort) dest
        LEFT JOIN baseline d ON d.cohort=dest.cohort AND d.profiler=s.profiler
                            AND d.feature=s.feature
        LEFT JOIN certificate c ON c.holdout_cohort=dest.cohort
                               AND c.profiler=s.profiler AND c.feature=s.feature
        WHERE s.q <= 0.05 AND dest.cohort != s.cohort
        ORDER BY s.profiler, s.feature, s.cohort, dest.cohort
        """
        cursor = con.execute(sql)
        names = [x[0] for x in cursor.description]
        raw = [dict(zip(names, values)) for values in cursor.fetchall()]
        if not raw:
            raise ValueError("no eligible baseline CRC calls")
        for row in raw:
            if row.pop("distinct_fits") != 1 or row.pop("destination_distinct_fits") != 1:
                raise ValueError("contradictory repeated baseline fits")
            if row["source_effect"] is None or row["source_q"] is None or row["source_n"] is None:
                raise ValueError("invalid source model fit")
            if row["destination_q"] is None and row["destination_effect"] is not None:
                raise ValueError("invalid destination model fit")
            if row["destination_q"] is not None and not 0 <= row["destination_q"] <= 1:
                raise ValueError("invalid destination q")
            score = row["destination_reliability"]
            if score is not None and not 0 <= score <= 1:
                raise ValueError("destination reliability outside [0,1]")
        return raw
    finally:
        con.close()


def write_tsv(path: Path, fields, rows) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_summary(rows: list[dict]) -> list[dict]:
    result = []
    for profiler in PROFILERS:
        for cohort in COHORTS:
            subset = [r for r in rows if r["profiler"] == profiler and r["destination_cohort"] == cohort]
            scored = [r for r in subset if r["destination_reliability"] is not None
                      and r["replication_status"] != "NOT_EVALUABLE"]
            scored.sort(key=lambda r: r["destination_reliability"])
            # Empirical tertiles are descriptive and defined within each destination stratum.
            for index, label in enumerate(("lower", "middle", "upper")):
                group = scored[index * len(scored) // 3:(index + 1) * len(scored) // 3]
                n = len(group)
                result.append(dict(profiler=profiler, destination_cohort=cohort,
                                   reliability_tertile=label, n=n,
                                   replicated=sum(r["replicated"] for r in group),
                                   replication_fraction=(sum(r["replicated"] for r in group) / n if n else ""),
                                   median_reliability=(group[n // 2]["destination_reliability"] if n else ""),
                                   missing_destination_fit=sum(r["replication_status"] == "NOT_EVALUABLE" for r in subset),
                                   missing_spike_score=sum(r["destination_reliability"] is None for r in subset)))
    return result


def plot(summary: list[dict], path: Path) -> None:
    width, height = 1440, 790
    panel_width, panel_height = 420, 300
    left, top = 65, 105
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#18232c}.title{font-size:27px;font-weight:700}'
        '.subtitle{font-size:15px}.panel{font-size:18px;font-weight:700}'
        '.tick{font-size:13px}.count{font-size:15px;font-weight:700}'
        '.note{font-size:12px;fill:#53616c}</style>',
        '<text class="title" x="65" y="42">Does spike-measured reliability track CRC biomarker replication?</text>',
        '<text class="subtitle" x="65" y="72">Same-direction BH q ≤ 0.05 in destination cohort; each bar shows replicated/evaluable calls</text>',
    ]
    for i, profiler in enumerate(PROFILERS):
        for j, cohort in enumerate(COHORTS):
            group = [r for r in summary if r["profiler"] == profiler and r["destination_cohort"] == cohort]
            x0, y0 = left + j * 455, top + i * 335
            plot_left, plot_top = x0 + 48, y0 + 42
            plot_width, plot_height = 345, 175
            title = f'{cohort.title()} · {"Kraken2 + Bracken" if i == 0 else "MetaPhlAn 4"}'
            elements.append(f'<rect x="{x0}" y="{y0}" width="{panel_width}" height="{panel_height}" '
                            'fill="#fff" stroke="#cad3da"/>')
            elements.append(f'<text class="panel" x="{x0 + 15}" y="{y0 + 28}">{html.escape(title)}</text>')
            for fraction in (0, .5, 1):
                y = plot_top + plot_height * (1 - fraction)
                elements.append(f'<line x1="{plot_left}" y1="{y:.1f}" x2="{plot_left + plot_width}" '
                                f'y2="{y:.1f}" stroke="#e0e5e8"/>')
                elements.append(f'<text class="tick" x="{plot_left - 10}" y="{y + 4:.1f}" '
                                f'text-anchor="end">{fraction:.0%}</text>')
            for k, row in enumerate(group):
                center = plot_left + 62 + k * 112
                if row["n"]:
                    bar_height = plot_height * row["replication_fraction"]
                    y = plot_top + plot_height - bar_height
                    color = "#168a78" if i == 0 else "#3b7097"
                    elements.append(f'<rect x="{center - 32}" y="{y:.1f}" width="64" '
                                    f'height="{bar_height:.1f}" fill="{color}"/>')
                    label = f'{row["replicated"]}/{row["n"]}'
                else:
                    label = "no data"
                elements.append(f'<text class="count" x="{center}" y="{max(plot_top + 15, y - 8) if row["n"] else plot_top + plot_height - 12:.1f}" '
                                f'text-anchor="middle">{label}</text>')
                elements.append(f'<text class="tick" x="{center}" y="{plot_top + plot_height + 20}" '
                                f'text-anchor="middle">{row["reliability_tertile"].title()}</text>')
            elements.append(f'<text class="note" x="{x0 + 15}" y="{y0 + 277}">'
                            f'Missing fit: {group[0]["missing_destination_fit"]}; '
                            f'no local spike score: {group[0]["missing_spike_score"]}</text>')
    elements.extend([
        '<text class="subtitle" x="65" y="785">Destination-cohort spike reliability tertile · DEVELOPMENT ONLY · descriptive, not causal</text>',
        '</svg>',
    ])
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calls", type=Path, required=True)
    parser.add_argument("--certificates", type=Path, required=True)
    parser.add_argument("--aliases", type=Path, default=Path("examples/spike_taxon_aliases.csv"))
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    for path in (args.calls, args.certificates, args.aliases):
        if not path.is_file() or not path.stat().st_size:
            parser.error(f"missing or empty input: {path}")
    if args.outdir.exists():
        parser.error(f"output already exists: {args.outdir}")
    rows = make_rows(args.calls, args.certificates, args.aliases)
    summary = make_summary(rows)
    if not any(row["n"] for row in summary):
        parser.error("no candidate pairs have both a destination fit and a destination spike score")
    args.outdir.mkdir(parents=True)
    write_tsv(args.outdir / "candidate_destination_bridge.tsv", FIELDS, rows)
    write_tsv(args.outdir / "replication_by_reliability_tertile.tsv", tuple(summary[0]), summary)
    plot(summary, args.outdir / "replication_by_reliability_tertile.svg")
    with (args.outdir / "source_sha256.tsv").open("w", encoding="utf-8") as handle:
        handle.write("source\tsha256\n")
        for path in (args.calls, args.certificates, args.aliases):
            handle.write(f"{path.resolve()}\t{sha256(path)}\n")
    (args.outdir / "DEVELOPMENT_ONLY.txt").write_text(
        "Descriptive source-to-destination association; not a causal explanation, "
        "prospective prediction, or final manuscript result.\n", encoding="utf-8")
    print(f"[PASS] {len(rows)} source-to-destination CRC candidate evaluations: {args.outdir}")


if __name__ == "__main__":
    main()
