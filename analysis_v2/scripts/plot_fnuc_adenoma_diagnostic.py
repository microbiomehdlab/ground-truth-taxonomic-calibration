#!/usr/bin/env python3
"""Focused, multi-metric F. nucleatum adenoma measurement audit."""
from __future__ import annotations

import argparse
import csv
import html
import math
import statistics
from collections import defaultdict
from pathlib import Path

COHORTS = ("feng", "yachida", "zeller")
PROFILERS = ("kraken2_bracken", "metaphlan4")
FEATURE = "Fusobacterium nucleatum"
FIELDS = ("cohort", "profiler", "feature", "adenoma_effect", "adenoma_lower_95",
          "adenoma_upper_95", "adenoma_q", "crc_effect", "crc_q", "adenoma_positive",
          "adenoma_n", "adenoma_positive_median_percent", "adenoma_all_iqr_percent",
          "community_n", "community_baseline_zero", "community_zero_rescued",
          "community_post_positive", "community_recovery_q1", "community_recovery_median",
          "community_recovery_q3", "community_below_half", "independent_n",
          "independent_recovery_q1", "independent_recovery_median",
          "independent_recovery_q3", "independent_below_half", "fuso_bystander_feature",
          "fuso_bystander_slope", "fuso_bystander_contexts")


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def quartiles(values):
    if not values:
        return (None, None, None)
    if len(values) == 1:
        return (values[0],) * 3
    return tuple(statistics.quantiles(values, n=4, method="inclusive"))


def get_unique(rows, label):
    if len(rows) != 1:
        raise ValueError(f"expected one {label}, found {len(rows)}")
    return rows[0]


def summarize(models, baselines, independent, endpoints, operators=None):
    model_rows = [r for r in read_tsv(models) if r["feature"] == FEATURE]
    base_rows = [r for r in read_tsv(baselines) if r["feature"] == FEATURE and
                 r["condition"] == "Adenoma"]
    ind_rows = [r for r in read_tsv(independent) if r["feature"] == FEATURE and
                r["condition"] == "Adenoma" and
                math.isclose(float(r["dose_percent"]), .01)]
    community = defaultdict(list)
    for r in read_tsv(endpoints):
        if (r["analysis_population"] != "community" or r["assembly_arm"] != "original" or
            r["target_label"] != "Fnuc" or r["condition"] != "Adenoma" or
            not math.isclose(float(r["spike_fraction_target"]), 1e-5, rel_tol=.05)):
            continue
        if r["profiler"] not in PROFILERS or r["cohort"] not in COHORTS:
            continue
        expected_ref = "genome_equivalent" if r["profiler"] == "metaphlan4" else "read_proportional"
        if r["reference_type"] != expected_ref:
            raise ValueError("wrong reference for " + repr((r["cohort"], r["profiler"])))
        implanted = float(r["implanted_signal_profiler_scale"])
        if implanted <= 0:
            raise ValueError("nonpositive implanted signal")
        community[r["cohort"], r["profiler"]].append((
            r["sample_id"], float(r["baseline_abundance_fraction"]),
            float(r["observed_abundance_fraction"]),
            float(r["recovered_spike_signal_profiler_scale"]) / implanted))
    bystanders = defaultdict(list)
    if operators is not None:
        for r in read_tsv(operators):
            if (r["holdout_cohort"] != "NONE" or r["analysis_population"] != "independent" or
                r["target_label"] != "Fnuc" or r["feature_role"] != "bystander" or
                not r["feature"].startswith("Fusobacterium ") or
                r["feature"] == FEATURE or not r["operator_slope"]):
                continue
            if r["cohort"] in COHORTS and r["profiler"] in PROFILERS:
                bystanders[r["cohort"], r["profiler"]].append(r)
    result = []
    for cohort in COHORTS:
        for profiler in PROFILERS:
            key = (cohort, profiler)
            aden = get_unique([r for r in model_rows if (r["cohort"], r["profiler"], r["contrast"]) ==
                               (cohort, profiler, "Adenoma_vs_Control")], "adenoma model " + repr(key))
            crc = get_unique([r for r in model_rows if (r["cohort"], r["profiler"], r["contrast"]) ==
                              (cohort, profiler, "CRC_vs_Control")], "CRC model " + repr(key))
            base = get_unique([r for r in base_rows if (r["cohort"], r["profiler"]) == key],
                              "baseline " + repr(key))
            ind = get_unique([r for r in ind_rows if (r["cohort"], r["profiler"]) == key],
                             "independent spike " + repr(key))
            samples = community[key]
            if not samples or len({r[0] for r in samples}) != len(samples):
                raise ValueError("missing/duplicated 0.001% community samples in " + repr(key))
            zeros = [v for _, v, _, _ in samples if v == 0]
            rescue = sum(b == 0 and o > 0 for _, b, o, _ in samples)
            ratios = [v for _, _, _, v in samples]
            q1, med, q3 = quartiles(ratios)
            edges = bystanders[key]
            strongest = max(edges, key=lambda r: abs(float(r["operator_slope"]))) if edges else None
            result.append(dict(cohort=cohort, profiler=profiler, feature=FEATURE,
                adenoma_effect=aden["effect"], adenoma_lower_95=aden["lower_95"],
                adenoma_upper_95=aden["upper_95"], adenoma_q=aden["q_value"],
                crc_effect=crc["effect"], crc_q=crc["q_value"],
                adenoma_positive=base["positive"], adenoma_n=base["n"],
                adenoma_positive_median_percent=base["positive_abundance_median_percent"],
                adenoma_all_iqr_percent=base["all_sample_iqr_percent"],
                community_n=len(samples), community_baseline_zero=len(zeros),
                community_zero_rescued=rescue,
                community_post_positive=sum(o > 0 for _, _, o, _ in samples),
                community_recovery_q1=q1, community_recovery_median=med,
                community_recovery_q3=q3,
                community_below_half=sum(v < .5 for v in ratios),
                independent_n=ind["n"], independent_recovery_q1=ind["recovery_q1"],
                independent_recovery_median=ind["recovery_median"],
                independent_recovery_q3=ind["recovery_q3"],
                independent_below_half=round(float(ind["fraction_below_half"])*float(ind["n"])),
                fuso_bystander_feature="" if strongest is None else strongest["feature"],
                fuso_bystander_slope="" if strongest is None else strongest["operator_slope"],
                fuso_bystander_contexts="" if strongest is None else strongest["eligible_contexts"]))
    return result


def text(parts, x, y, value, size=15, color="#263645", weight="normal", anchor=None):
    a = f' text-anchor="{anchor}"' if anchor else ""
    parts.append(f'<text x="{x}" y="{y}" font-family="Arial,sans-serif" font-size="{size}" '
                 f'font-weight="{weight}" fill="{color}"{a}>{html.escape(str(value))}</text>')


def scale_bar(parts, x, y, width, fraction, color):
    fraction = min(1, max(0, fraction))
    parts.append(f'<rect x="{x}" y="{y}" width="{width}" height="11" fill="#e3e9ed" rx="4"/>')
    parts.append(f'<rect x="{x}" y="{y}" width="{width*fraction:.1f}" height="11" fill="{color}" rx="4"/>')


def draw(path, rows):
    width, height = 1700, 1050
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>']
    text(parts, 32, 47, "Can measurement explain missing F. nucleatum adenoma calls?", 30, weight="bold")
    text(parts, 32, 78, "Six cohort–profiler audits; community member dose 0.001% (ten-species mixture), independent dose 0.01% · DEVELOPMENT ONLY", 16)
    text(parts, 32, 104, "Each metric retains its own scale. Recovery 1 is ideal; zero-rescue uses paired baseline-negative samples. These are diagnostics, not causal verdicts.", 14)
    colors = {"kraken2_bracken": "#bf5700", "metaphlan4": "#006da2"}
    lookup = {(r["cohort"], r["profiler"]): r for r in rows}
    for ci, cohort in enumerate(COHORTS):
        text(parts, 340 + ci*555, 151, cohort.title(), 22, weight="bold", anchor="middle")
    for pi, profiler in enumerate(PROFILERS):
        y0 = 170 + pi*390
        text(parts, 30, y0+28, "Kraken2 + Bracken" if pi == 0 else "MetaPhlAn 4", 19,
             color=colors[profiler], weight="bold")
        for ci, cohort in enumerate(COHORTS):
            row = lookup[cohort, profiler]
            x = 190+ci*555
            parts.append(f'<rect x="{x}" y="{y0}" width="530" height="367" rx="6" fill="#f4f7f9"/>')
            effect = float(row["adenoma_effect"])
            text(parts, x+13, y0+28,
                 f'Adenoma effect {effect:+.2f} [{float(row["adenoma_lower_95"]):+.2f},'
                 f'{float(row["adenoma_upper_95"]):+.2f}], q={float(row["adenoma_q"]):.2g}', 16, weight="bold")
            text(parts, x+13, y0+52,
                 f'CRC effect {float(row["crc_effect"]):+.2f}, q={float(row["crc_q"]):.2g}', 13)
            bpos, bn = int(row["adenoma_positive"]), int(row["adenoma_n"])
            text(parts, x+13, y0+88, f'Unspiked adenoma detection: {bpos}/{bn}', 15)
            scale_bar(parts, x+15, y0+98, 490, bpos/bn, colors[profiler])
            median = row["adenoma_positive_median_percent"]
            median_txt = "none positive" if median == "" else f'{float(median):.3g}% among positives'
            text(parts, x+13, y0+132,
                 f'Native abundance: {median_txt}; all-sample IQR {float(row["adenoma_all_iqr_percent"]):.3g}%', 13)
            rescued, zero = int(row["community_zero_rescued"]), int(row["community_baseline_zero"])
            text(parts, x+13, y0+169, f'0.001% mixture: zero→positive {rescued}/{zero}', 15)
            scale_bar(parts, x+15, y0+179, 490, rescued/zero if zero else 0, colors[profiler])
            text(parts, x+13, y0+213,
                 f'0.001% recovery: {float(row["community_recovery_median"]):.2f} '
                 f'[{float(row["community_recovery_q1"]):.2f}, {float(row["community_recovery_q3"]):.2f}]; '
                 f'<0.5: {row["community_below_half"]}/{row["community_n"]}', 14)
            text(parts, x+13, y0+247,
                 f'0.01% single-species recovery: {float(row["independent_recovery_median"]):.2f} '
                 f'[{float(row["independent_recovery_q1"]):.2f}, {float(row["independent_recovery_q3"]):.2f}]; '
                 f'<0.5: {row["independent_below_half"]}/{row["independent_n"]}', 14)
            parts.append(f'<line x1="{x+13}" y1="{y0+267}" x2="{x+517}" y2="{y0+267}" stroke="#cad6dc"/>')
            if row["fuso_bystander_feature"]:
                label = row["fuso_bystander_feature"].replace("Fusobacterium ", "F. ")
                text(parts, x+13, y0+292, f'Strongest Fuso bystander: {label}', 13)
                text(parts, x+13, y0+313,
                     f'Independent-spike response slope {float(row["fuso_bystander_slope"]):+.2g} '
                     f'(n={row["fuso_bystander_contexts"]}); not a misassignment fraction', 12)
            else:
                text(parts, x+13, y0+292, "Fusobacterium bystander: no evaluable operator row", 13)
            text(parts, x+13, y0+345, "A good spike response weakens only this measured failure mode.", 12,
                 color="#566672")
    text(parts, 32, height-70, "Community and independent doses are different experiments; do not join their points as one dose curve.", 15)
    text(parts, 32, height-42, "Nonsignificant adenoma q is not biological absence. Read-level QC, power and a direct CRC–Adenoma contrast remain separate checks.", 15)
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matched-dir", required=True, type=Path)
    p.add_argument("--endpoints", required=True, type=Path)
    p.add_argument("--operator", type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    a = p.parse_args()
    files = (a.matched_dir / "matched_species_models.tsv",
             a.matched_dir / "matched_species_baseline_summary.tsv",
             a.matched_dir / "matched_species_spike_summary.tsv", a.endpoints)
    if a.operator:
        files += (a.operator,)
    for path in files:
        if not path.is_file() or not path.stat().st_size:
            p.error("missing or empty input: " + str(path))
    if a.outdir.exists():
        p.error("output exists: " + str(a.outdir))
    try:
        rows = summarize(*files[:3], a.endpoints, a.operator)
        a.outdir.mkdir(parents=True)
        with (a.outdir / "fnuc_adenoma_diagnostic.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        draw(a.outdir / "fnuc_adenoma_diagnostic.svg", rows)
        (a.outdir / "DEVELOPMENT_ONLY.txt").write_text(
            "Descriptive technical diagnostics, not causal attribution or evidence of biological absence.\n",
            encoding="utf-8")
        (a.outdir / "SUCCESS").write_text("status\tPASS\n", encoding="utf-8")
    except (OSError, KeyError, ValueError) as error:
        p.error(str(error))
    print(f"[PASS] F. nucleatum six-context diagnostic: {a.outdir}")


if __name__ == "__main__":
    main()
