#!/usr/bin/env python3
"""Full-baseline and paired low-dose spike atlas for two discordant CRC taxa."""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import math
import statistics
from collections import defaultdict
from pathlib import Path

from diagnose_crc_candidate_recovery import (
    CASES, CONDITIONS, DOSES, PROFILERS, aliases, endpoint_rows, model_rows,
    quartiles, summarize,
)

MANIFEST_FIELDS = {"cohort", "sample_id", "condition", "analysis_population",
                   "assembly_arm", "profiler", "dose_level", "include", "source_profile"}
ABUNDANCE_FIELDS = {"profiler", "source_profile", "feature", "abundance_fraction"}
BASELINE_FIELDS = ("cohort", "target_label", "canonical_taxon", "condition",
                   "sample_id", "profiler", "reported_feature", "abundance_fraction")
SUMMARY_FIELDS = ("cohort", "target_label", "canonical_taxon", "condition", "profiler",
                  "n_full_baseline", "baseline_positive", "baseline_prevalence",
                  "baseline_zero", "positive_abundance_q1_percent",
                  "positive_abundance_median_percent", "positive_abundance_q3_percent",
                  "all_sample_median_percent", "all_sample_iqr_percent",
                  "paired_spike_n_at_0p01", "paired_spike_recovery_q1_at_0p01",
                  "paired_spike_recovery_median_at_0p01", "paired_spike_recovery_q3_at_0p01",
                  "paired_spike_fraction_below_half_at_0p01")
COLORS = {"kraken2_bracken": "#d55e00", "metaphlan4": "#0072b2"}


def write_tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def full_baselines(manifest_path, abundance_path, alias_map):
    wanted = {(cohort, profiler): (target, canonical, alias_map[canonical, profiler])
              for cohort, target, canonical in CASES for profiler in PROFILERS}
    profiles = {}
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not MANIFEST_FIELDS <= set(reader.fieldnames or ()):
            raise ValueError("native profile manifest lacks baseline context columns")
        for row in reader:
            context = row["cohort"], row["profiler"]
            if context not in wanted or row["analysis_population"] != "community" or row["assembly_arm"] != "original" or row["dose_level"] != "baseline" or row["include"] != "1":
                continue
            if row["condition"] not in CONDITIONS:
                continue
            key = row["cohort"], row["condition"], row["sample_id"], row["profiler"]
            profile = str(Path(row["source_profile"]).resolve())
            if key in profiles and profiles[key] != profile:
                raise ValueError("conflicting baseline profiles for " + repr(key))
            profiles[key] = profile
    for cohort, target, _ in CASES:
        for condition in CONDITIONS:
            sample_sets = [{key[2] for key in profiles if key[0] == cohort and key[1] == condition and key[3] == profiler}
                           for profiler in PROFILERS]
            if sample_sets[0] != sample_sets[1] or len(sample_sets[0]) < 2:
                raise ValueError("full baseline profiler samples differ or are too few for {} {}".format(cohort, condition))
    keys = {(profiler, path, wanted[cohort, profiler][2])
            for (cohort, _, _, profiler), path in profiles.items()}
    found = {}
    with abundance_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not ABUNDANCE_FIELDS <= set(reader.fieldnames or ()):
            raise ValueError("native abundance long table lacks required columns")
        for row in reader:
            key = row["profiler"], str(Path(row["source_profile"]).resolve()), row["feature"]
            if key not in keys:
                continue
            value = float(row["abundance_fraction"])
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError("invalid native baseline abundance")
            if key in found:
                raise ValueError("duplicate native baseline feature " + repr(key))
            found[key] = value
    result = []
    for (cohort, condition, sample, profiler), profile in sorted(profiles.items()):
        target, canonical, feature = wanted[cohort, profiler]
        result.append(dict(cohort=cohort, target_label=target, canonical_taxon=canonical,
                           condition=condition, sample_id=sample, profiler=profiler,
                           reported_feature=feature,
                           abundance_fraction=found.get((profiler, profile, feature), 0.0)))
    return result


def summary_rows(baselines, spikes, spike_samples):
    result = []
    for cohort, target, canonical in CASES:
        for condition in CONDITIONS:
            for profiler in PROFILERS:
                values = sorted(r["abundance_fraction"] for r in baselines if
                                r["cohort"] == cohort and r["condition"] == condition and
                                r["profiler"] == profiler)
                if not values:
                    raise ValueError("empty full baseline group")
                positive = [v for v in values if v > 0]
                p1, pm, p3 = quartiles(positive) if positive else (0, 0, 0)
                a1, am, a3 = quartiles(values)
                low = next(r for r in spikes if r["cohort"] == cohort and r["condition"] == condition and
                           r["profiler"] == profiler and r["dose_fraction"] == DOSES[0])
                low_samples = [r for r in spike_samples if r["cohort"] == cohort and
                               r["condition"] == condition and r["profiler"] == profiler and
                               r["dose_fraction"] == DOSES[0]]
                if len(low_samples) != low["n_paired"]:
                    raise ValueError("low-dose sample count disagrees with summary")
                result.append(dict(cohort=cohort, target_label=target, canonical_taxon=canonical,
                                   condition=condition, profiler=profiler,
                                   n_full_baseline=len(values), baseline_positive=len(positive),
                                   baseline_prevalence=len(positive) / len(values),
                                   baseline_zero=len(values) - len(positive),
                                   positive_abundance_q1_percent=100 * p1,
                                   positive_abundance_median_percent=100 * pm,
                                   positive_abundance_q3_percent=100 * p3,
                                   all_sample_median_percent=100 * am,
                                   all_sample_iqr_percent=100 * (a3 - a1),
                                   paired_spike_n_at_0p01=low["n_paired"],
                                   paired_spike_recovery_q1_at_0p01=low["recovery_ratio_q1"],
                                   paired_spike_recovery_median_at_0p01=low["recovery_ratio_median"],
                                   paired_spike_recovery_q3_at_0p01=low["recovery_ratio_q3"],
                                   paired_spike_fraction_below_half_at_0p01=sum(
                                       r["recovery_ratio"] < .5 for r in low_samples) / len(low_samples)))
    return result


def svg_text(parts, x, y, label, size=15, weight="normal", color="#1b2933", anchor=None):
    anchor_attr = f' text-anchor="{anchor}"' if anchor else ""
    parts.append(f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial,sans-serif" '
                 f'font-size="{size}" font-weight="{weight}" fill="{color}"{anchor_attr}>'
                 f'{html.escape(str(label))}</text>')


def draw(path, cohort, target, canonical, model, baselines, baseline_summary, spikes, spike_samples):
    width, height = 1600, 1120
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>']
    svg_text(parts, 42, 52, f'Why is the {canonical} CRC call profiler-dependent in {cohort.title()}?', 29, "bold")
    svg_text(parts, 42, 82, "Full unspiked community baselines above; same-sample independent spikes below · DEVELOPMENT ONLY", 16)
    for i, profiler in enumerate(PROFILERS):
        x = 42 + i * 785
        m = next(r for r in model if r["profiler"] == profiler)
        name = "Kraken2 + Bracken" if i == 0 else "MetaPhlAn 4"
        parts.append(f'<rect x="{x}" y="105" width="745" height="105" rx="5" fill="#f3f6f8"/>')
        svg_text(parts, x + 15, 137, name, 21, "bold", COLORS[profiler])
        svg_text(parts, x + 15, 165, f'CRC effect {m["effect"]:+.2f} [95% CI {m["lower_95"]:+.2f}, {m["upper_95"]:+.2f}]', 17)
        svg_text(parts, x + 15, 190, f'BH q={m["q_value"]:.3g} · model n={m["n_samples"]} · feature: {m["reported_feature"]}', 15)
    svg_text(parts, 42, 248, "A. Full-cohort unspiked baseline: detectability and abundance variability", 22, "bold")
    svg_text(parts, 42, 273, "Each dot is one sample. The zero column is separate; positive abundances use a shared log10-percent axis.", 15)
    positive = [100 * r["abundance_fraction"] for r in baselines if r["abundance_fraction"] > 0]
    low_log = math.floor(math.log10(min(positive))) if positive else -5
    high_log = math.ceil(math.log10(max(positive))) if positive else 1
    if high_log <= low_log:
        high_log = low_log + 1
    axis_x0, axis_x1 = 205, 1160
    def abundance_x(value):
        return axis_x0 + (math.log10(value) - low_log) / (high_log - low_log) * (axis_x1 - axis_x0)
    for tick in range(low_log, high_log + 1):
        x = abundance_x(10 ** tick)
        parts.append(f'<line x1="{x:.1f}" y1="302" x2="{x:.1f}" y2="575" stroke="#e5eaed"/>')
        svg_text(parts, x, 295, f'10^{tick}%', 13, anchor="middle")
    svg_text(parts, 92, 295, "zero", 13, anchor="middle")
    for i, profiler in enumerate(PROFILERS):
        for j, condition in enumerate(CONDITIONS):
            y = 335 + i * 128 + j * 58
            group = sorted((r for r in baselines if r["profiler"] == profiler and r["condition"] == condition),
                           key=lambda r: r["sample_id"])
            s = next(r for r in baseline_summary if r["profiler"] == profiler and r["condition"] == condition)
            svg_text(parts, 42, y + 5, ("K+B" if i == 0 else "Meta") + " " + condition, 16, "bold", COLORS[profiler])
            parts.append(f'<line x1="45" y1="{y+21}" x2="1550" y2="{y+21}" stroke="#eef0f2"/>')
            for k, row in enumerate(group):
                value = 100 * row["abundance_fraction"]
                x = 92 if value == 0 else abundance_x(value)
                jitter = ((k * 17) % 11 - 5) * 1.8
                parts.append(f'<circle cx="{x:.1f}" cy="{y+jitter:.1f}" r="2.5" '
                             f'fill="{COLORS[profiler]}" fill-opacity="0.34"/>')
            if s["baseline_positive"]:
                x1 = abundance_x(s["positive_abundance_q1_percent"])
                xm = abundance_x(s["positive_abundance_median_percent"])
                x3 = abundance_x(s["positive_abundance_q3_percent"])
                parts.append(f'<line x1="{x1:.1f}" y1="{y+12}" x2="{x3:.1f}" y2="{y+12}" '
                             f'stroke="{COLORS[profiler]}" stroke-width="5"/>')
                parts.append(f'<circle cx="{xm:.1f}" cy="{y+12}" r="5" fill="{COLORS[profiler]}"/>')
            svg_text(parts, 1190, y - 3,
                     f'{s["baseline_positive"]}/{s["n_full_baseline"]} positive ({100*s["baseline_prevalence"]:.0f}%)', 15)
            med = s["positive_abundance_median_percent"]
            iqr = s["positive_abundance_q1_percent"], s["positive_abundance_q3_percent"]
            svg_text(parts, 1190, y + 17,
                     f'positive median {med:.3g}% [IQR {iqr[0]:.3g}, {iqr[1]:.3g}]', 13)
    svg_text(parts, 42, 633, "B. Direct spike recovery in paired samples (0.01%, 0.05%, 0.10%)", 22, "bold")
    svg_text(parts, 42, 659, "Recovery ratio = measured added signal / expected added signal; 1 is ideal. Points are samples; thick marks are median and Q1–Q3.", 15)
    chart_x0, chart_x1 = 120, 1180
    chart_y0, chart_y1 = 698, 955
    ymin, ymax = -.5, 2.0
    def recovery_y(value):
        return chart_y1 - (min(max(value, ymin), ymax) - ymin) / (ymax - ymin) * (chart_y1 - chart_y0)
    for tick in (0, .5, 1, 1.5, 2):
        y = recovery_y(tick)
        parts.append(f'<line x1="{chart_x0}" y1="{y:.1f}" x2="{chart_x1}" y2="{y:.1f}" '
                     f'stroke="{("#75828c" if tick == 1 else "#e5eaed")}" '
                     f'stroke-dasharray="{("5 4" if tick == 1 else "none")}"/>')
        svg_text(parts, 105, y + 5, str(tick), 13, anchor="end")
    clipped = 0
    for i, profiler in enumerate(PROFILERS):
        for j, condition in enumerate(CONDITIONS):
            for k, dose in enumerate(DOSES):
                x = 215 + k * 325 + i * 135 + j * 55
                rows = sorted((r for r in spike_samples if r["profiler"] == profiler and
                               r["condition"] == condition and r["dose_fraction"] == dose),
                              key=lambda r: r["sample_id"])
                s = next(r for r in spikes if r["profiler"] == profiler and
                         r["condition"] == condition and r["dose_fraction"] == dose)
                for n, row in enumerate(rows):
                    value = row["recovery_ratio"]
                    clipped += value < ymin or value > ymax
                    dx = ((n * 7) % 9 - 4) * 2
                    parts.append(f'<circle cx="{x+dx:.1f}" cy="{recovery_y(value):.1f}" r="3" '
                                 f'fill="{COLORS[profiler]}" fill-opacity="0.35"/>')
                q1, median, q3 = (s[field] for field in ("recovery_ratio_q1", "recovery_ratio_median", "recovery_ratio_q3"))
                parts.append(f'<line x1="{x}" y1="{recovery_y(q1):.1f}" x2="{x}" y2="{recovery_y(q3):.1f}" '
                             f'stroke="{COLORS[profiler]}" stroke-width="5"/>')
                parts.append(f'<line x1="{x-10}" y1="{recovery_y(median):.1f}" x2="{x+10}" y2="{recovery_y(median):.1f}" '
                             f'stroke="#17242e" stroke-width="3"/>')
    for k, dose in enumerate(DOSES):
        for i, profiler in enumerate(PROFILERS):
            for j, condition in enumerate(CONDITIONS):
                x = 215 + k * 325 + i * 135 + j * 55
                svg_text(parts, x, 974,
                         ("K" if i == 0 else "M") + (" Ctrl" if j == 0 else " CRC"),
                         11, anchor="middle")
        svg_text(parts, 315 + k * 325, 1000, f'{dose*100:.2f}% implanted', 16, "bold", anchor="middle")
    svg_text(parts, 1200, 715, "Key", 16, "bold")
    for i, profiler in enumerate(PROFILERS):
        svg_text(parts, 1200, 745 + i * 30,
                 "● " + ("Kraken2 + Bracken" if i == 0 else "MetaPhlAn 4"), 15, color=COLORS[profiler])
    svg_text(parts, 1200, 825, "Each dose: Control then CRC", 14)
    svg_text(parts, 1200, 848, f'Points beyond plotted y range: {clipped}', 13)
    svg_text(parts, 1200, 880, "At 0.01% (paired subset):", 15, "bold")
    for i, profiler in enumerate(PROFILERS):
        for j, condition in enumerate(CONDITIONS):
            s = next(r for r in spikes if r["profiler"] == profiler and
                     r["condition"] == condition and r["dose_fraction"] == DOSES[0])
            samples = [r for r in spike_samples if r["profiler"] == profiler and
                       r["condition"] == condition and r["dose_fraction"] == DOSES[0]]
            under = sum(r["recovery_ratio"] < .5 for r in samples)
            svg_text(parts, 1200, 904 + (i * 2 + j) * 23,
                     f'{"K+B" if i == 0 else "Meta"} {condition}: detected '
                     f'{s["observed_positive"]}/{s["n_paired"]}; '
                     f'<0.5 recovery {under}/{len(samples)}', 13, color=COLORS[profiler])
    svg_text(parts, 42, 1030, "Caution: full-cohort baselines and the paired spike subset have different n. Abundance fractions are profiler-specific estimates.", 15)
    svg_text(parts, 42, 1057, "Good spike recovery does not prove a biological CRC association; a non-significant q does not prove biological absence.", 15)
    svg_text(parts, 42, 1085, "Spike reference: Kraken read-proportional; MetaPhlAn genome-equivalent. Original assembly; no pseudocounts.", 15)
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def build(calls, endpoints, manifest, abundance, alias_path, outdir):
    if outdir.exists():
        raise ValueError("output exists: " + str(outdir))
    alias_map = aliases(alias_path)
    models = model_rows(calls, alias_map)
    spikes, spike_samples = summarize(endpoint_rows(endpoints))
    baselines = full_baselines(manifest, abundance, alias_map)
    full_summary = summary_rows(baselines, spikes, spike_samples)
    outdir.mkdir(parents=True)
    write_tsv(outdir / "full_baseline_samples.tsv", BASELINE_FIELDS, baselines)
    write_tsv(outdir / "full_baseline_and_low_dose_summary.tsv", SUMMARY_FIELDS, full_summary)
    write_tsv(outdir / "paired_spike_samples.tsv", tuple(spike_samples[0]), spike_samples)
    write_tsv(outdir / "paired_spike_summary.tsv", tuple(spikes[0]), spikes)
    write_tsv(outdir / "baseline_crc_model.tsv", tuple(models[0]), models)
    for cohort, target, canonical in CASES:
        draw(outdir / f"{cohort}_{target}_mechanism_atlas.svg", cohort, target, canonical,
             [r for r in models if r["cohort"] == cohort],
             [r for r in baselines if r["cohort"] == cohort],
             [r for r in full_summary if r["cohort"] == cohort],
             [r for r in spikes if r["cohort"] == cohort],
             [r for r in spike_samples if r["cohort"] == cohort])
    with (outdir / "source_sha256.tsv").open("w", encoding="utf-8") as handle:
        handle.write("source\tsha256\n")
        for path in (calls, endpoints, manifest, abundance, alias_path):
            handle.write(f"{path.resolve()}\t{sha256(path)}\n")
    (outdir / "DEVELOPMENT_ONLY.txt").write_text(
        "Provisional full-baseline and paired-spike descriptive diagnostic; not a causal test.\n",
        encoding="utf-8")
    (outdir / "SUCCESS").write_text("status\tPASS\n", encoding="utf-8")
    print("[PASS] two taxon-specific baseline and low-dose recovery atlases: " + str(outdir))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calls", required=True, type=Path)
    parser.add_argument("--endpoints", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--abundance", required=True, type=Path)
    parser.add_argument("--aliases", default=Path("examples/spike_taxon_aliases.csv"), type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    for path in (args.calls, args.endpoints, args.manifest, args.abundance, args.aliases):
        if not path.is_file() or not path.stat().st_size:
            parser.error("missing or empty input: " + str(path))
    try:
        build(args.calls, args.endpoints, args.manifest, args.abundance, args.aliases, args.outdir)
    except (OSError, ValueError, KeyError, statistics.StatisticsError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
