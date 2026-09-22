#!/usr/bin/env python3
"""Paired-profiler diagnostic for two discordant baseline CRC candidates."""
import argparse
import csv
import hashlib
import html
import math
import statistics
from pathlib import Path

PROFILERS = ("kraken2_bracken", "metaphlan4")
CONDITIONS = ("Control", "CRC")
DOSES = (.0001, .0005, .001)
CASES = (("zeller", "Pmic", "Parvimonas micra"),
         ("yachida", "Dpne", "Dialister pneumosintes"))
CALL_FIELDS = {"cohort", "analysis_population", "assembly_arm", "profiler",
               "contrast", "include", "model_spec", "spike_fraction_target",
               "feature", "effect", "standard_error", "lower_95", "upper_95",
               "q_value", "n_samples"}
ENDPOINT_FIELDS = {"cohort", "sample_id", "condition", "analysis_population",
                   "assembly_arm", "profiler", "target_label", "spike_fraction_target",
                   "reference_type", "baseline_abundance_fraction",
                   "observed_abundance_fraction", "implanted_signal_profiler_scale",
                   "recovered_spike_signal_profiler_scale", "implanted_read_pairs_target"}


def number(value, label):
    try:
        result = float(value)
    except ValueError as error:
        raise ValueError("invalid {}: {}".format(label, value)) from error
    if not math.isfinite(result):
        raise ValueError("nonfinite {}".format(label))
    return result


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def aliases(path):
    found = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not {"canonical", "alias", "tool"} <= set(reader.fieldnames or ()):
            raise ValueError("alias map lacks canonical/alias/tool columns")
        for row in reader:
            key = row["canonical"], row["tool"]
            if key in found and found[key] != row["alias"]:
                raise ValueError("ambiguous alias for {}".format(key))
            found[key] = row["alias"]
    for _, _, canonical in CASES:
        for profiler in PROFILERS:
            if not found.get((canonical, profiler)):
                raise ValueError("missing alias for {} {}".format(canonical, profiler))
    if found["Dialister pneumosintes", "kraken2_bracken"] != "Allisonella pneumosintes":
        raise ValueError("Kraken Dialister/Allisonella alias contract changed")
    return found


def model_rows(path, alias_map):
    wanted = {(cohort, profiler, alias_map[canonical, profiler]): (label, canonical)
              for cohort, label, canonical in CASES for profiler in PROFILERS}
    found = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not CALL_FIELDS <= set(reader.fieldnames or ()):
            raise ValueError("disease calls lack model-uncertainty columns")
        for row in reader:
            key = row["cohort"], row["profiler"], row["feature"]
            if key not in wanted or row["analysis_population"] != "community" or row["assembly_arm"] != "original" or row["contrast"] != "CRC_vs_Control" or row["include"] != "1" or row["model_spec"] != "primary_age_sex":
                continue
            if number(row["spike_fraction_target"], "model dose") != 0:
                continue
            values = {field: number(row[field], field) for field in
                      ("effect", "standard_error", "lower_95", "upper_95", "q_value", "n_samples")}
            if not 0 <= values["q_value"] <= 1 or values["standard_error"] < 0 or values["n_samples"] <= 0 or values["lower_95"] > values["upper_95"]:
                raise ValueError("invalid baseline model uncertainty for {}".format(key))
            if key in found and found[key] != values:
                raise ValueError("contradictory repeated baseline model for {}".format(key))
            found[key] = values
    if set(found) != set(wanted):
        raise ValueError("missing baseline model fit(s): {}".format(sorted(set(wanted) - set(found))))
    rows = []
    for cohort, target, canonical in CASES:
        for profiler in PROFILERS:
            feature = alias_map[canonical, profiler]
            values = found[cohort, profiler, feature]
            rows.append(dict(cohort=cohort, target_label=target, canonical_taxon=canonical,
                             profiler=profiler, reported_feature=feature,
                             effect=values["effect"], standard_error=values["standard_error"],
                             lower_95=values["lower_95"], upper_95=values["upper_95"],
                             q_value=values["q_value"], n_samples=int(values["n_samples"])))
    return rows


def dose_name(value):
    dose = min(DOSES, key=lambda d: abs(value - d) / d)
    if abs(value - dose) / dose > .05:
        return None
    return dose


def endpoint_rows(path):
    selected = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not ENDPOINT_FIELDS <= set(reader.fieldnames or ()):
            raise ValueError("paired endpoints lack required fields")
        for row in reader:
            if (row["cohort"], row["target_label"]) not in {(c, t) for c, t, _ in CASES} or row["analysis_population"] != "independent" or row["assembly_arm"] != "original" or row["condition"] not in CONDITIONS or row["profiler"] not in PROFILERS:
                continue
            dose = dose_name(number(row["spike_fraction_target"], "spike dose"))
            if dose is None:
                continue
            expected_ref = "genome_equivalent" if row["profiler"] == "metaphlan4" else "read_proportional"
            if row["reference_type"] != expected_ref:
                raise ValueError("wrong reference type in paired endpoint")
            baseline = number(row["baseline_abundance_fraction"], "baseline abundance")
            observed = number(row["observed_abundance_fraction"], "observed abundance")
            implanted = number(row["implanted_signal_profiler_scale"], "implanted signal")
            recovered = number(row["recovered_spike_signal_profiler_scale"], "recovered signal")
            reads = number(row["implanted_read_pairs_target"], "implanted read pairs")
            if not (0 <= baseline <= 1 and 0 <= observed <= 1 and implanted > 0 and reads >= 0):
                raise ValueError("invalid paired endpoint abundance or implanted reads")
            key = row["cohort"], row["target_label"], row["condition"], dose, row["sample_id"], row["profiler"]
            if key in selected:
                raise ValueError("duplicate endpoint {}".format(key))
            selected[key] = (baseline, observed, implanted, recovered, reads)
    return selected


def quartiles(values):
    if len(values) == 1:
        return values[0], values[0], values[0]
    q1, q2, q3 = statistics.quantiles(values, n=4, method="inclusive")
    return q1, q2, q3


def summarize(endpoints):
    summaries, sample_rows = [], []
    for cohort, target, canonical in CASES:
        for condition in CONDITIONS:
            for dose in DOSES:
                samples = [{key[4] for key in endpoints if key[:4] == (cohort, target, condition, dose) and key[5] == profiler}
                           for profiler in PROFILERS]
                common = sorted(samples[0] & samples[1])
                if len(common) < 2:
                    raise ValueError("fewer than two paired profiler samples for {} {} {}".format(cohort, condition, dose))
                for profiler, available in zip(PROFILERS, samples):
                    values = [endpoints[cohort, target, condition, dose, sample, profiler]
                              for sample in common]
                    baseline = [v[0] for v in values]
                    observed = [v[1] for v in values]
                    ratios = [v[3] / v[2] for v in values]
                    errors = [abs(v[3] - v[2]) / v[2] for v in values]
                    q1, median, q3 = quartiles(ratios)
                    summaries.append(dict(cohort=cohort, target_label=target,
                        canonical_taxon=canonical, condition=condition, dose_fraction=dose,
                        profiler=profiler, n_paired=len(common), n_profiler_available=len(available),
                        n_unpaired_excluded=len(available - set(common)),
                        baseline_positive=sum(v > 0 for v in baseline),
                        baseline_prevalence=sum(v > 0 for v in baseline) / len(common),
                        baseline_median=statistics.median(baseline),
                        observed_positive=sum(v > 0 for v in observed),
                        observed_prevalence=sum(v > 0 for v in observed) / len(common),
                        recovery_ratio_q1=q1, recovery_ratio_median=median,
                        recovery_ratio_q3=q3, median_absolute_relative_error=statistics.median(errors),
                        median_implanted_read_pairs=statistics.median(v[4] for v in values)))
                    for sample, value in zip(common, values):
                        sample_rows.append(dict(cohort=cohort, target_label=target,
                            canonical_taxon=canonical, condition=condition, dose_fraction=dose,
                            sample_id=sample, profiler=profiler,
                            baseline_abundance_fraction=value[0],
                            observed_abundance_fraction=value[1],
                            implanted_signal_profiler_scale=value[2],
                            recovered_spike_signal_profiler_scale=value[3],
                            recovery_ratio=value[3] / value[2],
                            implanted_read_pairs_target=value[4]))
    # A paired sample's unspiked abundance cannot change with spike dose.
    baselines = {}
    for row in sample_rows:
        key = row["cohort"], row["target_label"], row["condition"], row["sample_id"], row["profiler"]
        value = row["baseline_abundance_fraction"]
        if key in baselines and not math.isclose(baselines[key], value, abs_tol=1e-10, rel_tol=1e-9):
            raise ValueError("baseline abundance varies across doses for {}".format(key))
        baselines[key] = value
    return summaries, sample_rows


def draw(path, model, summary, cohort, target, canonical):
    width, height = 1550, 680
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}">'.format(width, height, width, height),
             '<rect width="100%" height="100%" fill="white"/>',
             '<text x="35" y="48" font-family="sans-serif" font-size="31" font-weight="bold">Can known spikes clarify the {} CRC call in {}?</text>'.format(html.escape(canonical), cohort.title()),
             '<text x="35" y="80" font-family="sans-serif" font-size="17">Original community CRC model above; paired independent single-species spikes below · DEVELOPMENT ONLY</text>',
             '<text x="35" y="114" font-family="sans-serif" font-size="17">Kraken feature: {} · MetaPhlAn feature: {}</text>'.format(
                 html.escape(next(r["reported_feature"] for r in model if r["profiler"] == PROFILERS[0])),
                 html.escape(next(r["reported_feature"] for r in model if r["profiler"] == PROFILERS[1])))]
    colors = {"kraken2_bracken": "#d55e00", "metaphlan4": "#0072b2"}
    for i, profiler in enumerate(PROFILERS):
        m = next(r for r in model if r["profiler"] == profiler)
        x = 35 + i * 750
        name = "Kraken2 + Bracken" if i == 0 else "MetaPhlAn 4"
        parts.append('<rect x="{}" y="139" width="705" height="112" rx="6" fill="#f3f5f7"/>'.format(x))
        parts.append('<text x="{}" y="171" font-family="sans-serif" font-size="21" font-weight="bold" fill="{}">{}</text>'.format(x + 16, colors[profiler], name))
        parts.append('<text x="{}" y="199" font-family="sans-serif" font-size="18">CRC effect {:+.2f} [95% CI {:+.2f}, {:+.2f}] · BH q={:.3g} · model n={}</text>'.format(
            x + 16, m["effect"], m["lower_95"], m["upper_95"], m["q_value"], m["n_samples"]))
        parts.append('<text x="{}" y="228" font-family="sans-serif" font-size="15">Original unspiked abundance, not the spike-in recovery estimate</text>'.format(x + 16))
    parts.append('<text x="35" y="297" font-family="sans-serif" font-size="22" font-weight="bold">Known low-dose target additions: paired-sample response</text>')
    for j, condition in enumerate(CONDITIONS):
        y = 326 + j * 150
        parts.append('<text x="35" y="{}" font-family="sans-serif" font-size="19" font-weight="bold">{}</text>'.format(y + 40, condition))
        for k, dose in enumerate(DOSES):
            x = 195 + k * 440
            parts.append('<text x="{}" y="{}" font-family="sans-serif" font-size="18" font-weight="bold">{:.2f}% implanted</text>'.format(x, y, 100 * dose))
            for pi, profiler in enumerate(PROFILERS):
                r = next(r for r in summary if r["condition"] == condition and r["dose_fraction"] == dose and r["profiler"] == profiler)
                yy = y + 26 + pi * 48
                label = "K+B" if pi == 0 else "MetaPhlAn"
                parts.append('<rect x="{}" y="{}" width="422" height="44" rx="3" fill="#f3f5f7"/>'.format(x, yy))
                parts.append('<text x="{}" y="{}" font-family="sans-serif" font-size="14" fill="{}">{}: n={} · baseline {:.0f}% · observed {:.0f}% positive</text>'.format(
                    x + 6, yy + 17, colors[profiler], label, r["n_paired"],
                    100 * r["baseline_prevalence"], 100 * r["observed_prevalence"]))
                parts.append('<text x="{}" y="{}" font-family="sans-serif" font-size="14" fill="{}">Recovery median {:.2f} [Q1 {:.2f}, Q3 {:.2f}]</text>'.format(
                    x + 6, yy + 35, colors[profiler], r["recovery_ratio_median"],
                    r["recovery_ratio_q1"], r["recovery_ratio_q3"]))
    parts.append('<text x="35" y="635" font-family="sans-serif" font-size="15">Recovery = (observed - diluted baseline) / implanted signal on each profiler’s scale; 1 is ideal. Brackets are sample Q1–Q3, not CIs.</text>')
    parts.append('<text x="35" y="659" font-family="sans-serif" font-size="15">A recovery deficit supports measurement sensitivity as a contributor; it does not prove why an unspiked CRC association differs across cohorts.</text>')
    parts.append('</svg>\n')
    path.write_text("".join(parts), encoding="utf-8")


def build(calls, endpoints, alias_path, outdir):
    if outdir.exists():
        raise ValueError("output exists: {}".format(outdir))
    alias_map = aliases(alias_path)
    model = model_rows(calls, alias_map)
    summary, samples = summarize(endpoint_rows(endpoints))
    outdir.mkdir(parents=True)
    write_tsv(outdir / "baseline_crc_model.tsv", model)
    write_tsv(outdir / "paired_recovery_summary.tsv", summary)
    write_tsv(outdir / "paired_recovery_samples.tsv", samples)
    for cohort, target, canonical in CASES:
        draw(outdir / (cohort + "_" + target + "_diagnostic.svg"),
             [r for r in model if r["cohort"] == cohort and r["target_label"] == target],
             [r for r in summary if r["cohort"] == cohort and r["target_label"] == target],
             cohort, target, canonical)
    (outdir / "source.sha256").write_text("".join("{}  {}\n".format(digest(p), p.resolve())
                                                 for p in (calls, endpoints, alias_path)), encoding="utf-8")
    (outdir / "DEVELOPMENT_ONLY.txt").write_text("status\tDEVELOPMENT_ONLY\n", encoding="utf-8")
    (outdir / "SUCCESS").write_text("status\tPASS\n", encoding="utf-8")
    print("[PASS] paired-profiler low-dose recovery diagnostic for Zeller Pmic and Yachida Dpne")
    print(outdir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calls", type=Path, required=True)
    parser.add_argument("--endpoints", type=Path, required=True)
    parser.add_argument("--aliases", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    try:
        build(args.calls, args.endpoints, args.aliases, args.outdir)
    except (OSError, ValueError, KeyError, statistics.StatisticsError) as error:
        raise SystemExit("[ERROR] {}".format(error)) from error
