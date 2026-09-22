#!/usr/bin/env python3
"""CRC-anchored Adenoma contrast and direct-spike evidence in three cohorts."""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import math
import statistics
from collections import defaultdict
from pathlib import Path

COHORTS = ("feng", "yachida", "zeller")
PROFILERS = ("kraken2_bracken", "metaphlan4")
CONDITIONS = ("Control", "Adenoma", "CRC")
CONTRASTS = ("CRC_vs_Control", "Adenoma_vs_Control")
DOSES = (.0001, .0005, .001)
CALL_FIELDS = {"cohort", "analysis_population", "assembly_arm", "profiler", "contrast",
               "include", "model_spec", "spike_fraction_target", "feature", "effect",
               "lower_95", "upper_95", "q_value", "n_samples"}
MANIFEST_FIELDS = {"cohort", "analysis_population", "assembly_arm", "dose_level",
                   "include", "sample_id", "condition", "profiler", "source_profile"}
ABUNDANCE_FIELDS = {"profiler", "source_profile", "feature", "abundance_fraction"}
ENDPOINT_FIELDS = {"cohort", "sample_id", "condition", "analysis_population",
                   "assembly_arm", "profiler", "target_label", "spike_fraction_target",
                   "reference_type", "baseline_abundance_fraction", "observed_abundance_fraction",
                   "implanted_signal_profiler_scale", "recovered_spike_signal_profiler_scale"}
CALL_OUT = ("profiler", "feature", "display_feature", "cohort", "contrast", "evaluable",
            "effect", "lower_95", "upper_95", "q_value", "n_samples", "significant",
            "same_direction_as_crc", "direct_spike_target")
BASE_OUT = ("profiler", "feature", "cohort", "condition", "n", "positive", "prevalence",
            "zero", "positive_abundance_q1_percent", "positive_abundance_median_percent",
            "positive_abundance_q3_percent", "all_sample_median_percent",
            "all_sample_iqr_percent")
SPIKE_OUT = ("profiler", "feature", "cohort", "condition", "dose_percent", "n",
             "baseline_positive", "observed_positive", "recovery_q1", "recovery_median",
             "recovery_q3", "fraction_below_half")


def number(raw, field):
    try:
        value = float(raw)
    except (TypeError, ValueError) as error:
        raise ValueError(f"invalid {field}: {raw!r}") from error
    if not math.isfinite(value):
        raise ValueError(f"nonfinite {field}")
    return value


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def quartiles(values):
    if not values:
        return None, None, None
    if len(values) == 1:
        return values[0], values[0], values[0]
    return statistics.quantiles(values, n=4, method="inclusive")


def panel_aliases(panel_path, alias_path):
    with panel_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"label", "taxon_name"} <= set(reader.fieldnames or ()):
            raise ValueError("panel needs label and taxon_name")
        names = {}
        for row in reader:
            if row["taxon_name"] in names:
                raise ValueError("duplicate panel taxon")
            names[row["taxon_name"]] = row["label"]
    mapping = {}
    with alias_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not {"canonical", "alias", "tool"} <= set(reader.fieldnames or ()):
            raise ValueError("aliases need canonical, alias, tool")
        for row in reader:
            if row["canonical"] not in names:
                continue
            key = row["tool"], row["alias"]
            value = names[row["canonical"]], row["canonical"]
            if key in mapping and mapping[key] != value:
                raise ValueError("ambiguous panel feature alias")
            mapping[key] = value
    return mapping


def read_calls(path, panel_map):
    fits = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not CALL_FIELDS <= set(reader.fieldnames or ()):
            raise ValueError("disease calls lack required model columns")
        for row in reader:
            if (row["cohort"] not in COHORTS or row["profiler"] not in PROFILERS or
                row["contrast"] not in CONTRASTS or row["analysis_population"] != "community" or
                row["assembly_arm"] != "original" or row["include"] != "1" or
                row["model_spec"] != "primary_age_sex" or
                number(row["spike_fraction_target"], "model dose") != 0):
                continue
            key = row["cohort"], row["profiler"], row["feature"], row["contrast"]
            value = tuple(number(row[field], field) for field in
                          ("effect", "lower_95", "upper_95", "q_value", "n_samples"))
            if not (value[1] <= value[2] and 0 <= value[3] <= 1 and value[4] >= 2 and
                    value[4].is_integer()):
                raise ValueError("invalid disease model fit for " + repr(key))
            if key in fits and fits[key] != value:
                raise ValueError("contradictory repeated disease model fit for " + repr(key))
            fits[key] = value
    shared = {}
    for profiler in PROFILERS:
        features = {f for c, p, f, contrast in fits if p == profiler and contrast == CONTRASTS[0]}
        candidates = []
        for feature in features:
            values = [fits.get((cohort, profiler, feature, CONTRASTS[0])) for cohort in COHORTS]
            if all(v is not None and v[3] <= .05 for v in values) and (
                all(v[0] > 0 for v in values) or all(v[0] < 0 for v in values)
            ):
                candidates.append(feature)
        if not candidates:
            raise ValueError("no three-cohort shared CRC candidates for " + profiler)
        shared[profiler] = sorted(candidates, key=lambda f: (max(
            fits[cohort, profiler, f, CONTRASTS[0]][3] for cohort in COHORTS), f))
    rows = []
    for profiler in PROFILERS:
        for feature in shared[profiler]:
            display = ("Dialister pneumosintes [Allisonella in Kraken]"
                       if profiler == "kraken2_bracken" and feature == "Allisonella pneumosintes"
                       else feature)
            direct = panel_map.get((profiler, feature))
            for cohort in COHORTS:
                crc = fits[cohort, profiler, feature, CONTRASTS[0]]
                for contrast in CONTRASTS:
                    value = fits.get((cohort, profiler, feature, contrast))
                    rows.append(dict(profiler=profiler, feature=feature, display_feature=display,
                                     cohort=cohort, contrast=contrast, evaluable=int(value is not None),
                                     effect="" if value is None else value[0],
                                     lower_95="" if value is None else value[1],
                                     upper_95="" if value is None else value[2],
                                     q_value="" if value is None else value[3],
                                     n_samples="" if value is None else int(value[4]),
                                     significant=int(value is not None and value[3] <= .05),
                                     same_direction_as_crc="" if value is None else int(
                                         value[0] * crc[0] > 0),
                                     direct_spike_target="" if direct is None else direct[0]))
    return shared, rows


def full_baselines(manifest_path, abundance_path, shared):
    profiles = {}
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not MANIFEST_FIELDS <= set(reader.fieldnames or ()):
            raise ValueError("native manifest lacks required baseline columns")
        for row in reader:
            if (row["cohort"] not in COHORTS or row["profiler"] not in PROFILERS or
                row["condition"] not in CONDITIONS or row["analysis_population"] != "community" or
                row["assembly_arm"] != "original" or row["dose_level"] != "baseline" or
                row["include"] != "1"):
                continue
            key = row["cohort"], row["condition"], row["sample_id"], row["profiler"]
            source = str(Path(row["source_profile"]).resolve())
            if key in profiles and profiles[key] != source:
                raise ValueError("conflicting physical baseline profile for " + repr(key))
            profiles[key] = source
    for cohort in COHORTS:
        for condition in CONDITIONS:
            sample_sets = [{k[2] for k in profiles if k[0] == cohort and k[1] == condition and k[3] == p}
                           for p in PROFILERS]
            if sample_sets[0] != sample_sets[1] or len(sample_sets[0]) < 2:
                raise ValueError("unmatched or insufficient full baseline samples in " + cohort + " " + condition)
    wanted = {(p, source, f) for (c, _, _, p), source in profiles.items() for f in shared[p]}
    found = {}
    with abundance_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not ABUNDANCE_FIELDS <= set(reader.fieldnames or ()):
            raise ValueError("native abundance table lacks required columns")
        for row in reader:
            key = row["profiler"], str(Path(row["source_profile"]).resolve()), row["feature"]
            if key not in wanted:
                continue
            value = number(row["abundance_fraction"], "native abundance")
            if not 0 <= value <= 1 or key in found:
                raise ValueError("invalid or duplicate native abundance for " + repr(key))
            found[key] = value
    raw, summary = [], []
    for (cohort, condition, sample, profiler), source in sorted(profiles.items()):
        for feature in shared[profiler]:
            raw.append(dict(profiler=profiler, feature=feature, cohort=cohort,
                            condition=condition, sample_id=sample,
                            abundance_fraction=found.get((profiler, source, feature), 0.0)))
    for profiler in PROFILERS:
        for feature in shared[profiler]:
            for cohort in COHORTS:
                for condition in CONDITIONS:
                    values = sorted(r["abundance_fraction"] for r in raw if
                                    r["profiler"] == profiler and r["feature"] == feature and
                                    r["cohort"] == cohort and r["condition"] == condition)
                    if not values:
                        raise ValueError("missing full baseline stratum")
                    positive = [v for v in values if v > 0]
                    p1, pm, p3 = quartiles(positive)
                    a1, am, a3 = quartiles(values)
                    summary.append(dict(profiler=profiler, feature=feature, cohort=cohort,
                                        condition=condition, n=len(values), positive=len(positive),
                                        prevalence=len(positive) / len(values), zero=len(values)-len(positive),
                                        positive_abundance_q1_percent="" if p1 is None else 100*p1,
                                        positive_abundance_median_percent="" if pm is None else 100*pm,
                                        positive_abundance_q3_percent="" if p3 is None else 100*p3,
                                        all_sample_median_percent=100*am,
                                        all_sample_iqr_percent=100*(a3-a1)))
    return raw, summary


def direct_spikes(path, shared, panel_map):
    selected = {(p, panel_map[p, f][0]): f for p in PROFILERS for f in shared[p]
                if (p, f) in panel_map}
    paired = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not ENDPOINT_FIELDS <= set(reader.fieldnames or ()):
            raise ValueError("paired endpoints lack required response columns")
        for row in reader:
            key = row["profiler"], row["target_label"]
            if (key not in selected or row["cohort"] not in COHORTS or
                row["condition"] not in CONDITIONS or row["analysis_population"] != "independent" or
                row["assembly_arm"] != "original"):
                continue
            dose = number(row["spike_fraction_target"], "endpoint dose")
            nearest = min(DOSES, key=lambda d: abs(dose-d)/d)
            if abs(dose-nearest)/nearest > .05:
                continue
            expected_ref = "genome_equivalent" if row["profiler"] == "metaphlan4" else "read_proportional"
            if row["reference_type"] != expected_ref:
                raise ValueError("wrong paired endpoint reference type")
            base = number(row["baseline_abundance_fraction"], "endpoint baseline")
            observed = number(row["observed_abundance_fraction"], "endpoint observed")
            implanted = number(row["implanted_signal_profiler_scale"], "implanted signal")
            recovered = number(row["recovered_spike_signal_profiler_scale"], "recovered signal")
            if not (0 <= base <= 1 and 0 <= observed <= 1 and implanted > 0):
                raise ValueError("invalid paired endpoint values")
            identity = row["cohort"], row["condition"], nearest, row["sample_id"], row["profiler"], selected[key]
            if identity in paired:
                raise ValueError("duplicate direct spike endpoint")
            paired[identity] = (base, observed, recovered/implanted)
    raw, summary = [], []
    for cohort in COHORTS:
        for profiler in PROFILERS:
            for feature in shared[profiler]:
                if (profiler, feature) not in panel_map:
                    continue
                for condition in CONDITIONS:
                    for dose in DOSES:
                        records = [(sample, values) for (c, cond, d, sample, p, f), values in paired.items()
                                   if (c, cond, d, p, f) == (cohort, condition, dose, profiler, feature)]
                        if not records:
                            continue
                        ratios = [v[2] for _, v in records]
                        q1, med, q3 = quartiles(ratios)
                        summary.append(dict(profiler=profiler, feature=feature, cohort=cohort,
                                            condition=condition, dose_percent=100*dose,
                                            n=len(records), baseline_positive=sum(v[0]>0 for _,v in records),
                                            observed_positive=sum(v[1]>0 for _,v in records),
                                            recovery_q1=q1, recovery_median=med, recovery_q3=q3,
                                            fraction_below_half=sum(v[2]<.5 for _,v in records)/len(records)))
                        for sample, values in records:
                            raw.append(dict(profiler=profiler, feature=feature, cohort=cohort,
                                            condition=condition, dose_percent=100*dose, sample_id=sample,
                                            baseline_abundance_fraction=values[0],
                                            observed_abundance_fraction=values[1], recovery_ratio=values[2]))
    return raw, summary


def txt(parts, x, y, value, size=15, color="#202d37", weight="normal", anchor=None):
    extra = f' text-anchor="{anchor}"' if anchor else ""
    parts.append(f'<text x="{x}" y="{y}" font-family="Arial,sans-serif" font-size="{size}" '
                 f'font-weight="{weight}" fill="{color}"{extra}>{html.escape(str(value))}</text>')


def draw_associations(path, profiler, features, calls):
    width, row_h = 1610, 98
    height = 215 + len(features)*row_h + 95
    title = "Kraken2 + Bracken" if profiler == "kraken2_bracken" else "MetaPhlAn 4"
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>']
    txt(parts, 30, 45, f'Do shared {title} CRC candidates appear in adenomas?', 29, weight="bold")
    txt(parts, 30, 76, "Rows selected by same-direction CRC q ≤ 0.05 in all three cohorts; adenoma results shown regardless of significance · DEVELOPMENT ONLY", 15)
    txt(parts, 30, 106, "Each cell: age/sex-adjusted effect [95% CI], BH q, model n. Colour describes Adenoma vs Control only.", 15)
    for i, cohort in enumerate(COHORTS):
        x = 510 + i*365
        txt(parts, x+169, 165, cohort.title(), 21, weight="bold", anchor="middle")
    for row_index, feature in enumerate(features):
        y = 188 + row_index*row_h
        display = next(r["display_feature"] for r in calls if r["profiler"]==profiler and r["feature"]==feature)
        txt(parts, 30, y+35, display if len(display)<54 else display[:51]+"...", 18)
        for i, cohort in enumerate(COHORTS):
            x = 510 + i*365
            crc = next(r for r in calls if (r["profiler"],r["feature"],r["cohort"],r["contrast"]) ==
                       (profiler,feature,cohort,"CRC_vs_Control"))
            aden = next(r for r in calls if (r["profiler"],r["feature"],r["cohort"],r["contrast"]) ==
                        (profiler,feature,cohort,"Adenoma_vs_Control"))
            if not aden["evaluable"]:
                fill, label = "#e5e9ed", "Adenoma: not evaluable"
            elif aden["significant"] and aden["same_direction_as_crc"]:
                fill, label = "#c6eadc", "Adenoma: same-direction call"
            elif aden["same_direction_as_crc"]:
                fill, label = "#f8e7bd", "Adenoma: same direction, not called"
            elif aden["effect"] == 0:
                fill, label = "#e5e9ed", "Adenoma: zero estimated effect"
            else:
                fill, label = "#f1cbd5", "Adenoma: opposite direction"
            parts.append(f'<rect x="{x}" y="{y}" width="345" height="80" rx="4" fill="{fill}"/>')
            txt(parts, x+9, y+19, label, 13, weight="bold")
            txt(parts, x+9, y+41, f'CRC {crc["effect"]:+.2f} [{crc["lower_95"]:+.2f},{crc["upper_95"]:+.2f}] q={crc["q_value"]:.2g} n={crc["n_samples"]}', 13)
            if aden["evaluable"]:
                txt(parts, x+9, y+64,
                    f'Aden {aden["effect"]:+.2f} [{aden["lower_95"]:+.2f},{aden["upper_95"]:+.2f}] q={aden["q_value"]:.2g} n={aden["n_samples"]}', 13)
    txt(parts, 30, height-65, "Adenoma q > 0.05 is not evidence of absence. Compare effect size and CI before invoking analytical detection limits.", 15)
    txt(parts, 30, height-38, "CRC anchors are provisional partial-snapshot results; directly spiked species are examined separately for recovery.", 15)
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def draw_measurement(path, profiler, features, calls, baselines, spikes):
    width, row_h = 1710, 196
    height = 220 + len(features)*row_h + 95
    title = "Kraken2 + Bracken" if profiler == "kraken2_bracken" else "MetaPhlAn 4"
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>']
    txt(parts, 28, 45, f'{title}: are shared CRC markers measurable in adenomas?', 28, weight="bold")
    txt(parts, 28, 75, "Full unspiked baselines: detected n/N and positive-only median [Q1,Q3] abundance (%). Direct 0.01% spikes: recovery median [Q1,Q3].", 14)
    txt(parts, 28, 100, "A blank direct-spike result means this species was not implanted individually; it is not evidence of good recovery. DEVELOPMENT ONLY", 14)
    for i, cohort in enumerate(COHORTS):
        txt(parts, 630+i*355, 167, cohort.title(), 20, weight="bold", anchor="middle")
    for ri, feature in enumerate(features):
        y = 188+ri*row_h
        display = next(r["display_feature"] for r in calls if r["profiler"]==profiler and r["feature"]==feature)
        txt(parts, 28, y+28, display if len(display)<49 else display[:46]+"...", 17)
        direct = next(r["direct_spike_target"] for r in calls if r["profiler"]==profiler and r["feature"]==feature)
        txt(parts, 28, y+53, "direct spike: " + (direct if direct else "not in panel"), 13, color="#5a6975")
        for ci, cohort in enumerate(COHORTS):
            x = 455+ci*355
            parts.append(f'<rect x="{x}" y="{y}" width="340" height="180" rx="4" fill="#f2f6f8"/>')
            for j, condition in enumerate(CONDITIONS):
                b = next(r for r in baselines if (r["profiler"],r["feature"],r["cohort"],r["condition"]) ==
                         (profiler,feature,cohort,condition))
                med = b["positive_abundance_median_percent"]
                med_text = ("no positive reports" if med == "" else
                            f'pos {med:.2g} [{b["positive_abundance_q1_percent"]:.2g},'
                            f'{b["positive_abundance_q3_percent"]:.2g}]%')
                txt(parts, x+8, y+20+j*22,
                    f'{condition[:3]}: {b["positive"]}/{b["n"]}; {med_text}', 12)
            parts.append(f'<line x1="{x+8}" y1="{y+77}" x2="{x+330}" y2="{y+77}" stroke="#cbd7de"/>')
            for j, condition in enumerate(CONDITIONS):
                low = next((r for r in spikes if (r["profiler"],r["feature"],r["cohort"],r["condition"]) ==
                            (profiler,feature,cohort,condition) and math.isclose(r["dose_percent"],.01)), None)
                if direct and low:
                    label = (f'{condition[:3]} spike n={low["n"]}, det={low["observed_positive"]}; '
                             f'rec {low["recovery_median"]:.2f} '
                             f'[{low["recovery_q1"]:.2f},{low["recovery_q3"]:.2f}]')
                elif direct:
                    label = f'{condition[:3]} spike: not available'
                else:
                    label = "Direct spike: not in panel" if j == 0 else ""
                txt(parts, x+8, y+98+j*22, label, 12, color="#0a6079")
            txt(parts, x+8, y+169, "Control → Adenoma → CRC", 11, color="#5a6975")
    txt(parts, 28, height-64, "Full baseline n differs from paired spike n. Positive-only abundance summaries exclude zeros; detection counts include them.", 14)
    txt(parts, 28, height-39, "Known-spike recovery can reveal an analytical floor but cannot prove whether an unspiked Adenoma association is biological.", 14)
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def build(calls, manifest, abundance, endpoints, panel, alias_path, outdir):
    if outdir.exists():
        raise ValueError("output exists: " + str(outdir))
    panel_map = panel_aliases(panel, alias_path)
    shared, call_rows = read_calls(calls, panel_map)
    baseline_samples, baseline_summary = full_baselines(manifest, abundance, shared)
    spike_samples, spike_summary = direct_spikes(endpoints, shared, panel_map)
    missing_direct = [(p, f) for p in PROFILERS for f in shared[p]
                      if (p, f) in panel_map and not any(
                          r["profiler"] == p and r["feature"] == f for r in spike_summary)]
    if missing_direct:
        raise ValueError("panel candidates have no direct-spike endpoints: " + repr(missing_direct))
    outdir.mkdir(parents=True)
    write_tsv(outdir / "shared_crc_adenoma_model.tsv", CALL_OUT, call_rows)
    write_tsv(outdir / "shared_crc_baseline_samples.tsv",
              ("profiler","feature","cohort","condition","sample_id","abundance_fraction"), baseline_samples)
    write_tsv(outdir / "shared_crc_baseline_summary.tsv", BASE_OUT, baseline_summary)
    write_tsv(outdir / "shared_crc_direct_spike_samples.tsv",
              ("profiler","feature","cohort","condition","dose_percent","sample_id",
               "baseline_abundance_fraction","observed_abundance_fraction","recovery_ratio"), spike_samples)
    write_tsv(outdir / "shared_crc_direct_spike_summary.tsv", SPIKE_OUT, spike_summary)
    overview = []
    call_lookup = {(r["profiler"], r["feature"], r["cohort"], r["contrast"]): r
                   for r in call_rows}
    for profiler in PROFILERS:
        features = shared[profiler]
        adenoma_all = sum(all(
            call_lookup[profiler, feature, cohort, "Adenoma_vs_Control"]["significant"] and
            call_lookup[profiler, feature, cohort, "Adenoma_vs_Control"]["same_direction_as_crc"]
            for cohort in COHORTS) for feature in features)
        overview.extend([
            dict(profiler=profiler, metric="shared_crc_candidates", value=len(features)),
            dict(profiler=profiler, metric="same_direction_adenoma_calls_all_three", value=adenoma_all),
            dict(profiler=profiler, metric="directly_spiked_shared_candidates",
                 value=sum((profiler, f) in panel_map for f in features)),
        ])
    write_tsv(outdir / "candidate_summary.tsv", ("profiler","metric","value"), overview)
    for profiler in PROFILERS:
        draw_associations(outdir / (profiler+"_crc_anchored_adenoma.svg"), profiler,
                          shared[profiler], call_rows)
        draw_measurement(outdir / (profiler+"_adenoma_measurement.svg"), profiler,
                         shared[profiler], call_rows, baseline_summary, spike_summary)
    with (outdir / "source_sha256.tsv").open("w", encoding="utf-8") as handle:
        handle.write("source\tsha256\n")
        for path in (calls, manifest, abundance, endpoints, panel, alias_path):
            handle.write(f"{path.resolve()}\t{sha256(path)}\n")
    (outdir / "DEVELOPMENT_ONLY.txt").write_text(
        "Provisional CRC-anchored Adenoma comparison and direct-spike descriptive evidence.\n",
        encoding="utf-8")
    (outdir / "SUCCESS").write_text("status\tPASS\n", encoding="utf-8")
    print(f"[PASS] shared CRC-to-Adenoma bridge: Kraken {len(shared[PROFILERS[0]])}; "
          f"MetaPhlAn {len(shared[PROFILERS[1]])}; {outdir}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("calls", "manifest", "abundance", "endpoints", "panel", "aliases", "outdir"):
        parser.add_argument("--"+flag, required=True, type=Path)
    args = parser.parse_args()
    inputs = (args.calls, args.manifest, args.abundance, args.endpoints, args.panel, args.aliases)
    for path in inputs:
        if not path.is_file() or not path.stat().st_size:
            parser.error("missing or empty input: " + str(path))
    try:
        build(*inputs, args.outdir)
    except (OSError, ValueError, KeyError, statistics.StatisticsError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
