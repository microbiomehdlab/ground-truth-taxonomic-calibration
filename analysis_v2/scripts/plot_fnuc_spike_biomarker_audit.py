#!/usr/bin/env python3
"""Audit F. nucleatum disease calls and related-taxon response after a shared community spike."""
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
CONTRASTS = ("CRC_vs_Control", "Adenoma_vs_Control")
CONDITIONS = ("Control", "Adenoma", "CRC")
TARGET = "Fusobacterium nucleatum"
MEMBER_DOSE = 1e-5  # fraction: 0.001% per member; ten-member mixture totals about 0.01%
LEDGER_REQUIRED = {"cohort", "analysis_population", "assembly_arm", "profiler", "contrast",
                   "target_label", "spike_fraction_target", "q_threshold", "feature", "feature_role",
                   "baseline_called", "dose_called", "baseline_effect", "dose_effect",
                   "baseline_q_value", "dose_q_value", "transition"}
ENDPOINT_REQUIRED = {"cohort", "analysis_population", "assembly_arm", "profiler", "condition",
                     "target_label", "spike_fraction_target", "sample_id",
                     "baseline_abundance_fraction", "observed_abundance_fraction",
                     "implanted_signal_profiler_scale", "recovered_spike_signal_profiler_scale",
                     "reference_type"}
CALL_FIELDS = ("cohort", "profiler", "contrast", "feature", "feature_role",
               "baseline_effect", "dose_effect", "baseline_q_value", "dose_q_value",
               "baseline_called", "dose_called", "transition", "category")
SUMMARY_FIELDS = ("cohort", "profiler", "contrast", "target_baseline_effect",
                  "target_dose_effect", "target_baseline_q", "target_dose_q",
                  "target_baseline_called", "target_dose_called", "related_gained",
                  "related_lost", "related_retained", "other_bystander_gained",
                  "other_bystander_lost", "other_bystander_retained")
RESPONSE_FIELDS = ("cohort", "profiler", "condition", "n", "baseline_zero",
                   "zero_rescued", "post_positive", "recovery_q1", "recovery_median",
                   "recovery_q3", "below_half")
BYSTANDER_FIELDS = ("cohort", "profiler", "condition", "feature", "n",
                    "response_q1", "response_median", "response_q3",
                    "positive_response_fraction")


def rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("empty table: " + str(path))
        yield set(reader.fieldnames)
        yield from reader


def write(path, fields, records):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def quartiles(values):
    if not values:
        return None, None, None
    if len(values) == 1:
        return (values[0],) * 3
    return tuple(statistics.quantiles(values, n=4, method="inclusive"))


def parse_ledger(path, member_dose=MEMBER_DOSE):
    source = rows(path)
    if not LEDGER_REQUIRED <= next(source):
        raise ValueError("transition ledger lacks required fields")
    selected = {}
    for r in source:
        if (r["analysis_population"] != "community" or r["assembly_arm"] != "original" or
            r["target_label"] != "CRCpanel" or r["cohort"] not in COHORTS or
            r["profiler"] not in PROFILERS or r["contrast"] not in CONTRASTS or
            not math.isclose(float(r["spike_fraction_target"]), member_dose, rel_tol=.05) or
            not math.isclose(float(r["q_threshold"]), .05, abs_tol=1e-10)):
            continue
        category = ("target" if r["feature"] == TARGET else
                    "related_Fusobacterium" if r["feature"].startswith("Fusobacterium ") else
                    "other_bystander" if r["feature_role"] == "bystander" else
                    "other_implanted_target")
        record = {field: r[field] for field in CALL_FIELDS if field != "category"}
        record["category"] = category
        key = r["cohort"], r["profiler"], r["contrast"], r["feature"]
        if key in selected:
            raise ValueError("duplicate transition " + repr(key))
        for field in ("baseline_effect", "dose_effect", "baseline_q_value", "dose_q_value"):
            if not math.isfinite(float(record[field])):
                raise ValueError("nonfinite transition: " + repr(key))
        selected[key] = record
    detail = list(selected.values())
    summary = []
    for cohort in COHORTS:
        for profiler in PROFILERS:
            for contrast in CONTRASTS:
                key = cohort, profiler, contrast, TARGET
                if key not in selected:
                    raise ValueError("missing F. nucleatum target transition " + repr(key))
                target = selected[key]
                peers = [r for r in detail if (r["cohort"], r["profiler"], r["contrast"]) ==
                         (cohort, profiler, contrast)]
                counts = {}
                for category in ("related_Fusobacterium", "other_bystander"):
                    group = [r for r in peers if r["category"] == category]
                    counts[category + "_gained"] = sum(r["transition"] == "gained" for r in group)
                    counts[category + "_lost"] = sum(r["transition"] == "lost" for r in group)
                    counts[category + "_retained"] = sum(r["transition"].startswith("retained") for r in group)
                summary.append(dict(cohort=cohort, profiler=profiler, contrast=contrast,
                    target_baseline_effect=target["baseline_effect"],
                    target_dose_effect=target["dose_effect"],
                    target_baseline_q=target["baseline_q_value"],
                    target_dose_q=target["dose_q_value"],
                    target_baseline_called=target["baseline_called"],
                    target_dose_called=target["dose_called"],
                    related_gained=counts["related_Fusobacterium_gained"],
                    related_lost=counts["related_Fusobacterium_lost"],
                    related_retained=counts["related_Fusobacterium_retained"],
                    other_bystander_gained=counts["other_bystander_gained"],
                    other_bystander_lost=counts["other_bystander_lost"],
                    other_bystander_retained=counts["other_bystander_retained"]))
    return detail, summary


def parse_endpoints(path, member_dose=MEMBER_DOSE):
    source = rows(path)
    if not ENDPOINT_REQUIRED <= next(source):
        raise ValueError("paired endpoints lack required fields")
    groups = defaultdict(list)
    for r in source:
        if (r["analysis_population"] != "community" or r["assembly_arm"] != "original" or
            r["target_label"] != "Fnuc" or r["cohort"] not in COHORTS or
            r["profiler"] not in PROFILERS or r["condition"] not in CONDITIONS or
            not math.isclose(float(r["spike_fraction_target"]), member_dose, rel_tol=.05)):
            continue
        expected = "genome_equivalent" if r["profiler"] == "metaphlan4" else "read_proportional"
        if r["reference_type"] != expected:
            raise ValueError("wrong profiler-scale spike reference")
        signal = float(r["implanted_signal_profiler_scale"])
        if signal <= 0:
            raise ValueError("nonpositive implanted signal")
        groups[r["cohort"], r["profiler"], r["condition"]].append((
            r["sample_id"], float(r["baseline_abundance_fraction"]),
            float(r["observed_abundance_fraction"]),
            float(r["recovered_spike_signal_profiler_scale"]) / signal))
    result = []
    for cohort in COHORTS:
        for profiler in PROFILERS:
            for condition in CONDITIONS:
                key = cohort, profiler, condition
                values = groups[key]
                if not values or len({v[0] for v in values}) != len(values):
                    raise ValueError("missing/duplicate paired target samples " + repr(key))
                ratio = [v[3] for v in values]
                q1, med, q3 = quartiles(ratio)
                result.append(dict(cohort=cohort, profiler=profiler, condition=condition,
                    n=len(values), baseline_zero=sum(v[1] == 0 for v in values),
                    zero_rescued=sum(v[1] == 0 and v[2] > 0 for v in values),
                    post_positive=sum(v[2] > 0 for v in values),
                    recovery_q1=q1, recovery_median=med, recovery_q3=q3,
                    below_half=sum(v < .5 for v in ratio)))
    return result


def parse_bystanders(path, dose_index=1):
    try:
        import duckdb
    except ImportError as error:
        raise ValueError("DuckDB is needed for --paired-features; use the analysis Apptainer image") from error
    con = duckdb.connect()
    try:
        columns = {r[0] for r in con.execute("SELECT * FROM read_parquet(?) LIMIT 0", [str(path)]).description}
        required = {"cohort", "profiler", "condition", "analysis_population", "dose_index", "feature",
                    "is_direct_target", "response_delta"}
        if not required <= columns:
            raise ValueError("paired-feature parquet lacks condition-specific response fields")
        query = """
          SELECT cohort, profiler, condition, feature,
                 count(*) n, quantile_cont(response_delta, 0.25) q1,
                 quantile_cont(response_delta, 0.5) med,
                 quantile_cont(response_delta, 0.75) q3,
                 avg(CASE WHEN response_delta > 0 THEN 1.0 ELSE 0.0 END) positive_fraction
          FROM read_parquet(?)
          WHERE analysis_population = 'community' AND dose_index = ?
            AND NOT is_direct_target AND feature LIKE 'Fusobacterium %'
          GROUP BY cohort, profiler, condition, feature
        """
        raw = con.execute(query, [str(path), dose_index]).fetchall()
    finally:
        con.close()
    result = []
    for cohort, profiler, condition, feature, n, q1, med, q3, positive in raw:
        if cohort not in COHORTS or profiler not in PROFILERS or condition not in CONDITIONS:
            continue
        result.append(dict(cohort=cohort, profiler=profiler, condition=condition,
                           feature=feature, n=n, response_q1=q1, response_median=med,
                           response_q3=q3, positive_response_fraction=positive))
    return sorted(result, key=lambda r: (r["cohort"], r["profiler"], r["condition"], r["feature"]))


def txt(parts, x, y, value, size=15, color="#22313e", weight="normal", anchor=None):
    a = f' text-anchor="{anchor}"' if anchor else ""
    parts.append(f'<text x="{x}" y="{y}" font-family="Arial,sans-serif" font-size="{size}" '
                 f'font-weight="{weight}" fill="{color}"{a}>{html.escape(str(value))}</text>')


def draw(path, summary, response, bystanders, member_dose_percent=.001):
    width, height = 1900, 1100
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>']
    txt(parts, 28, 46, "Does equal F. nucleatum spiking change disease-biomarker calls?", 29, weight="bold")
    txt(parts, 28, 76, f"{member_dose_percent:g}% per member in a ten-species community mixture, added across clinical groups · matched baseline versus spiked model · DEVELOPMENT ONLY", 15)
    txt(parts, 28, 101, "Open dot = baseline effect; filled dot = post-spike effect. A gained call is not a proven biological false positive; q is BH-adjusted.", 14)
    lookup = {(r["cohort"], r["profiler"], r["contrast"]): r for r in summary}
    rec = {(r["cohort"], r["profiler"], r["condition"]): r for r in response}
    all_effects = [float(r[field]) for r in summary
                   for field in ("target_baseline_effect", "target_dose_effect")]
    effect_min = min(-.5, min(all_effects)-.5)
    effect_max = max(.5, max(all_effects)+.5)
    for ci, cohort in enumerate(COHORTS):
        txt(parts, 451+ci*560, 153, cohort.title(), 22, weight="bold", anchor="middle")
    for pi, profiler in enumerate(PROFILERS):
        top = 174 + pi*420
        name = "Kraken2 + Bracken" if pi == 0 else "MetaPhlAn 4"
        color = "#bd5600" if pi == 0 else "#006e9e"
        txt(parts, 25, top+27, name, 17, color=color, weight="bold")
        for ci, cohort in enumerate(COHORTS):
            x = 184 + ci*560
            parts.append(f'<rect x="{x}" y="{top}" width="535" height="395" rx="5" fill="#f3f7f9"/>')
            for ri, (contrast, label) in enumerate((("CRC_vs_Control", "CRC vs Control"),
                                                     ("Adenoma_vs_Control", "Adenoma vs Control"))):
                row = lookup[cohort, profiler, contrast]
                y = top + 30 + ri*91
                txt(parts, x+12, y, label, 16, weight="bold")
                old = "called" if row["target_baseline_called"] == "1" else "not called"
                new = "called" if row["target_dose_called"] == "1" else "not called"
                txt(parts, x+12, y+22,
                    f'F. nucleatum: {old} → {new}; effect {float(row["target_baseline_effect"]):+.2f} → '
                    f'{float(row["target_dose_effect"]):+.2f}', 14, color=color)
                txt(parts, x+12, y+42,
                    f'q {float(row["target_baseline_q"]):.2g} → {float(row["target_dose_q"]):.2g}', 13)
                txt(parts, x+12, y+62,
                    f'Related Fuso calls: +{row["related_gained"]} / −{row["related_lost"]}; '
                    f'other bystander calls: +{row["other_bystander_gained"]} / −{row["other_bystander_lost"]}', 13)
                def effect_x(value):
                    return x+316 + 190*(value-effect_min)/(effect_max-effect_min)
                old_x = effect_x(float(row["target_baseline_effect"]))
                new_x = effect_x(float(row["target_dose_effect"]))
                axis_y = y+77
                parts.append(f'<line x1="{x+316}" y1="{axis_y}" x2="{x+506}" y2="{axis_y}" stroke="#adbac3"/>')
                zero_x = effect_x(0)
                parts.append(f'<line x1="{zero_x:.1f}" y1="{axis_y-8}" x2="{zero_x:.1f}" y2="{axis_y+8}" stroke="#7d8992"/>')
                parts.append(f'<line x1="{old_x:.1f}" y1="{axis_y}" x2="{new_x:.1f}" y2="{axis_y}" stroke="{color}" stroke-width="2"/>')
                parts.append(f'<circle cx="{old_x:.1f}" cy="{axis_y}" r="4" fill="white" stroke="{color}" stroke-width="2"/>')
                parts.append(f'<circle cx="{new_x:.1f}" cy="{axis_y}" r="4" fill="{color}"/>')
            parts.append(f'<line x1="{x+12}" y1="{top+213}" x2="{x+523}" y2="{top+213}" stroke="#cbd9de"/>')
            txt(parts, x+12, top+237, "Same known addition, different clinical backgrounds:", 14, weight="bold")
            for gi, condition in enumerate(CONDITIONS):
                r = rec[cohort, profiler, condition]
                txt(parts, x+12, top+260+gi*26,
                    f'{condition}: zero→positive {r["zero_rescued"]}/{r["baseline_zero"]}; '
                    f'recovery {float(r["recovery_median"]):.2f} '
                    f'[{float(r["recovery_q1"]):.2f},{float(r["recovery_q3"]):.2f}]', 13)
            if bystanders:
                relevant = [r for r in bystanders if (r["cohort"], r["profiler"]) == (cohort, profiler)]
                if relevant:
                    strongest = max(relevant, key=lambda r: abs(float(r["response_median"])))
                    label = strongest["feature"].replace("Fusobacterium ", "F. ")
                    feature_rows = {r["condition"]: r for r in relevant
                                    if r["feature"] == strongest["feature"]}
                    txt(parts, x+12, top+355, f'Largest related-taxon shift: {label}', 12)
                    if all(c in feature_rows for c in CONDITIONS):
                        values = " / ".join(
                            f'{c[:3]} {100*float(feature_rows[c]["response_median"]):+.2g} pp'
                            for c in CONDITIONS)
                        txt(parts, x+12, top+374, values, 12)
            else:
                txt(parts, x+12, top+355, "Related-species paired response not supplied; call changes shown above.", 12)
    txt(parts, 28, height-70,
        "Equal spiking cannot create a true CRC or adenoma association; post-spike calls test how the measurement/model responds to the known perturbation.", 14)
    txt(parts, 28, height-42,
        "The community mix also adds nine other taxa. Related-species calls and response require cautious attribution; baseline biological truth is not known.", 14)
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ledger", type=Path, required=True)
    p.add_argument("--endpoints", type=Path, required=True)
    p.add_argument("--paired-features", type=Path)
    p.add_argument("--member-dose-percent", type=float, choices=(.001, .1), default=.001,
                   help="community spike dose per species; 0.1%% is approximately 1%% total mixture")
    p.add_argument("--outdir", type=Path, required=True)
    a = p.parse_args()
    for path in (a.ledger, a.endpoints, a.paired_features):
        if path is not None and (not path.is_file() or not path.stat().st_size):
            p.error("missing or empty input: " + str(path))
    if a.outdir.exists():
        p.error("output exists: " + str(a.outdir))
    try:
        dose_index = 1 if a.member_dose_percent == .001 else 5
        detail, summary = parse_ledger(a.ledger, a.member_dose_percent / 100)
        response = parse_endpoints(a.endpoints, a.member_dose_percent / 100)
        bystanders = parse_bystanders(a.paired_features, dose_index) if a.paired_features else []
        a.outdir.mkdir(parents=True)
        write(a.outdir / "fnuc_call_transition_detail.tsv", CALL_FIELDS, detail)
        write(a.outdir / "fnuc_call_transition_summary.tsv", SUMMARY_FIELDS, summary)
        write(a.outdir / "fnuc_target_condition_response.tsv", RESPONSE_FIELDS, response)
        write(a.outdir / "fnuc_related_species_response.tsv", BYSTANDER_FIELDS, bystanders)
        draw(a.outdir / "fnuc_spike_biomarker_audit.svg", summary, response, bystanders,
             a.member_dose_percent)
        (a.outdir / "dose_metadata.tsv").write_text(
            f"member_dose_percent\tapprox_total_mix_percent\n{a.member_dose_percent:g}\t{10*a.member_dose_percent:g}\n",
            encoding="utf-8")
        (a.outdir / "DEVELOPMENT_ONLY.txt").write_text(
            "Spike-induced model-call changes are not proven biological false positives.\n", encoding="utf-8")
        (a.outdir / "SUCCESS").write_text("status\tPASS\n", encoding="utf-8")
    except (OSError, ValueError, KeyError) as error:
        p.error(str(error))
    print(f"[PASS] F. nucleatum spike-to-biomarker audit: {a.outdir}")


if __name__ == "__main__":
    main()
