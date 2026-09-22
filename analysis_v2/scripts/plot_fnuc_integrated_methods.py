#!/usr/bin/env python3
"""One-page, evidence-graded F. nucleatum adenoma measurement audit."""
from __future__ import annotations

import argparse
import csv
import html
import math
from pathlib import Path

COHORTS = ("feng", "yachida", "zeller")
PROFILERS = ("kraken2_bracken", "metaphlan4")
CONDITIONS = ("Control", "Adenoma", "CRC")
TARGET = "Fusobacterium nucleatum"
INDEPENDENT_DOSES = (1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2)


def read(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def unique(rows, key, label):
    found = [row for row in rows if key(row)]
    if len(found) != 1:
        raise ValueError(f"expected one {label}, found {len(found)}")
    return found[0]


def num(row, field):
    value = float(row[field])
    if not math.isfinite(value):
        raise ValueError(f"nonfinite {field}")
    return value


def paired_calls(path):
    """First significant tested independent single-species dose by background."""
    required = {"cohort", "analysis_population", "target_label", "assembly_arm",
                "profiler", "contrast", "spike_fraction_target", "q_threshold",
                "target_called", "target_effect", "target_q_value"}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not required <= set(reader.fieldnames or ()):
            raise ValueError("paired metrics lack required target-call columns")
        records = list(reader)
    selected = {}
    for cohort in COHORTS:
        for profiler in PROFILERS:
            for condition in CONDITIONS:
                contrast = "spiked_vs_matched_baseline__background_" + condition
                matches = [r for r in records if r["cohort"] == cohort and r["profiler"] == profiler and
                    r["contrast"] == contrast and r["target_label"] == "Fnuc" and
                    r["analysis_population"] == "independent" and r["assembly_arm"] == "original" and
                    math.isclose(num(r, "q_threshold"), .05, abs_tol=1e-10)]
                doses = {}
                for row in matches:
                    achieved = num(row, "spike_fraction_target")
                    dose = min(INDEPENDENT_DOSES, key=lambda d: abs(achieved-d)/d)
                    if abs(achieved-dose)/dose > .05 or dose in doses:
                        raise ValueError("invalid/duplicate independent dose " + repr((cohort, profiler, condition, achieved)))
                    if row["target_called"] not in ("0", "1"):
                        raise ValueError("invalid paired target_called")
                    effect, q = num(row, "target_effect"), num(row, "target_q_value")
                    if not 0 <= q <= 1 or (row["target_called"] == "1") != (effect > 0 and q <= .05):
                        raise ValueError("paired target call disagrees with effect/q")
                    doses[dose] = row
                if set(doses) != set(INDEPENDENT_DOSES):
                    raise ValueError("incomplete independent dose grid " + repr((cohort, profiler, condition)))
                first = next(((dose, doses[dose]) for dose in INDEPENDENT_DOSES
                              if doses[dose]["target_called"] == "1"), None)
                selected[cohort, profiler, condition] = first
    return selected


def build(matched, audit, high_audit, paired_metrics):
    files = {
        "models": matched / "matched_species_models.tsv",
        "baseline": matched / "matched_species_baseline_summary.tsv",
        "transitions": audit / "fnuc_call_transition_summary.tsv",
        "response": audit / "fnuc_target_condition_response.tsv",
    }
    for name, path in files.items():
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"missing {name}: {path}")
    tables = {name: read(path) for name, path in files.items()}
    metadata = read(high_audit / "dose_metadata.tsv")
    if len(metadata) != 1 or not math.isclose(num(metadata[0], "member_dose_percent"), .1):
        raise ValueError("high audit must be the 0.1%-per-member community spike")
    high_transitions = read(high_audit / "fnuc_call_transition_summary.tsv")
    high_response = read(high_audit / "fnuc_target_condition_response.tsv")
    paired = paired_calls(paired_metrics)
    if not tables["models"] or not tables["baseline"]:
        raise ValueError("matched-species tables are empty")
    result = []
    for profiler in PROFILERS:
        for cohort in COHORTS:
            same = lambda r: r["profiler"] == profiler and r["cohort"] == cohort
            model = {contrast: unique(tables["models"],
                lambda r, contrast=contrast: same(r) and r["feature"] == TARGET and
                r["contrast"] == contrast, f"{cohort}/{profiler}/{contrast}")
                for contrast in ("CRC_vs_Control", "Adenoma_vs_Control")}
            baseline = {condition: unique(tables["baseline"],
                lambda r, condition=condition: same(r) and r["feature"] == TARGET and
                r["condition"] == condition, f"baseline {cohort}/{profiler}/{condition}")
                for condition in CONDITIONS}
            transition = {contrast: unique(tables["transitions"],
                lambda r, contrast=contrast: same(r) and r["contrast"] == contrast,
                f"transition {cohort}/{profiler}/{contrast}")
                for contrast in ("CRC_vs_Control", "Adenoma_vs_Control")}
            community = {condition: unique(tables["response"],
                lambda r, condition=condition: same(r) and r["condition"] == condition,
                f"community {cohort}/{profiler}/{condition}") for condition in CONDITIONS}
            high_t = {contrast: unique(high_transitions,
                lambda r, contrast=contrast: same(r) and r["contrast"] == contrast,
                f"high-dose transition {cohort}/{profiler}/{contrast}")
                for contrast in ("CRC_vs_Control", "Adenoma_vs_Control")}
            high_c = {condition: unique(high_response,
                lambda r, condition=condition: same(r) and r["condition"] == condition,
                f"high-dose response {cohort}/{profiler}/{condition}") for condition in CONDITIONS}
            result.append(dict(profiler=profiler, cohort=cohort, model=model,
                               baseline=baseline, transition=transition,
                               community=community,
                               high_transition=high_t, high_community=high_c,
                               paired={condition: paired[cohort, profiler, condition]
                                       for condition in CONDITIONS}))
    return result


def text(parts, x, y, value, size=15, color="#243240", weight="normal", anchor=None):
    align = f' text-anchor="{anchor}"' if anchor else ""
    parts.append(f'<text x="{x}" y="{y}" font-family="Arial,sans-serif" font-size="{size}" '
                 f'font-weight="{weight}" fill="{color}"{align}>{html.escape(str(value))}</text>')


def line(parts, x1, y1, x2, y2, color="#dae3e9", width=1):
    parts.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" '
                 f'stroke="{color}" stroke-width="{width}"/>')


def ratio(row, numerator, denominator):
    n, d = int(row[numerator]), int(row[denominator])
    return f"{n}/{d}" if d else "0/0"


def pct(row):
    return 100 * int(row["positive"]) / int(row["n"])


def draw(path, contexts):
    width, height = 2000, 1600
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>']
    text(parts, 40, 50, "Can measurement explain the missing F. nucleatum adenoma call?", 32, weight="bold")
    text(parts, 40, 80, "Six contexts · community mixes at 0.001% and 0.1% per member; independent single-species paired-dose test · DEVELOPMENT ONLY", 17)
    text(parts, 40, 106, "Each row tests a different mechanism; these are not interchangeable scales or a causal score.", 15)
    x0, cw, gap = 270, 545, 20
    colors = {"kraken2_bracken": "#b75000", "metaphlan4": "#006b9b"}
    for ci, cohort in enumerate(COHORTS):
        text(parts, x0 + ci*(cw+gap)+cw/2, 135, cohort.title(), 22, weight="bold", anchor="middle")
    for pi, profiler in enumerate(PROFILERS):
        block = 182 + pi*590
        color = colors[profiler]
        text(parts, 40, block+18, "Kraken2 + Bracken" if pi == 0 else "MetaPhlAn 4", 17,
             color=color, weight="bold")
        labels = ((27, "Native detection C/A/CRC"), (57, "Adenoma positive abundance"),
                  (106, "Adenoma effect [CI], q"), (137, "CRC effect [CI], q"),
                  (187, "0.001% zero→positive"), (217, "0.001% adenoma recovery"),
                  (247, "0.001% recovery <0.5"),
                  (280, "Independent first paired call"),
                  (320, "0.001% adenoma effect"),
                  (350, "0.001% adenoma q / call"), (380, "0.001% CRC / other calls"),
                  (420, "0.1% zero→positive"), (450, "0.1% adenoma recovery"),
                  (480, "0.1% adenoma call"), (510, "0.1% CRC call"),
                  (540, "0.1% other calls"))
        for offset, label in labels:
            text(parts, 40, block+offset, label, 13)
        for offset in (75, 150, 260, 393):
            line(parts, 40, block+offset, 1968, block+offset)
        for ci, cohort in enumerate(COHORTS):
            r = next(z for z in contexts if z["cohort"] == cohort and z["profiler"] == profiler)
            x = x0 + ci*(cw+gap)
            b, c, t, m = r["baseline"], r["community"], r["transition"], r["model"]
            def put(offset, value, size=14, emphasis=False):
                text(parts, x+10, block+offset, value, size, color if emphasis else "#243240",
                     "bold" if emphasis else "normal")
            parts.append(f'<rect x="{x}" y="{block}" width="{cw}" height="555" rx="6" fill="#f5f8fa"/>')
            put(27, " / ".join(f'{k[:3]} {b[k]["positive"]}/{b[k]["n"]} ({pct(b[k]):.0f}%)'
                              for k in CONDITIONS), 13)
            med = b["Adenoma"]["positive_abundance_median_percent"]
            q1 = b["Adenoma"].get("positive_abundance_q1_percent", "")
            q3 = b["Adenoma"].get("positive_abundance_q3_percent", "")
            native = (f'{float(med):.3g}% [{float(q1):.3g},{float(q3):.3g}] positive-only'
                      if med and q1 and q3 else
                      f'{float(med):.3g}% positive-only' if med else "no positives")
            put(57, native + f'; all-IQR {num(b["Adenoma"], "all_sample_iqr_percent"):.3g}%', 13)
            am, cm = m["Adenoma_vs_Control"], m["CRC_vs_Control"]
            put(106, f'{num(am,"effect"):+.2f} [{num(am,"lower_95"):+.2f},{num(am,"upper_95"):+.2f}]  q={num(am,"q_value"):.2g}', 14, True)
            put(137, f'{num(cm,"effect"):+.2f} [{num(cm,"lower_95"):+.2f},{num(cm,"upper_95"):+.2f}]  q={num(cm,"q_value"):.2g}')
            put(187, " / ".join(f'{k[:3]} {ratio(c[k],"zero_rescued","baseline_zero")}' for k in CONDITIONS))
            ac = c["Adenoma"]
            put(217, f'{num(ac,"recovery_median"):.2f} [{num(ac,"recovery_q1"):.2f},{num(ac,"recovery_q3"):.2f}]')
            put(247, f'<0.5: {ac["below_half"]}/{ac["n"]}')
            paired = r["paired"]
            put(280, " / ".join(
                f'{k[:3]} {100*paired[k][0]:g}% q={num(paired[k][1],"target_q_value"):.2g}'
                if paired[k] else f'{k[:3]} NR' for k in CONDITIONS), 12)
            at, ct = t["Adenoma_vs_Control"], t["CRC_vs_Control"]
            put(320, f'{num(at,"target_baseline_effect"):+.2f} → {num(at,"target_dose_effect"):+.2f}', 14, True)
            put(350, f'{num(at,"target_baseline_q"):.2g} → {num(at,"target_dose_q"):.2g} '
                     f'({"called" if at["target_dose_called"] == "1" else "not called"})')
            status = lambda v: "call" if v == "1" else "no call"
            put(380, f'{status(ct["target_baseline_called"])} → {status(ct["target_dose_called"])}; '
                     f'Fuso +{ct["related_gained"]}/−{ct["related_lost"]}, '
                     f'other +{ct["other_bystander_gained"]}/−{ct["other_bystander_lost"]}', 13)
            hc, ht = r["high_community"], r["high_transition"]
            put(420, " / ".join(f'{k[:3]} {ratio(hc[k],"zero_rescued","baseline_zero")}'
                                for k in CONDITIONS))
            ha = hc["Adenoma"]
            put(450, f'{num(ha,"recovery_median"):.2f} [{num(ha,"recovery_q1"):.2f},{num(ha,"recovery_q3"):.2f}]; '
                     f'<0.5 {ha["below_half"]}/{ha["n"]}')
            for offset, contrast in ((480, "Adenoma_vs_Control"), (510, "CRC_vs_Control")):
                h = ht[contrast]
                put(offset, f'{num(h,"target_baseline_effect"):+.2f} → {num(h,"target_dose_effect"):+.2f}; '
                    f'q {num(h,"target_baseline_q"):.2g} → {num(h,"target_dose_q"):.2g}; '
                    f'{"called" if h["target_dose_called"] == "1" else "not called"}', 12, True)
            high_crc = ht["CRC_vs_Control"]
            put(540, f'CRC: Fuso +{high_crc["related_gained"]}/−{high_crc["related_lost"]}; '
                     f'other +{high_crc["other_bystander_gained"]}/−{high_crc["other_bystander_lost"]}', 12)
    text(parts, 40, 1428, "Reading the figure", 20, weight="bold")
    text(parts, 40, 1458, "Independent first paired call = lowest tested single-species dose with F. nucleatum increased vs own unspiked profile, BH q≤0.05; NR = none.", 16)
    text(parts, 40, 1486, "0.1% per member is ~1% total ten-species mixture: a strong perturbation, not an endogenous-equivalent detection limit.", 16)
    text(parts, 40, 1514, "A recovered spike weakens one measured failure mode; it does not establish native accuracy or biological absence.", 16)
    text(parts, 40, 1542, "Call changes cannot be attributed solely to F. nucleatum; direct CRC–adenoma contrast/power, read-level QC and native validation remain unresolved.", 16)
    text(parts, 40, 1570, "C/A/CRC = Control/Adenoma/CRC; recovery 1 is ideal. Independent paired and community disease-call tests are distinct.", 15)
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matched-dir", type=Path, required=True)
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--high-audit-dir", type=Path, required=True)
    parser.add_argument("--paired-metrics", type=Path, required=True,
                        help="phenotype-stratified biomarker_propagation_metrics.tsv with all independent doses")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if args.outdir.exists():
        parser.error("output exists: " + str(args.outdir))
    try:
        contexts = build(args.matched_dir, args.audit_dir, args.high_audit_dir,
                         args.paired_metrics)
        args.outdir.mkdir(parents=True)
        draw(args.outdir / "fnuc_integrated_methods.svg", contexts)
        (args.outdir / "SUCCESS").write_text("status\tPASS\n", encoding="utf-8")
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    print("[PASS] integrated F. nucleatum methods audit: " + str(args.outdir))


if __name__ == "__main__":
    main()
