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


def build(matched, audit):
    files = {
        "models": matched / "matched_species_models.tsv",
        "baseline": matched / "matched_species_baseline_summary.tsv",
        "spikes": matched / "matched_species_spike_summary.tsv",
        "transitions": audit / "fnuc_call_transition_summary.tsv",
        "response": audit / "fnuc_target_condition_response.tsv",
        "related": audit / "fnuc_related_species_response.tsv",
    }
    for name, path in files.items():
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"missing {name}: {path}")
    tables = {name: read(path) for name, path in files.items()}
    if not tables["models"] or not tables["baseline"] or not tables["spikes"]:
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
            direct = {condition: unique(tables["spikes"],
                lambda r, condition=condition: same(r) and r["feature"] == TARGET and
                r["condition"] == condition and math.isclose(num(r, "dose_percent"), .01),
                f"direct spike {cohort}/{profiler}/{condition}") for condition in CONDITIONS}
            transition = {contrast: unique(tables["transitions"],
                lambda r, contrast=contrast: same(r) and r["contrast"] == contrast,
                f"transition {cohort}/{profiler}/{contrast}")
                for contrast in ("CRC_vs_Control", "Adenoma_vs_Control")}
            community = {condition: unique(tables["response"],
                lambda r, condition=condition: same(r) and r["condition"] == condition,
                f"community {cohort}/{profiler}/{condition}") for condition in CONDITIONS}
            related = [r for r in tables["related"] if same(r)]
            result.append(dict(profiler=profiler, cohort=cohort, model=model,
                               baseline=baseline, direct=direct, transition=transition,
                               community=community, related=related))
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
    width, height = 2000, 1370
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>']
    text(parts, 40, 50, "Can measurement explain the missing F. nucleatum adenoma call?", 32, weight="bold")
    text(parts, 40, 80, "Six cohort–profiler contexts · 0.001% per member community mix; 0.01% independent F. nucleatum spike · DEVELOPMENT ONLY", 17)
    text(parts, 40, 106, "Each row tests a different mechanism; these are not interchangeable scales or a causal score.", 15)
    x0, cw, gap = 270, 545, 20
    colors = {"kraken2_bracken": "#b75000", "metaphlan4": "#006b9b"}
    for ci, cohort in enumerate(COHORTS):
        text(parts, x0 + ci*(cw+gap)+cw/2, 135, cohort.title(), 22, weight="bold", anchor="middle")
    for pi, profiler in enumerate(PROFILERS):
        block = 182 + pi*485
        color = colors[profiler]
        text(parts, 40, block+18, "Kraken2 + Bracken" if pi == 0 else "MetaPhlAn 4", 17,
             color=color, weight="bold")
        labels = ((27, "Native detection C/A/CRC"), (57, "Adenoma positive abundance"),
                  (106, "Adenoma effect [CI], q"), (137, "CRC effect [CI], q"),
                  (187, "0.001% zero→positive"), (217, "Adenoma recovery [Q1,Q3]"),
                  (247, "Adenoma recovery <0.5"), (295, "0.01% independent recovery"),
                  (325, "0.01% detected"), (375, "Adenoma effect before→after"),
                  (405, "Adenoma q before→after"), (435, "CRC call; other calls"),
                  (464, "Related Fuso paired shift"))
        for offset, label in labels:
            text(parts, 40, block+offset, label, 13)
        for offset in (75, 150, 260, 337):
            line(parts, 40, block+offset, 1968, block+offset)
        for ci, cohort in enumerate(COHORTS):
            r = next(z for z in contexts if z["cohort"] == cohort and z["profiler"] == profiler)
            x = x0 + ci*(cw+gap)
            b, d, c, t, m = r["baseline"], r["direct"], r["community"], r["transition"], r["model"]
            def put(offset, value, size=14, emphasis=False):
                text(parts, x+10, block+offset, value, size, color if emphasis else "#243240",
                     "bold" if emphasis else "normal")
            parts.append(f'<rect x="{x}" y="{block}" width="{cw}" height="470" rx="6" fill="#f5f8fa"/>')
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
            ad = d["Adenoma"]
            put(295, f'{num(ad,"recovery_median"):.2f} [{num(ad,"recovery_q1"):.2f},{num(ad,"recovery_q3"):.2f}]')
            put(325, f'{ad["observed_positive"]}/{ad["n"]}')
            at, ct = t["Adenoma_vs_Control"], t["CRC_vs_Control"]
            put(375, f'{num(at,"target_baseline_effect"):+.2f} → {num(at,"target_dose_effect"):+.2f}', 14, True)
            put(405, f'{num(at,"target_baseline_q"):.2g} → {num(at,"target_dose_q"):.2g} '
                     f'({"called" if at["target_dose_called"] == "1" else "not called"})')
            status = lambda v: "call" if v == "1" else "no call"
            put(435, f'{status(ct["target_baseline_called"])} → {status(ct["target_dose_called"])}; '
                     f'Fuso +{ct["related_gained"]}/−{ct["related_lost"]}, '
                     f'other +{ct["other_bystander_gained"]}/−{ct["other_bystander_lost"]}', 13)
            # A related-species response is optional; the model-call count is shown independently.
            if r["related"]:
                strongest = max(r["related"], key=lambda z: abs(num(z,"response_median")))
                label = strongest["feature"].replace("Fusobacterium ", "F. ")
                put(464, f'{label}: {100*num(strongest,"response_median"):+.2g} percentage points', 12)
            else:
                put(464, "paired response not supplied", 12)
    text(parts, 40, 1191, "Reading the figure", 20, weight="bold")
    text(parts, 40, 1221, "A recovered spike weakens a specific detection-failure explanation; it does not establish native measurement accuracy or biological absence.", 16)
    text(parts, 40, 1249, "The equal spike is a ten-species mixture: changes in calls or related taxa cannot be attributed solely to F. nucleatum.", 16)
    text(parts, 40, 1277, "Still unresolved: direct CRC–adenoma contrast and power; read depth/batch/marker QC; native-signal validation and species-specific misassignment.", 16)
    text(parts, 40, 1310, "C/A/CRC = Control/Adenoma/CRC; recovery 1 is ideal. Profiler-scale recoveries are not cross-platform abundance comparisons.", 15)
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matched-dir", type=Path, required=True)
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if args.outdir.exists():
        parser.error("output exists: " + str(args.outdir))
    try:
        contexts = build(args.matched_dir, args.audit_dir)
        args.outdir.mkdir(parents=True)
        draw(args.outdir / "fnuc_integrated_methods.svg", contexts)
        (args.outdir / "SUCCESS").write_text("status\tPASS\n", encoding="utf-8")
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    print("[PASS] integrated F. nucleatum methods audit: " + str(args.outdir))


if __name__ == "__main__":
    main()
