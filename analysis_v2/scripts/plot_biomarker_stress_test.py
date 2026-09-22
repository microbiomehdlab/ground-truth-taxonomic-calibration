#!/usr/bin/env python3
"""Slide-sized, development-only CRC biomarker stress test from corrected ledger."""
import argparse
import csv
import hashlib
import html
import math
from pathlib import Path

PROFILERS = ("kraken2_bracken", "metaphlan4")
COLORS = {
    "retained": "#168a71", "lost": "#e2a317",
    "reversed": "#b72878", "implanted": "#dce1e7",
}


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("empty table: {}".format(path))
        return reader.fieldnames, list(reader)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def build(ledger_path, panel_path, outdir, cohort, dose_percent, top):
    if outdir.exists():
        raise ValueError("output exists: {}".format(outdir))
    panel_fields, panel = read_tsv(panel_path)
    if "label" not in panel_fields or not panel:
        raise ValueError("invalid spike panel")
    targets = [row["label"] for row in panel]
    if len(targets) != len(set(targets)):
        raise ValueError("duplicate target label")
    fields, ledger = read_tsv(ledger_path)
    required = {"cohort", "analysis_population", "assembly_arm", "profiler",
                "contrast", "target_label", "spike_fraction_target", "q_threshold",
                "feature", "feature_role", "baseline_called", "dose_called",
                "baseline_q_value", "effect_sign_changed"}
    if not required <= set(fields):
        raise ValueError("transition ledger lacks required columns")
    selected = {}
    baseline_q = {}
    for row in ledger:
        if (row["cohort"] != cohort or row["analysis_population"] != "independent"
                or row["assembly_arm"] != "original" or row["contrast"] != "CRC_vs_Control"):
            continue
        try:
            dose = float(row["spike_fraction_target"])
            q_cut = float(row["q_threshold"])
            q = float(row["baseline_q_value"])
        except ValueError as error:
            raise ValueError("nonnumeric ledger value") from error
        if not math.isclose(dose, dose_percent / 100, rel_tol=.05) or not math.isclose(q_cut, .05):
            continue
        if row["profiler"] not in PROFILERS or row["target_label"] not in targets:
            raise ValueError("unexpected profiler or target in selected ledger")
        if row["baseline_called"] != "1":
            continue
        if not (math.isfinite(q) and 0 <= q <= .05):
            raise ValueError("invalid baseline q")
        feature = row["feature"]
        base_key = row["profiler"], feature
        if base_key in baseline_q and not math.isclose(baseline_q[base_key], q, abs_tol=1e-12):
            raise ValueError("inconsistent baseline q for {}".format(base_key))
        baseline_q[base_key] = q
        key = row["profiler"], feature, row["target_label"]
        if key in selected:
            raise ValueError("duplicate stress-test cell: {}".format(key))
        if row["feature_role"] == "implanted_target":
            status = "implanted"
        elif row["feature_role"] == "bystander":
            if row["dose_called"] not in ("0", "1") or row["effect_sign_changed"] not in ("0", "1"):
                raise ValueError("invalid call/sign flag")
            status = ("reversed" if row["effect_sign_changed"] == "1" else
                      "retained" if row["dose_called"] == "1" else "lost")
        else:
            raise ValueError("invalid feature role")
        selected[key] = status
    chosen = {}
    for profiler in PROFILERS:
        candidates = sorted(((q, feature) for (tool, feature), q in baseline_q.items()
                             if tool == profiler))
        if not candidates:
            raise ValueError("no baseline CRC biomarkers for {}".format(profiler))
        chosen[profiler] = [feature for _, feature in candidates[:top]]
    cells = []
    for profiler in PROFILERS:
        for feature in chosen[profiler]:
            for target in targets:
                key = profiler, feature, target
                if key not in selected:
                    raise ValueError("incomplete challenge grid: {}".format(key))
                cells.append(dict(cohort=cohort, profiler=profiler, feature=feature,
                                  target_label=target, baseline_q=baseline_q[(profiler, feature)],
                                  status=selected[key]))
    outdir.mkdir(parents=True)
    source = outdir / "biomarker_stress_test_cells.tsv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cells[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(cells)
    width = 1530
    row_height = 31
    left, top_y, cell_w = 440, 143, 79
    height = top_y + sum(len(chosen[p]) * row_height + 70 for p in PROFILERS) + 100
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}">'.format(width,height,width,height),
             '<rect width="100%" height="100%" fill="white"/>',
             '<text x="35" y="40" font-family="sans-serif" font-size="29" font-weight="bold">Which CRC biomarker calls survive a known challenge?</text>',
             '<text x="35" y="70" font-family="sans-serif" font-size="17">{} · independent spikes at {:.3f}% · baseline CRC vs Control, BH q ≤ 0.05 · DEVELOPMENT ONLY</text>'.format(
                 html.escape(cohort.title()), dose_percent)]
    for i, target in enumerate(targets):
        parts.append('<text x="{}" y="110" font-family="sans-serif" font-size="16" text-anchor="middle">{}</text>'.format(
            left + i * cell_w + cell_w / 2, html.escape(target)))
    y = top_y
    for profiler in PROFILERS:
        label = "Kraken2 + Bracken" if profiler == "kraken2_bracken" else "MetaPhlAn 4"
        parts.append('<text x="35" y="{}" font-family="sans-serif" font-size="21" font-weight="bold">{}</text>'.format(y-10,label))
        for feature in chosen[profiler]:
            relevant = [r for r in cells if r["profiler"] == profiler and r["feature"] == feature]
            eligible = sum(r["status"] != "implanted" for r in relevant)
            retained = sum(r["status"] == "retained" for r in relevant)
            parts.append('<text x="35" y="{}" font-family="sans-serif" font-size="15">{}</text>'.format(
                y + 21, html.escape(feature)))
            for i, target in enumerate(targets):
                status = selected[(profiler, feature, target)]
                parts.append('<rect x="{}" y="{}" width="{}" height="25" rx="3" fill="{}"/>'.format(
                    left + i * cell_w + 2, y, cell_w - 4, COLORS[status]))
            parts.append('<text x="{}" y="{}" font-family="sans-serif" font-size="15">{}/{}</text>'.format(
                left + len(targets) * cell_w + 12, y + 20, retained, eligible))
            y += row_height
        y += 70
    legend = (("retained", "Retained, same direction"), ("lost", "Lost significance"),
              ("reversed", "Direction reversed"), ("implanted", "Direct target excluded"))
    for i, (status, label) in enumerate(legend):
        x = 38 + i * 365
        parts.append('<rect x="{}" y="{}" width="20" height="20" fill="{}"/>'.format(x,height-66,COLORS[status]))
        parts.append('<text x="{}" y="{}" font-family="sans-serif" font-size="16">{}</text>'.format(x+28,height-50,label))
    parts.append('<text x="38" y="{}" font-family="sans-serif" font-size="14">Rows selected by unspiked q only; right counts omit direct-target cells. Technical robustness is not biological validation.</text>'.format(height-16))
    parts.append('</svg>\n')
    (outdir / "biomarker_stress_test.svg").write_text("".join(parts), encoding="utf-8")
    (outdir / "source.sha256").write_text("{}  {}\n{}  {}\n".format(
        digest(ledger_path), ledger_path.resolve(), digest(panel_path), panel_path.resolve()), encoding="utf-8")
    (outdir / "DEVELOPMENT_ONLY.txt").write_text("status\tDEVELOPMENT_ONLY\n", encoding="utf-8")
    (outdir / "SUCCESS").write_text("status\tPASS\n", encoding="utf-8")
    print("[PASS] biomarker stress-test slide: {}".format(outdir))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", required=True, type=Path)
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--cohort", default="yachida")
    parser.add_argument("--dose-percent", type=float, default=.1)
    parser.add_argument("--top", type=int, default=10)
    args = parser.parse_args()
    if not (math.isfinite(args.dose_percent) and args.dose_percent > 0 and args.top > 0):
        parser.error("dose percent and top must be positive")
    try:
        build(args.ledger, args.panel, args.outdir, args.cohort, args.dose_percent, args.top)
    except (OSError, ValueError) as error:
        raise SystemExit("[ERROR] {}".format(error)) from error
