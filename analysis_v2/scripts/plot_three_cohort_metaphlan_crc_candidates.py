#!/usr/bin/env python3
"""Provisional MetaPhlAn CRC candidates from unspiked community fits."""
import argparse
import csv
import hashlib
import html
import math
from pathlib import Path

COHORTS = ("feng", "yachida", "zeller")
REQUIRED = {"cohort", "analysis_population", "assembly_arm", "profiler",
            "contrast", "spike_fraction_target", "include", "model_spec",
            "feature", "effect", "q_value", "n_samples"}


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


def build(calls, outdir, top):
    if outdir.exists():
        raise ValueError("output exists: {}".format(outdir))
    fits = {}
    with calls.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not REQUIRED <= set(reader.fieldnames or ()):
            raise ValueError("disease calls lack required columns")
        for row in reader:
            if (row["analysis_population"] != "community" or
                row["assembly_arm"] != "original" or
                row["profiler"] != "metaphlan4" or
                row["contrast"] != "CRC_vs_Control" or
                row["include"] != "1" or
                row["model_spec"] != "primary_age_sex"):
                continue
            if number(row["spike_fraction_target"], "dose") != 0:
                continue
            cohort, feature = row["cohort"], row["feature"].strip()
            if cohort not in COHORTS or not feature:
                raise ValueError("unexpected cohort or empty feature")
            effect = number(row["effect"], "effect")
            q = number(row["q_value"], "q_value")
            n = number(row["n_samples"], "n_samples")
            if not 0 <= q <= 1 or n <= 0:
                raise ValueError("invalid q-value or sample count")
            key = cohort, feature
            value = effect, q, int(n)
            if key in fits and fits[key] != value:
                raise ValueError("contradictory baseline fit: {}".format(key))
            fits[key] = value
    if any(not any(c == cohort for c, _ in fits) for cohort in COHORTS):
        raise ValueError("one or more cohorts have no eligible baseline fits")
    features = sorted({feature for _, feature in fits})
    candidates = [f for f in features if any((c, f) in fits and fits[c, f][1] <= .05
                                               for c in COHORTS)]
    if not candidates:
        raise ValueError("no MetaPhlAn CRC calls at BH q <= 0.05")

    def classify(feature):
        observed = [(c, fits[c, feature]) for c in COHORTS if (c, feature) in fits]
        sig = [(c, v) for c, v in observed if v[1] <= .05]
        positive = sum(v[0] > 0 for _, v in sig)
        negative = sum(v[0] < 0 for _, v in sig)
        replicated = max(positive, negative)
        all_same_direction = (len(observed) == 3 and
                              all(v[0] > 0 for _, v in observed) or
                              len(observed) == 3 and
                              all(v[0] < 0 for _, v in observed))
        return len(observed), len(sig), replicated, all_same_direction, min(v[1] for _, v in sig)

    ranked = sorted(candidates, key=lambda f: (-classify(f)[2],
                                                -int(classify(f)[3]),
                                                -classify(f)[0], classify(f)[4], f))
    rows = []
    for feature in ranked:
        evaluable, significant, replicated, same_direction, best_q = classify(feature)
        for cohort in COHORTS:
            value = fits.get((cohort, feature))
            rows.append(dict(feature=feature, cohort=cohort,
                             evaluable=int(value is not None),
                             effect="" if value is None else "{:.9g}".format(value[0]),
                             q_value="" if value is None else "{:.9g}".format(value[1]),
                             n_samples="" if value is None else value[2],
                             significant=int(value is not None and value[1] <= .05),
                             evaluable_cohorts=evaluable, significant_cohorts=significant,
                             same_direction_significant_cohorts=replicated,
                             same_direction_all_three=int(same_direction),
                             best_q="{:.9g}".format(best_q)))
    outdir.mkdir(parents=True)
    table = outdir / "metaphlan_crc_candidates_all.tsv"
    write_tsv(table, rows)
    summary = [dict(metric="eligible_baseline_features", value=len(features)),
               dict(metric="discovered_in_any_cohort", value=len(candidates)),
               dict(metric="same_direction_significant_in_two_or_more", value=sum(classify(f)[2] >= 2 for f in candidates)),
               dict(metric="same_direction_significant_in_all_three", value=sum(classify(f)[2] == 3 for f in candidates)),
               dict(metric="same_effect_direction_all_three_evaluable", value=sum(classify(f)[3] for f in candidates)),
               dict(metric="plotted_features", value=min(top, len(candidates)))]
    write_tsv(outdir / "candidate_summary.tsv", summary)

    shown = ranked[:top]
    width, left, cell_w, row_h = 1320, 530, 215, 47
    height = 195 + row_h * len(shown) + 90
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}">'.format(width, height, width, height),
             '<rect width="100%" height="100%" fill="white"/>',
             '<text x="30" y="42" font-family="sans-serif" font-size="28" font-weight="bold">Which MetaPhlAn CRC candidates recur across cohorts?</text>',
             '<text x="30" y="72" font-family="sans-serif" font-size="17">Original unspiked community abundance · CRC vs Control · BH q ≤ 0.05 · DEVELOPMENT ONLY</text>',
             '<text x="30" y="105" font-family="sans-serif" font-size="15">Rows ranked by same-direction significant cohorts, then cross-cohort direction and best q; no pooled effect.</text>']
    for i, cohort in enumerate(COHORTS):
        parts.append('<text x="{}" y="152" font-family="sans-serif" font-size="20" text-anchor="middle" font-weight="bold">{}</text>'.format(left + i * cell_w + 96, cohort.title()))
    for ri, feature in enumerate(shown):
        y = 170 + ri * row_h
        short = feature if len(feature) <= 48 else feature[:45] + "..."
        parts.append('<title>{}</title><text x="30" y="{}" font-family="sans-serif" font-size="17">{}</text>'.format(
            html.escape(feature), y + 27, html.escape(short)))
        for ci, cohort in enumerate(COHORTS):
            value = fits.get((cohort, feature))
            if value is None:
                fill, label = "#e6e9ed", "not evaluable"
            else:
                effect, q, _ = value
                significant = q <= .05
                fill = ("#087e75" if effect > 0 else "#ad5862") if significant else "#e7f2f0" if effect > 0 else "#f7edef"
                label = "{:+.2f}{}  q={:.2g}".format(effect, "*" if significant else "", q)
            x = left + ci * cell_w
            parts.append('<rect x="{}" y="{}" width="205" height="40" rx="3" fill="{}"/>'.format(x, y, fill))
            parts.append('<text x="{}" y="{}" font-family="sans-serif" font-size="16" fill="{}">{}</text>'.format(
                x + 10, y + 26, "white" if value is not None and value[1] <= .05 else "#26313b", html.escape(label)))
    footer = height - 49
    parts.append('<text x="30" y="{}" font-family="sans-serif" font-size="15">Cell: signed model effect and BH q; * q ≤ 0.05. Dark teal/red = significant positive/negative; pale = not significant; grey = not evaluable.</text>'.format(footer))
    parts.append('<text x="30" y="{}" font-family="sans-serif" font-size="15">Shared direction is not proof of biological causality or spike-in robustness. All candidates and counts are in the source TSV.</text>'.format(footer + 24))
    parts.append('</svg>\n')
    (outdir / "metaphlan_crc_candidates.svg").write_text("".join(parts), encoding="utf-8")
    (outdir / "source.sha256").write_text("{}  {}\n{}  {}\n".format(
        digest(calls), calls.resolve(), digest(table), table.resolve()), encoding="utf-8")
    (outdir / "DEVELOPMENT_ONLY.txt").write_text("status\tDEVELOPMENT_ONLY\n", encoding="utf-8")
    (outdir / "SUCCESS").write_text("status\tPASS\n", encoding="utf-8")
    print("[PASS] {} MetaPhlAn candidates; {} shared significant direction in >=2 cohorts; {} in all 3".format(
        len(candidates), summary[2]["value"], summary[3]["value"]))
    print(outdir / "metaphlan_crc_candidates.svg")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calls", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--top", type=int, default=12)
    args = parser.parse_args()
    if args.top < 1:
        parser.error("--top must be positive")
    try:
        build(args.calls, args.outdir, args.top)
    except (OSError, ValueError) as error:
        raise SystemExit("[ERROR] {}".format(error)) from error
