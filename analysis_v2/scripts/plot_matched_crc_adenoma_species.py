#!/usr/bin/env python3
"""Compare three explicitly matched spike-panel species across profilers and disease stages."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from plot_shared_crc_adenoma_bridge import (
    COHORTS, CONDITIONS, CONTRASTS, PROFILERS, CALL_FIELDS, CALL_OUT,
    BASE_OUT, SPIKE_OUT, direct_spikes, full_baselines, number,
    panel_aliases, sha256, txt, write_tsv,
)

LABELS = ("Fnuc", "Pmic", "Dpne")
NAMES = {"Fnuc": "Fusobacterium nucleatum", "Pmic": "Parvimonas micra",
         "Dpne": "Dialister pneumosintes"}


def matched_features(panel_map):
    selected = {}
    for profiler in PROFILERS:
        selected[profiler] = []
        for label in LABELS:
            matches = [feature for (tool, feature), (lab, _) in panel_map.items()
                       if tool == profiler and lab == label]
            if len(matches) != 1:
                raise ValueError(f"expected one {label} alias for {profiler}; found {matches}")
            selected[profiler].append(matches[0])
    return selected


def matched_calls(path, selected):
    wanted = {(p, f) for p in PROFILERS for f in selected[p]}
    fits = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not CALL_FIELDS <= set(reader.fieldnames or ()):
            raise ValueError("disease calls lack required model columns")
        for row in reader:
            if ((row["profiler"], row["feature"]) not in wanted or
                row["cohort"] not in COHORTS or row["contrast"] not in CONTRASTS or
                row["analysis_population"] != "community" or
                row["assembly_arm"] != "original" or row["include"] != "1" or
                row["model_spec"] != "primary_age_sex" or
                number(row["spike_fraction_target"], "model dose") != 0):
                continue
            key = row["profiler"], row["feature"], row["cohort"], row["contrast"]
            value = tuple(number(row[field], field) for field in
                          ("effect", "lower_95", "upper_95", "q_value", "n_samples"))
            if key in fits and fits[key] != value:
                raise ValueError("conflicting disease model fit: " + repr(key))
            fits[key] = value
    rows = []
    for profiler in PROFILERS:
        for feature in selected[profiler]:
            for cohort in COHORTS:
                crc = fits.get((profiler, feature, cohort, "CRC_vs_Control"))
                for contrast in CONTRASTS:
                    value = fits.get((profiler, feature, cohort, contrast))
                    rows.append(dict(profiler=profiler, feature=feature,
                                     display_feature=feature, cohort=cohort, contrast=contrast,
                                     evaluable=int(value is not None),
                                     effect="" if value is None else value[0],
                                     lower_95="" if value is None else value[1],
                                     upper_95="" if value is None else value[2],
                                     q_value="" if value is None else value[3],
                                     n_samples="" if value is None else int(value[4]),
                                     significant=int(value is not None and value[3] <= .05),
                                     same_direction_as_crc="" if value is None or crc is None
                                     else int(value[0] * crc[0] > 0),
                                     direct_spike_target=next(label for label in LABELS if
                                         selected[profiler][LABELS.index(label)] == feature)))
    return rows


def draw(path, selected, calls, baselines, spikes):
    width, height = 1710, 1580
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>']
    txt(parts, 28, 43, "Matched species: CRC, adenoma, baseline and low-dose spike response", 28, weight="bold")
    txt(parts, 28, 72, "Same biological spike-panel species in both profilers; rows are NOT selected by CRC significance. DEVELOPMENT ONLY", 16)
    txt(parts, 28, 99, "C and A: age/sex-adjusted effect vs Control [95% CI], BH q. Baseline: detected/total. Spike: 0.01% direct recovery median [Q1,Q3].", 14)
    lookup = {(r["profiler"], r["feature"], r["cohort"], r["contrast"]): r for r in calls}
    base = {(r["profiler"], r["feature"], r["cohort"], r["condition"]): r for r in baselines}
    spike = {(r["profiler"], r["feature"], r["cohort"], r["condition"]): r for r in spikes
             if math.isclose(r["dose_percent"], .01)}
    for ci, cohort in enumerate(COHORTS):
        txt(parts, 472 + ci*410, 147, cohort.title(), 21, anchor="middle", weight="bold")
    for li, label in enumerate(LABELS):
        section_y = 178 + li*450
        txt(parts, 28, section_y, NAMES[label], 22, weight="bold")
        for pi, profiler in enumerate(PROFILERS):
            feature = selected[profiler][li]
            y = section_y + 17 + pi*198
            profiler_name = "Kraken2 + Bracken" if pi == 0 else "MetaPhlAn 4"
            txt(parts, 28, y+30, profiler_name, 17, color="#c05600" if pi == 0 else "#006a9c", weight="bold")
            if feature != NAMES[label]:
                txt(parts, 28, y+54, "Reported: " + feature, 13)
            for ci, cohort in enumerate(COHORTS):
                x = 267 + ci*410
                parts.append(f'<rect x="{x}" y="{y}" width="391" height="183" rx="5" fill="#f1f5f7"/>')
                crc = lookup[profiler, feature, cohort, "CRC_vs_Control"]
                aden = lookup[profiler, feature, cohort, "Adenoma_vs_Control"]
                for j, (tag, row) in enumerate((("C", crc), ("A", aden))):
                    label_text = (f'{tag}: {row["effect"]:+.2f} [{row["lower_95"]:+.2f},{row["upper_95"]:+.2f}] q={row["q_value"]:.2g}'
                                  if row["evaluable"] else f"{tag}: not evaluable")
                    color = "#117866" if row["significant"] else "#2b3945"
                    txt(parts, x+8, y+22+j*23, label_text, 13, color=color)
                for j, condition in enumerate(CONDITIONS):
                    b = base[profiler, feature, cohort, condition]
                    s = spike.get((profiler, feature, cohort, condition))
                    med = b["positive_abundance_median_percent"]
                    baseline = "none positive" if med == "" else f"positive median {med:.2g}%"
                    response = (f'rec {s["recovery_median"]:.2f} [{s["recovery_q1"]:.2f},{s["recovery_q3"]:.2f}], '
                                f'low {s["fraction_below_half"]:.0%}, det {s["observed_positive"]}/{s["n"]}'
                                if s else "spike unavailable")
                    txt(parts, x+8, y+79+j*33, f'{condition[:3]} {b["positive"]}/{b["n"]}; {baseline}', 12)
                    txt(parts, x+8, y+94+j*33, response, 12, color="#006a9c")
    txt(parts, 28, height-69, "Spike recovery is on each profiler's own scale; compare within profiler. Positive-only medians exclude zeros.", 14)
    txt(parts, 28, height-42, "Good spike recovery does not prove biological absence; a nonsignificant adenoma q does not establish absence.", 14)
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def build(calls, manifest, abundance, endpoints, panel, aliases, outdir):
    if outdir.exists():
        raise ValueError("output exists: " + str(outdir))
    panel_map = panel_aliases(panel, aliases)
    selected = matched_features(panel_map)
    model = matched_calls(calls, selected)
    baseline_samples, baseline = full_baselines(manifest, abundance, selected)
    spike_samples, spike = direct_spikes(endpoints, selected, panel_map)
    missing = [(p, f) for p in PROFILERS for f in selected[p]
               if not any(r["profiler"] == p and r["feature"] == f for r in spike)]
    if missing:
        raise ValueError("matched species lack direct-spike endpoints: " + repr(missing))
    outdir.mkdir(parents=True)
    write_tsv(outdir / "matched_species_models.tsv", CALL_OUT, model)
    write_tsv(outdir / "matched_species_baseline_summary.tsv", BASE_OUT, baseline)
    write_tsv(outdir / "matched_species_spike_summary.tsv", SPIKE_OUT, spike)
    write_tsv(outdir / "matched_species_baseline_samples.tsv",
              ("profiler","feature","cohort","condition","sample_id","abundance_fraction"), baseline_samples)
    write_tsv(outdir / "matched_species_spike_samples.tsv",
              ("profiler","feature","cohort","condition","dose_percent","sample_id",
               "baseline_abundance_fraction","observed_abundance_fraction","recovery_ratio"), spike_samples)
    draw(outdir / "matched_crc_adenoma_species.svg", selected, model, baseline, spike)
    write_tsv(outdir / "source_sha256.tsv", ("source", "sha256"),
              [dict(source=str(p.resolve()), sha256=sha256(p)) for p in
               (calls, manifest, abundance, endpoints, panel, aliases)])
    (outdir / "DEVELOPMENT_ONLY.txt").write_text(
        "Descriptive matched-species comparison; no direct CRC-vs-Adenoma contrast is inferred.\n",
        encoding="utf-8")
    (outdir / "SUCCESS").write_text("status\tPASS\n", encoding="utf-8")
    print(f"[PASS] three matched species across both profilers and three cohorts: {outdir}")


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
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
