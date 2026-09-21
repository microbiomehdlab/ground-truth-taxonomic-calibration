#!/usr/bin/env python3
"""Descriptive directed CRC biomarker replication from frozen baseline calls.

Only unspiked original-assembly community CRC-vs-Control primary fits enter.
Missing validation features are excluded from the evaluable denominator rather
than silently counted as nonsignificant. This is DEVELOPMENT_ONLY.
"""
import argparse
import csv
import hashlib
import html
import math
from pathlib import Path

COHORTS = ("feng", "yachida", "zeller")
PROFILERS = ("kraken2_bracken", "metaphlan4")
REQUIRED = {"cohort", "analysis_population", "assembly_arm", "profiler", "feature",
            "contrast", "spike_fraction_target", "effect", "q_value",
            "include", "model_spec", "n_samples"}


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def build(calls, outdir):
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
                row["contrast"] != "CRC_vs_Control" or
                row["include"] != "1" or
                row["model_spec"] != "primary_age_sex"):
                continue
            try:
                dose = float(row["spike_fraction_target"])
            except ValueError:
                raise ValueError("nonnumeric dose in eligible disease row")
            if dose != 0:
                continue
            cohort, profiler, feature = row["cohort"], row["profiler"], row["feature"].strip()
            if cohort not in COHORTS or profiler not in PROFILERS or not feature:
                raise ValueError("unexpected cohort/profiler/feature in baseline calls")
            try:
                q, effect, n = float(row["q_value"]), float(row["effect"]), int(row["n_samples"])
            except ValueError:
                raise ValueError("nonnumeric baseline q/effect/sample count")
            if not (math.isfinite(q) and 0 <= q <= 1 and math.isfinite(effect) and n > 0):
                raise ValueError("invalid baseline q/effect/sample count")
            key = cohort, profiler, feature
            value = q, effect, n
            if key in fits and fits[key] != value:
                raise ValueError("contradictory repeated baseline fit: {}".format(key))
            fits[key] = value
    for cohort in COHORTS:
        for profiler in PROFILERS:
            if not any(k[0] == cohort and k[1] == profiler for k in fits):
                raise ValueError("missing baseline fits for {}/{}".format(cohort, profiler))
    rows = []
    for profiler in PROFILERS:
        for discovery in COHORTS:
            found = {feature: value for (cohort, tool, feature), value in fits.items()
                     if cohort == discovery and tool == profiler and value[0] <= .05}
            for validation in COHORTS:
                if validation == discovery:
                    rows.append(dict(profiler=profiler, discovery=discovery,
                                     validation=validation, discovered=len(found),
                                     evaluable=len(found), replicated=len(found),
                                     directionally_discordant=0, rate=1.0 if found else "NA"))
                    continue
                valid = {feature: value for (cohort, tool, feature), value in fits.items()
                         if cohort == validation and tool == profiler}
                common = found.keys() & valid.keys()
                replicated = sum(valid[f][0] <= .05 and found[f][1] * valid[f][1] > 0
                                 for f in common)
                discordant = sum(valid[f][0] <= .05 and found[f][1] * valid[f][1] < 0
                                 for f in common)
                rows.append(dict(profiler=profiler, discovery=discovery,
                                 validation=validation, discovered=len(found),
                                 evaluable=len(common), replicated=replicated,
                                 directionally_discordant=discordant,
                                 rate=replicated / len(common) if common else "NA"))
    outdir.mkdir(parents=True)
    table = outdir / "directional_replication_matrix.tsv"
    with table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="1180" height="650" viewBox="0 0 1180 650">',
             '<rect width="1180" height="650" fill="white"/>',
             '<text x="40" y="44" font-size="26" font-family="sans-serif">Directional CRC biomarker replication</text>',
             '<text x="40" y="70" font-size="15" font-family="sans-serif">Unspiked community baselines · partial cohorts · BH q ≤ 0.05 · DEVELOPMENT ONLY</text>']
    for pi, profiler in enumerate(PROFILERS):
        left = 65 + pi * 565
        label = "Kraken2 + Bracken" if pi == 0 else "MetaPhlAn 4"
        parts.append('<text x="{}" y="118" font-size="20" font-family="sans-serif">{}</text>'.format(left, label))
        for ci, cohort in enumerate(COHORTS):
            parts.append('<text x="{}" y="147" font-size="14" font-family="sans-serif">{}</text>'.format(
                left + ci * 160 + 48, cohort.title()))
            parts.append('<text x="{}" y="{}" font-size="14" font-family="sans-serif">{}</text>'.format(
                left - 55, 240 + ci * 145, cohort.title()))
        for di, discovery in enumerate(COHORTS):
            for vi, validation in enumerate(COHORTS):
                row = next(r for r in rows if r["profiler"] == profiler and
                           r["discovery"] == discovery and r["validation"] == validation)
                x, y = left + vi * 160, 165 + di * 145
                rate = row["rate"]
                if discovery == validation:
                    fill = "#e7ebee"
                    main = "{} found".format(row["discovered"])
                    detail = "discovery"
                else:
                    p = rate if isinstance(rate, float) else 0
                    fill = "#{:02x}{:02x}{:02x}".format(int(244 - 110*p), int(247 - 45*p), int(224 - 100*p))
                    main = "{:.1f}%".format(100*p) if isinstance(rate, float) else "NA"
                    detail = "{}/{} evaluable".format(row["replicated"], row["evaluable"])
                parts.append('<rect x="{}" y="{}" width="150" height="132" fill="{}" stroke="white"/>'.format(x,y,fill))
                parts.append('<text x="{}" y="{}" font-size="20" font-family="sans-serif">{}</text>'.format(x+18,y+58,html.escape(main)))
                parts.append('<text x="{}" y="{}" font-size="14" font-family="sans-serif">{}</text>'.format(x+18,y+84,html.escape(detail)))
    parts.append('<text x="50" y="630" font-size="14" font-family="sans-serif">Rows: discovery cohort; columns: independent validation cohort. Missing features are not evaluable.</text></svg>\n')
    svg = outdir / "directional_replication_matrix.svg"
    svg.write_text("".join(parts), encoding="utf-8")
    (outdir / "source.sha256").write_text(
        "{}  {}\n{}  {}\n".format(digest(calls), calls.resolve(), digest(table), table.resolve()),
        encoding="utf-8")
    (outdir / "DEVELOPMENT_ONLY.txt").write_text("status\tDEVELOPMENT_ONLY\n")
    (outdir / "SUCCESS").write_text("status\tPASS\n", encoding="utf-8")
    print("[PASS] directional replication matrix: {}".format(outdir))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calls", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    try:
        build(args.calls, args.outdir)
    except (OSError, ValueError) as error:
        raise SystemExit("[ERROR] {}".format(error)) from error
