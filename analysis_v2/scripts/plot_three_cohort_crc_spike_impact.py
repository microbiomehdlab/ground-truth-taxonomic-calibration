#!/usr/bin/env python3
"""Baseline CRC calls and their community-mixture spike fate, for both profilers."""
import argparse
import csv
import hashlib
import html
import math
from pathlib import Path

COHORTS = ("feng", "yachida", "zeller")
PROFILERS = ("kraken2_bracken", "metaphlan4")
CALL_FIELDS = {"cohort", "analysis_population", "assembly_arm", "profiler", "contrast",
               "include", "model_spec", "spike_fraction_target", "feature", "effect",
               "q_value", "n_samples"}
LEDGER_FIELDS = {"cohort", "analysis_population", "assembly_arm", "profiler", "contrast",
                 "target_label", "spike_fraction_target", "spike_fraction_total",
                 "q_threshold", "feature", "feature_role", "baseline_called",
                 "dose_called", "baseline_q_value", "effect_sign_changed"}


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


def read_calls(path):
    fits = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not CALL_FIELDS <= set(reader.fieldnames or ()):
            raise ValueError("disease calls lack required columns")
        for row in reader:
            if (row["analysis_population"] != "community" or
                row["assembly_arm"] != "original" or
                row["profiler"] not in PROFILERS or
                row["contrast"] != "CRC_vs_Control" or row["include"] != "1" or
                row["model_spec"] != "primary_age_sex"):
                continue
            if number(row["spike_fraction_target"], "baseline dose") != 0:
                continue
            cohort, profiler, feature = row["cohort"], row["profiler"], row["feature"].strip()
            if cohort not in COHORTS or not feature:
                raise ValueError("unexpected cohort or empty feature")
            effect, q = number(row["effect"], "effect"), number(row["q_value"], "q")
            n = number(row["n_samples"], "sample count")
            if not 0 <= q <= 1 or n <= 0 or not n.is_integer():
                raise ValueError("invalid q or sample count")
            key, value = (cohort, profiler, feature), (effect, q, int(n))
            if key in fits and fits[key] != value:
                raise ValueError("contradictory repeated baseline fit: {}".format(key))
            fits[key] = value
    if any(not any(k[0] == c and k[1] == p for k in fits)
           for c in COHORTS for p in PROFILERS):
        raise ValueError("one or more cohort/profiler baseline strata missing")
    return fits


def read_ledger(path, dose_percent):
    chosen = {}
    totals = set()
    context_seen = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not LEDGER_FIELDS <= set(reader.fieldnames or ()):
            raise ValueError("transition ledger lacks required columns")
        for row in reader:
            if (row["analysis_population"] != "community" or
                row["assembly_arm"] != "original" or
                row["profiler"] not in PROFILERS or
                row["contrast"] != "CRC_vs_Control" or
                row["target_label"] != "CRCpanel"):
                continue
            dose = number(row["spike_fraction_target"], "member dose")
            if not math.isclose(dose, dose_percent / 100, rel_tol=.05):
                continue
            if not math.isclose(number(row["q_threshold"], "q threshold"), .05,
                                abs_tol=1e-9):
                continue
            cohort, profiler, feature = row["cohort"], row["profiler"], row["feature"].strip()
            if cohort not in COHORTS or not feature:
                raise ValueError("unexpected cohort or empty feature in ledger")
            total = number(row["spike_fraction_total"], "total mixture dose")
            if total <= 0:
                raise ValueError("invalid total mixture dose")
            totals.add(round(total, 12))
            context_seen.add((cohort, profiler))
            key = cohort, profiler, feature
            if key in chosen:
                raise ValueError("duplicate selected community challenge row: {}".format(key))
            chosen[key] = row
    if context_seen != {(c, p) for c in COHORTS for p in PROFILERS}:
        raise ValueError("selected community challenge missing cohort/profiler contexts")
    if len(totals) != 1:
        raise ValueError("selected community mixtures have different total fractions")
    return chosen, totals.pop()


def rank(fits, profiler):
    features = {f for c, p, f in fits if p == profiler}
    candidates = [f for f in features if any((c, profiler, f) in fits and
                                            fits[c, profiler, f][1] <= .05 for c in COHORTS)]
    def score(feature):
        observed = [fits[c, profiler, feature] for c in COHORTS if (c, profiler, feature) in fits]
        sig = [v for v in observed if v[1] <= .05]
        shared = max(sum(v[0] > 0 for v in sig), sum(v[0] < 0 for v in sig))
        same_direction = len(observed) == 3 and (all(v[0] > 0 for v in observed) or
                                                  all(v[0] < 0 for v in observed))
        return shared, same_direction, len(observed), min(v[1] for v in sig)
    return sorted(candidates, key=lambda f: (-score(f)[0], -int(score(f)[1]),
                                             -score(f)[2], score(f)[3], f)), score


def fate(fit, row):
    if fit is None:
        return "not_evaluable"
    effect, q, _ = fit
    if row is not None:
        ledger_q = number(row["baseline_q_value"], "ledger baseline q")
        if not math.isclose(ledger_q, q, rel_tol=1e-8, abs_tol=1e-10):
            raise ValueError("ledger and baseline q disagree for {}".format(row["feature"]))
        if row["baseline_called"] != str(int(q <= .05)):
            raise ValueError("ledger and baseline call disagree for {}".format(row["feature"]))
    if row is None:
        if q <= .05:
            raise ValueError("missing spike fate for significant baseline feature")
        return "not_called"
    if row["feature_role"] == "implanted_target":
        return "direct_target_excluded"
    if row["feature_role"] != "bystander" or row["dose_called"] not in ("0", "1") or row["effect_sign_changed"] not in ("0", "1"):
        raise ValueError("invalid bystander spike fate")
    if q > .05:
        return "gained_after_spike" if row["dose_called"] == "1" else "not_called"
    if row["dose_called"] == "1" and row["effect_sign_changed"] == "1":
        return "direction_reversed"
    return "retained_same_direction" if row["dose_called"] == "1" else "lost_significance"


COLORS = {
    "retained_same_direction": "#168a71", "lost_significance": "#dfa51f",
    "direction_reversed": "#b72878", "gained_after_spike": "#5377b5",
    "not_called": "#eef1f3", "not_evaluable": "#dce1e7",
    "direct_target_excluded": "#a9adb4",
}


def draw(outdir, profiler, rows, total_fraction, dose_percent):
    width, height = 1500, 200 + 63 * (len(rows) // 3) + 100
    label = "Kraken2 + Bracken" if profiler == "kraken2_bracken" else "MetaPhlAn 4"
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}">'.format(width, height, width, height),
             '<rect width="100%" height="100%" fill="white"/>',
             '<text x="30" y="43" font-family="sans-serif" font-size="29" font-weight="bold">Which {} CRC calls survive a known spike?</text>'.format(html.escape(label)),
             '<text x="30" y="74" font-family="sans-serif" font-size="17">Unspiked community CRC calls; {:.3f}% per target, {:.3f}% total community mixture · BH q ≤ 0.05 · DEVELOPMENT ONLY</text>'.format(dose_percent, total_fraction * 100),
             '<text x="30" y="106" font-family="sans-serif" font-size="16">Each cell: original signed effect and q; colour: fate after the controlled community spike.</text>']
    for ci, cohort in enumerate(COHORTS):
        parts.append('<text x="{}" y="159" font-family="sans-serif" font-size="21" text-anchor="middle" font-weight="bold">{}</text>'.format(650 + ci * 265, cohort.title()))
    for ri in range(len(rows) // 3):
        feature = rows[ri * 3]["feature"]
        y = 181 + ri * 63
        short = feature if len(feature) <= 54 else feature[:51] + "..."
        parts.append('<title>{}</title><text x="30" y="{}" font-family="sans-serif" font-size="17">{}</text>'.format(
            html.escape(feature), y + 28, html.escape(short)))
        for ci, row in enumerate(rows[ri * 3:ri * 3 + 3]):
            x = 530 + ci * 265
            status = row["spike_fate"]
            effect = row["baseline_effect"]
            q = row["baseline_q"]
            text = "not evaluable" if not effect else "{:+.2f}{}  q={:.2g}".format(
                float(effect), "*" if row["baseline_called"] == 1 else "", float(q))
            parts.append('<rect x="{}" y="{}" width="250" height="49" rx="4" fill="{}"/>'.format(x, y, COLORS[status]))
            parts.append('<text x="{}" y="{}" font-family="sans-serif" font-size="17" fill="{}">{}</text>'.format(
                x + 10, y + 30, "white" if status in ("retained_same_direction", "direction_reversed", "gained_after_spike") else "#26313b", html.escape(text)))
    footer = height - 61
    parts.append('<text x="30" y="{}" font-family="sans-serif" font-size="15">Green retained · gold lost q≤0.05 · magenta reversed · blue gained · grey direct implanted target excluded · pale not called/evaluable.</text>'.format(footer))
    parts.append('<text x="30" y="{}" font-family="sans-serif" font-size="15">* baseline BH q≤0.05. A changed call indicates technical perturbation sensitivity, not proof of a biological false positive.</text>'.format(footer + 26))
    parts.append('</svg>\n')
    (outdir / (profiler + "_top_candidates.svg")).write_text("".join(parts), encoding="utf-8")


def build(calls_path, ledger_path, outdir, dose_percent, top):
    if outdir.exists():
        raise ValueError("output exists: {}".format(outdir))
    fits = read_calls(calls_path)
    ledger, total = read_ledger(ledger_path, dose_percent)
    rows, summary = [], []
    for profiler in PROFILERS:
        ranked, score = rank(fits, profiler)
        if not ranked:
            raise ValueError("no baseline CRC candidates for {}".format(profiler))
        shown = ranked[:top]
        profiler_rows = []
        for feature in shown:
            shared, same_direction, _, _ = score(feature)
            for cohort in COHORTS:
                fit = fits.get((cohort, profiler, feature))
                status = fate(fit, ledger.get((cohort, profiler, feature)))
                record = dict(profiler=profiler, rank=len(profiler_rows) // 3 + 1,
                              feature=feature, cohort=cohort,
                              baseline_effect="" if fit is None else "{:.9g}".format(fit[0]),
                              baseline_q="" if fit is None else "{:.9g}".format(fit[1]),
                              baseline_called=int(fit is not None and fit[1] <= .05),
                              spike_fate=status, same_direction_significant_cohorts=shared,
                              same_direction_all_three=int(same_direction),
                              member_dose_percent=dose_percent,
                              total_mixture_percent=100 * total)
                profiler_rows.append(record)
        rows.extend(profiler_rows)
        summary.extend([dict(profiler=profiler, metric="baseline_candidates_any_cohort", value=len(ranked)),
                        dict(profiler=profiler, metric="shared_significant_direction_all_three", value=sum(score(f)[0] == 3 for f in ranked)),
                        dict(profiler=profiler, metric="plotted_candidates", value=len(shown)),
                        dict(profiler=profiler, metric="plotted_baseline_calls_retained", value=sum(r["spike_fate"] == "retained_same_direction" for r in profiler_rows)),
                        dict(profiler=profiler, metric="plotted_baseline_calls_lost", value=sum(r["spike_fate"] == "lost_significance" for r in profiler_rows)),
                        dict(profiler=profiler, metric="plotted_baseline_calls_direct_target_excluded", value=sum(r["spike_fate"] == "direct_target_excluded" for r in profiler_rows))])
    outdir.mkdir(parents=True)
    write_tsv(outdir / "top_candidate_spike_fates.tsv", rows)
    write_tsv(outdir / "candidate_summary.tsv", summary)
    for profiler in PROFILERS:
        draw(outdir, profiler, [r for r in rows if r["profiler"] == profiler], total, dose_percent)
    (outdir / "source.sha256").write_text("{}  {}\n{}  {}\n".format(
        digest(calls_path), calls_path.resolve(), digest(ledger_path), ledger_path.resolve()), encoding="utf-8")
    (outdir / "DEVELOPMENT_ONLY.txt").write_text("status\tDEVELOPMENT_ONLY\n", encoding="utf-8")
    (outdir / "SUCCESS").write_text("status\tPASS\n", encoding="utf-8")
    print("[PASS] top {} baseline candidates per profiler, three cohorts, community spike fate".format(top))
    print(outdir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calls", type=Path, required=True)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--dose-percent", type=float, default=.1)
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.dose_percent) or args.dose_percent <= 0 or args.top < 1:
        parser.error("dose percent and top must be positive")
    try:
        build(args.calls, args.ledger, args.outdir, args.dose_percent, args.top)
    except (OSError, ValueError) as error:
        raise SystemExit("[ERROR] {}".format(error)) from error
