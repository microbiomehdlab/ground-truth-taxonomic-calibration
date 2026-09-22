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
                 "dose_called", "baseline_q_value", "dose_q_value",
                 "effect_sign_changed"}


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
    totals = {}
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
            context = cohort, profiler
            if context in totals and not math.isclose(totals[context], total, rel_tol=1e-9,
                                                       abs_tol=1e-12):
                raise ValueError("selected community mixture has inconsistent total fraction within {}".format(context))
            totals[context] = total
            context_seen.add(context)
            key = cohort, profiler, feature
            if key in chosen:
                raise ValueError("duplicate selected community challenge row: {}".format(key))
            dose_q = number(row["dose_q_value"], "post-spike q")
            if not 0 <= dose_q <= 1:
                raise ValueError("post-spike q outside [0,1]")
            if row["dose_called"] not in ("0", "1") or row["dose_called"] != str(int(dose_q <= .05)):
                raise ValueError("post-spike q and call flag disagree for {}".format(key))
            chosen[key] = row
    if context_seen != {(c, p) for c in COHORTS for p in PROFILERS}:
        raise ValueError("selected community challenge missing cohort/profiler contexts")
    return chosen, totals


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


def draw(outdir, profiler, rows, totals, dose_percent):
    width, height = 1500, 200 + 63 * (len(rows) // 3) + 100
    label = "Kraken2 + Bracken" if profiler == "kraken2_bracken" else "MetaPhlAn 4"
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}">'.format(width, height, width, height),
             '<rect width="100%" height="100%" fill="white"/>',
             '<text x="30" y="43" font-family="sans-serif" font-size="29" font-weight="bold">Which {} CRC calls survive a known spike?</text>'.format(html.escape(label)),
             '<text x="30" y="74" font-family="sans-serif" font-size="17">Unspiked community CRC calls; {:.3f}% median member dose; cohort-specific total mixture shown below · BH q ≤ 0.05 · DEVELOPMENT ONLY</text>'.format(dose_percent),
             '<text x="30" y="106" font-family="sans-serif" font-size="16">Each cell: original signed effect and q; colour: fate after the controlled community spike.</text>']
    for ci, cohort in enumerate(COHORTS):
        parts.append('<text x="{}" y="159" font-family="sans-serif" font-size="20" text-anchor="middle" font-weight="bold">{} (total {:.3f}%)</text>'.format(
            650 + ci * 265, cohort.title(), 100 * totals[cohort, profiler]))
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


def draw_shared(outdir, profiler, rows, totals, dose_percent):
    width, row_h = 1650, 79
    height = 218 + row_h * (len(rows) // 3) + 104
    label = "Kraken2 + Bracken" if profiler == "kraken2_bracken" else "MetaPhlAn 4"
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}">'.format(width, height, width, height),
             '<rect width="100%" height="100%" fill="white"/>',
             '<text x="30" y="44" font-family="sans-serif" font-size="30" font-weight="bold">Three-cohort {} CRC candidates: before and after a known spike</text>'.format(html.escape(label)),
             '<text x="30" y="77" font-family="sans-serif" font-size="17">Original unspiked community calls; {:.3f}% median member dose; one physical community mixture per cohort · DEVELOPMENT ONLY</text>'.format(dose_percent),
             '<text x="30" y="108" font-family="sans-serif" font-size="16">Each cohort: original CRC effect and q (left) | fate and post-spike q (right). Only same-direction BH q ≤ 0.05 in all three shown.</text>']
    for ci, cohort in enumerate(COHORTS):
        x = 500 + ci * 378
        parts.append('<text x="{}" y="159" font-family="sans-serif" font-size="20" text-anchor="middle" font-weight="bold">{} · total {:.3f}%</text>'.format(
            x + 174, cohort.title(), 100 * totals[cohort, profiler]))
        parts.append('<text x="{}" y="184" font-family="sans-serif" font-size="15" text-anchor="middle">Original CRC call</text>'.format(x + 91))
        parts.append('<text x="{}" y="184" font-family="sans-serif" font-size="15" text-anchor="middle">After spike</text>'.format(x + 268))
    status_label = {"retained_same_direction": "retained", "lost_significance": "lost",
                    "direction_reversed": "direction reversed",
                    "direct_target_excluded": "direct target excluded"}
    for ri in range(len(rows) // 3):
        feature = rows[3 * ri]["feature"]
        y = 202 + ri * row_h
        short = feature if len(feature) <= 53 else feature[:50] + "..."
        parts.append('<title>{}</title><text x="30" y="{}" font-family="sans-serif" font-size="19">{}</text>'.format(
            html.escape(feature), y + 37, html.escape(short)))
        for ci, row in enumerate(rows[3 * ri:3 * ri + 3]):
            x = 500 + ci * 378
            effect, q = float(row["baseline_effect"]), float(row["baseline_q"])
            before_fill = "#176c70" if effect > 0 else "#a84c62"
            status = row["spike_fate"]
            after_q = row["post_spike_q"]
            after_text = status_label.get(status, status.replace("_", " "))
            if status != "direct_target_excluded" and after_q:
                after_text += " · q={:.2g}".format(float(after_q))
            parts.append('<rect x="{}" y="{}" width="177" height="55" rx="4" fill="{}"/>'.format(x, y, before_fill))
            parts.append('<rect x="{}" y="{}" width="173" height="55" rx="4" fill="{}"/>'.format(x + 181, y, COLORS[status]))
            parts.append('<text x="{}" y="{}" font-family="sans-serif" font-size="15" fill="white">{:+.2f} · q={:.2g}</text>'.format(x + 7, y + 34, effect, q))
            parts.append('<text x="{}" y="{}" font-family="sans-serif" font-size="14" fill="{}">{}</text>'.format(
                x + 187, y + 34, "white" if status in ("retained_same_direction", "direction_reversed") else "#25303a", html.escape(after_text)))
    footer = height - 65
    parts.append('<text x="30" y="{}" font-family="sans-serif" font-size="15">Original: teal positive effect, rose negative effect. After spike: green retained, gold lost significance, magenta reversed, grey direct target excluded.</text>'.format(footer))
    parts.append('<text x="30" y="{}" font-family="sans-serif" font-size="15">Directly implanted species cannot be evaluated as bystanders. Loss of significance is technical sensitivity, not proof of a biological false positive.</text>'.format(footer + 26))
    parts.append('</svg>\n')
    (outdir / (profiler + "_shared_candidates.svg")).write_text("".join(parts), encoding="utf-8")


def build(calls_path, ledger_path, outdir, dose_percent, top):
    if outdir.exists():
        raise ValueError("output exists: {}".format(outdir))
    fits = read_calls(calls_path)
    ledger, totals = read_ledger(ledger_path, dose_percent)
    rows, shared_rows, summary = [], [], []
    for profiler in PROFILERS:
        ranked, score = rank(fits, profiler)
        if not ranked:
            raise ValueError("no baseline CRC candidates for {}".format(profiler))
        shown = ranked[:top]
        shared_features = [feature for feature in ranked if score(feature)[0] == 3]
        profiler_rows, profiler_shared_rows = [], []
        for feature in dict.fromkeys(shown + shared_features):
            shared_count, same_direction, _, _ = score(feature)
            feature_rows = []
            for cohort in COHORTS:
                fit = fits.get((cohort, profiler, feature))
                challenge = ledger.get((cohort, profiler, feature))
                status = fate(fit, challenge)
                record = dict(profiler=profiler, rank=ranked.index(feature) + 1,
                              feature=feature, cohort=cohort,
                              baseline_effect="" if fit is None else "{:.9g}".format(fit[0]),
                              baseline_q="" if fit is None else "{:.9g}".format(fit[1]),
                              baseline_called=int(fit is not None and fit[1] <= .05),
                              spike_fate=status,
                              post_spike_q="" if challenge is None or status == "direct_target_excluded" else challenge["dose_q_value"],
                              same_direction_significant_cohorts=shared_count,
                              same_direction_all_three=int(same_direction),
                              member_dose_percent=dose_percent,
                              total_mixture_percent=100 * totals[cohort, profiler])
                feature_rows.append(record)
            if feature in shown:
                profiler_rows.extend(feature_rows)
            if shared_count == 3:
                profiler_shared_rows.extend(feature_rows)
        rows.extend(profiler_rows)
        shared_rows.extend(profiler_shared_rows)
        summary.extend([dict(profiler=profiler, metric="baseline_candidates_any_cohort", value=len(ranked)),
                        dict(profiler=profiler, metric="shared_significant_direction_all_three", value=sum(score(f)[0] == 3 for f in ranked)),
                        dict(profiler=profiler, metric="plotted_candidates", value=len(shown)),
                        dict(profiler=profiler, metric="plotted_baseline_calls_retained", value=sum(r["spike_fate"] == "retained_same_direction" for r in profiler_rows)),
                        dict(profiler=profiler, metric="plotted_baseline_calls_lost", value=sum(r["spike_fate"] == "lost_significance" for r in profiler_rows)),
                        dict(profiler=profiler, metric="plotted_baseline_calls_direct_target_excluded", value=sum(r["spike_fate"] == "direct_target_excluded" for r in profiler_rows))])
    outdir.mkdir(parents=True)
    write_tsv(outdir / "top_candidate_spike_fates.tsv", rows)
    write_tsv(outdir / "shared_candidate_spike_fates.tsv", shared_rows)
    write_tsv(outdir / "candidate_summary.tsv", summary)
    for profiler in PROFILERS:
        draw(outdir, profiler, [r for r in rows if r["profiler"] == profiler], totals, dose_percent)
        draw_shared(outdir, profiler, [r for r in shared_rows if r["profiler"] == profiler], totals, dose_percent)
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
