#!/usr/bin/env python3
"""Measure changes in observed disease-biomarker call sets across spike doses."""
from __future__ import annotations
import argparse, csv, hashlib, math, statistics
from collections import defaultdict
from pathlib import Path

KEY = ["cohort", "study", "analysis_population", "target_label", "assembly_arm", "profiler", "contrast"]
REQUIRED = KEY + ["spike_fraction_target", "feature", "effect", "p_value", "q_value", "include", "exclusion_reason"]
OUTPUT = KEY + ["spike_fraction_target", "q_threshold", "target_alias", "baseline_biomarkers",
 "dose_biomarkers", "retained_biomarkers", "lost_biomarkers", "gained_biomarkers",
 "baseline_retention_rate", "dose_overlap_fraction", "biomarker_set_jaccard_vs_baseline",
 "direction_flips_among_retained", "median_abs_effect_change_baseline_biomarkers",
 "max_abs_effect_change_baseline_biomarkers", "target_significant", "target_effect",
 "target_q_value", "target_effect_change_from_baseline"]

def read(path, delimiter="\t"):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None: raise ValueError("empty table: {}".format(path))
        return list(reader), set(reader.fieldnames)

def number(row, field):
    try: value = float(row[field])
    except ValueError as error: raise ValueError("invalid {} for {}".format(field, row.get("feature"))) from error
    if not math.isfinite(value): raise ValueError("non-finite {}".format(field))
    return value

def render(value): return "NA" if value is None else format(value, ".17g")
def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""): value.update(block)
    return value.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calls", type=Path, required=True); parser.add_argument("--aliases", type=Path, required=True)
    parser.add_argument("--spike-panel", type=Path, required=True); parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--q-thresholds", default="0.05,0.10"); args = parser.parse_args()
    try:
        rows, fields = read(args.calls)
        if not set(REQUIRED) <= fields: raise ValueError("disease results lack required columns")
        alias_rows, alias_fields = read(args.aliases, ","); panel, panel_fields = read(args.spike_panel)
        if not {"canonical", "alias", "tool"} <= alias_fields or not {"label", "taxon_name"} <= panel_fields:
            raise ValueError("invalid alias or spike-panel table")
        thresholds = [float(x) for x in args.q_thresholds.split(",")]
        if any(not 0 < x < 1 for x in thresholds) or len(set(thresholds)) != len(thresholds): raise ValueError("invalid q thresholds")
        taxa = {x["label"]: x["taxon_name"] for x in panel}
        alias_lookup = {(x["canonical"], x["tool"]): x["alias"] for x in alias_rows}
        target_alias = {(label, profiler): alias_lookup[(taxon, profiler)] for label, taxon in taxa.items()
                        for profiler in ("kraken2_bracken", "metaphlan4") if (taxon, profiler) in alias_lookup}
        included=[]; identities=set()
        for row in rows:
            if row["include"] not in {"0", "1"} or ((row["include"] == "0") != bool(row["exclusion_reason"].strip())): raise ValueError("invalid include/exclusion encoding")
            dose=number(row,"spike_fraction_target"); effect=number(row,"effect"); q=number(row,"q_value")
            if dose < 0 or not 0 <= q <= 1: raise ValueError("dose or q-value outside range")
            identity=tuple(row[x] for x in KEY)+(render(dose),row["feature"])
            if identity in identities: raise ValueError("duplicate disease feature row")
            identities.add(identity)
            if row["include"] == "1": included.append(dict(row,_dose=dose,_effect=effect,_q=q))
        groups=defaultdict(list)
        for row in included: groups[tuple(row[x] for x in KEY)+(row["_dose"],)].append(row)
        baselines={key[:-1]:value for key,value in groups.items() if key[-1] == 0}; output=[]
        for group_key, dose_rows in sorted(groups.items(), key=lambda x: tuple(map(str,x[0]))):
            dose=group_key[-1]
            if dose == 0: continue
            context=group_key[:-1]; baseline=baselines.get(context)
            if baseline is None: raise ValueError("positive-dose context lacks baseline")
            label=context[KEY.index("target_label")]; profiler=context[KEY.index("profiler")]; target=target_alias.get((label,profiler))
            if target is None: raise ValueError("missing target alias for {} {}".format(label,profiler))
            dose_by={x["feature"]:x for x in dose_rows}; base_by={x["feature"]:x for x in baseline}
            if target not in dose_by or target not in base_by: raise ValueError("target absent from disease result context")
            if set(dose_by) != set(base_by): raise ValueError("feature universe differs between baseline and dose")
            for threshold in thresholds:
                base_calls={f for f,x in base_by.items() if x["_q"] <= threshold}; dose_calls={f for f,x in dose_by.items() if x["_q"] <= threshold}
                retained=base_calls & dose_calls; lost=base_calls-dose_calls; gained=dose_calls-base_calls; union=base_calls|dose_calls
                changes=[abs(dose_by[f]["_effect"]-base_by[f]["_effect"]) for f in base_calls]
                flips=sum(1 for f in retained if dose_by[f]["_effect"]*base_by[f]["_effect"] < 0)
                target_row=dose_by[target]; base_target=base_by[target]
                record={field:value for field,value in zip(KEY,context)}
                record.update(spike_fraction_target=render(dose),q_threshold=render(threshold),target_alias=target,
                  baseline_biomarkers=str(len(base_calls)),dose_biomarkers=str(len(dose_calls)),retained_biomarkers=str(len(retained)),
                  lost_biomarkers=str(len(lost)),gained_biomarkers=str(len(gained)),
                  baseline_retention_rate=render(len(retained)/len(base_calls) if base_calls else None),
                  dose_overlap_fraction=render(len(retained)/len(dose_calls) if dose_calls else None),
                  biomarker_set_jaccard_vs_baseline=render(len(retained)/len(union) if union else 1.0),direction_flips_among_retained=str(flips),
                  median_abs_effect_change_baseline_biomarkers=render(statistics.median(changes) if changes else None),
                  max_abs_effect_change_baseline_biomarkers=render(max(changes) if changes else None),
                  target_significant=str(int(target_row["_q"] <= threshold)),target_effect=render(target_row["_effect"]),
                  target_q_value=render(target_row["_q"]),target_effect_change_from_baseline=render(target_row["_effect"]-base_target["_effect"]))
                output.append(record)
        if not output: raise ValueError("no positive-dose disease contexts")
        args.outdir.mkdir(parents=True,exist_ok=True); result=args.outdir/"disease_biomarker_propagation_metrics.tsv"
        with result.open("w",newline="",encoding="utf-8") as handle:
            writer=csv.DictWriter(handle,fieldnames=OUTPUT,delimiter="\t",lineterminator="\n"); writer.writeheader(); writer.writerows(output)
        summary=args.outdir/"disease_biomarker_propagation_summary.tsv"
        summary.write_text("metric\tvalue\ninput_rows\t{}\nevaluated_rows\t{}\nstatus\tPASS\n".format(len(rows),len(output)),encoding="utf-8")
        (args.outdir/"disease_biomarker_propagation.sha256").write_text("".join("{}  {}\n".format(digest(p),p.resolve()) for p in (args.calls,args.aliases,args.spike_panel,result,summary)),encoding="utf-8")
        (args.outdir/"SUCCESS").write_text("evaluated_rows\t{}\nstatus\tPASS\n".format(len(output)),encoding="utf-8")
        print("[PASS] Evaluated {} disease-biomarker propagation rows".format(len(output)))
    except (ValueError,FileNotFoundError,KeyError) as error: raise SystemExit("[ERROR] {}".format(error)) from error
if __name__ == "__main__": main()
