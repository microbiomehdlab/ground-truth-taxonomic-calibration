#!/usr/bin/env python3
"""Measure changes in observed disease-biomarker call sets across spike doses."""
from __future__ import annotations
import argparse, csv, hashlib, math, statistics
from collections import defaultdict
from pathlib import Path

KEY = ["cohort", "study", "analysis_population", "target_label", "assembly_arm", "profiler", "contrast"]
REQUIRED = KEY + ["dose_level", "spike_fraction_target", "feature", "effect", "p_value", "q_value", "include", "exclusion_reason"]
OUTPUT = KEY + ["spike_fraction_target", "spike_fraction_total", "q_threshold", "target_alias", "baseline_biomarkers",
 "dose_biomarkers", "retained_biomarkers", "lost_biomarkers", "gained_biomarkers",
 "baseline_retention_rate", "dose_overlap_fraction", "biomarker_set_jaccard_vs_baseline",
 "baseline_bystander_biomarkers", "dose_bystander_biomarkers", "retained_bystanders",
 "lost_bystanders", "gained_bystanders", "bystander_retention_rate", "bystander_jaccard",
 "bystander_induced_call_rate",
 "bystander_direction_flips_among_retained",
 "median_abs_effect_change_baseline_bystanders", "max_abs_effect_change_baseline_bystanders",
 "direction_flips_among_retained", "median_abs_effect_change_baseline_biomarkers",
 "max_abs_effect_change_baseline_biomarkers", "target_significant", "target_effect",
 "target_q_value", "target_baseline_significant", "target_baseline_effect",
 "target_baseline_q_value", "target_effect_change_from_baseline"]
LEDGER = KEY + ["spike_fraction_target", "spike_fraction_total", "q_threshold", "target_alias", "feature",
 "feature_role", "baseline_called", "dose_called", "transition", "baseline_effect",
 "dose_effect", "effect_change", "baseline_q_value", "dose_q_value", "effect_sign_changed"]

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

def collapse_community(rows, panel_labels):
    """Collapse ten repeated panel-member fits to one physical mixture.

    The disease model historically emits the same community profile once per
    panel member.  A community stress test is nevertheless one physical
    perturbation, not ten independent observations.  Repeated model results
    must agree exactly on all fitted quantities. Preserve the per-member dose
    used on the figure's x axis, and also record the total mixture fraction.
    """
    independent=[]; buckets=defaultdict(list)
    physical_fields=["cohort","study","analysis_population","assembly_arm",
                     "profiler","contrast","dose_level","feature"]
    for row in rows:
        if row["analysis_population"] != "community":
            independent.append(row); continue
        buckets[tuple(row[x] for x in physical_fields)].append(row)
    collapsed=[]
    compare_fields=["effect","p_value","q_value","include","exclusion_reason"]
    for key, members in buckets.items():
        labels={row["target_label"] for row in members}
        if labels != panel_labels:
            raise ValueError("community context does not contain every implanted target: "
                             "missing={} unexpected={}".format(
                                 sorted(panel_labels-labels), sorted(labels-panel_labels)))
        for field in compare_fields:
            if len({row[field] for row in members}) != 1:
                raise ValueError("community repeated fits disagree on {}".format(field))
        by_label=defaultdict(list)
        for row in members: by_label[row["target_label"]].append(row)
        if any(len(values) != 1 for values in by_label.values()):
            raise ValueError("community context repeats a target label")
        representative=dict(members[0])
        representative["target_label"]="CRCpanel"
        fractions=[number(row,"spike_fraction_target") for row in members]
        if any(f < 0 for f in fractions):
            raise ValueError("negative community member fraction")
        representative["spike_fraction_target"]=render(statistics.median(fractions))
        representative["spike_fraction_total"]=render(sum(fractions))
        collapsed.append(representative)
    return independent+collapsed

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calls", type=Path, required=True); parser.add_argument("--aliases", type=Path, required=True)
    parser.add_argument("--spike-panel", type=Path, required=True); parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--q-thresholds", default="0.05,0.10"); args = parser.parse_args()
    try:
        rows, fields = read(args.calls); input_row_count=len(rows)
        if not set(REQUIRED) <= fields: raise ValueError("disease results lack required columns")
        alias_rows, alias_fields = read(args.aliases, ","); panel, panel_fields = read(args.spike_panel)
        if not {"canonical", "alias", "tool"} <= alias_fields or not {"label", "taxon_name"} <= panel_fields:
            raise ValueError("invalid alias or spike-panel table")
        thresholds = [float(x) for x in args.q_thresholds.split(",")]
        if any(not 0 < x < 1 for x in thresholds) or len(set(thresholds)) != len(thresholds): raise ValueError("invalid q thresholds")
        taxa = {x["label"]: x["taxon_name"] for x in panel}
        if len(taxa) != len(panel) or not taxa: raise ValueError("duplicate or empty spike panel")
        alias_lookup = {(x["canonical"], x["tool"]): x["alias"] for x in alias_rows}
        profilers_present={row["profiler"] for row in rows}
        target_alias = {(label, profiler): alias_lookup[(taxon, profiler)] for label, taxon in taxa.items()
                        for profiler in profilers_present if (taxon, profiler) in alias_lookup}
        canonical_lookup={}
        for label,taxon in taxa.items():
            for profiler in profilers_present:
                alias=target_alias.get((label,profiler))
                if alias is None: raise ValueError("missing target alias for {} {}".format(label,profiler))
                for feature in (taxon,alias):
                    old=canonical_lookup.setdefault((profiler,feature),taxon)
                    if old != taxon: raise ValueError("ambiguous feature alias {} {}".format(profiler,feature))
        # Canonicalization precedes exclusion.  This is essential when a panel
        # member appears under a profiler alias rather than its panel taxon name.
        implanted_canonical=set(taxa.values())
        rows=collapse_community(rows,set(taxa)); collapsed_row_count=len(rows)
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
        baselines={key[:-1]:value for key,value in groups.items() if key[-1] == 0}; output=[]; ledger=[]
        for group_key, dose_rows in sorted(groups.items(), key=lambda x: tuple(map(str,x[0]))):
            dose=group_key[-1]
            if dose == 0: continue
            context=group_key[:-1]; baseline=baselines.get(context)
            if baseline is None: raise ValueError("positive-dose context lacks baseline")
            label=context[KEY.index("target_label")]; profiler=context[KEY.index("profiler")]
            community=context[KEY.index("analysis_population")] == "community"
            target=target_alias.get((label,profiler)) if not community else None
            if not community and target is None: raise ValueError("missing target alias for {} {}".format(label,profiler))
            excluded_canonical=(implanted_canonical if community else {taxa[label]})
            def is_implanted(feature):
                return canonical_lookup.get((profiler,feature),feature) in excluded_canonical
            dose_by={x["feature"]:x for x in dose_rows}; base_by={x["feature"]:x for x in baseline}
            if not community and (target not in dose_by or target not in base_by): raise ValueError("target absent from disease result context")
            if set(dose_by) != set(base_by): raise ValueError("feature universe differs between baseline and dose")
            if community and not all(any(canonical_lookup.get((profiler,f),f) == taxon for f in base_by)
                                     for taxon in implanted_canonical):
                raise ValueError("community feature universe lacks an implanted panel member")
            for threshold in thresholds:
                base_calls={f for f,x in base_by.items() if x["_q"] <= threshold}; dose_calls={f for f,x in dose_by.items() if x["_q"] <= threshold}
                retained=base_calls & dose_calls; lost=base_calls-dose_calls; gained=dose_calls-base_calls; union=base_calls|dose_calls
                base_bystanders={f for f in base_calls if not is_implanted(f)}
                dose_bystanders={f for f in dose_calls if not is_implanted(f)}
                retained_bystanders=base_bystanders & dose_bystanders
                lost_bystanders=base_bystanders-dose_bystanders; gained_bystanders=dose_bystanders-base_bystanders
                bystander_union=base_bystanders|dose_bystanders
                eligible_nonbaseline_bystanders={f for f in base_by if not is_implanted(f)}-base_bystanders
                changes=[abs(dose_by[f]["_effect"]-base_by[f]["_effect"]) for f in base_calls]
                flips=sum(1 for f in retained if dose_by[f]["_effect"]*base_by[f]["_effect"] < 0)
                bystander_changes=[abs(dose_by[f]["_effect"]-base_by[f]["_effect"]) for f in base_bystanders]
                bystander_flips=sum(1 for f in retained_bystanders
                                    if dose_by[f]["_effect"]*base_by[f]["_effect"] < 0)
                target_row=dose_by[target] if not community else None
                base_target=base_by[target] if not community else None
                record={field:value for field,value in zip(KEY,context)}
                record.update(spike_fraction_target=render(dose),
                  spike_fraction_total=dose_rows[0].get("spike_fraction_total", render(dose)),
                  q_threshold=render(threshold),
                  target_alias=target if not community else ";".join(sorted(
                      target_alias[(member,profiler)] for member in taxa)),
                  baseline_biomarkers=str(len(base_calls)),dose_biomarkers=str(len(dose_calls)),retained_biomarkers=str(len(retained)),
                  lost_biomarkers=str(len(lost)),gained_biomarkers=str(len(gained)),
                  baseline_retention_rate=render(len(retained)/len(base_calls) if base_calls else None),
                  dose_overlap_fraction=render(len(retained)/len(dose_calls) if dose_calls else None),
                  biomarker_set_jaccard_vs_baseline=render(len(retained)/len(union) if union else None),direction_flips_among_retained=str(flips),
                  baseline_bystander_biomarkers=str(len(base_bystanders)),dose_bystander_biomarkers=str(len(dose_bystanders)),
                  retained_bystanders=str(len(retained_bystanders)),lost_bystanders=str(len(lost_bystanders)),
                  gained_bystanders=str(len(gained_bystanders)),
                  bystander_retention_rate=render(len(retained_bystanders)/len(base_bystanders) if base_bystanders else None),
                  bystander_jaccard=render(len(retained_bystanders)/len(bystander_union) if bystander_union else None),
                  bystander_induced_call_rate=render(len(gained_bystanders)/len(eligible_nonbaseline_bystanders) if eligible_nonbaseline_bystanders else None),
                  bystander_direction_flips_among_retained=str(bystander_flips),
                  median_abs_effect_change_baseline_bystanders=render(statistics.median(bystander_changes) if bystander_changes else None),
                  max_abs_effect_change_baseline_bystanders=render(max(bystander_changes) if bystander_changes else None),
                  median_abs_effect_change_baseline_biomarkers=render(statistics.median(changes) if changes else None),
                  max_abs_effect_change_baseline_biomarkers=render(max(changes) if changes else None),
                  target_significant="NA" if community else str(int(target_row["_q"] <= threshold)),
                  target_effect=render(None if community else target_row["_effect"]),
                  target_q_value=render(None if community else target_row["_q"]),
                  target_baseline_significant="NA" if community else str(int(base_target["_q"] <= threshold)),
                  target_baseline_effect=render(None if community else base_target["_effect"]),
                  target_baseline_q_value=render(None if community else base_target["_q"]),
                  target_effect_change_from_baseline=render(None if community else target_row["_effect"]-base_target["_effect"]))
                output.append(record)
                # Sparse feature-level ledger: the union of disease calls plus the
                # implanted target. Uncalled bystanders are deliberately omitted.
                required_features={f for f in base_by if is_implanted(f)}
                for feature in sorted(union | required_features):
                    baseline_called=feature in base_calls; dose_called=feature in dose_calls
                    sign_changed=(dose_by[feature]["_effect"]*base_by[feature]["_effect"] < 0)
                    if baseline_called and dose_called:
                        transition="retained_sign_flipped" if sign_changed else "retained"
                    elif baseline_called:
                        transition="lost"
                    elif dose_called:
                        transition="gained"
                    else:
                        transition="direct_target_not_called"
                    entry={field:value for field,value in zip(KEY,context)}
                    entry.update(spike_fraction_target=render(dose),
                      spike_fraction_total=record["spike_fraction_total"],
                      q_threshold=render(threshold),
                      target_alias=record["target_alias"],feature=feature,
                      feature_role="implanted_target" if is_implanted(feature) else "bystander",
                      baseline_called=str(int(baseline_called)),dose_called=str(int(dose_called)),
                      transition=transition,baseline_effect=render(base_by[feature]["_effect"]),
                      dose_effect=render(dose_by[feature]["_effect"]),
                      effect_change=render(dose_by[feature]["_effect"]-base_by[feature]["_effect"]),
                      baseline_q_value=render(base_by[feature]["_q"]),
                      dose_q_value=render(dose_by[feature]["_q"]),
                      effect_sign_changed=str(int(sign_changed)))
                    ledger.append(entry)
        if not output: raise ValueError("no positive-dose disease contexts")
        args.outdir.mkdir(parents=True,exist_ok=True); result=args.outdir/"disease_biomarker_propagation_metrics.tsv"
        with result.open("w",newline="",encoding="utf-8") as handle:
            writer=csv.DictWriter(handle,fieldnames=OUTPUT,delimiter="\t",lineterminator="\n"); writer.writeheader(); writer.writerows(output)
        ledger_path=args.outdir/"disease_biomarker_transition_ledger.tsv"
        with ledger_path.open("w",newline="",encoding="utf-8") as handle:
            writer=csv.DictWriter(handle,fieldnames=LEDGER,delimiter="\t",lineterminator="\n"); writer.writeheader(); writer.writerows(ledger)
        summary=args.outdir/"disease_biomarker_propagation_summary.tsv"
        summary.write_text(
            "metric\tvalue\ninput_rows\t{}\nphysical_rows_after_community_collapse\t{}\n"
            "community_repeated_rows_collapsed\t{}\nevaluated_rows\t{}\n"
            "transition_rows\t{}\ncommunity_target_exclusion\tall_panel_members\n"
            "feature_identity_policy\tcanonicalize_before_exclusion\nstatus\tPASS\n".format(
                input_row_count,collapsed_row_count,input_row_count-collapsed_row_count,
                len(output),len(ledger)),encoding="utf-8")
        (args.outdir/"disease_biomarker_propagation.sha256").write_text("".join("{}  {}\n".format(digest(p),p.resolve()) for p in (args.calls,args.aliases,args.spike_panel,result,ledger_path,summary)),encoding="utf-8")
        (args.outdir/"SUCCESS").write_text("evaluated_rows\t{}\nstatus\tPASS\n".format(len(output)),encoding="utf-8")
        print("[PASS] Evaluated {} disease-biomarker propagation rows".format(len(output)))
    except (ValueError,FileNotFoundError,KeyError) as error: raise SystemExit("[ERROR] {}".format(error)) from error
if __name__ == "__main__": main()
