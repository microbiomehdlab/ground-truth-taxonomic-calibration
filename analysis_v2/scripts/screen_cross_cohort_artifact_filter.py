#!/usr/bin/env python3
"""Screen a blind cross-cohort artificial-call filter without oracle target protection.

This is a DEVELOPMENT_ONLY post-fit screen. A manuscript claim requires applying
the frozen rule before model fitting and recomputing the complete BH family.
"""
from __future__ import annotations
import argparse, csv, hashlib, math
from collections import defaultdict
from pathlib import Path

GROUP = ("cohort", "analysis_population", "profiler")
CONTEXT = ("cohort", "study", "analysis_population", "target_label",
           "assembly_arm", "profiler", "contrast", "spike_fraction_target")
REQUIRED = set(CONTEXT) | {"feature", "effect", "q_value", "include"}

def num(x: str) -> float:
    z = float(x)
    if not math.isfinite(z): raise ValueError(f"non-finite number: {x}")
    return z

def render(x): return format(x, ".17g")
def key(row, fields):
    # R and Python serialize the same IEEE-754 dose differently (for example,
    # 1e-04 versus 0.0001). Match contexts numerically, never by float text.
    return tuple(render(num(row[x])) if x == "spike_fraction_target" else row[x]
                 for x in fields)
def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""): h.update(block)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--paired-results",type=Path,required=True)
    p.add_argument("--metrics",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    p.add_argument("--thresholds",default="0,0.1,0.2,0.3,0.5,0.75,1")
    p.add_argument("--q-threshold",type=float,default=.05)
    p.add_argument("--minimum-relative-target-recall",type=float,default=.95)
    a=p.parse_args()
    if a.outdir.exists() and any(a.outdir.iterdir()): raise SystemExit("[ERROR] OUTDIR must be new or empty")
    thresholds=sorted(set(num(x) for x in a.thresholds.split(",")))
    if (any(x<0 or x>1 for x in thresholds) or not 0<a.q_threshold<1 or
            not 0 <= a.minimum_relative_target_recall <= 1):
        raise SystemExit("[ERROR] invalid threshold")

    # Resolve exact target aliases from the sealed evaluator rather than fuzzy names.
    targets={}; cohorts=set()
    with a.metrics.open(newline="",encoding="utf-8") as f:
        r=csv.DictReader(f,delimiter="\t")
        for row in r:
            if "analysis_scope" in row and row["analysis_scope"]!="pooled_primary": continue
            if abs(num(row["q_threshold"])-a.q_threshold)>1e-12: continue
            ck=key(row,CONTEXT)
            old=targets.setdefault(ck,row["target_alias"])
            if old!=row["target_alias"]: raise SystemExit("[ERROR] conflicting target aliases")
            cohorts.add(row["cohort"])
    if len(cohorts)<2: raise SystemExit("[ERROR] at least two cohorts are required")

    # Score recurrence across implanted targets, not correlated dose cells.
    eligible=defaultdict(set); called=defaultdict(set); fields=None; rows_seen=0
    with a.paired_results.open(newline="",encoding="utf-8") as f:
        r=csv.DictReader(f,delimiter="\t"); fields=set(r.fieldnames or [])
        if not REQUIRED<=fields: raise SystemExit("[ERROR] paired result columns are incomplete")
        for row in r:
            if row["include"]!="1" or num(row["spike_fraction_target"])<=0: continue
            ck=key(row,CONTEXT); target=targets.get(ck)
            if target is None: continue
            rows_seen+=1
            if row["feature"]==target: continue
            sk=key(row,GROUP)+(row["feature"],)
            eligible[sk].add(row["target_label"])
            if num(row["effect"])>0 and num(row["q_value"])<=a.q_threshold:
                called[sk].add(row["target_label"])
    scores={k:len(called[k])/len(v) for k,v in eligible.items() if v}

    # Apply each other-cohort score blindly. The held-out target can be removed.
    stats=defaultdict(lambda:dict(eligible_non_target=0,tp_before=0,fp_before=0,tp_after=0,fp_after=0))
    with a.paired_results.open(newline="",encoding="utf-8") as f:
        for row in csv.DictReader(f,delimiter="\t"):
            if row["include"]!="1" or num(row["spike_fraction_target"])<=0: continue
            ck=key(row,CONTEXT); target=targets.get(ck)
            if target is None: continue
            test=row["cohort"]
            for train in cohorts-{test}:
                score=scores.get((train,row["analysis_population"],row["profiler"],row["feature"]),0.0)
                for threshold in thresholds:
                    z=stats[(train,test)+ck[1:]+(threshold,)]
                    is_target=row["feature"]==target
                    if not is_target: z["eligible_non_target"]+=1
                    significant=num(row["effect"])>0 and num(row["q_value"])<=a.q_threshold
                    if not significant: continue
                    flagged=score>=threshold if threshold>0 else False
                    if is_target:
                        z["tp_before"]+=1
                        if not flagged: z["tp_after"]+=1
                    else:
                        z["fp_before"]+=1
                        if not flagged: z["fp_after"]+=1

    a.outdir.mkdir(parents=True,exist_ok=True)
    context_fields=("train_cohort","test_cohort")+CONTEXT[1:]+("score_threshold","eligible_non_target_tests",
      "target_called_before","target_called_after","off_target_calls_before","off_target_calls_after")
    context_path=a.outdir/"heldout_context_metrics.tsv"
    with context_path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=context_fields,delimiter="\t",lineterminator="\n");w.writeheader()
        for k,z in sorted(stats.items(),key=lambda q:tuple(map(str,q[0]))):
            row=dict(zip(("train_cohort","test_cohort")+CONTEXT[1:]+("score_threshold",),k))
            row.update(eligible_non_target_tests=z["eligible_non_target"],target_called_before=z["tp_before"],
              target_called_after=z["tp_after"],off_target_calls_before=z["fp_before"],off_target_calls_after=z["fp_after"]);w.writerow(row)

    # Aggregate micro-averaged held-out performance; threshold 0 is unfiltered.
    sums=defaultdict(lambda:defaultdict(int))
    for k,z in stats.items():
        g=(k[0],k[1],k[3],k[6],k[-1]) # train,test,population,profiler,threshold
        for name,value in z.items(): sums[g][name]+=value
    summary_fields=("train_cohort","test_cohort","analysis_population","profiler","score_threshold",
      "contexts","features_scored_in_training","eligible_non_target_tests","off_target_calls_before",
      "off_target_calls_after","off_target_reduction","off_target_rate_before","off_target_rate_after",
      "target_recall_before","target_recall_after","precision_before","precision_after")
    summary_path=a.outdir/"heldout_filter_summary.tsv"
    summary_rows=[]
    with summary_path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=summary_fields,delimiter="\t",lineterminator="\n");w.writeheader()
        for g,z in sorted(sums.items(),key=lambda q:tuple(map(str,q[0]))):
            train,test,pop,prof,threshold=g; nctx=sum(1 for k in stats if (k[0],k[1],k[3],k[6],k[-1])==g)
            btp,bfp,atp,afp=z["tp_before"],z["fp_before"],z["tp_after"],z["fp_after"]
            prec=lambda tp,fp:tp/(tp+fp) if tp+fp else 0
            outrow=dict(train_cohort=train,test_cohort=test,analysis_population=pop,profiler=prof,
              score_threshold=render(threshold),contexts=nctx,
              features_scored_in_training=sum(1 for k in scores if k[:3]==(train,pop,prof)),
              eligible_non_target_tests=z["eligible_non_target"],off_target_calls_before=bfp,
              off_target_calls_after=afp,off_target_reduction=render(1-afp/bfp if bfp else 0),
              off_target_rate_before=render(bfp/z["eligible_non_target"]),off_target_rate_after=render(afp/z["eligible_non_target"]),
              target_recall_before=render(btp/nctx),target_recall_after=render(atp/nctx),
              precision_before=render(prec(btp,bfp)),precision_after=render(prec(atp,afp)))
            summary_rows.append(outrow); w.writerow(outrow)

    # Select the best *screening* point under a prespecified recall constraint.
    # Threshold zero is the unfiltered comparator and guarantees a valid choice.
    choice_fields=summary_fields+("relative_target_recall","minimum_relative_target_recall",
      "screen_decision")
    choice_path=a.outdir/"recall_preserving_choices.tsv"
    grouped=defaultdict(list)
    for row in summary_rows:
        grouped[tuple(row[x] for x in summary_fields[:4])].append(row)
    with choice_path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=choice_fields,delimiter="\t",lineterminator="\n");w.writeheader()
        for group,rows in sorted(grouped.items()):
            candidate_points=[]
            for row in rows:
                before=num(row["target_recall_before"]); after=num(row["target_recall_after"])
                relative=after/before if before else 0.0
                if relative+1e-12>=a.minimum_relative_target_recall:
                    candidate_points.append((num(row["off_target_reduction"]),relative,
                                             num(row["score_threshold"]),row))
            if not candidate_points: raise SystemExit(f"[ERROR] no recall-preserving comparator for {group}")
            best_reduction=max(x[0] for x in candidate_points)
            if best_reduction<=0:
                reduction,relative,threshold,row=next(x for x in candidate_points if x[2]==0)
            else:
                reduction,relative,threshold,row=max(candidate_points,key=lambda x:(x[0],x[1],x[2]))
            selected=dict(row,relative_target_recall=render(relative),
              minimum_relative_target_recall=render(a.minimum_relative_target_recall),
              screen_decision=("NO_FILTER" if threshold==0 or reduction<=0 else "CANDIDATE_FILTER"))
            w.writerow(selected)

    score_path=a.outdir/"training_artifact_scores.tsv"
    with score_path.open("w",newline="",encoding="utf-8") as f:
        names=GROUP+("feature","eligible_targets","off_target_called_targets","artifact_score")
        w=csv.DictWriter(f,fieldnames=names,delimiter="\t",lineterminator="\n");w.writeheader()
        for k,s in sorted(scores.items()):
            w.writerow(dict(zip(GROUP+("feature",),k),eligible_targets=len(eligible[k]),
              off_target_called_targets=len(called[k]),artifact_score=render(s)))
    note=a.outdir/"DEVELOPMENT_ONLY.txt"
    note.write_text("status=DEVELOPMENT_ONLY\npostfit_screen_only=YES\nrefit_required_for_manuscript=YES\ntarget_oracle_protection=NO\n",encoding="utf-8")
    prov=a.outdir/"artifact_filter_screen.sha256"
    prov.write_text("".join(f"{digest(x)}  {x.resolve()}\n" for x in (a.paired_results,a.metrics,context_path,summary_path,choice_path,score_path,note)),encoding="utf-8")
    (a.outdir/"SUCCESS").write_text(f"paired_rows={rows_seen}\ntransfers={len(cohorts)*(len(cohorts)-1)}\nstatus=PASS\n",encoding="utf-8")
    print(f"[PASS] blind cross-cohort artifact-filter screen: {a.outdir}")

if __name__=="__main__": main()
