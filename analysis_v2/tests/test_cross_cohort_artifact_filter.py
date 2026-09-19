#!/usr/bin/env python3
import csv,subprocess,tempfile
from pathlib import Path
repo=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory() as t:
 r=Path(t); paired=r/"paired.tsv"; metrics=r/"metrics.tsv"; out=r/"out"
 pf=["cohort","study","analysis_population","target_label","assembly_arm","profiler","contrast","spike_fraction_target","feature","effect","q_value","include"]
 mf=list(pf[:8])+["q_threshold","target_alias","analysis_scope"]
 prows=[];mrows=[]
 for cohort in ("A","B"):
  for target,alias in (("T1","one"),("T2","two")):
   heldout_alias="bad" if cohort=="B" and target=="T1" else alias
   # Deliberately use a different valid float serialization from paired.tsv.
   mrows.append([cohort,"S","community",target,"original","P", "paired","0.10000000000000001",.05,heldout_alias,"pooled_primary"])
   # 'bad' is recurrent off-target; in B/T1 it is also the true target alias.
   for feature in ("one","two","bad","clean"):
    q=.01 if feature in (alias,"bad") else .5
    prows.append([cohort,"S","community",target,"original","P","paired",.1,feature,1,q,1])
 with paired.open("w",newline="") as h:w=csv.writer(h,delimiter="\t");w.writerow(pf);w.writerows(prows)
 with metrics.open("w",newline="") as h:w=csv.writer(h,delimiter="\t");w.writerow(mf);w.writerows(mrows)
 subprocess.run(["python3",str(repo/"analysis_v2/scripts/screen_cross_cohort_artifact_filter.py"),"--paired-results",str(paired),"--metrics",str(metrics),"--outdir",str(out),"--thresholds","0,0.5"],check=True,capture_output=True,text=True)
 rows=list(csv.DictReader((out/"heldout_filter_summary.tsv").open(),delimiter="\t")); filt=[x for x in rows if x["score_threshold"]=="0.5"]
 assert filt and all(int(x["off_target_calls_after"])<int(x["off_target_calls_before"]) for x in filt)
 assert any(float(x["target_recall_after"])<float(x["target_recall_before"]) for x in filt)
 choices=list(csv.DictReader((out/"recall_preserving_choices.tsv").open(),delimiter="\t"))
 assert len(choices)==2
 assert all(float(x["relative_target_recall"])>=.95 for x in choices)
 assert any(x["screen_decision"]=="NO_FILTER" for x in choices)
 assert (out/"SUCCESS").is_file() and (out/"DEVELOPMENT_ONLY.txt").is_file()
print("[PASS] blind cross-cohort artifact-filter fixture")
