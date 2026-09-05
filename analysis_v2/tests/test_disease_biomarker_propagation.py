#!/usr/bin/env python3
import csv, subprocess, tempfile
from pathlib import Path
repo=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp); calls=root/"calls.tsv"; aliases=root/"aliases.csv"; panel=root/"panel.tsv"; out=root/"out"
    fields=["cohort","study","analysis_population","target_label","assembly_arm","profiler","contrast","spike_fraction_target","feature","effect","p_value","q_value","include","exclusion_reason"]
    rows=[]
    for dose in (0,.01):
      values={"Target":(.2,.01 if dose else .2),"Retained":(.5,.01),"Lost":(.4,.2 if dose else .01),"Gained":(-.3,.01 if dose else .2)}
      for feature,(effect,q) in values.items(): rows.append(["yachida","S","independent","Pana","clean","metaphlan4","CRC_vs_Control",dose,feature,effect,q*q,q,1,""])
    with calls.open("w",newline="") as h: w=csv.writer(h,delimiter="\t"); w.writerow(fields); w.writerows(rows)
    aliases.write_text("canonical,alias,tool\nSpecies,Target,metaphlan4\n"); panel.write_text("label\ttaxon_name\nPana\tSpecies\n")
    subprocess.run(["python3",str(repo/"analysis_v2/scripts/evaluate_disease_biomarker_propagation.py"),"--calls",str(calls),"--aliases",str(aliases),"--spike-panel",str(panel),"--outdir",str(out)],check=True,capture_output=True,text=True)
    data=list(csv.DictReader((out/"disease_biomarker_propagation_metrics.tsv").open(),delimiter="\t")); q05=[x for x in data if abs(float(x["q_threshold"])-.05)<1e-9][0]
    assert (q05["baseline_biomarkers"],q05["retained_biomarkers"],q05["lost_biomarkers"],q05["gained_biomarkers"]) == ("2","1","1","2")
    assert abs(float(q05["biomarker_set_jaccard_vs_baseline"])-.25)<1e-12 and q05["target_significant"]=="1"
print("[PASS] disease-biomarker propagation fixture")
