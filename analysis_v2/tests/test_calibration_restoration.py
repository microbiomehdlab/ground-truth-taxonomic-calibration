#!/usr/bin/env python3
import csv, subprocess, tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
H=["cohort","study","analysis_population","target_label","assembly_arm","profiler","contrast",
   "dose_level","feature","effect","q_value","include","spike_fraction_target"]
M=["cohort","study","analysis_population","target_label","assembly_arm","profiler","target_feature"]
with tempfile.TemporaryDirectory() as t:
 t=Path(t);u=t/'u.tsv';c=t/'c.tsv';m=t/'m.tsv';out=t/'out'
 base=["zeller","study","independent","Target","original","metaphlan4","CRC_vs_Control"]
 ur=[base+["baseline","Bystander","1","0.01","1","0"],base+["dose_01","Bystander","0","0.2","1","0.1"],
     base+["baseline","Target species","0.5","0.02","1","0"],base+["dose_01","Target species","1.5","0.01","1","0.1"]]
 cr=[base+["baseline","Bystander","1","0.01","1","0"],base+["dose_01","Bystander","0.9","0.01","1","0.1"],
     base+["baseline","Target species","0.5","0.02","1","0"],base+["dose_01","Target species","1.5","0.01","1","0.1"]]
 for p,h,rows in ((u,H,ur),(c,H,cr),(m,M,[["zeller","study","independent","Target","original","metaphlan4","Target species"]])):
  with p.open('w',newline='',encoding='utf-8') as f:w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(h);w.writerows(rows)
 subprocess.run(['python3',str(ROOT/'analysis_v2/scripts/evaluate_calibration_restoration.py'),
   '--uncorrected-results',str(u),'--corrected-results',str(c),'--corrected-manifest',str(m),'--outdir',str(out)],check=True)
 with (out/'calibration_restoration_summary.tsv').open(newline='',encoding='utf-8') as f:r=list(csv.DictReader(f,delimiter='\t'))
 assert len(r)==1 and r[0]['baseline_biomarkers_rescued']=='1'
 assert float(r[0]['relative_mean_effect_restoration'])>.89
 print('[PASS] calibration-restoration fixture')
