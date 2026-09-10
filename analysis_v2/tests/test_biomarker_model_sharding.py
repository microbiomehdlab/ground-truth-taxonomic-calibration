#!/usr/bin/env python3
import csv, os, pathlib, subprocess, tempfile
root=pathlib.Path(__file__).resolve().parents[2]; script=root/'analysis_v2/scripts/shard_biomarker_model_input.py'
with tempfile.TemporaryDirectory() as t:
 d=pathlib.Path(t); manifest=d/'manifest.tsv'; abundance=d/'abundance.tsv'; out=d/'shards'
 fields=['cohort','analysis_population','target_label','profiler','source_profile','include']
 rows=[]
 for cohort in ['feng','zeller']:
  for pop in ['independent','community']:
   for target in ['Pana','Pint']:
    path=str((d/f'{cohort}_{pop}_{target}.tsv').resolve()); rows.append(dict(cohort=cohort,analysis_population=pop,target_label=target,profiler='metaphlan4',source_profile=path,include='1'))
 with manifest.open('w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
 with abundance.open('w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=['profiler','source_profile','feature','abundance_fraction'],delimiter='\t',lineterminator='\n'); w.writeheader()
  for r in rows: w.writerow(dict(profiler=r['profiler'],source_profile=r['source_profile'],feature='Species x',abundance_fraction='0.1'))
 counts=subprocess.check_output(['python3',str(script),'--profile-manifest',str(manifest),'--counts-only'],text=True).strip()
 assert counts=='8 4',counts
 subprocess.run(['python3',str(script),'--profile-manifest',str(manifest),'--abundance-long',str(abundance),'--outdir',str(out)],check=True)
 assert len(list((out/'paired').glob('*/profile_manifest.tsv')))==8
 assert len(list((out/'disease').glob('*/profile_manifest.tsv')))==4
 assert (out/'SUCCESS').is_file()
 fake=d/'sbatch'; counter=d/'counter'; calls=d/'calls'
 fake.write_text(f'''#!/usr/bin/env bash
n=$(cat "{counter}" 2>/dev/null || echo 2000)
n=$((n+1)); echo "$n" > "{counter}"; echo "$*" >> "{calls}"; echo "$n"
'''); fake.chmod(0o755)
 input_root=d/'input'; input_root.mkdir(); (input_root/'canonical_input.tsv').write_bytes(manifest.read_bytes())
 env=os.environ|{'DEV_INPUT_ROOT':str(input_root),'ANALYSIS_SIF':str(d/'image.sif'),
                 'OUTDIR':str(d/'run'),'SBATCH_BIN':str(fake)}
 submitted=subprocess.check_output(['bash',str(root/'analysis_v2/submit_legacy_biomarker_mapreduce.sh')],cwd=root,env=env,text=True)
 lines=calls.read_text().splitlines(); assert len(lines)==9
 assert '--array=1-8%12' in lines[1] and '--array=1-4%4' in lines[2]
 assert '--dependency=afterok:2006:2007:2008' in lines[8]
 assert 'seal=2009' in submitted
print('[PASS] biomarker model-sharding fixture')
