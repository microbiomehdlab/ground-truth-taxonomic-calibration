#!/usr/bin/env python3
import csv,pathlib,subprocess,tempfile
root=pathlib.Path(__file__).resolve().parents[2]; script=root/'analysis_v2/scripts/merge_biomarker_model_shards.py'
with tempfile.TemporaryDirectory() as temporary:
 t=pathlib.Path(temporary); records=[]
 for index,(cohort,target) in enumerate([('feng','Pana'),('zeller','Pint')],1):
  d=t/'paired'/f'{index:04d}'/'pooled'; d.mkdir(parents=True)
  common='cohort\tanalysis_population\ttarget_label\tfeature\teffect\n'
  (d/'paired_da_results.tsv').write_text(common+f'{cohort}\tindependent\t{target}\tX\t{index}\n')
  (d/'sample_feature_log2_changes.tsv').write_text(common+f'{cohort}\tindependent\t{target}\tX\t{index}\n')
  (d/'paired_da_settings.tsv').write_text('setting\tvalue\ntest\tfixture\n'); (d/'SUCCESS').write_text('status\tPASS\n')
  records.append(dict(kind='paired',index=index,cohort=cohort,population='independent',target_label=target,rows=1,directory=str(d.parent),manifest_sha256='x',abundance_sha256='y'))
 manifest=t/'shard_manifest.tsv'
 with manifest.open('w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=list(records[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(records)
 subprocess.run(['python3',str(script),'--shard-manifest',str(manifest),'--kind','paired','--model-subdir','pooled','--outdir',str(t/'merged')],check=True)
 assert (t/'merged/SUCCESS').is_file()
 assert len((t/'merged/paired_da_results.tsv').read_text().splitlines())==3
print('[PASS] bounded-memory biomarker shard-merge fixture')
