#!/usr/bin/env python3
import argparse, csv, hashlib
from collections import defaultdict
from pathlib import Path

p=argparse.ArgumentParser(); p.add_argument('--profile-manifest',required=True); p.add_argument('--abundance-long'); p.add_argument('--outdir'); p.add_argument('--counts-only',action='store_true'); a=p.parse_args()
if a.counts_only:
    with open(a.profile_manifest,newline='',encoding='utf-8') as h: rs=list(csv.DictReader(h,delimiter='\t'))
    rs=[r for r in rs if r.get('include')=='1']
    print(len({(r['cohort'],r['analysis_population'],r['target_label']) for r in rs}),
          len({(r['cohort'],r['analysis_population']) for r in rs}))
    raise SystemExit
if not a.abundance_long or not a.outdir: p.error('--abundance-long and --outdir are required unless --counts-only is used')
out=Path(a.outdir); out.mkdir(parents=True,exist_ok=False)
with open(a.profile_manifest,newline='',encoding='utf-8') as h: rows=list(csv.DictReader(h,delimiter='\t'))
pair=defaultdict(list); disease=defaultdict(list)
for r in rows:
    if r['include']!='1': continue
    pair[(r['cohort'],r['analysis_population'],r['target_label'])].append(r)
    disease[(r['cohort'],r['analysis_population'])].append(r)
fields=list(rows[0]); profile_keys=defaultdict(set); records=[]
def write(kind,key,subset,index):
    d=out/kind/f'{index:04d}'; d.mkdir(parents=True)
    path=d/'profile_manifest.tsv'
    with path.open('w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(subset)
    for r in subset: profile_keys[(r['profiler'],str(Path(r['source_profile']).resolve()))].add((kind,index))
    records.append({'kind':kind,'index':index,'cohort':key[0],'population':key[1],
                    'target_label':key[2] if len(key)>2 else 'ALL','rows':len(subset),'directory':str(d.resolve())})
for i,key in enumerate(sorted(pair),1): write('paired',key,pair[key],i)
for i,key in enumerate(sorted(disease),1): write('disease',key,disease[key],i)
handles={}; writers={}
try:
    with open(a.abundance_long,newline='',encoding='utf-8') as h:
        reader=csv.DictReader(h,delimiter='\t'); afields=list(reader.fieldnames or [])
        for r in reader:
            key=(r['profiler'],str(Path(r['source_profile']).resolve()))
            for kind,index in profile_keys.get(key,()):
                token=(kind,index)
                if token not in handles:
                    fh=(out/kind/f'{index:04d}'/'abundance_long.tsv').open('w',newline='',encoding='utf-8'); handles[token]=fh
                    writers[token]=csv.DictWriter(fh,fieldnames=afields,delimiter='\t',lineterminator='\n'); writers[token].writeheader()
                writers[token].writerow(r)
finally:
    for h in handles.values(): h.close()
for r in records:
    d=Path(r['directory']); ap=d/'abundance_long.tsv'
    if not ap.is_file(): raise SystemExit(f'[ERROR] no abundance rows for {d}')
    r['manifest_sha256']=hashlib.sha256((d/'profile_manifest.tsv').read_bytes()).hexdigest()
    r['abundance_sha256']=hashlib.sha256(ap.read_bytes()).hexdigest()
with (out/'shard_manifest.tsv').open('w',newline='',encoding='utf-8') as h:
    w=csv.DictWriter(h,fieldnames=list(records[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(records)
(out/'SUCCESS').write_text(f'paired_shards\t{len(pair)}\ndisease_shards\t{len(disease)}\nstatus\tPASS\n')
print(f'[PASS] Created {len(pair)} paired and {len(disease)} disease model shards')
