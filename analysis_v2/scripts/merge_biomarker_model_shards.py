#!/usr/bin/env python3
import argparse,csv,hashlib
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--shard-manifest',required=True); p.add_argument('--kind',choices=['paired','disease'],required=True); p.add_argument('--model-subdir'); p.add_argument('--outdir',required=True); a=p.parse_args()
root=Path(a.shard_manifest).parent; out=Path(a.outdir); out.mkdir(parents=True,exist_ok=False)
with open(a.shard_manifest,newline='',encoding='utf-8') as h: shards=[r for r in csv.DictReader(h,delimiter='\t') if r['kind']==a.kind]
if not shards: raise SystemExit('[ERROR] no expected shards')
files = (['paired_da_results.tsv','sample_feature_log2_changes.tsv'] if a.kind=='paired' else
         ['primary_disease_da_results.tsv','sensitivity_bmi_disease_da_results.tsv','disease_da_exclusions.tsv','disease_sample_panel_audit.tsv'])
counts={}
for name in files:
    output=out/name; header=None; total=0
    with output.open('w',newline='',encoding='utf-8') as oh:
        writer=None
        for shard in shards:
            # Shards are disjoint by construction. Validate that every output row
            # belongs to its declared shard, then retain duplicate signatures only
            # for the current shard. This bounds RAM by the largest shard instead of
            # the full merged result while still rejecting within- and cross-shard
            # identity violations.
            shard_seen=set()
            d=Path(shard['directory'])/(a.model_subdir or '')
            if not (d/'SUCCESS').is_file(): raise SystemExit(f"[ERROR] incomplete shard: {d}")
            path=d/name
            if not path.is_file(): raise SystemExit(f'[ERROR] missing {path}')
            with path.open(newline='',encoding='utf-8') as ih:
                reader=csv.DictReader(ih,delimiter='\t')
                if header is None: header=reader.fieldnames; writer=csv.DictWriter(oh,fieldnames=header,delimiter='\t',lineterminator='\n'); writer.writeheader()
                elif reader.fieldnames != header: raise SystemExit(f'[ERROR] schema mismatch: {path}')
                for row in reader:
                    if row.get('cohort') != shard['cohort'] or row.get('analysis_population') != shard['population']:
                        raise SystemExit(f'[ERROR] row outside declared shard in {path}')
                    if shard['target_label'] != 'ALL' and row.get('target_label') != shard['target_label']:
                        raise SystemExit(f'[ERROR] target outside declared shard in {path}')
                    signature=tuple(row[x] for x in header)
                    if signature in shard_seen: raise SystemExit(f'[ERROR] duplicate row within shard while merging {name}')
                    shard_seen.add(signature); writer.writerow(row); total+=1
    counts[name]=total
settings='paired_da_settings.tsv' if a.kind=='paired' else 'disease_da_settings.tsv'
first=Path(shards[0]['directory'])/(a.model_subdir or '')/settings
(out/settings).write_bytes(first.read_bytes())
with (out/('paired_da_summary.tsv' if a.kind=='paired' else 'disease_da_summary.tsv')).open('w',encoding='utf-8') as h:
    h.write('metric\tvalue\n'); [h.write(f'{k}\t{v}\n') for k,v in counts.items()]; h.write(f'shards\t{len(shards)}\nstatus\tPASS\n')
with (out/'model_shards.sha256').open('w',encoding='utf-8') as h:
    for name in files+[settings]: h.write(f'{hashlib.sha256((out/name).read_bytes()).hexdigest()}  {out/name}\n')
(out/'SUCCESS').write_text(f'analysis\t{a.kind}_model_shard_merge\nshards\t{len(shards)}\nstatus\tPASS\n')
print(f'[PASS] Merged {len(shards)} {a.kind} model shards into {out}')
