#!/usr/bin/env python3
"""Combine validated canonical inputs into one development-only input."""
import argparse, csv, hashlib, shutil, subprocess, sys
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--input', action='append', required=True, type=Path)
p.add_argument('--outdir', required=True, type=Path)
a=p.parse_args()
if len(a.input) < 2: p.error('--input must be supplied at least twice')
if a.outdir.exists() and any(a.outdir.iterdir()): raise SystemExit('[ERROR] outdir must be new or empty')
a.outdir.mkdir(parents=True, exist_ok=True)
rows=[]; fields=None; seen=set()
for source in a.input:
    marker=source.parent/'validation'/'SUCCESS'
    if not source.is_file() or not marker.is_file():
        raise SystemExit(f'[ERROR] input is not validated: {source}')
    with source.open(newline='',encoding='utf-8') as h:
        reader=csv.DictReader(h,delimiter='\t')
        if fields is None: fields=reader.fieldnames
        elif reader.fieldnames != fields: raise SystemExit(f'[ERROR] schema mismatch: {source}')
        for row in reader:
            key=(row['cohort'],row['sample_id'],row['analysis_population'],row['target_label'],
                 row['assembly_arm'],row['profiler'],row['profile_id'])
            if key in seen: raise SystemExit(f'[ERROR] duplicate canonical key: {key}')
            seen.add(key); rows.append(row)
out=a.outdir/'canonical_input.tsv'
with out.open('w',newline='',encoding='utf-8') as h:
    w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
validator=Path(__file__).with_name('validate_canonical_input.py')
subprocess.run([sys.executable,str(validator),'--input',str(out),'--outdir',str(a.outdir/'validation')],check=True)
(a.outdir/'DEVELOPMENT_ONLY.txt').write_text(
    'status\tDEVELOPMENT_ONLY\nuse_for_manuscript\tNO\nsource_policy\tlegacy_feng_zeller_plus_definitive_yachida\n')
with (a.outdir/'source_inputs.tsv').open('w',encoding='utf-8') as h:
    h.write('canonical_input\tsha256\n')
    for source in a.input: h.write(f'{source.resolve()}\t{hashlib.sha256(source.read_bytes()).hexdigest()}\n')
print(f'[PASS] Combined development canonical input: {len(rows)} rows from {len(a.input)} inputs')
