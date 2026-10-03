#!/usr/bin/env python3
"""Immutable DA3 plan, actual MaAsLin2, atomic resumable batches and collection."""
from __future__ import annotations
import argparse
import csv
import fcntl
import gzip
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import time
from collections import Counter
from decimal import Decimal
from pathlib import Path

from audit_bracken_denominators import digest, require, table, write_table
from build_biomarker_abundance_input import parse_profile
from prepare_da_pilot import verify
from plan_da3_canary import finish
from maaslin_timing import checked, fresh

CODE = ('run_maaslin_production.sbatch','submit_maaslin_production.sh',
        'scripts/maaslin_production.py','scripts/fit_maaslin_batch.R',
        'scripts/time_maaslin_context.R','scripts/maaslin_timing.py',
        'scripts/build_biomarker_abundance_input.py','scripts/audit_bracken_denominators.py',
        'scripts/prepare_da_pilot.py','scripts/plan_da3_canary.py','scripts/da_allocations.py',
        'lib/maaslin_context.R','lib/maaslin_contract.R','lib/unpaired_null.R')
FEATURE_FIELDS = ['context_id','feature','model_feature','beta','stderr','raw_p','native_q',
                  'estimable','status','p_for_BH','wrapper_q','positive_discovery']


def identity(root, repo):
    require(digest(root/'identity.json')==(root/'identity.sha256').read_text().split()[0], 'Run identity changed')
    state = json.loads((root/'identity.json').read_text())
    image=root/'image.sha256'
    if 'image_sha256' in state:
        sha,name=image.read_text().strip().split(None,1)
        require(sha==state['image_sha256'] and os.path.abspath(name)==state['image_path'], 'Image identity changed')
    require(digest(root/'plan/SHA256SUMS') == state['plan_sha256'], 'Plan identity changed')
    for name, sha in state['code'].items():
        require(digest(repo/'analysis_v2'/name) == sha, 'Source snapshot changed: '+name)
    return state


def prepare(source, root, repo, image):
    """Copy the already approved plan; never construct a new allocation."""
    require(not (root/'identity.json').exists(), 'Preparation already exists; use resume')
    fresh(root/'plan', [source])
    checked(source, ('tasks.tsv','catalog.json'))
    require('status\tPASS_DRAFT_BATCH_PLAN_NOT_PRODUCTION_AUTHORIZATION' in
            (source/'SUCCESS').read_text().splitlines(), 'Wrong source plan')
    tasks = table(source/'tasks.tsv')
    require(len(tasks)==4800 and [int(t['index']) for t in tasks]==list(range(4800)) and
            all(int(t['contexts'])==25 for t in tasks), 'Expected 4800 batches of 25')
    hashes = {n:s for s,n in (l.split(None,1) for l in (source/'SHA256SUMS').read_text().splitlines())}
    catalog = json.loads((source/'catalog.json').read_text())
    seen = set()
    for t in tasks:
        require(t['file'] in hashes, 'Unsealed batch file')
        contexts = json.loads((source/t['file']).read_text())
        require(len(contexts)==25, 'Wrong batch size')
        for c in contexts:
            require(c['context_id'] not in seen, 'Duplicate context')
            seen.add(c['context_id'])
            context_features(c,catalog)
    shutil.copytree(source,root/'plan')
    verify(root/'plan')
    require(digest(root/'plan/SHA256SUMS')==digest(source/'SHA256SUMS'), 'Source plan changed')
    image_hash,image_name=(root/'image.sha256').read_text().strip().split(None,1)
    require(re.fullmatch(r'[0-9a-f]{64}',image_hash) is not None and
            os.path.abspath(image_name)==os.path.abspath(image), 'Invalid host image fingerprint')
    state = dict(contract='DIRECT_MAASLIN2_DA3_V1', contexts=len(seen), batches=len(tasks),
        plan_sha256=digest(root/'plan/SHA256SUMS'), image_path=os.path.abspath(image),
        image_sha256=image_hash, code={n:digest(repo/'analysis_v2'/n) for n in CODE},
        source_commit=(root/'source_commit.txt').read_text().strip(),
        model='MaAsLin2_1.18.0_LM_group_only', transform='log2(1+a/1e-8)',
        normalization='NONE', multiplicity='BH_full_frozen_family_nonestimable_p1',
        retention='full_feature_stats_targets_compressed_console_logs_no_persistent_models')
    (root/'identity.json').write_text(json.dumps(state,sort_keys=True,indent=2)+'\n')
    (root/'identity.sha256').write_text(digest(root/'identity.json')+'  identity.json\n')
    print('[PASS] Exact existing plan copied:',len(seen),'contexts')


def context_features(c,catalog):
    require(re.fullmatch(r'da3_[0-9]{7}',c['context_id']) is not None, 'Unsafe context ID')
    obs=c['observations']; n=int(c['n'])
    require(n in (5,10,15,20) and len(obs)==2*n and
            len({r['sample_id'] for r in obs})==2*n and
            Counter(r['group'] for r in obs)=={0:n,1:n}, 'Invalid disjoint-person design')
    features=catalog['families']['|'.join((c['cohort'],c['profiler'],c['background']))]
    require(features and len(features)==len(set(features)), 'Invalid family')
    targets=catalog['targets'][c['profiler']]
    require(len(targets)==10 and len(set(targets.values()))==10 and set(targets)<=set(features),
            'Expected exact ten targets in frozen family')
    for r in obs:
        key='|'.join((c['cohort'],c['profiler'],r['sample_id'],str(Decimal(r['dose']))))
        require(key==r['key'], 'Observation key differs')
        p=catalog['profiles'][key]
        require(all(p[k]==v for k,v in dict(cohort=c['cohort'],profiler=c['profiler'],
                    sample_id=r['sample_id'],condition=c['background'],analysis_population='community').items())
                and Decimal(p['nominal_total_dose'])==Decimal(r['dose']), 'Profile metadata differs')
        require(r['group']==1 or Decimal(r['dose'])==0, 'Control is not a baseline')
    return features


def load_batch(root,index):
    plan=root/'plan'
    hashes={n:s for s,n in (l.split(None,1) for l in (plan/'SHA256SUMS').read_text().splitlines())}
    for name in ('tasks.tsv','catalog.json','SUCCESS'):
        require(digest(plan/name)==hashes[name], 'Plan file changed')
    tasks=table(plan/'tasks.tsv')
    require(0<=index<len(tasks) and int(tasks[index]['index'])==index, 'Invalid batch index')
    name=tasks[index]['file']
    require(digest(plan/name)==hashes[name], 'Batch definition changed')
    contexts=json.loads((plan/name).read_text())
    require(len(contexts)==int(tasks[index]['contexts']), 'Batch length differs')
    return contexts,json.loads((plan/'catalog.json').read_text())


def verify_batch(folder,state,contexts,catalog):
    checked(folder, ('features.tsv.gz','targets.tsv.gz','summary.tsv','input_hashes.tsv','identity.json'))
    require(json.loads((folder/'identity.json').read_text())==state, 'Completed batch identity differs')
    rows=table(folder/'summary.tsv')
    require([r['context_id'] for r in rows]==[c['context_id'] for c in contexts], 'Context coverage differs')
    for row,c in zip(rows,contexts):
        require(all(row[k]==str(c[k]) for k in ('cohort','background','profiler','n','arm','anchor','allocation_id')),
                'Context metadata differs')
        require(int(row['family_n'])==len(context_features(c,catalog)) and
                int(row['target_rows'])==10 and row['backend_version']=='1.18.0', 'Feature/backend count differs')
    with gzip.open(folder/'targets.tsv.gz','rt',newline='') as h:
        targets=list(csv.DictReader(h,delimiter='\t'))
    expected={(c['context_id'],f,label) for c in contexts for f,label in catalog['targets'][c['profiler']].items()}
    require(len(targets)==len(expected) and
            {(r['context_id'],r['feature'],r['target_label']) for r in targets}==expected, 'Target identity differs')
    return rows,targets


def validate_feature_rows(folder,contexts,catalog):
    """Prove every family member was retained once before sealing the batch."""
    with gzip.open(folder/'features.tsv.gz','rt',newline='') as h:
        reader=csv.DictReader(h,delimiter='\t')
        require(reader.fieldnames==FEATURE_FIELDS, 'Unexpected full-feature schema')
        for c in contexts:
            for feature in context_features(c,catalog):
                row=next(reader,None)
                require(row is not None and None not in row and None not in row.values() and
                        row['context_id']==c['context_id'] and row['feature']==feature,
                        'Missing/duplicate/misordered full-family result')
        require(next(reader,None) is None, 'Extra full-family results')


def run(root,index,repo):
    state=identity(root,repo); contexts,catalog=load_batch(root,index)
    results=root/'results'; results.mkdir(exist_ok=True)
    locks=root/'locks'; locks.mkdir(exist_ok=True)
    name='batch_%05d'%index; final=results/name
    with (locks/(name+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if final.exists():
            verify_batch(final,state,contexts,catalog)
            for r in table(final/'input_hashes.tsv'):
                require(digest(Path(r['path']))==r['sha256'], 'Completed input changed')
            print('[PASS] Already complete:',name); return
        attempts=root/'attempts'; attempts.mkdir(exist_ok=True)
        attempt=Path(tempfile.mkdtemp(prefix=name+'_',dir=str(attempts)))
        stage=attempt/'bundle'; stage.mkdir()
        started=time.monotonic(); hashes={}; cache={}
        with tempfile.TemporaryDirectory(prefix='maaslin_da3_') as temporary:
            scratch=Path(temporary); inp=scratch/'input'; inp.mkdir()
            try:
                meta=[]
                for c in contexts:
                    features=context_features(c,catalog); folder=inp/c['context_id']; folder.mkdir()
                    native=[]; people=[]
                    for r in c['observations']:
                        p=catalog['profiles'][r['key']]; path=Path(p['source_profile'])
                        require(str(path) not in hashes or hashes[str(path)]==p['sha256'], 'Conflicting hashes')
                        if str(path) not in cache:
                            require(digest(path)==p['sha256'], 'Profile checksum changed: '+str(path))
                            hashes[str(path)]=p['sha256']; cache[str(path)]=parse_profile(path,c['profiler'])
                        values=cache[str(path)]
                        require(values and all(math.isfinite(v) and 0<=v<=1.00001 for v in values.values()),'Invalid native data')
                        native.append(dict(observation_id=r['sample_id'],**{f:values.get(f,0) for f in features}))
                        people.append(dict(observation_id=r['sample_id'],group=r['group']))
                    write_table(folder/'abundance.tsv',native); write_table(folder/'metadata.tsv',people)
                    write_table(folder/'targets.tsv',[dict(feature=f,target_label=l) for f,l in catalog['targets'][c['profiler']].items()])
                    meta.append({k:c[k] for k in ('context_id','cohort','background','profiler','n','arm','anchor','allocation_id')})
                write_table(inp/'contexts.tsv',meta)
                with (attempt/'R.stdout').open('w') as stdout,(attempt/'R.stderr').open('w') as stderr:
                    subprocess.run(['Rscript',str(repo/'analysis_v2/scripts/fit_maaslin_batch.R'),
                                    str(inp),str(stage),str(repo)],stdout=stdout,stderr=stderr,check=True)
                validate_feature_rows(stage,contexts,catalog)
                require(all(digest(Path(p))==h for p,h in hashes.items()), 'Profile changed during batch')
                require(identity(root,repo)==state, 'Identity changed during batch')
            except Exception:
                # A failed batch never becomes final; retain its scratch evidence for review.
                shutil.copytree(inp,attempt/'failed_inputs')
                raise
        for name_log in ('R.stdout','R.stderr'):
            with (attempt/name_log).open('rb') as src,gzip.open(stage/(name_log+'.gz'),'wb') as dst:
                shutil.copyfileobj(src,dst)
        write_table(stage/'input_hashes.tsv',[dict(path=p,sha256=h) for p,h in sorted(hashes.items())])
        (stage/'identity.json').write_text(json.dumps(state,sort_keys=True)+'\n')
        write_table(stage/'batch_timing.tsv',[dict(contexts=len(contexts),seconds=time.monotonic()-started,
            output_bytes=sum(p.stat().st_size for p in stage.rglob('*') if p.is_file()))])
        finish(stage,'PASS_DIRECT_MAASLIN_BATCH')
        verify_batch(stage,state,contexts,catalog)
        os.rename(stage,final)
        # Only transient files created in this successful attempt are removed.
        shutil.rmtree(attempt)
        print('[PASS] Completed:',final.name)


def collect(root,repo,out):
    fresh(out,[root/'plan',root/'results']); state=identity(root,repo)
    verify(root/'plan'); tasks=table(root/'plan/tasks.tsv')
    catalog=json.loads((root/'plan/catalog.json').read_text())
    failures=[]; completed=0; total=0; mismatches=0; ledger=[]; out.mkdir(parents=True)
    with gzip.open(out/'targets.tsv.gz','wt',newline='') as th,(out/'context_summary.tsv').open('w',newline='') as sh:
        tw=sw=None
        for t in tasks:
            index=int(t['index']); folder=root/'results'/('batch_%05d'%index)
            try:
                contexts=json.loads((root/'plan'/t['file']).read_text())
                rows,targets=verify_batch(folder,state,contexts,catalog)
                if sw is None:
                    sw=csv.DictWriter(sh,fieldnames=list(rows[0]),delimiter='\t'); sw.writeheader()
                    tw=csv.DictWriter(th,fieldnames=list(targets[0]),delimiter='\t'); tw.writeheader()
                sw.writerows(rows); tw.writerows(targets)
                members={n:s for s,n in (l.split(None,1) for l in (folder/'SHA256SUMS').read_text().splitlines())}
                ledger.append(dict(index=index,path=str(folder/'features.tsv.gz'),sha256=members['features.tsv.gz'],
                    contexts=len(rows),feature_rows=sum(int(r['family_n']) for r in rows),
                    bundle_sha256=digest(folder/'SHA256SUMS')))
                completed+=1; total+=len(rows); mismatches+=sum(int(r['mismatched_features']) for r in rows)
            except (OSError,ValueError,KeyError) as exc:
                failures.append(dict(index=index,error=str(exc)))
    require(identity(root,repo)==state, 'Run identity changed during collection')
    status=dict(status='INCOMPLETE' if failures else ('COMPLETE_REVIEW_DIFFERENCES' if mismatches else
        'COMPLETE_DIRECT_MAASLIN_PENDING_SCIENTIFIC_REVIEW'), expected_batches=len(tasks),completed_batches=completed,
        expected_contexts=state['contexts'],completed_contexts=total,mismatched_features=mismatches,
        failures=failures,production_authorized=False)
    (out/'status.json').write_text(json.dumps(status,indent=2)+'\n')
    (out/'identity.json').write_text(json.dumps(state,sort_keys=True,indent=2)+'\n')
    if ledger: write_table(out/'feature_results_inventory.tsv',ledger)
    finish(out,'INCOMPLETE' if failures else 'PASS_DIRECT_MAASLIN_COLLECTION')
    print(json.dumps(status,indent=2)); return int(bool(failures))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); s=p.add_subparsers(dest='command',required=True)
    q=s.add_parser('prepare')
    for k in ('source','root','repo','image'): q.add_argument('--'+k,type=Path,required=True)
    q=s.add_parser('run'); q.add_argument('--index',type=int,required=True)
    for k in ('root','repo'): q.add_argument('--'+k,type=Path,required=True)
    q=s.add_parser('collect')
    for k in ('root','repo','out'): q.add_argument('--'+k,type=Path,required=True)
    args=vars(p.parse_args()); command=args.pop('command')
    raise SystemExit(globals()[command](**args) or 0)
