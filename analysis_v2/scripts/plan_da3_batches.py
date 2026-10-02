#!/usr/bin/env python3
"""Compact 100-allocation design and batched workers; never auto-launch."""
from __future__ import annotations
import argparse,csv,gzip,json,math,subprocess,tempfile,time,fcntl,os
from decimal import Decimal
from pathlib import Path
from audit_bracken_denominators import table,require,digest,write_table
from prepare_da_pilot import verify
from da_allocations import allocations
from plan_da3_canary import assign,DOSES,TOOLS,finish,score
from build_biomarker_abundance_input import parse_profile

def plan(inventory,out,repetitions=100,batch_size=25):
    inventory,out=inventory.resolve(),out.resolve()
    require(repetitions==100 and batch_size>0,'Approved draft is 100 allocations; positive batch size required')
    require(not out.exists() and inventory not in out.parents and out not in inventory.parents,'Fresh nonoverlapping plan required')
    catalog={};family={};hashes={};pools={}
    for cohort in ('yachida','feng','zeller'):
        folder=inventory/cohort;verify(folder)
        require('status\tPASS_DESIGN_INVENTORY' in (folder/'SUCCESS').read_text().splitlines(),'Wrong inventory')
        hashes[str(folder/'SHA256SUMS')]=digest(folder/'SHA256SUMS')
        pools[cohort]={}
        for r in table(folder/'profile_inventory.tsv'):
            if r['analysis_population']!='community':continue
            key='|'.join((cohort,r['profiler'],r['sample_id'],str(Decimal(r['nominal_total_dose']))))
            require(key not in catalog,'Duplicate physical context');catalog[key]=r
            if Decimal(r['nominal_total_dose'])==0:pools[cohort][r['sample_id']]=r['condition']
        for r in table(folder/'feature_families.tsv'):
            key='|'.join((cohort,r['profiler'],r['condition']));family.setdefault(key,[]).append(r['feature'])
    repo=Path(__file__).resolve().parents[2]
    panel=table(repo/'spikes/spike_panel.tsv');canonical={r['taxon_name']:r['label'] for r in panel};targets={}
    with (repo/'examples/spike_taxon_aliases.csv').open(newline='') as h:
        for r in csv.DictReader(h):
            if r['canonical'] in canonical:targets.setdefault(r['tool'],{})[r['alias']]=canonical[r['canonical']]
    require(all(len(targets[t])==10 for t in TOOLS),'Exact ten-target coverage required')
    for cohort in pools:
        for tool in TOOLS:
            for condition in ('Adenoma','CRC'):
                features=family.get('|'.join((cohort,tool,condition)),[])
                require(len(features)==len(set(features)) and set(targets[tool])<=set(features),'Incomplete frozen target family')
    for path in (repo/'spikes/spike_panel.tsv',repo/'examples/spike_taxon_aliases.csv'):hashes[str(path)]=digest(path)
    out.mkdir(parents=True);(out/'batches').mkdir();buffer=[];ledger=[];count=0
    def flush():
        if not buffer:return
        name='batch_%05d.json'%len(ledger);(out/'batches'/name).write_text(json.dumps(buffer))
        ledger.append(dict(index=len(ledger),file='batches/'+name,contexts=len(buffer)));buffer.clear()
    scenarios=[('U','1',d) for d in DOSES+('0.05','0.1')]+[(a,p,d) for a,p in (('P25','.25'),('P50','.5'),('P75','.75')) for d in DOSES]+[('V','1','0.001'),('PV','.5','0.001'),('N','0','0.001')]
    for cohort in pools:
        for background in ('Adenoma','CRC'):
            ids=sorted(s for s,c in pools[cohort].items() if c==background);require(len(ids)>=40,'Insufficient pool')
            for draw in allocations(ids,cohort,background,repetitions=repetitions):
                namespace='|'.join(('da3_heterogeneity_v1',cohort,background,draw['allocation_id']))
                for arm,p,d in scenarios:
                    doses=assign(draw['cases'],namespace,arm,p,d)
                    seed=int(score(namespace+'|permutation|'+str(draw['n']),arm+'|'+d),16)%2147483646+1
                    for tool in TOOLS:
                        observations=[]
                        for group,samples in ((1,draw['cases']),(0,draw['controls'])):
                            for sid in samples:
                                dose=doses[sid] if group else '0';key='|'.join((cohort,tool,sid,str(Decimal(dose))))
                                require(key in catalog,'Missing planned profile')
                                observations.append(dict(sample_id=sid,group=group,dose=dose,key=key))
                        buffer.append(dict(context_id='da3_%07d'%count,cohort=cohort,background=background,profiler=tool,n=draw['n'],arm=arm,requested_exposure=p,anchor=d if arm not in ('N','V','PV') else '',allocation_id=draw['allocation_id'],permutation_seed=seed,observations=observations));count+=1
                        if len(buffer)==batch_size:flush()
    flush();require(count==120000,'Expected complete 100-allocation grid')
    (out/'catalog.json').write_text(json.dumps(dict(profiles=catalog,families=family,targets=targets)))
    write_table(out/'tasks.tsv',ledger);write_table(out/'input_hashes.tsv',[dict(path=p,sha256=h) for p,h in sorted(hashes.items())])
    finish(out,'PASS_DRAFT_BATCH_PLAN_NOT_PRODUCTION_AUTHORIZATION');print('[PASS]',count,'contexts;',len(ledger),'batches')

def run(plan_root,results,index,repo):
    plan_root,results,repo=plan_root.resolve(),results.resolve(),repo.resolve()
    require(results!=plan_root and plan_root not in results.parents and results not in plan_root.parents,'Overlapping results')
    require(results!=repo and results not in repo.parents,'Results contain source repository')
    files=[repo/'analysis_v2'/n for n in ('scripts/plan_da3_batches.py','scripts/fit_da3_batch_context.R','scripts/build_biomarker_abundance_input.py','scripts/audit_bracken_denominators.py','scripts/prepare_da_pilot.py','scripts/plan_da3_canary.py','lib/exact_da3_candidate.R','lib/monte_carlo_da3_candidate.R','lib/unpaired_null.R')]
    identity=dict(plan_sha256=digest(plan_root/'SHA256SUMS'),code={str(p):digest(p) for p in files})
    results.mkdir(parents=True,exist_ok=True);locks=results/'locks';locks.mkdir(exist_ok=True)
    name='batch_%05d'%index;final=results/name
    with (locks/(name+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if final.exists():
            verify(final)
            require(json.loads((final/'resume_identity.json').read_text())==identity,'Resume identity differs')
            for r in table(final/'input_hashes.tsv'):require(digest(Path(r['path']))==r['sha256'],'Completed batch source changed')
            print('[PASS] Already completed:',name);return
        attempts=results/'attempts';attempts.mkdir(exist_ok=True)
        attempt=Path(tempfile.mkdtemp(prefix=name+'_',dir=str(attempts)))
        _run_new(plan_root,attempt,index,repo)
        staged=attempt/name
        require(identity['code']=={str(p):digest(p) for p in files},'Code changed during batch')
        (staged/'resume_identity.json').write_text(json.dumps(identity,sort_keys=True)+'\n')
        finish(staged,'PASS_BATCH_COMPUTATION_NOT_PRODUCTION_AUTHORIZATION')
        os.rename(staged,final)
        attempt.rmdir()
        print('[PASS] Completed:',name)

def _run_new(plan_root,results,index,repo):
    # Verify only this task and shared catalog against the plan checksum ledger;
    # do not rehash all 4,800 immutable batch files in every worker.
    hashes={name:h for h,name in (line.split(None,1) for line in (plan_root/'SHA256SUMS').read_text().splitlines())}
    for name in ('SUCCESS','tasks.tsv','catalog.json'):require(digest(plan_root/name)==hashes[name],'Plan checksum mismatch')
    tasks=table(plan_root/'tasks.tsv');require(0<=index<len(tasks),'Invalid batch index');name=tasks[index]['file']
    require(digest(plan_root/name)==hashes[name],'Batch checksum mismatch')
    catalog=json.loads((plan_root/'catalog.json').read_text());contexts=json.loads((plan_root/name).read_text())
    out=results/('batch_%05d'%index);require(not out.exists(),'Fresh batch required')
    require(plan_root.resolve() not in out.resolve().parents and out.resolve() not in plan_root.resolve().parents,'Overlapping output')
    cache={};source_hashes={};summaries=[];targets=[];start=time.monotonic()
    codefiles=[repo/'analysis_v2'/n for n in ('scripts/plan_da3_batches.py','scripts/fit_da3_batch_context.R','lib/exact_da3_candidate.R','lib/monte_carlo_da3_candidate.R','lib/unpaired_null.R')]
    codehash={str(p):digest(p) for p in codefiles}
    for context in contexts:
        tool=context['profiler'];features=catalog['families']['|'.join((context['cohort'],tool,context['background']))]
        obs=context['observations'];require(len({r['sample_id'] for r in obs})==2*context['n'],'Biological overlap')
        native=[];metadata=[]
        for r in obs:
            profile=catalog['profiles'][r['key']];path=Path(profile['source_profile'])
            require(str(path) not in source_hashes or source_hashes[str(path)]==profile['sha256'],'Conflicting profile hashes')
            if str(path) not in cache:
                require(digest(path)==profile['sha256'],'Source profile changed');source_hashes[str(path)]=profile['sha256'];cache[str(path)]=parse_profile(path,tool)
            values=cache[str(path)];require(values and all(math.isfinite(v) and 0<=v<=1.00001 for v in values.values()),'Invalid native profile')
            native.append(dict(observation_id=r['sample_id'],**{f:values.get(f,0) for f in features}));metadata.append(dict(observation_id=r['sample_id'],group=r['group']))
        with tempfile.TemporaryDirectory(prefix='da3_batch_') as d:
            temp=Path(d);write_table(temp/'abundance.tsv',native);write_table(temp/'metadata.tsv',metadata)
            write_table(temp/'context.tsv',[{k:v for k,v in context.items() if k!='observations'}])
            subprocess.run(['Rscript',str(repo/'analysis_v2/scripts/fit_da3_batch_context.R'),str(temp),str(repo)],check=True)
            fit=table(temp/'results.tsv');summaries.extend(table(temp/'summary.tsv'))
            for row in fit:
                if row['feature'] in catalog['targets'][tool]:targets.append(dict(context_id=context['context_id'],target_label=catalog['targets'][tool][row['feature']],**row))
    require(all(digest(Path(p))==h for p,h in dict(source_hashes,**codehash).items()),'Source/code changed during batch')
    out.mkdir(parents=True);write_table(out/'summary.tsv',summaries)
    with gzip.open(out/'targets.tsv.gz','wt',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(targets[0]),delimiter='\t');w.writeheader();w.writerows(targets)
    write_table(out/'input_hashes.tsv',[dict(path=p,sha256=h) for p,h in sorted(dict(source_hashes,**codehash).items())]);write_table(out/'timing.tsv',[dict(contexts=len(contexts),elapsed_seconds=time.monotonic()-start)])
    finish(out,'PASS_BATCH_COMPUTATION_NOT_PRODUCTION_AUTHORIZATION')

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='command',required=True)
    q=s.add_parser('plan');q.add_argument('--inventory',type=Path,required=True);q.add_argument('--out',type=Path,required=True)
    q=s.add_parser('run');q.add_argument('--plan',type=Path,required=True);q.add_argument('--results',type=Path,required=True);q.add_argument('--repo',type=Path,required=True);q.add_argument('--index',type=int,required=True)
    a=p.parse_args()
    if a.command=='plan':plan(a.inventory,a.out)
    else:run(a.plan,a.results,a.index,a.repo)
