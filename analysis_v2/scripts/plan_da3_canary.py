#!/usr/bin/env python3
"""Plan and assemble every proposed DA3 arm; no inference or production launch."""
from __future__ import annotations
import argparse
import hashlib
import math
import time
import subprocess
import json
from decimal import Decimal,ROUND_HALF_UP
from pathlib import Path
from audit_bracken_denominators import table,digest,require,write_table
from build_biomarker_abundance_input import parse_profile
from prepare_da_pilot import verify
from da_allocations import allocations

DOSES=('0.0001','0.0005','0.001','0.005','0.01')
TOOLS=('kraken2_bracken','metaphlan4')
ARMS=(('U','1'),('P25','.25'),('P50','.5'),('P75','.75'),('V','1'),('PV','.5'),('N','0'))

def score(namespace,sid):
    return hashlib.sha256((namespace+'|'+sid).encode()).hexdigest()

def assign(cases,namespace,arm,exposure,anchor):
    require(arm in dict(ARMS) and Decimal(exposure)==Decimal(dict(ARMS)[arm]),'Invalid arm/exposure')
    require(len(cases)==len(set(cases)) and bool(cases),'Unique cases required')
    count=int((Decimal(exposure)*len(cases)).quantize(Decimal('1'),rounding=ROUND_HALF_UP))
    exposed=set(sorted(cases,key=lambda s:(score(namespace+'|exposure',s),s))[:count])
    return {sid:('0' if sid not in exposed else
        DOSES[int(score(namespace+'|dose',sid),16)%len(DOSES)] if arm in ('V','PV') else anchor)
        for sid in sorted(cases)}

def finish(folder,status):
    (folder/'SUCCESS').write_text('status\t'+status+'\nproduction_authorized\t0\n')
    (folder/'SHA256SUMS').write_text(''.join(digest(p)+'  '+str(p.relative_to(folder))+'\n' for p in sorted(folder.rglob('*')) if p.is_file() and p.name!='SHA256SUMS'))

def plan(inventory,out):
    inventory,out=inventory.resolve(),out.resolve()
    require(not out.exists() and inventory not in out.parents and out not in inventory.parents,'Fresh nonoverlapping output required')
    contexts=[];people=[];families=[];hashes={}
    repo=Path(__file__).resolve().parents[2]
    for rel in ('scripts/plan_da3_canary.py','scripts/fit_da3_canary.R','scripts/da_allocations.py',
                'scripts/build_biomarker_abundance_input.py','scripts/audit_bracken_denominators.py',
                'scripts/prepare_da_pilot.py','lib/exact_da3_candidate.R',
                'lib/monte_carlo_da3_candidate.R','lib/unpaired_null.R'):
        path=repo/'analysis_v2'/rel;hashes[path]=digest(path)
    for cohort in ('yachida','feng','zeller'):
        source=inventory/cohort;verify(source)
        require('status\tPASS_DESIGN_INVENTORY' in (source/'SUCCESS').read_text().splitlines(),'Wrong inventory status')
        for p in source.iterdir():
            if p.is_file():hashes[p]=digest(p)
        profiles=table(source/'profile_inventory.tsv');index={};pool={}
        for r in profiles:
            if r['analysis_population']!='community':continue
            require(r['cohort']==cohort and r['profiler'] in TOOLS,'Invalid inventory identity')
            key=(r['sample_id'],r['profiler'],Decimal(r['nominal_total_dose']))
            require(key not in index,'Duplicate profile key');index[key]=r
            if Decimal(r['nominal_total_dose'])==0:
                require(r['sample_id'] not in pool or pool[r['sample_id']]==r['condition'],'Condition conflict')
                pool[r['sample_id']]=r['condition']
        frozen=table(source/'feature_families.tsv')
        for background in ('Adenoma','CRC'):
            ids=sorted(s for s,c in pool.items() if c==background)
            require(len(ids)>=40,'Canary requires all four sample sizes')
            for draw in allocations(ids,cohort,background,repetitions=1):
                namespace='|'.join(('da3_heterogeneity_v1',cohort,background,draw['allocation_id']))
                # Two uniform anchors, one anchor for each partial arm; variable
                # arms draw from the fixed grid. This tests assembly, not all doses.
                scenarios=[(a,e,d) for a,e in ARMS for d in (('0.0001','0.005') if a=='U' else ('0.001',))]
                for arm,exposure,anchor in scenarios:
                    doses=assign(draw['cases'],namespace,arm,exposure,anchor)
                    for tool in TOOLS:
                        cid='canary_%04d'%len(contexts)
                        ff=[r for r in frozen if r['cohort']==cohort and r['condition']==background and r['profiler']==tool]
                        require(ff and len({r['feature'] for r in ff})==len(ff),'Invalid frozen family')
                        contexts.append(dict(context_id=cid,cohort=cohort,background=background,profiler=tool,n=draw['n'],arm=arm,requested_exposure=exposure,exposed_n=sum(v!='0' for v in doses.values()),anchor_total_dose=anchor if arm not in ('V','PV','N') else '',allocation_id=draw['allocation_id'],seed_namespace=namespace,family_n=len(ff),stage='ASSEMBLY_CANARY_NO_INFERENCE'))
                        families.extend(dict(context_id=cid,feature=r['feature']) for r in ff)
                        for group,samples in (('1',draw['cases']),('0',draw['controls'])):
                            for sid in samples:
                                dose=doses[sid] if group=='1' else '0';key=(sid,tool,Decimal(dose))
                                require(key in index,'Missing planned native profile')
                                r=index[key];path=Path(r['source_profile'])
                                require(path.is_absolute(),'Nonabsolute source profile')
                                people.append(dict(context_id=cid,observation_id=sid,biological_sample_id=sid,group=group,exposed=int(dose!='0'),nominal_total_dose=dose,nominal_target_dose=str(Decimal(dose)/10),source_profile=str(path),sha256=r['sha256']))
    require(len(contexts)==384,'Expected 384 canary contexts')
    require(all(digest(p)==h for p,h in hashes.items()),'Inventory changed')
    out.mkdir(parents=True)
    for name,rows in [('contexts.tsv',contexts),('observations.tsv',people),('families.tsv',families),('input_hashes.tsv',[dict(path=str(p),sha256=h) for p,h in sorted(hashes.items())])]:write_table(out/name,rows)
    (out/'workload.tsv').write_text('stage\tcontexts\tinference\nassembly_canary\t384\t0\nproposed_full_da3\t1200000\tPENDING_BUDGET_AND_METHOD_APPROVAL\n')
    finish(out,'PASS_DA3_ASSEMBLY_CANARY_PLAN');print('[PASS] Planned 384 assembly-only contexts:',out)

def assemble(plan_root,out,index,fit=False):
    plan_root,out=plan_root.resolve(),out.resolve();verify(plan_root)
    if (plan_root/'input_hashes.tsv').exists():
        for r in table(plan_root/'input_hashes.tsv'):
            require(digest(Path(r['path']))==r['sha256'],'Planned code/inventory changed')
    require(not out.exists() and plan_root not in out.parents and out not in plan_root.parents,'Fresh nonoverlapping task directory required')
    contexts=table(plan_root/'contexts.tsv');require(0<=index<len(contexts),'Invalid task index');context=contexts[index]
    cid=context['context_id'];obs=[r for r in table(plan_root/'observations.tsv') if r['context_id']==cid]
    features=[r['feature'] for r in table(plan_root/'families.tsv') if r['context_id']==cid]
    require(len(obs)==2*int(context['n']) and len({r['biological_sample_id'] for r in obs})==len(obs),'Invalid biological allocation')
    require(all(sum(r['group']==g for r in obs)==int(context['n']) for g in ('0','1')),'Unbalanced group')
    require(sum(r['exposed']=='1' for r in obs)==int(context['exposed_n']),'Exposure mismatch')
    started=time.monotonic();native=[];hashes={}
    for r in obs:
        path=Path(r['source_profile']).resolve()
        require(out!=path and out not in path.parents,'Output overlaps physical input')
        require(digest(path)==r['sha256'],'Native profile checksum mismatch');hashes[path]=r['sha256']
        values=parse_profile(path,context['profiler'])
        require(values and all(math.isfinite(v) and 0<=v<=1.00001 for v in values.values()),'Invalid native profile')
        native.append(dict(observation_id=r['observation_id'],**{f:values.get(f,0.0) for f in features}))
    require(all(digest(p)==h for p,h in hashes.items()),'Profile changed')
    out.mkdir(parents=True)
    write_table(out/'abundance.tsv',native);write_table(out/'metadata.tsv',obs);write_table(out/'context.tsv',[context])
    write_table(out/'timing.tsv',[dict(context_id=cid,assembly_seconds=time.monotonic()-started,n_observations=len(obs),family_n=len(features),inference_seconds='NOT_RUN')])
    write_table(out/'input_hashes.tsv',[dict(path=str(p),sha256=h) for p,h in sorted(hashes.items())])
    if fit:
        repo=Path(__file__).resolve().parents[2]
        subprocess.run(['Rscript',str(repo/'analysis_v2/scripts/fit_da3_canary.R'),str(out),str(repo)],check=True)
    finish(out,'PASS_ENGINEERING_CANARY_NOT_PRODUCTION' if fit else 'PASS_ASSEMBLY_ONLY_NOT_STATISTICAL_VALIDATION');print('[PASS]',cid,'engineering canary' if fit else 'native assembly; no inference')

def collect(plan_root,results,out):
    verify(plan_root);require(not out.exists(),'Fresh report required')
    require(all(out.resolve()!=p.resolve() and out.resolve() not in p.resolve().parents and p.resolve() not in out.resolve().parents for p in (plan_root,results)),'Report overlaps inputs')
    rows=[];failures=[]
    for context in table(plan_root/'contexts.tsv'):
        try:
            folder=results/context['context_id'];verify(folder)
            require('status\tPASS_ENGINEERING_CANARY_NOT_PRODUCTION' in (folder/'SUCCESS').read_text().splitlines(),'Inference incomplete')
            timing=table(folder/'fit_timing.tsv')[0]
            require(timing['context_id']==context['context_id'],'Context mismatch')
            rows.append(dict(context,**{k:v for k,v in timing.items() if k not in context}))
        except (ValueError,OSError,KeyError) as e:failures.append(dict(context_id=context['context_id'],error=str(e)))
    out.mkdir(parents=True)
    if rows:write_table(out/'timing_summary.tsv',rows)
    (out/'status.json').write_text(json.dumps(dict(status='INCOMPLETE' if failures else 'COMPLETE_PENDING_COST_AND_METHOD_REVIEW',completed=len(rows),expected=len(table(plan_root/'contexts.tsv')),failures=failures,production_authorized=False),indent=2)+'\n')
    finish(out,'INCOMPLETE' if failures else 'PASS_ENGINEERING_CANARY_REPORT');print((out/'status.json').read_text())

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('plan');q.add_argument('--inventory-root',type=Path,required=True);q.add_argument('--outdir',type=Path,required=True)
    q=sub.add_parser('assemble');q.add_argument('--plan',type=Path,required=True);q.add_argument('--outdir',type=Path,required=True);q.add_argument('--index',type=int,required=True)
    q.add_argument('--fit',action='store_true',help='Run candidate inference and record measured runtime')
    q=sub.add_parser('collect');q.add_argument('--plan',type=Path,required=True);q.add_argument('--results',type=Path,required=True);q.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args()
    if a.command=='plan':plan(a.inventory_root,a.outdir)
    elif a.command=='assemble':assemble(a.plan,a.outdir,a.index,a.fit)
    else:collect(a.plan,a.results,a.outdir)
