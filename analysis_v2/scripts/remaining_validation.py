#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,subprocess,shutil
from pathlib import Path
from audit_bracken_denominators import table,require,write_table,digest
from prepare_da_pilot import verify
from plan_da3_canary import finish

def plan(pilot,out):
    pilot,out=pilot.resolve(),out.resolve();require(not out.exists() and pilot not in out.parents and out not in pilot.parents,'Fresh nonoverlapping output required')
    tasks=[];identities=set()
    for source in sorted(pilot.glob('pilot_*')):
        if not source.is_dir():continue
        context=table(source/'context.tsv')[0]
        if context['analysis']!='DA1':continue
        key=(context['cohort'],context['profiler'],context['background'])
        require(key not in identities,'Duplicate DA1 context');identities.add(key)
        verify(source)
        for mode,scenario in [('clinical_backend','actual'),('clinical_null','gaussian'),('clinical_null','skewed'),('clinical_null','heteroskedastic')]:
            tasks.append(dict(index=len(tasks),mode=mode,scenario=scenario,source=str(source),n='',family_n='',source_sha256=digest(source/'SHA256SUMS')))
    require(len(tasks)==48,'Expected all twelve completed DA1 pilot contexts')
    require(identities=={(c,t,b) for c in ('yachida','feng','zeller') for t in ('kraken2_bracken','metaphlan4') for b in ('Adenoma','CRC')},'Incomplete clinical context grid')
    for n in (5,10,15,20):
        for m in (474,3471):
            for scenario in ('gaussian','correlated','skewed'):
                tasks.append(dict(index=len(tasks),mode='partial_null',scenario=scenario,source='',n=n,family_n=m,source_sha256=''))
    out.mkdir(parents=True);write_table(out/'tasks.tsv',tasks);finish(out,'PASS_REMAINING_VALIDATION_PLAN');print('[PASS] Planned',len(tasks),'validation tasks')

def run(plan_root,results,index,repo):
    verify(plan_root);tasks=table(plan_root/'tasks.tsv');require(0<=index<len(tasks),'Invalid index');r=tasks[index]
    task=results/('validation_%03d'%index);require(not task.exists(),'Fresh task required')
    require(results.resolve()!=plan_root.resolve() and plan_root.resolve() not in results.resolve().parents and results.resolve() not in plan_root.resolve().parents,'Overlapping results')
    task.mkdir(parents=True);write_table(task/'settings.tsv',[r]);hashes={}
    if r['source']:
        source=Path(r['source']);require(digest(source/'SHA256SUMS')==r['source_sha256'],'Source checkpoint changed');verify(source)
        for name in ('abundance.tsv','metadata.tsv','context.tsv'):shutil.copyfile(source/name,task/name);hashes[str(source/name)]=digest(source/name)
    for name in ('scripts/run_remaining_validation.R','lib/maaslin_context.R','lib/maaslin_contract.R','lib/exact_da3_candidate.R','lib/monte_carlo_da3_candidate.R'):
        path=repo/'analysis_v2'/name;hashes[str(path)]=digest(path)
    with (task/'model.out').open('w') as stdout,(task/'model.err').open('w') as stderr:
        subprocess.run(['Rscript',str(repo/'analysis_v2/scripts/run_remaining_validation.R'),str(task),r['mode'],str(repo)],stdout=stdout,stderr=stderr,check=True)
    require(all(digest(Path(p))==h for p,h in hashes.items()),'Source/code changed')
    write_table(task/'input_hashes.tsv',[dict(path=p,sha256=h) for p,h in sorted(hashes.items())]);finish(task,'PASS_VALIDATION_COMPUTATION')

def collect(plan_root,results,out):
    verify(plan_root);require(not out.exists(),'Fresh report required');rows=[];errors=[]
    for r in table(plan_root/'tasks.tsv'):
        try:
            task=results/('validation_%03d'%int(r['index']));verify(task)
            require('status\tPASS_VALIDATION_COMPUTATION' in (task/'SUCCESS').read_text().splitlines(),'Wrong validation status')
            for s in table(task/'summary.tsv'):rows.append(dict(r,**s))
        except (OSError,ValueError,KeyError) as e:errors.append(dict(index=r['index'],error=str(e)))
    out.mkdir(parents=True)
    if rows:
        fields=sorted({k for r in rows for k in r});write_table(out/'comparison.tsv',[{k:r.get(k,'') for k in fields} for r in rows])
    (out/'status.json').write_text(json.dumps(dict(status='INCOMPLETE' if errors else 'COMPLETE_PENDING_SCIENTIFIC_REVIEW',failures=errors,production_authorized=False),indent=2)+'\n');finish(out,'REPORT_ONLY_NOT_AUTHORIZATION')

if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('plan');q.add_argument('--pilot',type=Path,required=True);q.add_argument('--out',type=Path,required=True)
    q=sub.add_parser('run');q.add_argument('--plan',type=Path,required=True);q.add_argument('--results',type=Path,required=True);q.add_argument('--index',type=int,required=True);q.add_argument('--repo',type=Path,required=True)
    q=sub.add_parser('collect');q.add_argument('--plan',type=Path,required=True);q.add_argument('--results',type=Path,required=True);q.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.command=='plan':plan(a.pilot,a.out)
    elif a.command=='run':run(a.plan,a.results,a.index,a.repo)
    else:collect(a.plan,a.results,a.out)
