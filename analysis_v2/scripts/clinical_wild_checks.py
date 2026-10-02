#!/usr/bin/env python3
"""Separate engineering checks, never production authorization."""
from __future__ import annotations
import argparse, json, shutil, subprocess, math
from pathlib import Path
from audit_bracken_denominators import table, require, write_table, digest
from prepare_da_pilot import verify
from plan_da3_canary import finish

def plan(pilot, out):
    pilot, out = pilot.resolve(), out.resolve()
    require(not out.exists() and pilot != out and pilot not in out.parents and out not in pilot.parents, 'Fresh nonoverlapping output required')
    rows, identities = [], set()
    for source in sorted(pilot.glob('pilot_*')):
        if not source.is_dir(): continue
        context = table(source/'context.tsv')[0]
        if context['analysis'] != 'DA1': continue
        verify(source)
        key = (context['cohort'], context['profiler'], context['background'])
        require(key not in identities, 'Duplicate clinical context'); identities.add(key)
        for scenario in ('actual', 'gaussian', 'skewed', 'heteroskedastic'):
            for chunk in range(1 if scenario == 'actual' else 4):
                rows.append(dict(index=len(rows), source=str(source), source_sha256=digest(source/'SHA256SUMS'), scenario=scenario, chunk=chunk, draws=9999, simulations=0 if scenario=='actual' else 50, family_n=20))
    require(identities == {(c,t,b) for c in ('yachida','feng','zeller') for t in ('kraken2_bracken','metaphlan4') for b in ('Adenoma','CRC')}, 'Incomplete twelve-context clinical grid')
    require(len(rows)==156, 'Unexpected task count')
    out.mkdir(parents=True); write_table(out/'tasks.tsv', rows); finish(out, 'PASS_WILD_ENGINEERING_PLAN')

def run(plan_root, results, index, repo):
    verify(plan_root); rows=table(plan_root/'tasks.tsv')
    require(0<=index<len(rows), 'Invalid task index')
    require(results.resolve()!=plan_root.resolve() and plan_root.resolve() not in results.resolve().parents and results.resolve() not in plan_root.resolve().parents, 'Overlapping results')
    row=rows[index]; source=Path(row['source'])
    require(digest(source/'SHA256SUMS')==row['source_sha256'], 'Source checkpoint changed'); verify(source)
    task=results/('task_%03d'%index); require(not task.exists(), 'Fresh task required')
    task.mkdir(parents=True); write_table(task/'settings.tsv',[row])
    hashes={}
    for name in ('metadata.tsv','abundance.tsv','context.tsv'):
        p=source/name; hashes[str(p)]=digest(p); shutil.copyfile(p,task/name)
    for name in ('scripts/run_clinical_wild_checks.R','lib/clinical_hc3.R','lib/clinical_wild_bootstrap.R'):
        p=repo/'analysis_v2'/name; hashes[str(p)]=digest(p)
    with (task/'model.out').open('w') as stdout, (task/'model.err').open('w') as stderr:
        subprocess.run(['Rscript',str(repo/'analysis_v2/scripts/run_clinical_wild_checks.R'),str(task),str(repo)],stdout=stdout,stderr=stderr,check=True)
    require(all(digest(Path(p))==h for p,h in hashes.items()), 'Inputs/code changed')
    write_table(task/'input_hashes.tsv',[dict(path=p,sha256=h) for p,h in sorted(hashes.items())]); finish(task,'PASS_WILD_ENGINEERING_TASK')

def collect(plan_root, results, out):
    verify(plan_root); require(not out.exists(),'Fresh report required')
    draws, actual, failures = [], [], []
    for row in table(plan_root/'tasks.tsv'):
        try:
            task=results/('task_%03d'%int(row['index'])); verify(task)
            require('status\tPASS_WILD_ENGINEERING_TASK' in (task/'SUCCESS').read_text().splitlines(), 'Wrong task status')
            if row['scenario']=='actual':
                actual.extend(dict(context=Path(row['source']).name, **r) for r in table(task/'comparison.tsv'))
            else:
                d=table(task/'draws.tsv'); require(len(d)==150, 'Expected 50 simulations by three methods')
                draws.extend(dict(context=Path(row['source']).name, scenario=row['scenario'], **r) for r in d)
        except (OSError, ValueError, KeyError) as e: failures.append(dict(index=row['index'],error=str(e)))
    out.mkdir(parents=True)
    if draws: write_table(out/'draws.tsv',draws)
    if actual: write_table(out/'actual_engineering_subset.tsv',actual)
    summaries=[]
    for context, scenario, method in sorted({(r['context'],r['scenario'],r['method']) for r in draws}):
        selected=[r for r in draws if (r['context'],r['scenario'],r['method'])==(context,scenario,method)]
        require(len({r['simulation'] for r in selected})==len(selected),'Duplicate simulation')
        n=len(selected); hits=sum(int(r['discoveries'])>0 for r in selected)
        rate=hits/n; z=1.959963984540054; den=1+z*z/n
        center=(rate+z*z/(2*n))/den
        half=z*math.sqrt(rate*(1-rate)/n+z*z/(4*n*n))/den
        summaries.append(dict(context=context,scenario=scenario,method=method,simulations=n,any_discoveries=hits,any_discovery_rate=rate,ci_low=max(0,center-half),ci_high=min(1,center+half),ci_method='WILSON_95_PERCENT',pointwise_rejection_rate=sum(int(r['pointwise_rejections']) for r in selected)/(n*20),family_n=20,bootstrap_draws=9999,interpretation='ENGINEERING_SUBSET_NOT_FULL_FAMILY_VALIDATION'))
    if summaries: write_table(out/'summary.tsv',summaries)
    (out/'status.json').write_text(json.dumps(dict(status='INCOMPLETE' if failures else 'COMPLETE_PENDING_SCIENTIFIC_REVIEW',expected_tasks=156,failures=failures,production_authorized=False),indent=2)+'\n')
    finish(out,'REPORT_ONLY_NOT_AUTHORIZATION')
    if failures: raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser(); s=p.add_subparsers(dest='command',required=True)
    q=s.add_parser('plan'); q.add_argument('--pilot',type=Path,required=True);q.add_argument('--out',type=Path,required=True)
    q=s.add_parser('run');q.add_argument('--plan',type=Path,required=True);q.add_argument('--results',type=Path,required=True);q.add_argument('--index',type=int,required=True);q.add_argument('--repo',type=Path,required=True)
    q=s.add_parser('collect');q.add_argument('--plan',type=Path,required=True);q.add_argument('--results',type=Path,required=True);q.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.command=='plan':plan(a.pilot,a.out)
    elif a.command=='run':run(a.plan,a.results,a.index,a.repo)
    else:collect(a.plan,a.results,a.out)
