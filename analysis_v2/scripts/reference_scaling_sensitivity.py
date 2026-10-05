#!/usr/bin/env python3
"""Prespecified half/double MetaPhlAn conversion stress test; not biological truth."""
from __future__ import annotations
import argparse
import json
import shutil
import subprocess
import sys
from collections import Counter
from decimal import Decimal
from pathlib import Path
from audit_bracken_denominators import digest, require, table, write_table
from prepare_da_pilot import verify
from reference_response_pilot import CODE, fresh, seal, expected_profile, collect

VARIANTS = {'scale_050':.5, 'scale_200':2.0}
TOOL = 'metaphlan4'

def prepare(base, out, image=None):
    base,out=base.resolve(),out.resolve();fresh(out,[base]);verify(base);verify(base/'REPORT')
    report=json.loads((base/'REPORT/status.json').read_text())
    require(report['contexts']==720 and report['target_rows']==7200,'Require completed original pilot')
    require(report['status']=='COMPLETE_PENDING_SCIENTIFIC_REVIEW','Unexpected base report state')
    require(json.loads((base/'REPORT/collector_provenance.json').read_text())['plan_sha256']==digest(base/'SHA256SUMS'),
            'Base report identity differs')
    state=json.loads((base/'identity.json').read_text())
    require(digest(Path(state['image_path']))==state['image_sha256'],'Original image changed')
    if image is not None:require(image.resolve()==Path(state['image_path']).resolve(),'Runtime image path differs')
    all_contexts=json.loads((base/'contexts.json').read_text())
    require(len(all_contexts)==720 and len({c['context_id'] for c in all_contexts})==720,'Base context coverage differs')
    original_rows=table(base/'REPORT/target_comparisons.tsv')
    require(len(original_rows)==7200 and len({(r['context_id'],r['target_label']) for r in original_rows})==7200,
            'Base target report coverage differs')
    contexts=[c for c in all_contexts if c['profiler']==TOOL]
    require(len(contexts)==360 and len({c['context_id'] for c in contexts})==360,'Incomplete MetaPhlAn pilot')
    cells=Counter((c['cohort'],c['n'],c['anchor']) for c in contexts)
    require(len(cells)==18 and set(cells.values())=={20} and all(c['arm']=='U' for c in contexts),
            'Unexpected scope')
    catalog=json.loads((base/'catalog.json').read_text())
    native=json.loads((base/'native_profiles.json').read_text())
    audits={r['key']:r for r in table(base/'reference_audit.tsv') if r['profiler']==TOOL}
    observed=[r for r in table(base/'observed_targets.tsv') if r['context_id'] in {c['context_id'] for c in contexts}]
    require(len(observed)==3600 and len({(r['context_id'],r['target_label']) for r in observed})==3600,
            'Observed coverage differs')
    staged={};construction=[]
    # Compute and validate all scenarios before publishing any new plan.
    for scenario,multiplier in VARIANTS.items():
        references={}
        for c in contexts:
            targets=catalog['targets'][TOOL]
            require(len(targets)==10 and len(set(targets.values()))==10,'Target mapping differs')
            for obs in c['observations']:
                key=obs['key']
                if key in references:continue
                baseline_key='|'.join((c['cohort'],TOOL,obs['sample_id'],'0'))
                baseline=native[catalog['profiles'][baseline_key]['source_profile']]
                if Decimal(obs['dose'])==0:references[key]=baseline;continue
                a=audits[key]
                require(a['cohort']==c['cohort'] and a['sample_id']==obs['sample_id'],'Audit identity differs')
                counts=json.loads(a['inserted_counts_json']);sizes=json.loads(a['target_genome_sizes_json'])
                original,inserted=int(a['original_pairs']),int(a['inserted_pairs'])
                require(original>0 and inserted>0 and sum(counts.values())==inserted,'Inserted counts differ')
                fraction=inserted/(original+inserted)
                require(abs(fraction-float(a['achieved_total_fraction']))<=1e-12,'Achieved dose differs')
                geff=float(a['effective_genome_size_bp'])*multiplier
                expected,scale=expected_profile(baseline,TOOL,counts,fraction,targets,geff=geff,genome_sizes=sizes)
                references[key]=expected
                construction.append(dict(scenario=scenario,multiplier=multiplier,key=key,
                    baseline_geff_bp=a['effective_genome_size_bp'],stress_geff_bp=geff,
                    achieved_fraction=fraction,retained_scale=scale,reference_sum=sum(expected.values())))
        staged[scenario]=references
    verify(base);verify(base/'REPORT')
    out.mkdir(parents=True)
    # Freeze orchestration/collection as well as the original calculation worker.
    for name in CODE+('scripts/reference_scaling_sensitivity.py',):
        src=Path(__file__).resolve().parents[1]/name;dest=out/'source/analysis_v2'/name
        dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
    for scenario,references in staged.items():
        folder=out/scenario;folder.mkdir()
        for name in ('catalog.json','native_profiles.json'):
            shutil.copy2(base/name,folder/name)
        for name in CODE:
            source=base/'source/analysis_v2'/name;dest=folder/'source/analysis_v2'/name
            dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dest)
            require(digest(source)==digest(dest),'Calculation snapshot copy differs')
        (folder/'contexts.json').write_text(json.dumps(contexts))
        (folder/'reference_profiles.json').write_text(json.dumps(references,allow_nan=False))
        write_table(folder/'observed_targets.tsv',observed)
        tasks=[]
        for i in range(0,len(contexts),10):
            name='batches/batch_%04d.json'%(i//10);path=folder/name;path.parent.mkdir(exist_ok=True)
            (path).write_text(json.dumps(contexts[i:i+10]))
            tasks.append(dict(index=i//10,file=name,contexts=len(contexts[i:i+10])))
        write_table(folder/'tasks.tsv',tasks)
        (folder/'identity.json').write_text(json.dumps(dict(image_path=state['image_path'],image_sha256=state['image_sha256'],
            contexts=360,scenario=scenario,multiplier=VARIANTS[scenario],base_plan_sha256=digest(base/'SHA256SUMS'),
            scope='GLOBAL_CONVERSION_STRESS_NOT_ESTIMATED_UNCERTAINTY',production_authorized=False)))
        seal(folder,'PASS_SCALING_STRESS_PLAN_NOT_PRODUCTION_AUTHORIZATION')
    write_table(out/'construction.tsv',construction)
    (out/'identity.json').write_text(json.dumps(dict(base_root=str(base),base_plan_sha256=digest(base/'SHA256SUMS'),
        base_report_sha256=digest(base/'REPORT/SHA256SUMS'),image_sha256=state['image_sha256'],
        scenarios=VARIANTS,new_contexts=720,batches=72,
        preparation_script_sha256=digest(Path(__file__)),production_authorized=False),indent=2)+'\n')
    seal(out,'PASS_SENSITIVITY_PREPARATION_NOT_PRODUCTION_AUTHORIZATION')
    print('[PASS] 720 MetaPhlAn stress contexts / 72 batches; original outputs untouched')

def run(root,index,image=None):
    verify(root);require(0<=index<72,'Invalid sensitivity index')
    scenario=list(VARIANTS)[index//36]
    folder=root/scenario
    if image is not None:
        state=json.loads((folder/'identity.json').read_text())
        require(image.resolve()==Path(state['image_path']).resolve() and digest(image)==state['image_sha256'],
                'Runtime image differs')
    subprocess.run([sys.executable,'-B',str(folder/'source/analysis_v2/scripts/reference_response_pilot.py'),
                    'run','--root',str(folder),'--index',str(index%36)],check=True)

def report(root):
    verify(root)
    identity=json.loads((root/'identity.json').read_text());base=Path(identity['base_root'])
    verify(base/'REPORT');require(digest(base/'REPORT/SHA256SUMS')==identity['base_report_sha256'],'Base report changed')
    reports={}
    for scenario in VARIANTS:
        folder=root/scenario;path=folder/'REPORT'
        if not path.exists():collect(folder,path)
        verify(path)
        require(json.loads((path/'collector_provenance.json').read_text())['plan_sha256']==digest(folder/'SHA256SUMS'),
                'Variant report identity differs')
        reports[scenario]=path
    out=root/'REPORT';fresh(out,[])
    original=[r for r in table(base/'REPORT/target_comparisons.tsv') if r['profiler']==TOOL]
    require(len(original)==3600,'Original MetaPhlAn coverage differs')
    original_index={(r['context_id'],r['target_label']):r for r in original}
    merged=[]
    for scenario,path in [('scale_100',base/'REPORT')]+list(reports.items()):
        rows=[r for r in table(path/'target_comparisons.tsv') if r['profiler']==TOOL]
        require(len(rows)==3600 and {(r['context_id'],r['target_label']) for r in rows}==set(original_index),
                'Scenario target coverage differs')
        for r in rows:
            old=original_index[r['context_id'],r['target_label']]
            require(all(r[k]==old[k] for k in r if k.startswith('observed_')),'Observed result changed')
            merged.append(dict(r,scenario=scenario,multiplier=1 if scenario=='scale_100' else VARIANTS[scenario]))
    counts=Counter((r['scenario'],r['cohort'],r['n'],r['nominal_total_dose'],r['target_label'],r['category']) for r in merged)
    summary=[dict(zip(('scenario','cohort','n','nominal_total_dose','target_label','category'),key),count=count,allocations=20,
                  conditional_frequency=count/20) for key,count in sorted(counts.items())]
    out.mkdir();write_table(out/'target_comparisons.tsv',merged);write_table(out/'category_frequencies.tsv',summary)
    (out/'status.json').write_text(json.dumps(dict(status='COMPLETE_SCALING_STRESS_PENDING_REVIEW',
        new_reference_contexts=720,target_rows=len(merged),observed_refits=0,production_authorized=False,
        interpretation='HALF_DOUBLE_GLOBAL_CONVERSION_STRESS_NOT_CONFIDENCE_BOUNDS'),indent=2)+'\n')
    (out/'provenance.json').write_text(json.dumps(dict(identity,collector_script_sha256=digest(Path(__file__))),indent=2)+'\n')
    seal(out,'PASS_SCALING_SENSITIVITY_NOT_PRODUCTION_AUTHORIZATION')
    print('[PASS] All scenarios collected; review effects/raw p/full-family q and NE separately')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
    q=s.add_parser('prepare');q.add_argument('--base',type=Path,required=True);q.add_argument('--out',type=Path,required=True);q.add_argument('--image',type=Path,required=True)
    q=s.add_parser('run');q.add_argument('--root',type=Path,required=True);q.add_argument('--index',type=int,required=True);q.add_argument('--image',type=Path,required=True)
    q=s.add_parser('collect');q.add_argument('--root',type=Path,required=True)
    a=p.parse_args()
    if a.command=='prepare':prepare(a.base,a.out,a.image)
    elif a.command=='run':run(a.root,a.index,a.image)
    else:report(a.root)
