#!/usr/bin/env python3
"""One locked, resumable engineering model context per array task."""
from __future__ import annotations
import argparse
import csv
import fcntl
import json
import math
import os
import subprocess
import tempfile
import time
from pathlib import Path
from audit_bracken_denominators import table, digest, require, write_table
from build_biomarker_abundance_input import parse_profile
from prepare_da_pilot import verify


def run(plan,output,index,repo,method='maaslin2'):
    plan,output,repo=[p.resolve() for p in (plan,output,repo)]
    require(output!=plan and output not in plan.parents and plan not in output.parents,'Overlapping output')
    verify(plan)
    contexts=table(plan/'contexts.tsv')
    require(0<=index<len(contexts),'Invalid array index')
    context=contexts[index]; cid=context['context_id']
    require(method in ('maaslin2','paired_differences'),'Invalid inference method')
    if method=='paired_differences':
        require(context['analysis']=='DA2' and context['paired']=='1','Paired differences require paired DA2')
    runner=repo/('analysis_v2/scripts/run_paired_difference_context.R' if method=='paired_differences'
                 else 'analysis_v2/scripts/run_da_pilot_context.R')
    library=repo/('analysis_v2/lib/paired_difference_context.R' if method=='paired_differences'
                  else 'analysis_v2/lib/maaslin_context.R')
    require(Path(cid).name==cid and cid not in ('','.','..'),'Unsafe context ID')
    provenance={str(p):digest(p) for p in (plan/'SHA256SUMS',library,
                  repo/'analysis_v2/lib/maaslin_contract.R',Path(__file__).resolve(),
                  repo/'analysis_v2/scripts/build_biomarker_abundance_input.py',
                  runner)}
    output.mkdir(parents=True,exist_ok=True)
    (output/'locks').mkdir(exist_ok=True)
    with (output/'locks'/(cid+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        final=output/cid
        if final.exists():
            verify(final)
            require((final/'SUCCESS').is_file(),'Existing context has no SUCCESS')
            require(json.loads((final/'provenance.json').read_text())==provenance,'Existing result belongs to different inputs/code')
            print('[PASS] Already completed:',cid); return
        attempts=output/'attempts'; attempts.mkdir(exist_ok=True)
        task=Path(tempfile.mkdtemp(prefix=cid+'_',dir=str(attempts)))
        start=time.time()
        try:
            family=[r['feature'] for r in table(plan/'families.tsv') if r['context_id']==cid]
            rows=[r for r in table(plan/'observations.tsv') if r['context_id']==cid]
            require(len(family)==int(context['family_n']) and len(rows)==int(context['observations_n']),'Context size mismatch')
            require(len(family)==len(set(family)) and len(rows)==len({r['observation_id'] for r in rows}),'Duplicate observation/feature')
            write_table(task/'context.tsv',[dict(context,inference_method=method)])
            write_table(task/'metadata.tsv',[{k:r[k] for k in ('observation_id','biological_sample_id','group','spike_state','age','sex')} for r in rows])
            with (task/'abundance.tsv').open('w',newline='') as h:
                writer=csv.writer(h,delimiter='\t'); writer.writerow(['observation_id']+family)
                for r in rows:
                    path=Path(r['source_profile'])
                    require(digest(path)==r['sha256'],'Changed native profile')
                    values=parse_profile(path,context['profiler'])
                    require(values and all(math.isfinite(v) and 0<=v<=1 for v in values.values()),'Invalid native abundances')
                    writer.writerow([r['observation_id']]+[values.get(f,0) for f in family])
            with (task/'model.out').open('w') as stdout,(task/'model.err').open('w') as stderr:
                subprocess.run(['Rscript',str(runner),str(task),str(repo)],
                               stdout=stdout,stderr=stderr,check=True)
            require((task/'fit/SUCCESS').is_file(),'Model did not pass')
            require(all(digest(Path(r['source_profile']))==r['sha256'] for r in rows),
                    'Native profile changed during fit')
            require(all(digest(Path(p))==sha for p,sha in provenance.items()),'Code/plan changed during fit')
            (task/'provenance.json').write_text(json.dumps(provenance,sort_keys=True))
            (task/'runtime.json').write_text(json.dumps(dict(elapsed_seconds=time.time()-start)))
            (task/'SUCCESS').write_text('status\tPASS_ENGINEERING_CONTEXT\ndefinitive\t0\n')
            members=sorted(p for p in task.rglob('*') if p.is_file())
            (task/'SHA256SUMS').write_text(''.join(digest(p)+'  '+str(p.relative_to(task))+'\n' for p in members))
            os.rename(task,final)
            print('[PASS] Engineering context:',cid)
        except Exception as e:
            (task/'FAILED.json').write_text(json.dumps(dict(error=str(e),elapsed_seconds=time.time()-start)))
            print('[FAIL] Retained diagnostics:',task)
            raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('plan','output','repo'): p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--index',type=int,required=True)
    p.add_argument('--method',choices=('maaslin2','paired_differences'),default='maaslin2')
    a=p.parse_args(); run(a.plan,a.output,a.index,a.repo,a.method)
