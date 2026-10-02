#!/usr/bin/env python3
from __future__ import annotations
import argparse
import fcntl
import os
import subprocess
import tempfile
from pathlib import Path
from prepare_da_pilot import verify
from audit_bracken_denominators import table,digest,require,write_table

def run(repo,root,validation,pilot,plan,index):
    repo,root,validation,pilot,plan=[p.resolve() for p in (repo,root,validation,pilot,plan)]
    require(0<=index<52,'Invalid candidate index')
    for source in (validation,pilot,plan):
        require(root!=source and root not in source.parents and source not in root.parents,'Overlapping output')
    root.mkdir(parents=True,exist_ok=True);name='candidate_%03d'%index
    with (root/(name+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX);require(not(root/name).exists(),'Fresh output required')
        task=Path(tempfile.mkdtemp(prefix=name+'_attempt_',dir=str(root)));hashes={}
        def record(p): hashes[p.resolve()]=digest(p)
        for rel in ('analysis_v2/scripts/run_da3_candidate.py','analysis_v2/scripts/run_da3_candidate.R',
                    'analysis_v2/lib/exact_da3_candidate.R','analysis_v2/lib/unpaired_null.R',
                    'analysis_v2/scripts/prepare_da_pilot.py','analysis_v2/scripts/audit_bracken_denominators.py'):
            record(repo/rel)
        if index<12:
            mode='null';source=validation/('task_%03d'%(index+12));verify(source)
            require((source/'SUCCESS').is_file(),'Source validation did not complete')
            allocations=[r for r in table(source/'allocations.tsv') if r['n']=='5']
            keys={(r['pool'],r['allocation_id']) for r in allocations}
            require(sum(pool=='full' for pool,_ in keys)==1000 and
                    sum(pool=='independent_subset' for pool,_ in keys)==252,'Incomplete null allocations')
            write_table(task/'allocations.tsv',allocations)
            members=('abundance.tsv','context.tsv')
        else:
            mode='positive';verify(plan)
            contexts=[r for r in table(plan/'contexts.tsv') if r['analysis']=='DA3']
            require(len(contexts)==40 and len({r['context_id'] for r in contexts})==40,'Expected 40 unique DA3 pilot contexts')
            source=pilot/contexts[index-12]['context_id'];verify(source)
            require((source/'SUCCESS').is_file(),'Source positive pilot did not complete')
            context=table(source/'context.tsv')[0]
            require(context['analysis']=='DA3' and context['paired']=='0' and
                    context['context_id']==contexts[index-12]['context_id'],'Wrong positive context')
            members=('abundance.tsv','context.tsv','metadata.tsv')
            record(plan/'SHA256SUMS')
        for p in source.rglob('*'):
            if p.is_file():record(p)
        for member in members:(task/member).write_bytes((source/member).read_bytes())
        with (task/'model.out').open('w') as stdout,(task/'model.err').open('w') as stderr:
            subprocess.run(['Rscript',str(repo/'analysis_v2/scripts/run_da3_candidate.R'),str(task),mode,str(repo)],
                           stdout=stdout,stderr=stderr,check=True)
        require(all(digest(p)==sha for p,sha in hashes.items()),'Source/code changed during candidate test')
        write_table(task/'input_hashes.tsv',[dict(path=str(p),sha256=sha) for p,sha in sorted(hashes.items())])
        write_table(task/'task.tsv',[dict(index=index,mode=mode)])
        (task/'SUCCESS').write_text('status\tPASS_CANDIDATE_COMPUTATION\nproduction_authorized\t0\n')
        (task/'SHA256SUMS').write_text(''.join(digest(p)+'  '+str(p.relative_to(task))+'\n'
            for p in sorted(task.rglob('*')) if p.is_file()))
        os.rename(task,root/name);print('[PASS]',name,mode)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for flag in ('repo','root','validation','pilot','plan'):p.add_argument('--'+flag,type=Path,required=True)
    p.add_argument('--index',type=int,required=True);a=p.parse_args()
    run(a.repo,a.root,a.validation,a.pilot,a.plan,a.index)
