#!/usr/bin/env python3
from __future__ import annotations
import argparse
import fcntl
import hashlib
import os
import subprocess
import tempfile
from pathlib import Path
from prepare_da_pilot import verify
from audit_bracken_denominators import table,digest,require,write_table

def run(repo,root,validation,index):
    repo,root,validation=[p.resolve() for p in (repo,root,validation)]
    require(0<=index<388,'Invalid robustness index')
    require(root!=validation and root not in validation.parents and validation not in root.parents,'Overlapping output')
    root.mkdir(parents=True,exist_ok=True);name='robust_%03d'%index
    with (root/(name+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX);require(not(root/name).exists(),'Fresh output required')
        task=Path(tempfile.mkdtemp(prefix=name+'_attempt_',dir=str(root)));hashes={}
        def record(p):hashes[p.resolve()]=digest(p)
        for rel in ('analysis_v2/scripts/run_da_robustness.py','analysis_v2/scripts/run_da_robustness.R',
                    'analysis_v2/lib/monte_carlo_da3_candidate.R','analysis_v2/lib/paired_direction_candidate.R',
                    'analysis_v2/lib/unpaired_null.R','analysis_v2/lib/maaslin_contract.R',
                    'analysis_v2/lib/paired_difference_context.R','analysis_v2/scripts/prepare_da_pilot.py',
                    'analysis_v2/scripts/audit_bracken_denominators.py'):record(repo/rel)
        if index<360:
            mode='larger_null';stratum=index//10;block=index%10
            source=validation/('task_%03d'%(12+stratum//3));verify(source)
            context=table(source/'context.tsv')[0];n=(10,15,20)[stratum%3]
            all_rows=[r for r in table(source/'allocations.tsv') if r['pool']=='full' and int(r['n'])==n]
            keys=sorted({r['allocation_id'] for r in all_rows})
            require(len(keys)==1000,'Expected 1000 larger-n allocations')
            chosen=set(keys[block*100:(block+1)*100]);rows=[dict(r) for r in all_rows if r['allocation_id'] in chosen]
            for r in rows:
                namespace='|'.join(('da3_mc_20261002',context['cohort'],context['background'],str(n),r['allocation_id']))
                r['seed']=str(int.from_bytes(hashlib.sha256(namespace.encode()).digest()[:4],'big')%2147483646+1)
            write_table(task/'allocations.tsv',rows)
            context.update(n=n,block=block,resamples=9999);write_table(task/'context.tsv',[context])
            (task/'abundance.tsv').write_bytes((source/'abundance.tsv').read_bytes())
        elif index<372:
            mode='paired_direction';source=validation/('task_%03d'%(index-360));verify(source)
            context=table(source/'context.tsv')[0]
            require(context['analysis']=='DA2','Expected DA2 source')
            for member in ('abundance.tsv','metadata.tsv','context.tsv'):
                (task/member).write_bytes((source/member).read_bytes())
        else:
            mode='direction_synthetic';n=(10,42,47,67)[(index-372)//4]
            scenario=('gaussian','sparse_symmetric','skewed_sign_null','skewed_mean_null')[(index-372)%4]
            write_table(task/'settings.tsv',[dict(n=n,scenario=scenario)])
            write_table(task/'context.tsv',[dict(n=n,scenario=scenario)])
            source=None
        if source is not None:
            require((source/'SUCCESS').is_file(),'Missing source SUCCESS')
            for p in source.rglob('*'):
                if p.is_file():record(p)
        with (task/'model.out').open('w') as stdout,(task/'model.err').open('w') as stderr:
            subprocess.run(['Rscript',str(repo/'analysis_v2/scripts/run_da_robustness.R'),str(task),mode,str(repo)],
                           stdout=stdout,stderr=stderr,check=True)
        require(all(digest(p)==sha for p,sha in hashes.items()),'Source/code changed')
        write_table(task/'input_hashes.tsv',[dict(path=str(p),sha256=sha) for p,sha in sorted(hashes.items())])
        write_table(task/'task.tsv',[dict(index=index,mode=mode)])
        (task/'SUCCESS').write_text('status\tPASS_ROBUSTNESS_COMPUTATION\nproduction_authorized\t0\n')
        (task/'SHA256SUMS').write_text(''.join(digest(p)+'  '+str(p.relative_to(task))+'\n'
            for p in sorted(task.rglob('*')) if p.is_file()))
        os.rename(task,root/name);print('[PASS]',name,mode)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for flag in ('repo','root','validation'):p.add_argument('--'+flag,type=Path,required=True)
    p.add_argument('--index',type=int,required=True);a=p.parse_args();run(a.repo,a.root,a.validation,a.index)
