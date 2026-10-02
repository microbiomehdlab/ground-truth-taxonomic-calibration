#!/usr/bin/env python3
"""Independent validation tasks; no definitive analysis or policy changes."""
from __future__ import annotations
import argparse
import csv
import fcntl
import json
import math
import os
import subprocess
import tempfile
from decimal import Decimal
from pathlib import Path
from audit_bracken_denominators import digest, table, require, write_table
from prepare_da_pilot import verify
from build_biomarker_abundance_input import parse_profile

INDICES = (2,3,6,7,10,11,14,15,18,19,22,23)
STRATA = [(c,p,b) for c in ('yachida','feng','zeller')
          for p in ('kraken2_bracken','metaphlan4') for b in ('Adenoma','CRC')]
SCENARIOS = [(n,s) for n in (10,42,47,67) for s in
             ('gaussian','correlated_gaussian','sparse_symmetric','skewed_zero_mean')]

def write_matrix(path, rows, family, profiler):
    with path.open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter='\t'); writer.writerow(['observation_id']+family)
        for sid,source,sha in rows:
            require(digest(source)==sha,'Changed native profile')
            values=parse_profile(source,profiler)
            require(values and all(math.isfinite(v) and 0<=v<=1 for v in values.values()),'Invalid native fractions')
            writer.writerow([sid]+[values.get(f,0) for f in family])

def run(repo, output, paired, inventory, index):
    repo,output,paired,inventory=[p.resolve() for p in (repo,output,paired,inventory)]
    require(0<=index<40,'Invalid validation task index')
    for source in (paired,inventory,repo):
        # output may reside inside repository work/, but never inside data sources.
        if source!=repo:
            require(output!=source and output not in source.parents and source not in output.parents,'Overlapping validation/source roots')
    output.mkdir(parents=True,exist_ok=True)
    name='task_%03d'%index
    with (output/(name+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        final=output/name
        require(not final.exists(),'Fresh task output required')
        task=Path(tempfile.mkdtemp(prefix=name+'_attempt_',dir=str(output)))
        hashes={}
        def record(path):
            path=path.resolve(); hashes[path]=digest(path)
        for rel in ('analysis_v2/scripts/run_final_validation_task.py',
                    'analysis_v2/scripts/run_final_validation.R','analysis_v2/lib/unpaired_null.R',
                    'analysis_v2/lib/paired_null.R','analysis_v2/lib/paired_difference_context.R',
                    'analysis_v2/lib/maaslin_contract.R','analysis_v2/lib/maaslin_context.R',
                    'analysis_v2/scripts/build_biomarker_abundance_input.py',
                    'analysis_v2/scripts/prepare_da_pilot.py',
                    'analysis_v2/scripts/audit_bracken_denominators.py'):
            record(repo/rel)
        try:
            if index<12:
                mode='sensitivity'; source=paired/('pilot_%04d'%INDICES[index]); verify(source)
                require((source/'SUCCESS').is_file(),'Missing paired source SUCCESS')
                for p in source.rglob('*'):
                    if p.is_file(): record(p)
                context=table(source/'context.tsv')[0]
                require(context['analysis']=='DA2' and context['inference_method']=='paired_differences','Wrong paired source')
                for member in ('abundance.tsv','metadata.tsv','context.tsv'):
                    (task/member).write_bytes((source/member).read_bytes())
            elif index<24:
                mode='baseline_null'; cohort,profiler,background=STRATA[index-12]
                source=inventory/cohort; verify(source)
                for p in source.iterdir():
                    if p.is_file(): record(p)
                profiles=[r for r in table(source/'profile_inventory.tsv') if
                          r['condition']==background and r['profiler']==profiler and
                          r['analysis_population']=='community' and Decimal(r['nominal_total_dose'])==0]
                require(profiles and len(profiles)==len({r['sample_id'] for r in profiles}),'Invalid baseline pool')
                family=[r['feature'] for r in table(source/'feature_families.tsv') if
                        r['condition']==background and r['profiler']==profiler]
                require(family and len(family)==len(set(family)),'Invalid frozen family')
                allocations=[r for r in table(source/'allocations.tsv') if r['condition']==background]
                keys={(r['pool'],r['n'],r['allocation_id']) for r in allocations}
                for n in ('5','10','15','20'):
                    require(sum(pool=='full' and size==n for pool,size,_ in keys)==1000,'Expected 1000 full allocations per n')
                require(sum(pool=='independent_subset' for pool,_,_ in keys)==252,'Expected 252 subset allocations')
                require(all(r['group'] in ('cases','controls') for r in allocations),'Invalid allocation group')
                write_table(task/'allocations.tsv',allocations)
                matrix_rows=[]
                for r in profiles:
                    p=Path(r['source_profile']).resolve(); record(p)
                    require(hashes[p]==r['sha256'],'Baseline hash mismatch')
                    matrix_rows.append((r['sample_id'],p,r['sha256']))
                write_matrix(task/'abundance.tsv',sorted(matrix_rows),family,profiler)
                write_table(task/'context.tsv',[dict(cohort=cohort,profiler=profiler,background=background)])
            else:
                mode='synthetic'; n,scenario=SCENARIOS[index-24]
                # Maximal observed DA2 family: conservative multiplicity stress size.
                write_table(task/'settings.tsv',[dict(n_pairs=n,scenario=scenario,family_n=3471)])
            with (task/'model.out').open('w') as stdout,(task/'model.err').open('w') as stderr:
                subprocess.run(['Rscript',str(repo/'analysis_v2/scripts/run_final_validation.R'),str(task),mode,str(repo)],
                               stdout=stdout,stderr=stderr,check=True)
            require((task/'summary.tsv').is_file(),'Missing summary')
            require(all(digest(p)==sha for p,sha in hashes.items()),'Inputs/code changed during validation')
            (task/'input_hashes.tsv').write_text('path\tsha256\n'+''.join(str(p)+'\t'+sha+'\n' for p,sha in sorted(hashes.items())))
            write_table(task/'task.tsv',[dict(index=index,mode=mode)])
            (task/'SUCCESS').write_text('status\tPASS_VALIDATION_COMPUTATION\ndefinitive\t0\n')
            members=sorted(p for p in task.rglob('*') if p.is_file())
            (task/'SHA256SUMS').write_text(''.join(digest(p)+'  '+str(p.relative_to(task))+'\n' for p in members))
            os.rename(task,final)
            print('[PASS]',name,mode)
        except Exception as e:
            (task/'FAILED.json').write_text(json.dumps(dict(error=str(e))))
            raise

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('repo','output','paired','inventory'): p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--index',type=int,required=True)
    a=p.parse_args(); run(a.repo,a.output,a.paired,a.inventory,a.index)
