#!/usr/bin/env python3
"""Account for every batch; do not convert absent tasks into nondiscoveries."""
import argparse,csv,gzip,json
from pathlib import Path
from audit_bracken_denominators import table,require
from prepare_da_pilot import verify
from plan_da3_canary import finish

def collect(plan,results,out):
    plan,results,out=plan.resolve(),results.resolve(),out.resolve();verify(plan)
    require(not out.exists() and all(out!=p and out not in p.parents and p not in out.parents for p in (plan,results)),'Fresh nonoverlapping report required')
    tasks=table(plan/'tasks.tsv');errors=[];completed=0;contexts=0;out.mkdir(parents=True)
    with (out/'context_summary.tsv').open('w',newline='') as summary_handle,gzip.open(out/'targets.tsv.gz','wt',newline='') as target_handle:
        summary_writer=target_writer=None
        for task in tasks:
            try:
                folder=results/('batch_%05d'%int(task['index']));verify(folder)
                from audit_bracken_denominators import digest
                require(json.loads((folder/'resume_identity.json').read_text())['plan_sha256']==digest(plan/'SHA256SUMS'),'Batch belongs to different plan')
                rows=table(folder/'summary.tsv')
                expected=json.loads((plan/task['file']).read_text())
                require({r['context_id'] for r in rows}=={r['context_id'] for r in expected} and len(rows)==len(expected),'Context coverage mismatch')
                with gzip.open(folder/'targets.tsv.gz','rt',newline='') as handle:
                    targets=list(csv.DictReader(handle,delimiter='\t'))
                require(len(targets)==len(expected)*10 and len({(r['context_id'],r['target_label']) for r in targets})==len(targets),'Target count/identity mismatch')
                require({r['context_id'] for r in targets}=={r['context_id'] for r in expected},'Target context mismatch')
                if summary_writer is None:summary_writer=csv.DictWriter(summary_handle,fieldnames=list(rows[0]),delimiter='\t');summary_writer.writeheader()
                # Exact and Monte Carlo schemas differ: canonical union, explicit blanks.
                fields=['context_id','target_label','feature','beta','variable','exact_p','exact_bh_q','mc_p','mc_bh_q','max_statistic_fwer_p','raw_exceedances','max_exceedances','resamples','permutations','minimum_mc_p','seed','parametric_p','parametric_bh_q']
                if target_writer is None:target_writer=csv.DictWriter(target_handle,fieldnames=fields,delimiter='\t');target_writer.writeheader()
                require(all(set(r)<=set(fields) for r in targets),'Unexpected target schema')
                summary_writer.writerows(rows);target_writer.writerows([{k:r.get(k,'') for k in fields} for r in targets])
                completed+=1;contexts+=len(rows)
            except (OSError,ValueError,KeyError) as e:errors.append(dict(batch=task['index'],error=str(e)))
    (out/'status.json').write_text(json.dumps(dict(status='INCOMPLETE' if errors else 'COMPLETE_PENDING_SCIENTIFIC_REVIEW',expected_batches=len(tasks),completed_batches=completed,completed_contexts=contexts,failures=errors,production_authorized=False),indent=2)+'\n')
    finish(out,'INCOMPLETE' if errors else 'PASS_COMPLETE_DA3_BATCH_COLLECTION')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for flag in ('plan','results','out'):p.add_argument('--'+flag,type=Path,required=True)
    a=p.parse_args();collect(a.plan,a.results,a.out)
