#!/usr/bin/env python3
import argparse,json
from pathlib import Path
from remaining_validation import collect
from prepare_da_pilot import verify
from audit_bracken_denominators import table,write_table,require
from plan_da3_canary import finish

def report(root):
    errors=[];rows=[]
    try:collect(root/'validation_plan',root/'validation_results',root/'validation_report')
    except (OSError,ValueError,KeyError) as e:errors.append(dict(stage='validation_collector',error=str(e)))
    status_path=root/'validation_report/status.json'
    if status_path.exists():errors.extend(json.loads(status_path.read_text())['failures'])
    for i in range(48):
        index=(i//8)*800+i%8;folder=root/'batch_canary_results'/('batch_%05d'%index)
        try:
            verify(folder)
            require('status\tPASS_BATCH_COMPUTATION_NOT_PRODUCTION_AUTHORIZATION' in (folder/'SUCCESS').read_text().splitlines(),'Wrong batch status')
            data=table(folder/'summary.tsv');require(len(data)==25,'Incomplete batch')
            require(len({r['context_id'] for r in data})==25,'Duplicate context')
            rows.extend(data)
        except (OSError,ValueError,KeyError) as e:errors.append(dict(batch=index,error=str(e)))
    out=root/'REPORT';require(not out.exists(),'Fresh combined report required');out.mkdir()
    if rows:write_table(out/'batch_timing_and_discoveries.tsv',rows)
    (out/'status.json').write_text(json.dumps(dict(status='INCOMPLETE' if errors else 'COMPLETE_PENDING_METHOD_AND_RESOURCE_REVIEW',expected_validation_tasks=72,expected_batch_canary_tasks=48,completed_batch_contexts=len(rows),failures=errors,production_authorized=False),indent=2)+'\n')
    finish(out,'PREPRODUCTION_REPORT_NOT_AUTHORIZATION');print((out/'status.json').read_text())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);report(p.parse_args().root)
