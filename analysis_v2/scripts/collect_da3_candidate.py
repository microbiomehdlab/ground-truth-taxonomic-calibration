#!/usr/bin/env python3
from __future__ import annotations
import argparse
import json
from pathlib import Path
from prepare_da_pilot import verify
from audit_bracken_denominators import table,write_table,digest

def collect(root):
    rows=[];failures=[]
    for i in range(52):
        folder=root/('candidate_%03d'%i)
        try:
            verify(folder)
            context=table(folder/'context.tsv')[0]
            for summary in table(folder/'summary.tsv'):
                row={'task':i,'mode':'null' if i<12 else 'positive'}
                row.update(context);row.update(summary);rows.append(row)
        except (OSError,ValueError,KeyError) as e:failures.append(dict(task=i,error=str(e)))
    report=root/'REPORT';report.mkdir(exist_ok=False)
    fields=sorted({key for row in rows for key in row})
    if rows:write_table(report/'comparison_summary.tsv',[{k:r.get(k,'') for k in fields} for r in rows])
    status=dict(status='INCOMPLETE' if failures else 'COMPLETE_PENDING_METHOD_REVIEW',
                completed_tasks=52-len(failures),expected_tasks=52,failures=failures,production_authorized=False)
    (report/'status.json').write_text(json.dumps(status,indent=2)+'\n')
    (report/'SHA256SUMS').write_text(''.join(digest(p)+'  '+p.name+'\n' for p in sorted(report.iterdir()) if p.is_file()))
    print(json.dumps(status,indent=2));print('Report:',report)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);collect(p.parse_args().root)
