#!/usr/bin/env python3
"""Collect all validation outcomes, including failures; never launch production."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from prepare_da_pilot import verify
from audit_bracken_denominators import table, write_table, digest

def collect(root):
    rows=[]; missing=[]; flags=[]; changes=[]
    for i in range(40):
        folder=root/('task_%03d'%i)
        try:
            verify(folder)
            if not (folder/'SUCCESS').is_file(): raise ValueError('No SUCCESS')
            mode=table(folder/'task.tsv')[0]['mode']
            context=table(folder/'context.tsv')[0] if (folder/'context.tsv').is_file() else {}
            for summary in table(folder/'summary.tsv'):
                row=dict(task=i,mode=mode,**context,**summary); rows.append(row)
                if summary.get('review_flag')=='TRUE': flags.append(row)
            if mode=='sensitivity':
                results=table(folder/'sensitivity_results.tsv')
                primary={r['feature']:r for r in results if float(r['pseudocount'])==1e-8}
                for r in results:
                    base=primary[r['feature']]
                    if r['positive_discovery']!=base['positive_discovery'] or r['status']!=base['status']:
                        changes.append(dict(task=i,context_id=context['context_id'],cohort=context['cohort'],
                            profiler=context['profiler'],feature=r['feature'],pseudocount=r['pseudocount'],
                            primary_q=base['wrapper_q'],sensitivity_q=r['wrapper_q'],
                            primary_discovery=base['positive_discovery'],sensitivity_discovery=r['positive_discovery'],
                            primary_status=base['status'],sensitivity_status=r['status']))
        except (OSError,ValueError,KeyError) as e:
            missing.append(dict(task=i,error=str(e)))
    report=root/'REPORT'
    report.mkdir(exist_ok=False)
    if rows:
        fields=sorted({k for row in rows for k in row})
        write_table(report/'validation_summary.tsv',[{k:r.get(k,'') for k in fields} for r in rows])
    if changes: write_table(report/'pseudocount_changes.tsv',changes)
    status='INCOMPLETE' if missing else ('REVIEW_FLAGS' if flags else 'COMPLETE_PENDING_SCIENTIFIC_REVIEW')
    result=dict(status=status,expected_tasks=40,completed_tasks=40-len(missing),
                review_flags=len(flags),pseudocount_changed_rows=len(changes),failures=missing,production_authorized=False)
    (report/'status.json').write_text(json.dumps(result,indent=2)+'\n')
    (report/'README.txt').write_text(
        'No automatic production authorization. Review null flags, target pseudocount changes,\n'
        'non-estimability and retained backend diagnostics. Synthetic skewness is a stress\n'
        'scenario, not a claim that actual differences have that distribution.\n')
    (report/'SHA256SUMS').write_text(''.join(digest(p)+'  '+p.name+'\n' for p in sorted(report.iterdir()) if p.is_file()))
    print(json.dumps(result,indent=2))
    print('Combined report:',report)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('root',type=Path)
    collect(p.parse_args().root.resolve())
