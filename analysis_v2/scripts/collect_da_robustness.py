#!/usr/bin/env python3
from __future__ import annotations
import argparse
import json
import math
from collections import defaultdict
from pathlib import Path
from prepare_da_pilot import verify
from audit_bracken_denominators import table,write_table,digest

def collect(root):
    summaries=[];failures=[];null=defaultdict(list);completed=0
    for i in range(388):
        folder=root/('robust_%03d'%i)
        try:
            verify(folder);context=table(folder/'context.tsv')[0];mode=table(folder/'task.tsv')[0]['mode']
            if mode=='larger_null':
                key=(context['cohort'],context['profiler'],context['background'],context['n'])
                null[key].extend(table(folder/'draws.tsv'))
            else:
                for summary in table(folder/'summary.tsv'):
                    row=dict(task=i,mode=mode);row.update(context);row.update(summary);summaries.append(row)
            completed+=1
        except (OSError,ValueError,KeyError) as e:failures.append(dict(task=i,error=str(e)))
    for key,draws in sorted(null.items()):
        ids=[r['allocation_id'] for r in draws]
        if len(ids)!=len(set(ids)):raise ValueError('Duplicate larger-null allocation')
        for method in ('parametric_bh','mc_bh','mc_max_fwer'):
            count=sum(int(r[method])>0 for r in draws);total=len(draws);rate=count/total
            z=1.959963984540054;denominator=1+z*z/total
            center=(rate+z*z/(2*total))/denominator
            half=z*math.sqrt(rate*(1-rate)/total+z*z/(4*total*total))/denominator
            summaries.append(dict(mode='larger_null',cohort=key[0],profiler=key[1],background=key[2],n=key[3],
                method=method,repetitions=len(draws),expected_repetitions=1000,
                any_discovery_rate=rate,mc_wilson_low=max(0,center-half),mc_wilson_high=min(1,center+half),
                mean_discoveries=sum(int(r[method]) for r in draws)/len(draws),
                interpretation='MC_BH_RESOLUTION_LIMITED' if method=='mc_bh' else 'GLOBAL_EXCHANGEABLE_LABEL_NULL'))
        if len(draws)!=1000 and not failures:
            failures.append(dict(stratum=list(key),error='Expected 1000 complete null allocations'))
    report=root/'REPORT';report.mkdir(exist_ok=False)
    if summaries:
        fields=sorted({k for r in summaries for k in r});write_table(report/'comparison_summary.tsv',[{k:r.get(k,'') for k in fields} for r in summaries])
    status=dict(status='INCOMPLETE' if failures else 'COMPLETE_PENDING_ESTIMAND_AND_METHOD_REVIEW',
                expected_tasks=388,completed_tasks=completed,failures=failures,production_authorized=False)
    (report/'status.json').write_text(json.dumps(status,indent=2)+'\n')
    (report/'SHA256SUMS').write_text(''.join(digest(p)+'  '+p.name+'\n' for p in sorted(report.iterdir()) if p.is_file()))
    print(json.dumps(status,indent=2));print('Report:',report)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);collect(p.parse_args().root)
