#!/usr/bin/env python3
"""Inspect saved native inputs for NE targets; no profiling or model fitting."""
from __future__ import annotations
import argparse
import json
import math
import sys
from collections import Counter
from pathlib import Path

from audit_bracken_denominators import digest, require, table, write_table
from prepare_da_pilot import verify
from reference_response_pilot import fresh, seal

TOLERANCE = 100 * sys.float_info.epsilon

def describe(values, groups):
    require(len(values)==len(groups) and len(values)>2, 'Invalid vector length')
    require(set(groups)=={0,1}, 'Both artificial groups required')
    require(all(math.isfinite(v) and v>=0 for v in values), 'Invalid native abundance')
    transformed=[math.log2(1+v/1e-8) for v in values]
    tolerance=TOLERANCE*max(1,max(map(abs,transformed)))
    span=max(transformed)-min(transformed)
    means={g:math.fsum(v for v,k in zip(transformed,groups) if k==g)/groups.count(g) for g in (0,1)}
    residual=max(abs(v-means[g]) for v,g in zip(transformed,groups))
    status=('NON_ESTIMABLE_CONSTANT' if span<=tolerance else
            'NON_ESTIMABLE_PERFECT_FIXED_FIT' if residual<=tolerance else 'ESTIMABLE')
    result=dict(arithmetic_input_class=status,n_observations=len(values),native_zero_count=sum(v==0 for v in values),
                native_min=min(values),native_max=max(values),native_unique_values=len(set(values)),
                transformed_range=span,tolerance=tolerance,max_within_group_residual=residual,
                constant_kind=('ALL_ZERO' if max(values)==0 else 'NONZERO_CONSTANT_OR_NUMERICALLY_CONSTANT')
                    if status=='NON_ESTIMABLE_CONSTANT' else 'NOT_CONSTANT')
    for g in (0,1):
        x=[v for v,k in zip(values,groups) if k==g]
        result.update({f'group{g}_'+k:v for k,v in dict(n=len(x),zero_count=sum(v==0 for v in x),
            native_min=min(x),native_max=max(x),native_mean=math.fsum(x)/len(x),
            transformed_mean=means[g]).items()})
    return result

def build(root, report, out):
    root,report,out=map(lambda p:p.resolve(),(root,report,out))
    fresh(out,[root,report]);verify(root);verify(report)
    provenance=json.loads((report/'collector_provenance.json').read_text())
    require(provenance['plan_sha256']==digest(root/'SHA256SUMS'), 'Report/plan identity differs')
    state=json.loads((report/'status.json').read_text())
    require(state['status']=='COMPLETE_PENDING_SCIENTIFIC_REVIEW','Unexpected report state')
    rows=table(report/'target_comparisons.tsv')
    contexts={c['context_id']:c for c in json.loads((root/'contexts.json').read_text())}
    require(len(rows)==state['target_rows'] and len(contexts)==state['contexts'], 'Coverage differs')
    require(len({(r['context_id'],r['target_label']) for r in rows})==len(rows),'Duplicate target/context')
    selected=[r for r in rows if r['observed_estimable']=='FALSE']
    require(selected,'No observed NE rows to audit')
    catalog=json.loads((root/'catalog.json').read_text())
    native=json.loads((root/'native_profiles.json').read_text())
    summaries=[];vectors=[]
    for r in selected:
        c=contexts[r['context_id']]; feature=r['feature']
        require(catalog['targets'][c['profiler']][feature]==r['target_label'],'Target identity differs')
        values=[];groups=[]
        require(len({o['sample_id'] for o in c['observations']})==len(c['observations']), 'Repeated person')
        for o in c['observations']:
            profile=catalog['profiles'][o['key']]['source_profile']
            value=float(native[profile].get(feature,0));group=int(o['group'])
            values.append(value);groups.append(group)
            vectors.append(dict(context_id=c['context_id'],target_label=r['target_label'],
                cohort=c['cohort'],profiler=c['profiler'],sample_id=o['sample_id'],group=group,
                dose=o['dose'],feature=feature,native_fraction=value,
                transformed_value=math.log2(1+value/1e-8),source_profile=profile))
        detail=describe(values,groups)
        require(groups.count(0)==groups.count(1)==int(c['n']), 'Group sizes differ from context')
        summaries.append(dict(**{k:r[k] for k in ('context_id','cohort','profiler','n','nominal_total_dose',
            'allocation_id','target_label','observed_status')},**detail,
            saved_status_agrees=detail['arithmetic_input_class']==r['observed_status']))
    require(all(r['saved_status_agrees'] for r in summaries),
            'Saved NE label disagrees with arithmetic input check; review before interpreting')
    # Verify caches remained unchanged throughout the audit, without requiring
    # live source profiles or changing the original sealed calculation snapshot.
    verify(root);verify(report)
    out.mkdir(parents=True)
    write_table(out/'input_summary.tsv',summaries);write_table(out/'input_vectors.tsv',vectors)
    (out/'status.json').write_text(json.dumps(dict(status='PASS_SAVED_INPUT_NE_AUDIT_PENDING_REVIEW',
        target_contexts=len(summaries),person_target_rows=len(vectors),
        statuses=dict(Counter(r['observed_status'] for r in summaries)),
        constant_kinds=dict(Counter(r['constant_kind'] for r in summaries)),
        new_fits=0,production_authorized=False),indent=2)+'\n')
    (out/'provenance.json').write_text(json.dumps(dict(plan_sha256=digest(root/'SHA256SUMS'),
        report_sha256=digest(report/'SHA256SUMS'),script_sha256=digest(Path(__file__)),
        source='checksum-verified saved native cache; absent feature maps to zero as in worker',
        check='group means and residual arithmetic; not refitting or replacing MaAsLin2'),indent=2)+'\n')
    seal(out,'PASS_SAVED_INPUT_AUDIT_NOT_PRODUCTION_AUTHORIZATION')
    print('[PASS] Saved-input audit:',len(summaries),'NE target/context rows')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','report','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();build(a.root,a.report,a.out)
