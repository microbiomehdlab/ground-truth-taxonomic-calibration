#!/usr/bin/env python3
"""Conditional cohort consistency and exact-allocation profiler agreement."""
from __future__ import annotations
import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from audit_bracken_denominators import table, require, write_table, digest
from reference_response_pilot import COHORTS, TOOLS, read_gz, fresh, seal
from prepare_da_pilot import verify


def category(a,b):
    for r in (a,b):
        require(r['estimable'] in ('TRUE','FALSE') and r['positive_discovery'] in ('TRUE','FALSE'), 'Invalid flags')
        require(r['estimable']=='TRUE' or r['positive_discovery']=='FALSE', 'Nonestimable discovery')
    if a['estimable']!='TRUE' or b['estimable']!='TRUE': return 'AT_LEAST_ONE_NON_ESTIMABLE'
    return {('TRUE','TRUE'):'BOTH',('TRUE','FALSE'):'BRACKEN_ONLY',
            ('FALSE','TRUE'):'METAPHLAN_ONLY',('FALSE','FALSE'):'NEITHER'}[a['positive_discovery'],b['positive_discovery']]


def build(report,out):
    fresh(out,[report]); verify(report)
    status=json.loads((report/'status.json').read_text())
    require(status.get('completed_contexts')==120000 and not status.get('failures')
            and status.get('mismatched_features')==0, 'Incomplete/unvalidated direct MaAsLin2 report')
    contexts={r['context_id']:r for r in table(report/'context_summary.tsv')}
    require(len(contexts)==120000, 'Duplicate/missing contexts')
    groups=defaultdict(list); pairs=defaultdict(dict); seen=set()
    for r in read_gz(report/'targets.tsv.gz'):
        cid=r['context_id']; require(cid in contexts,'Unknown target context')
        require((cid,r['target_label']) not in seen,'Duplicate target'); seen.add((cid,r['target_label']))
        c=contexts[cid]; require(c['cohort'] in COHORTS and c['profiler'] in TOOLS,'Unknown cohort/tool')
        require(r['estimable'] in ('TRUE','FALSE') and r['positive_discovery'] in ('TRUE','FALSE'),'Invalid flags')
        require(r['estimable']=='TRUE' or r['positive_discovery']=='FALSE','Nonestimable discovery')
        setting=(c['background'],c['n'],c['arm'],c['anchor'],r['target_label'])
        groups[setting+(c['profiler'],c['cohort'])].append((r,int(c['family_n'])))
        key=setting+(c['cohort'],c['allocation_id'])
        require(c['profiler'] not in pairs[key], 'Duplicate profiler allocation')
        pairs[key][c['profiler']]=r
    counts=Counter(cid for cid,label in seen)
    require(len(seen)==1200000 and set(counts)==set(contexts) and set(counts.values())=={10},'Target coverage differs')
    summaries=[]; bysetting=defaultdict(dict)
    fields=('background','n','exposure_arm','nominal_total_dose','target_label','profiler')
    for key,values in sorted(groups.items()):
        rows=[r for r,family in values]; n=len(rows)
        require(n==100,'Expected complete 100-allocation cells')
        e=sum(r['estimable']=='TRUE' for r in rows); d=sum(r['positive_discovery']=='TRUE' for r in rows)
        record=dict(zip(fields+('cohort',),key),allocations=n,estimable=e,positive_discoveries=d,
            non_estimable=n-e,positive_frequency=d/n,estimable_frequency=e/n,
            family_n_min=min(f for r,f in values),family_n_max=max(f for r,f in values),
            median_estimable_beta=median(float(r['beta']) for r in rows if r['estimable']=='TRUE') if e else 'NA')
        summaries.append(record); bysetting[key[:-1]][key[-1]]=record
    across=[]
    for key,cohorts in sorted(bysetting.items()):
        require(set(cohorts)==set(COHORTS),'Missing cohort setting')
        for discovery,replication in (('yachida','feng'),('yachida','zeller'),('feng','zeller')):
            a,b=cohorts[discovery],cohorts[replication]
            across.append(dict(zip(fields,key),discovery_cohort=discovery,replication_cohort=replication,
                discovery_frequency=a['positive_frequency'],replication_frequency=b['positive_frequency'],
                independent_allocation_joint_frequency=a['positive_frequency']*b['positive_frequency'],
                replication_given_discovery_frequency=b['positive_frequency'] if a['positive_discoveries'] else 'NA',
                both_estimable_frequency=a['estimable_frequency']*b['estimable_frequency'],
                all_three_positive_frequency=cohorts['yachida']['positive_frequency']*
                    cohorts['feng']['positive_frequency']*cohorts['zeller']['positive_frequency'],
                interpretation='PRODUCT_OF_EMPIRICAL_CONDITIONAL_FREQUENCIES_INDEPENDENT_COHORT_DRAWS_NOT_POPULATION_POWER'))
    paired=defaultdict(Counter)
    for key,tools in pairs.items():
        require(set(tools)==set(TOOLS),'Missing matched profiler')
        paired[key[:-1]][category(tools[TOOLS[0]],tools[TOOLS[1]])]+=1
    agreement=[]
    for key,counts in sorted(paired.items()):
        n=sum(counts.values()); require(n==100,'Incomplete paired profiler cell')
        for name in ('BOTH','BRACKEN_ONLY','METAPHLAN_ONLY','NEITHER','AT_LEAST_ONE_NON_ESTIMABLE'):
            agreement.append(dict(zip(fields[:-1]+('cohort',),key),category=name,allocations=n,
                count=counts[name],conditional_frequency=counts[name]/n))
    out.mkdir(parents=True)
    for name,rows in (('cohort_frequencies.tsv',summaries),('cross_cohort_consistency.tsv',across),
                      ('matched_profiler_agreement.tsv',agreement)):
        write_table(out/name,rows)
    write_table(out/'input_hashes.tsv',[dict(path=str(report/n),sha256=digest(report/n))
        for n in ('SHA256SUMS','status.json','context_summary.tsv','targets.tsv.gz')])
    verify(report)
    (out/'status.json').write_text(json.dumps(dict(status='PASS_CONDITIONAL_CONSISTENCY_PENDING_REVIEW',
        new_model_fits=0,target_rows=len(seen),production_authorized=False),indent=2)+'\n')
    seal(out,'PASS_CONDITIONAL_CONSISTENCY_NOT_PRODUCTION_AUTHORIZATION')
    print('[PASS] Conditional consistency:',len(summaries),'cells')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--report',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();build(a.report.resolve(),a.out.resolve())
