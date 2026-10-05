#!/usr/bin/env python3
"""Descriptive reference-response summaries; no new fits or population inference."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

CATEGORIES = ('BOTH_POSITIVE', 'REFERENCE_ONLY', 'OBSERVED_ONLY', 'NEITHER_POSITIVE',
              'OBSERVED_NON_ESTIMABLE', 'REFERENCE_NON_ESTIMABLE', 'BOTH_NON_ESTIMABLE')
KEY = ('cohort', 'profiler', 'n', 'nominal_total_dose', 'target_label')

def quantile(values, p):
    values = sorted(values)
    if not values: return 'NA'
    x = (len(values)-1)*p; lo = int(x); hi = min(lo+1,len(values)-1)
    return values[lo] + (values[hi]-values[lo])*(x-lo)

def summarize(rows):
    groups = defaultdict(list); identities = set()
    for row in rows:
        identity = (row['context_id'], row['target_label'])
        if identity in identities: raise ValueError('Duplicate target/context')
        identities.add(identity)
        if row['background'] != 'Adenoma' or row['exposure_arm'] != 'U':
            raise ValueError('This figure scope requires Adenoma/U')
        for field in ('category', 'hc3_category'):
            if row[field] not in CATEGORIES: raise ValueError('Unknown category')
        groups[tuple(row[k] for k in KEY)].append(row)
    frequencies=[]; effects=[]; significance=[]
    for key, cell in sorted(groups.items()):
        common=dict(zip(KEY,key)); common['allocations']=len(cell)
        for arm in ('observed','reference'):
            count=sum(r[arm+'_positive_discovery']=='TRUE' for r in cell)
            significance.append(dict(common,arm=arm,count=count,frequency=count/len(cell)))
        for method,field in (('primary_MaAsLin2','category'),('sensitivity_HC3','hc3_category')):
            counts=Counter(r[field] for r in cell)
            for category in CATEGORIES:
                frequencies.append(dict(common,inference=method,category=category,
                    count=counts[category],frequency=counts[category]/len(cell)))
        matched=[r for r in cell if r['observed_estimable']=='TRUE' and r['reference_estimable']=='TRUE']
        for quantity in ('observed_beta','reference_beta','observed_stderr','reference_stderr','beta_difference'):
            values=[]
            for r in matched:
                value=(float(r['observed_beta'])-float(r['reference_beta']) if quantity=='beta_difference' else float(r[quantity]))
                if not math.isfinite(value): raise ValueError('Non-finite estimable value')
                values.append(value)
            effects.append(dict(common,quantity=quantity,matched_estimable=len(matched),
                lower=quantile(values,.25),median=quantile(values,.5),upper=quantile(values,.75)))
    return frequencies,effects,significance

def write(path, rows):
    with path.open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def build(report,out):
    if out.exists(): raise ValueError('Fresh output required')
    for line in (report/'SHA256SUMS').read_text().splitlines():
        expected,name=line.split('  ',1); path=(report/name).resolve()
        if report.resolve() not in path.parents: raise ValueError('Unsafe checksum member')
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected: raise ValueError('Checksum mismatch: '+name)
    status=json.loads((report/'status.json').read_text())
    if status['status']!='COMPLETE_PENDING_SCIENTIFIC_REVIEW': raise ValueError('Unexpected report status')
    with (report/'target_comparisons.tsv').open(newline='') as h: rows=list(csv.DictReader(h,delimiter='\t'))
    if len(rows)!=status['target_rows'] or len({r['context_id'] for r in rows})!=status['contexts']:
        raise ValueError('Report coverage differs')
    frequencies,effects,significance=summarize(rows)
    out.mkdir(parents=True)
    write(out/'frequencies.tsv',frequencies);write(out/'effects.tsv',effects)
    write(out/'significance.tsv',significance)
    reasons=Counter(tuple(r[k] for k in KEY)+(arm,r[arm+'_status'])
        for r in rows for arm in ('observed','reference') if r[arm+'_estimable']=='FALSE')
    write(out/'nonestimable_reasons.tsv',[dict(zip(KEY+('arm','status'),key),count=count)
        for key,count in sorted(reasons.items())] or [dict(zip(KEY+('arm','status'),['NA']*7),count=0)])
    (out/'provenance.json').write_text(json.dumps(dict(source_report=str(report.resolve()),
        source_sha256=hashlib.sha256((report/'SHA256SUMS').read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        rows=len(rows),new_fits=0,production_authorized=False),indent=2)+'\n')
    print('[PASS] Descriptive summaries:',len(frequencies)//14,'cells')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--report',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();build(a.report,a.out)
