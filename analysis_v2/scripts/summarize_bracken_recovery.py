#!/usr/bin/env python3
"""Describe paired legacy/all-input recovery, keeping all samples and signed ratios."""
from __future__ import annotations
import argparse
import json
import math
from collections import defaultdict, Counter
from decimal import Decimal
from pathlib import Path

from audit_bracken_denominators import digest, require, table, write_table
from derive_bracken_all_input_recovery import verify_checksums

COHORTS = ('yachida', 'feng', 'zeller')
GROUP = ('cohort', 'analysis_population', 'condition', 'target_label', 'nominal_total_dose')
REAL_COUNTS = {'yachida': (201, 15870), 'feng': (154, 12580), 'zeller': (156, 12720)}


def quantile(values, q):
    values = sorted(values)
    position = (len(values)-1)*q
    lower = int(position)
    fraction = position-lower
    return values[lower]*(1-fraction)+values[min(lower+1,len(values)-1)]*fraction


def numeric(row, name):
    value = float(row[name])
    require(math.isfinite(value), 'Nonfinite value: '+name)
    return value


def summarize_group(key, rows):
    common = dict(zip(GROUP, key))
    base = dict(common, n_samples=len(rows),
                achieved_target_fraction_min=min(numeric(r,'implanted_fraction_target_exact') for r in rows),
                achieved_target_fraction_max=max(numeric(r,'implanted_fraction_target_exact') for r in rows),
                achieved_total_fraction_min=min(numeric(r,'total_fraction_exact') for r in rows),
                achieved_total_fraction_max=max(numeric(r,'total_fraction_exact') for r in rows),
                native_zero_samples=sum(int(r['detected_native_nonzero']) == 0 for r in rows),
                estimated_count_zero_samples=sum(int(r['estimated_count_nonzero']) == 0 for r in rows))
    summaries = []
    for method, field in [('Native fraction / old reference','legacy_recovery_ratio'),
                          ('All input pairs / read reference','recovery_ratio_all_input')]:
        ratios = [numeric(r,field) for r in rows]
        errors = [abs(x-1) for x in ratios]
        summary = dict(base, method=method,
                       recovery_ratio_median=quantile(ratios,.5),
                       recovery_ratio_q1=quantile(ratios,.25),
                       recovery_ratio_q3=quantile(ratios,.75),
                       recovery_ratio_min=min(ratios), recovery_ratio_max=max(ratios),
                       absolute_relative_error_median=quantile(errors,.5),
                       absolute_relative_error_q1=quantile(errors,.25),
                       absolute_relative_error_q3=quantile(errors,.75),
                       negative_recovery_samples=sum(x<0 for x in ratios),
                       zero_recovery_samples=sum(x==0 for x in ratios))
        summaries.append(summary)
    deltas = [numeric(r,'recovery_ratio_all_input')-numeric(r,'legacy_recovery_ratio') for r in rows]
    error_delta = [abs(numeric(r,'recovery_ratio_all_input')-1)-abs(numeric(r,'legacy_recovery_ratio')-1) for r in rows]
    comparison = dict(base, paired_ratio_delta_median=quantile(deltas,.5),
                      paired_ratio_delta_q1=quantile(deltas,.25), paired_ratio_delta_q3=quantile(deltas,.75),
                      paired_error_delta_median=quantile(error_delta,.5),
                      paired_error_delta_q1=quantile(error_delta,.25), paired_error_delta_q3=quantile(error_delta,.75),
                      samples_error_lower=sum(x<0 for x in error_delta),
                      samples_error_equal=sum(x==0 for x in error_delta),
                      samples_error_higher=sum(x>0 for x in error_delta))
    return summaries, comparison


def build(root, out, fixture=False):
    root, out = root.resolve(), out.resolve()
    require(not out.exists(), 'Output exists; choose fresh directory')
    require(out != root and root not in out.parents and out not in root.parents, 'Output overlaps recovery inputs')
    rows, evidence = [], []
    for cohort in COHORTS:
        folder = root/cohort
        verify_checksums(folder,'SHA256SUMS',['summary.json','bracken_recovery_comparison.tsv','pair_reconciliation.tsv','input_hashes.tsv'])
        summary=json.loads((folder/'summary.json').read_text())
        require(summary['cohort'] == cohort and summary['status'] == 'PASS_ENDPOINT_CONSTRUCTION', 'Wrong source stage/cohort')
        part=table(folder/'bracken_recovery_comparison.tsv')
        require(len(part) == summary['target_endpoints'], 'Endpoint count mismatch')
        if not fixture:
            require((summary['samples'],len(part)) == REAL_COUNTS[cohort], 'Unexpected final cohort size')
        keys=set()
        for r in part:
            require(r['cohort'] == cohort and r['condition'] in ('Control','Adenoma','CRC') and
                    r['analysis_population'] in ('community','independent'), 'Unexpected context')
            key=(r['sample_id'],r['analysis_population'],r['profile_id'],r['target_label'])
            require(key not in keys, 'Duplicate target endpoint')
            keys.add(key)
            require(Decimal(r['nominal_total_dose']) > 0, 'Nonpositive nominal dose')
            for field in ('legacy_recovery_ratio','recovery_ratio_all_input',
                          'implanted_fraction_target_exact','total_fraction_exact'):
                numeric(r,field)
            require(0 < float(r['implanted_fraction_target_exact']) <= float(r['total_fraction_exact']) < 1,
                    'Invalid achieved fraction')
            require(r['detected_native_nonzero'] in ('0','1') and r['estimated_count_nonzero'] in ('0','1'), 'Invalid detection flag')
        require(len({r['sample_id'] for r in part}) == summary['samples'], 'Sample coverage mismatch')
        require(dict(Counter(r['analysis_population'] for r in part)) == summary['populations'], 'Population count mismatch')
        rows.extend(part)
        for name in ('SHA256SUMS','summary.json','bracken_recovery_comparison.tsv','pair_reconciliation.tsv','input_hashes.tsv'):
            path=folder/name
            evidence.append(dict(path=str(path),sha256=digest(path)))
    groups=defaultdict(list)
    for row in rows:
        key=tuple(str(Decimal(row[f])) if f == 'nominal_total_dose' else row[f] for f in GROUP)
        groups[key].append(row)
    summaries, comparisons=[],[]
    for key,group in sorted(groups.items()):
        require(len({r['sample_id'] for r in group}) == len(group), 'Repeated biological sample in summary cell')
        ss,comparison=summarize_group(key,group)
        summaries.extend(ss)
        comparisons.append(comparison)
    for e in evidence:
        require(digest(Path(e['path'])) == e['sha256'], 'Input changed during reporting')
    out.mkdir(parents=True,exist_ok=False)
    write_table(out/'recovery_summary.tsv',summaries)
    write_table(out/'paired_method_differences.tsv',comparisons)
    write_table(out/'input_hashes.tsv',evidence)
    result=dict(status='PASS_SUMMARY_TABLES',target_endpoints=len(rows),summary_cells=len(groups),
                source_cohorts=list(COHORTS),fixture=fixture,
                uncertainty='Sample interquartile range; no confidence intervals or independent-draw claims',
                script_sha256=digest(Path(__file__)))
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'SHA256SUMS').write_text(''.join(digest(p)+'  '+p.name+'\n' for p in sorted(out.iterdir()) if p.is_file()))
    print(json.dumps(result))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--recovery-root',required=True,type=Path)
    p.add_argument('--outdir',required=True,type=Path)
    p.add_argument('--scaled-fixture',action='store_true',help='Testing only; disable real production count checks')
    a=p.parse_args()
    try:
        build(a.recovery_root,a.outdir,a.scaled_fixture)
    except (ValueError,OSError,KeyError) as error:
        p.exit(1,'[FAIL] '+str(error)+'\n')
