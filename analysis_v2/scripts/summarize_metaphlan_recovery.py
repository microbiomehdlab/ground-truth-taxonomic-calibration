#!/usr/bin/env python3
"""Descriptive paired reference comparison; raw MetaPhlAn abundances unchanged."""
from __future__ import annotations
import argparse
import json
import math
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from summarize_bracken_recovery import COHORTS, GROUP, REAL_COUNTS, summarize_group
from audit_bracken_denominators import digest, require, table, write_table
from derive_bracken_all_input_recovery import verify_checksums

LABELS = ("Read-fraction reference / sensitivity", "Genome-equivalent reference / primary")


def convert(row):
    require(row['reference_type'] == 'genome_equivalent', 'Wrong primary reference')
    nominal = Decimal(row['profile_id'].rsplit('_f', 1)[1].replace('p', '.'))
    require(nominal > 0, 'Invalid nominal dose')
    old, new = float(row['response_ratio']), float(row['response_ratio_profiler_scale'])
    require(math.isfinite(old) and math.isfinite(new), 'Nonfinite recovery')
    require(row['observed_detected_native_nonzero'] in ('0','1'), 'Invalid detection flag')
    return dict(row, nominal_total_dose=str(nominal),
                legacy_recovery_ratio=old, recovery_ratio_all_input=new,
                implanted_fraction_target_exact=row['spike_fraction_target'],
                total_fraction_exact=row['spike_fraction_total'],
                detected_native_nonzero=row['observed_detected_native_nonzero'],
                estimated_count_nonzero=row['observed_detected_native_nonzero'])


def build(root, out, fixture=False):
    root, out = root.resolve(), out.resolve()
    require(not out.exists(), 'Output exists')
    require(out != root and root not in out.parents and out not in root.parents, 'Output overlaps input')
    required = ['SUCCESS'] + [c+'/endpoints/paired_endpoints.tsv' for c in COHORTS]
    # Root manifest uses ./relative paths; verify the entire bundle with safe paths.
    covered = set()
    for line in (root/'SHA256SUMS').read_text().splitlines():
        checksum, name = line.split(maxsplit=1)
        name = name.lstrip('*')
        path = (root/name).resolve()
        require(root in path.parents and path not in covered, 'Unsafe/duplicate checksum path')
        require(digest(path) == checksum, 'Checksum mismatch: '+name)
        covered.add(path)
    require(all((root/n).resolve() in covered for n in required), 'Uncovered report inputs')
    success = dict(line.split('\t',1) for line in (root/'SUCCESS').read_text().splitlines())
    require(success.get('status') == 'PASS_ENDPOINTS', 'Incomplete endpoint bundle')
    evidence = [dict(path=str(p),sha256=digest(p)) for p in sorted(covered)]
    groups = defaultdict(list)
    count = 0
    for cohort in COHORTS:
        rows = [r for r in table(root/cohort/'endpoints/paired_endpoints.tsv') if r['profiler']=='metaphlan4']
        if not fixture:
            require((len({r['sample_id'] for r in rows}),len(rows)) == REAL_COUNTS[cohort], 'Cohort count mismatch')
        keys=set()
        for row in rows:
            require(row['cohort']==cohort and row['condition'] in ('Control','Adenoma','CRC') and row['analysis_population'] in ('community','independent'), 'Wrong context')
            key=(row['sample_id'],row['analysis_population'],row['profile_id'],row['target_label'])
            require(key not in keys, 'Duplicate endpoint')
            keys.add(key)
            converted=convert(row)
            groups[tuple(converted[f] for f in GROUP)].append(converted)
        count+=len(rows)
    summaries, differences=[],[]
    for key, rows in sorted(groups.items()):
        require(len({r['sample_id'] for r in rows})==len(rows), 'Repeated sample in cell')
        ss, dd=summarize_group(key,rows)
        for i,s in enumerate(ss):
            s['method']=LABELS[i]
            del s['estimated_count_zero_samples']
        del dd['estimated_count_zero_samples']
        summaries.extend(ss)
        differences.append(dd)
    require(count>0, 'No endpoints')
    require(all(digest(Path(e['path']))==e['sha256'] for e in evidence), 'Inputs changed')
    out.mkdir(parents=True)
    write_table(out/'recovery_summary.tsv',summaries)
    write_table(out/'paired_method_differences.tsv',differences)
    write_table(out/'input_hashes.tsv',evidence)
    result=dict(status='PASS_SUMMARY_TABLES',profiler='metaphlan4',fixture=fixture,
                target_endpoints=count,summary_cells=len(groups),
                uncertainty='Sample IQR, not confidence intervals',
                reference_change='Expected reference only; observed abundances unchanged')
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'SHA256SUMS').write_text(''.join(digest(p)+'  '+p.name+'\n' for p in sorted(out.iterdir())))
    print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--recovery-root',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--scaled-fixture',action='store_true')
    a=p.parse_args()
    try: build(a.recovery_root,a.outdir,a.scaled_fixture)
    except (ValueError,OSError,KeyError,IndexError) as error: p.exit(1,'[ERROR] '+str(error)+'\n')
