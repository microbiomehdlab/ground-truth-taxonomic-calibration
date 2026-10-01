#!/usr/bin/env python3
"""Freeze DA1 complete-case eligibility and baseline families; no models."""
from __future__ import annotations
import argparse
import csv
import math
from collections import Counter
from decimal import Decimal
from pathlib import Path
from audit_bracken_denominators import table, digest, require, write_table
from build_biomarker_abundance_input import parse_profile


def design_rank(rows):
    # Centred age in years improves numerical conditioning without changing rank.
    mean = sum(float(r['age']) for r in rows) / len(rows)
    a = [[1., float(r['group']), float(r['age'])-mean,
          float(r['sex']=='Male')] for r in rows]
    rank = 0
    for col in range(4):
        pivot = max(range(rank, len(a)), key=lambda i: abs(a[i][col]))
        if abs(a[pivot][col]) < 1e-10:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][col]
        a[rank] = [v/scale for v in a[rank]]
        for i in range(rank+1, len(a)):
            scale = a[i][col]
            a[i] = [v-scale*w for v,w in zip(a[i], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def eligibility(row):
    reasons = []
    try:
        age = float(row.get('age', ''))
        if not math.isfinite(age) or not 0 < age <= 120:
            raise ValueError()
    except ValueError:
        reasons.append('missing_or_invalid_age')
    if row.get('sex', '') not in ('Female', 'Male'):
        reasons.append('missing_or_invalid_sex')
    return reasons


def build(inventory, out, cohort):
    inventory, out = inventory.resolve(), out.resolve()
    require(not out.exists(), 'Output exists')
    require(out not in inventory.parents and inventory not in out.parents and out != inventory,
            'Output overlaps inventory')
    covered = set()
    for line in (inventory/'SHA256SUMS').read_text().splitlines():
        sha, name = line.split(maxsplit=1)
        require(Path(name).name == name and name not in covered, 'Unsafe inventory member')
        require(digest(inventory/name) == sha, 'Inventory checksum mismatch')
        covered.add(name)
    require({'SUCCESS','profile_inventory.tsv','input_hashes.tsv'} <= covered, 'Incomplete inventory')
    require('status\tPASS_DESIGN_INVENTORY' in (inventory/'SUCCESS').read_text().splitlines(),
            'Inventory did not pass')
    hashes = {inventory/name: digest(inventory/name) for name in covered}
    hashes[inventory/'SHA256SUMS'] = digest(inventory/'SHA256SUMS')
    for row in table(inventory/'input_hashes.tsv'):
        path = Path(row['path']).resolve()
        require(path not in hashes or hashes[path] == row['sha256'], 'Conflicting source hash')
        require(digest(path) == row['sha256'], 'Changed inventory source: '+str(path))
        hashes[path] = row['sha256']
    manifests = [p for p in hashes if p.name == 'production_manifest.tsv']
    require(len(manifests) == 1, 'Ambiguous manifest')
    manifest = table(manifests[0])
    samples = {r['sample_id']: r for r in manifest}
    require(len(samples) == len(manifest), 'Duplicate sample')
    require(all(r['condition'] in ('Control','Adenoma','CRC') for r in manifest), 'Unexpected condition')
    require(all('age' in r and 'sex' in r for r in manifest), 'Missing primary covariate columns')
    profiles = table(inventory/'profile_inventory.tsv')
    require(all(r['cohort']==cohort and r['sample_id'] in samples for r in profiles), 'Wrong cohort/sample')
    aliases = Path(__file__).resolve().parents[2]/'examples/spike_taxon_aliases.csv'
    require(aliases in hashes, 'Aliases absent from verified sources')
    with aliases.open(newline='') as h:
        alias_rows = list(csv.DictReader(h))
    ledger, summaries, families = [], [], []
    for disease in ('Adenoma','CRC'):
        selected = []
        for sid, r in sorted(samples.items()):
            if r['condition'] not in ('Control', disease):
                continue
            reasons = eligibility(r)
            rec = dict(cohort=cohort,comparison=disease+'_vs_Control',sample_id=sid,
                       condition=r['condition'],group=int(r['condition']==disease),
                       age=r['age'],sex=r['sex'],eligible=int(not reasons),
                       exclusion_reason=';'.join(reasons),bmi=r.get('bmi',''))
            ledger.append(rec)
            if not reasons:
                selected.append(rec)
        require(selected, 'No eligible samples')
        counts = Counter(r['condition'] for r in selected)
        rank = design_rank(selected)
        require(counts['Control']>0 and counts[disease]>0 and rank==4 and len(selected)>4,
                'Unidentifiable age/sex-adjusted model: '+disease)
        summaries.append(dict(cohort=cohort,comparison=disease+'_vs_Control',
                              eligible_control=counts['Control'],eligible_disease=counts[disease],
                              excluded=sum(r['comparison']==disease+'_vs_Control' and not r['eligible'] for r in ledger),
                              design_rank=rank,design_columns=4,residual_df=len(selected)-rank,
                              sex_reference='Female',group_reference='Control',models_run=0))
        ids = {r['sample_id'] for r in selected}
        for profiler in ('kraken2_bracken','metaphlan4'):
            base = [r for r in profiles if r['sample_id'] in ids and r['profiler']==profiler
                    and r['analysis_population']=='community' and Decimal(r['nominal_total_dose'])==0]
            require(len(base)==len(ids) and {r['sample_id'] for r in base}==ids, 'Invalid baseline coverage')
            positives = Counter()
            for r in base:
                path = Path(r['source_profile']).resolve()
                require(path in hashes and hashes[path]==r['sha256'], 'Baseline absent from verified sources')
                values = parse_profile(path, profiler)
                require(values and all(math.isfinite(v) and 0<=v<=1 for v in values.values()), 'Invalid abundance')
                positives.update(f for f,v in values.items() if v>0)
            targets = {r['alias'] for r in alias_rows if r['tool']==profiler}
            require(len(targets)==10, 'Expected ten target aliases')
            for feature in sorted(set(positives)|targets):
                if positives[feature]*10>=len(ids) or feature in targets:
                    families.append(dict(cohort=cohort,comparison=disease+'_vs_Control',profiler=profiler,
                                         feature=feature,baseline_positive=positives[feature],baseline_n=len(ids),
                                         target=int(feature in targets)))
    require(all(digest(p)==sha for p,sha in hashes.items()), 'Sources changed during audit')
    require(all(out!=p and out not in p.parents for p in hashes), 'Output overlaps source')
    out.mkdir(parents=True)
    for name, rows in [('eligibility.tsv',ledger),('model_design.tsv',summaries),('feature_families.tsv',families),
                       ('input_hashes.tsv',[dict(path=str(p),sha256=s) for p,s in sorted(hashes.items())])]:
        write_table(out/name,rows)
    (out/'SUCCESS').write_text('status\tPASS_DA1_CHECKPOINT\nmodels_run\t0\n')
    (out/'SHA256SUMS').write_text(''.join(digest(p)+'  '+p.name+'\n' for p in sorted(out.iterdir())))
    print('[PASS] DA1 eligibility/design/families; no models:',cohort)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inventory',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--cohort',choices=('yachida','feng','zeller'),required=True)
    a=p.parse_args()
    try:
        build(a.inventory,a.outdir,a.cohort)
    except (ValueError,OSError,KeyError) as e:
        p.exit(1,'[ERROR] '+str(e)+'\n')
