#!/usr/bin/env python3
"""Matched DA3 reference extension. Never alters profiles or authorizes production."""
from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import itertools
import json
import math
import os
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from decimal import Decimal
from pathlib import Path

from audit_bracken_denominators import digest, require, table, write_table
from build_biomarker_abundance_input import parse_profile
from prepare_da_pilot import verify

COHORTS = ('yachida', 'feng', 'zeller')
TOOLS = ('kraken2_bracken', 'metaphlan4')
DOSES = ('0.0001', '0.001', '0.01')  # total community fraction, NOT per target
CODE = ('scripts/reference_response_pilot.py', 'scripts/fit_reference_response.R',
        'scripts/audit_bracken_denominators.py', 'scripts/build_biomarker_abundance_input.py',
        'scripts/prepare_da_pilot.py', 'lib/maaslin_contract.R', 'lib/maaslin_context.R',
        'lib/clinical_hc3.R')


def fresh(out, inputs):
    out = out.resolve()
    require(not out.exists(), 'Fresh output required: ' + str(out))
    for p in inputs:
        p = p.resolve()
        require(out != p and out not in p.parents and p not in out.parents, 'Overlapping output')


def seal(folder, status):
    (folder / 'SUCCESS').write_text('status\t' + status + '\nproduction_authorized\t0\n')
    files = sorted(p for p in folder.rglob('*') if p.is_file() and p.name != 'SHA256SUMS')
    (folder / 'SHA256SUMS').write_text(''.join(
        digest(p) + '  ' + p.relative_to(folder).as_posix() + '\n' for p in files))


def read_gz(path):
    with gzip.open(path, 'rt', newline='') as h:
        yield from csv.DictReader(h, delimiter='\t')


def expected_profile(baseline, tool, inserted, total_fraction, targets,
                     species_total=None, geff=None, genome_sizes=None):
    """Whole-profile, baseline-anchored expectation; never mix abundance scales."""
    require(baseline and all(math.isfinite(v) and 0 <= v <= 1.00001
                            for v in baseline.values()), 'Invalid baseline fractions')
    require(set(inserted) == set(targets.values()) and len(inserted) == 10,
            'Exact ten-member community required')
    require(all(isinstance(v, int) and v > 0 for v in inserted.values()), 'Invalid inserted counts')
    require(math.isfinite(total_fraction) and 0 < total_fraction < 1, 'Invalid achieved total')
    if tool == 'kraken2_bracken':
        require(species_total is not None and species_total > 0, 'Invalid estimated species total')
        denominator = species_total + sum(inserted.values())
        scale = species_total / denominator
        additions = {f: inserted[label] / denominator for f, label in targets.items()}
    else:
        require(tool == 'metaphlan4' and geff is not None and math.isfinite(geff) and geff > 0,
                'Invalid effective genome size')
        require(genome_sizes is not None and set(inserted) <= set(genome_sizes), 'Missing genome sizes')
        require(all(math.isfinite(genome_sizes[k]) and genome_sizes[k] > 0 for k in inserted),
                'Invalid target genome size')
        total = sum(inserted.values())
        q = {k: total_fraction * inserted[k] / total * geff / genome_sizes[k] for k in inserted}
        denominator = 1 - total_fraction + sum(q.values())
        scale = (1 - total_fraction) / denominator
        additions = {f: q[label] / denominator for f, label in targets.items()}
    result = {f: scale * v for f, v in baseline.items()}
    for f, addition in additions.items():
        result[f] = result.get(f, 0) + addition
    require(all(math.isfinite(v) and 0 <= v <= 1.00001 for v in result.values()),
            'Invalid reference profile')
    return result, scale


def reconcile_fractions(rows, design):
    """Honor the historical eight-decimal total, but use integer-count truth."""
    original, added = int(design['R']), int(design['N_total'])
    require(original > 0 and added > 0, 'Invalid design counts')
    inserted = {r['target_label']:int(r['implanted_read_pairs_target']) for r in rows}
    require(len(inserted) == len(rows) == 10 and all(n > 0 for n in inserted.values())
            and sum(inserted.values()) == added, 'Community inserted counts do not reconcile')
    exact = added/(original+added)
    recorded = float(design['f_hat'])
    require(math.isfinite(recorded) and abs(recorded-exact) <= 5.01e-9,
            f'Design f_hat differs beyond eight-decimal rounding: recorded={recorded:.17g}, '
            f'exact={exact:.17g}, R={original}, N_total={added}')
    for row in rows:
        label = row['target_label']
        total, target = float(row['spike_fraction_total']), float(row['spike_fraction_target'])
        require(math.isfinite(total) and abs(total-recorded) <= 1e-12,
                f'Canonical total differs from design f_hat: {label}, canonical={total:.17g}, '
                f'design={recorded:.17g}')
        require(math.isfinite(target) and abs(target-inserted[label]/(original+added)) <= 1e-12,
                'Canonical target differs from exact inserted count fraction: '+label)
    return original, added, inserted, exact


def selected_contexts(plan, allocations=20, arms=('U',)):
    require(allocations in (20, 100), 'Use 20 engineering or 100 extension allocations')
    require(arms and len(set(arms)) == len(arms) and set(arms) <= {'U', 'P25', 'P50', 'P75'},
            'Unsupported exposure arms')
    verify(plan)
    require('status\tPASS_DRAFT_BATCH_PLAN_NOT_PRODUCTION_AUTHORIZATION' in
            (plan / 'SUCCESS').read_text().splitlines(), 'Wrong source plan')
    candidates = []
    for t in table(plan / 'tasks.tsv'):
        for c in json.loads((plan / t['file']).read_text()):
            if (c['background'] == 'Adenoma' and c['n'] in (10, 20)
                    and c['arm'] in arms and c['anchor'] in DOSES):
                candidates.append(c)
    ids = sorted({c['allocation_id'] for c in candidates})[:allocations]
    require(len(ids) == allocations, 'Insufficient saved allocations')
    result = [c for c in candidates if c['allocation_id'] in ids]
    cells = Counter((c['cohort'], c['profiler'], c['n'], c['arm'], c['anchor']) for c in result)
    expected = set(itertools.product(COHORTS, TOOLS, (10, 20), arms, DOSES))
    require(set(cells) == expected and set(cells.values()) == {allocations}, 'Incomplete pilot grid')
    require(len({c['context_id'] for c in result}) == len(result), 'Duplicate context')
    for c in result:
        obs = c['observations']
        require(len(obs) == 2*c['n'] and len({r['sample_id'] for r in obs}) == len(obs)
                and Counter(r['group'] for r in obs) == {0:c['n'], 1:c['n']}, 'Invalid allocation')
        require(all(r['group'] == 1 or Decimal(r['dose']) == 0 for r in obs), 'Spiked comparison group')
    return result


def prepare(a):
    out, repo = a.out.resolve(), a.repo.resolve()
    fresh(out, [a.plan, a.inventory, a.geff_root, a.observed_report, Path(a.image)])
    require(out != repo and out not in repo.parents, 'Output contains repository')
    contexts = selected_contexts(a.plan, a.allocations, tuple(a.arms.split(',')))
    catalog = json.loads((a.plan / 'catalog.json').read_text())
    verify(a.observed_report)
    observed_status = json.loads((a.observed_report / 'status.json').read_text())
    require(not observed_status.get('failures') and observed_status.get('mismatched_features') == 0
            and observed_status.get('completed_contexts') == 120000, 'Observed DA3 incomplete/unvalidated')
    observed_identity = json.loads((a.observed_report/'identity.json').read_text())
    require(observed_identity['plan_sha256'] == digest(a.plan/'SHA256SUMS'),
            'Observed results used a different allocation plan')
    require(observed_identity['image_sha256'] == digest(Path(a.image)), 'Observed/reference image differs')
    for name in ('lib/maaslin_contract.R','lib/maaslin_context.R','scripts/build_biomarker_abundance_input.py'):
        require(observed_identity['code'][name] == digest(repo/'analysis_v2'/name),
                'Observed/reference model or parser differs: '+name)
    verify(a.geff_root)
    require('status\tPASS' in (a.geff_root/'SUCCESS').read_text().splitlines(), 'Wrong G_eff bundle')
    hashes = {}
    def track(p):
        p = p.resolve(); h = digest(p)
        require(str(p) not in hashes or hashes[str(p)] == h, 'Source changed during preparation')
        hashes[str(p)] = h
    for name in CODE:
        track(repo/'analysis_v2'/name)
    for p in (a.plan/'SHA256SUMS', a.observed_report/'SHA256SUMS', a.geff_root/'SHA256SUMS'):
        track(p)
    target_path = a.geff_root/'targets/target_genome_sizes.tsv'
    track(target_path)
    sizes = {r['target_label']:float(r['genome_size_bp']) for r in table(target_path)}
    require(len(sizes) == 10, 'Exact ten genome sizes required')
    geff, canonical = {}, {}
    for cohort in COHORTS:
        inv = a.inventory/cohort
        verify(inv); track(inv/'SHA256SUMS')
        require('status\tPASS_DESIGN_INVENTORY' in (inv/'SUCCESS').read_text().splitlines(), 'Wrong inventory')
        paths = [r for r in table(inv/'input_hashes.tsv') if Path(r['path']).name == 'canonical_input.tsv']
        require(len(paths) == 1, 'Ambiguous canonical source: '+cohort)
        path = Path(paths[0]['path']); require(digest(path) == paths[0]['sha256'], 'Canonical changed'); track(path)
        for r in table(path):
            if r['analysis_population'] == 'community' and float(r['spike_fraction_total']) > 0:
                canonical.setdefault((cohort,r['profiler'],r['sample_id'],r['source_profile']), []).append(r)
        path = a.geff_root/cohort/'geff_primary_cov95/effective_genome_size.tsv'; track(path)
        for r in table(path):
            key = (cohort,r['sample_id']); require(key not in geff, 'Duplicate G_eff')
            require(float(r['mapping_coverage']) >= .95, 'G_eff coverage below primary policy')
            geff[key] = float(r['effective_genome_size_bp'])
    wanted = {c['context_id']:c for c in contexts}
    for p in (a.observed_report/'targets.tsv.gz',a.observed_report/'context_summary.tsv',a.plan/'catalog.json'):
        track(p)
    summary = {r['context_id']:r for r in table(a.observed_report/'context_summary.tsv') if r['context_id'] in wanted}
    require(set(summary) == set(wanted), 'Observed context coverage differs')
    for cid,c in wanted.items():
        require(all(summary[cid][k] == str(c[k]) for k in
                    ('cohort','background','profiler','n','arm','anchor','allocation_id')), 'Observed metadata differs')
        require(int(summary[cid]['family_n']) == len(catalog['families']['|'.join((c['cohort'],c['profiler'],c['background']))]),
                'Observed frozen family differs')
    observed = [r for r in read_gz(a.observed_report/'targets.tsv.gz') if r['context_id'] in wanted]
    require(len(observed) == 10*len(contexts), 'Observed target count differs')
    require(len({(r['context_id'],r['target_label']) for r in observed}) == len(observed), 'Duplicate observed target')
    require({(r['context_id'],r['target_label']) for r in observed} ==
            {(c['context_id'],label) for c in contexts for label in catalog['targets'][c['profiler']].values()},
            'Observed target panel differs')
    profiles, references, audits = {}, {}, []
    for c in contexts:
        tool, cohort = c['profiler'], c['cohort']; targets = catalog['targets'][tool]
        require(len(targets) == 10 and set(targets.values()) == set(sizes), 'Invalid target aliases')
        for obs in c['observations']:
            key = obs['key']
            require(key == '|'.join((cohort,tool,obs['sample_id'],str(Decimal(obs['dose'])))), 'Observation identity differs')
            if key in references: continue
            p = catalog['profiles'][key]
            require(p['sample_id'] == obs['sample_id'] and p['cohort'] == cohort and
                    p['profiler'] == tool and p['condition'] == 'Adenoma' and
                    p['analysis_population'] == 'community' and
                    Decimal(p['nominal_total_dose']) == Decimal(obs['dose']), 'Catalog metadata differs')
            baseline_key = '|'.join((cohort,tool,obs['sample_id'],'0'))
            bp = catalog['profiles'][baseline_key]; baseline_path = Path(bp['source_profile'])
            require(all(bp[k] == p[k] for k in ('cohort','profiler','sample_id','condition','analysis_population'))
                    and Decimal(bp['nominal_total_dose']) == 0, 'Baseline identity differs')
            for entry in (p,bp):
                path = Path(entry['source_profile'])
                require(digest(path) == entry['sha256'], 'Profile changed: '+str(path)); track(path)
                if entry['source_profile'] not in profiles:
                    profiles[entry['source_profile']] = parse_profile(path, tool)
            b = profiles[str(baseline_path)]
            if Decimal(obs['dose']) == 0:
                references[key] = b; continue
            rows = canonical.get((cohort,tool,obs['sample_id'],p['source_profile']), [])
            require(len(rows) == 10 and {r['target_label'] for r in rows} == set(sizes), 'Incomplete canonical community')
            require(len({r['source_design'] for r in rows}) == 1, 'Conflicting community designs')
            design = Path(rows[0]['source_design']); track(design)
            dr = [r for r in table(design) if Decimal(r['fraction']) == Decimal(obs['dose'])]
            require(len(dr) == 1 and dr[0]['sample_id'] == obs['sample_id'], 'Ambiguous design dose')
            try:
                original, added, inserted, F = reconcile_fractions(rows,dr[0])
            except ValueError as error:
                raise ValueError(f'{key}; design={design}: {error}') from error
            S = None
            if tool == 'kraken2_bracken':
                S = sum(int(r['new_est_reads']) for r in table(baseline_path) if r['taxonomy_lvl'] == 'S')
            e, scale = expected_profile(b, tool, inserted, F, targets, S, geff.get((cohort,obs['sample_id'])), sizes)
            references[key] = e
            audits.append(dict(key=key,cohort=cohort,sample_id=obs['sample_id'],profiler=tool,
                nominal_total_dose=obs['dose'],original_pairs=original,inserted_pairs=added,
                achieved_total_fraction=F,baseline_species_counts=S or '',retained_scale=scale,
                recorded_design_total_fraction=dr[0]['f_hat'],
                exact_minus_recorded_total_fraction=F-float(dr[0]['f_hat']),
                inserted_counts_json=json.dumps(inserted,sort_keys=True),
                effective_genome_size_bp=geff[(cohort,obs['sample_id'])] if tool == 'metaphlan4' else '',
                target_genome_sizes_json=json.dumps(sizes,sort_keys=True) if tool == 'metaphlan4' else '',
                baseline_native_sum=sum(b.values()),reference_sum=sum(e.values()),
                baseline_anchored_reference='1',not_biological_truth='1'))
    for p in (a.observed_report/'targets.tsv.gz',a.observed_report/'context_summary.tsv',a.plan/'catalog.json'):
        track(p)
    require(all(digest(Path(p)) == h for p,h in hashes.items()), 'Source changed during preparation')
    out.mkdir(parents=True); (out/'source/analysis_v2').mkdir(parents=True)
    for name in CODE:
        src = repo/'analysis_v2'/name; dst = out/'source/analysis_v2'/name
        dst.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(src,dst)
        require(digest(dst) == hashes[str(src.resolve())], 'Calculation code changed during preparation')
    (out/'catalog.json').write_text(json.dumps(catalog))
    (out/'reference_profiles.json').write_text(json.dumps(references,allow_nan=False))
    (out/'native_profiles.json').write_text(json.dumps(profiles,allow_nan=False))
    (out/'contexts.json').write_text(json.dumps(contexts))
    (out/'batches').mkdir()
    tasks=[]
    for i in range(0,len(contexts),a.batch_size):
        name='batches/batch_%04d.json'%len(tasks)
        (out/name).write_text(json.dumps(contexts[i:i+a.batch_size]))
        tasks.append(dict(index=len(tasks),file=name,contexts=len(contexts[i:i+a.batch_size])))
    write_table(out/'tasks.tsv',tasks); write_table(out/'observed_targets.tsv',observed)
    write_table(out/'reference_audit.tsv',audits)
    write_table(out/'input_hashes.tsv',[dict(path=p,sha256=h) for p,h in sorted(hashes.items())])
    (out/'identity.json').write_text(json.dumps(dict(contract='REFERENCE_RESPONSE_PILOT_V1',
        contexts=len(contexts),allocations=a.allocations,arms=a.arms,image_path=str(Path(a.image).resolve()),
        image_sha256=digest(Path(a.image)),new_profiling=0,production_authorized=False,
        baseline_anchor='printed_native_abundance',bracken_denominator='estimated_species_counts_plus_inserted',
        interpretation='CONDITIONAL_REFERENCE_COMPARISON_NOT_CAUSAL_ATTRIBUTION'),indent=2)+'\n')
    seal(out,'PASS_REFERENCE_PLAN_NOT_PRODUCTION_AUTHORIZATION')
    print('[PASS] Reference contexts:',len(contexts),'batches:',len(tasks))


def run(root, index):
    verify(root)
    state=json.loads((root/'identity.json').read_text())
    require(digest(Path(state['image_path'])) == state['image_sha256'], 'Analysis image changed')
    tasks=table(root/'tasks.tsv'); require(0 <= index < len(tasks), 'Invalid task index')
    task=tasks[index]; contexts=json.loads((root/task['file']).read_text())
    require(len(contexts) == int(task['contexts']), 'Task count differs')
    output=root/'results'; output.mkdir(exist_ok=True)
    locks=root/'locks'; locks.mkdir(exist_ok=True)
    name='batch_%04d'%index; final=output/name
    with (locks/(name+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        identity=dict(plan_sha256=digest(root/'SHA256SUMS'),image_sha256=state['image_sha256'])
        if final.exists():
            verify(final); require(json.loads((final/'identity.json').read_text()) == identity, 'Resume identity changed')
            print('[PASS] Already completed:',name); return
        catalog=json.loads((root/'catalog.json').read_text())
        references=json.loads((root/'reference_profiles.json').read_text())
        natives=json.loads((root/'native_profiles.json').read_text())
        attempts=root/'attempts'; attempts.mkdir(exist_ok=True)
        attempt=Path(tempfile.mkdtemp(prefix=name+'_',dir=str(attempts)))
        inp=attempt/'input'; inp.mkdir(); meta=[]
        for c in contexts:
            folder=inp/c['context_id']; folder.mkdir()
            features=catalog['families']['|'.join((c['cohort'],c['profiler'],c['background']))]
            require(len(features) == len(set(features)), 'Duplicate frozen features')
            targets=catalog['targets'][c['profiler']]; require(set(targets) <= set(features), 'Missing targets in family')
            for arm,values in (('reference',references),('observed',natives)):
                rows=[]
                for r in c['observations']:
                    key=r['key'] if arm == 'reference' else catalog['profiles'][r['key']]['source_profile']
                    rows.append(dict(observation_id=r['sample_id'],**{f:values[key].get(f,0) for f in features}))
                write_table(folder/(arm+'.tsv'),rows)
            write_table(folder/'metadata.tsv',[dict(observation_id=r['sample_id'],group=r['group']) for r in c['observations']])
            write_table(folder/'targets.tsv',[dict(feature=f,target_label=l) for f,l in targets.items()])
            meta.append({k:v for k,v in c.items() if k != 'observations'})
        write_table(inp/'contexts.tsv',meta)
        subprocess.run(['Rscript',str(root/'source/analysis_v2/scripts/fit_reference_response.R'),
                        str(inp),str(attempt/'fit'),str(root/'source')],check=True)
        staged=attempt/'fit'
        rows=table(staged/'targets.tsv')
        require(len(rows) == len(contexts)*10*2, 'Missing target results')
        validate_features(staged,contexts,catalog)
        (staged/'identity.json').write_text(json.dumps(identity))
        seal(staged,'PASS_REFERENCE_BATCH_NOT_PRODUCTION_AUTHORIZATION')
        os.rename(staged,final)
        # Only our newly created successful scratch attempt; failures are retained.
        require(attempt.parent == attempts and attempt.name.startswith(name+'_'), 'Unsafe scratch cleanup')
        shutil.rmtree(attempt)
        print('[PASS]',name)


def classify(observed, reference):
    require(observed['estimable'] in ('TRUE','FALSE') and reference['estimable'] in ('TRUE','FALSE'), 'Invalid estimability')
    if observed['estimable'] != 'TRUE' or reference['estimable'] != 'TRUE':
        return ('BOTH_NON_ESTIMABLE' if observed['estimable'] == reference['estimable']
                else 'OBSERVED_NON_ESTIMABLE' if observed['estimable'] == 'FALSE' else 'REFERENCE_NON_ESTIMABLE')
    for r in (observed,reference):
        require(r['positive_discovery'] in ('TRUE','FALSE'), 'Invalid discovery flag')
    return {('TRUE','TRUE'):'BOTH_POSITIVE',('FALSE','TRUE'):'REFERENCE_ONLY',
            ('TRUE','FALSE'):'OBSERVED_ONLY',('FALSE','FALSE'):'NEITHER_POSITIVE'}[
                observed['positive_discovery'],reference['positive_discovery']]


def validate_features(folder, contexts, catalog):
    rows=iter(read_gz(folder/'features.tsv.gz'))
    for c in contexts:
        family=catalog['families']['|'.join((c['cohort'],c['profiler'],c['background']))]
        for arm in ('reference','observed'):
            for f in family:
                r=next(rows,None)
                require(r is not None and (r['context_id'],r['arm'],r['feature']) ==
                        (c['context_id'],arm,f), 'Full-family identity/order/coverage differs')
                require(r['estimable'] in ('TRUE','FALSE'), 'Invalid estimability')
                require(0 <= float(r['wrapper_q']) <= 1, 'Invalid full-family q')
    require(next(rows,None) is None, 'Extra full-family result')


def require_observed_equivalence(saved,current):
    require(saved['estimable'] == current['estimable'], 'Observed estimability changed')
    for field in ('beta','stderr','raw_p','wrapper_q'):
        a,b=saved[field],current[field]
        if a == b == 'NA': continue
        require(a != 'NA' and b != 'NA', 'Observed inference missingness differs')
        a,b=float(a),float(b)
        require(math.isfinite(a) and math.isfinite(b) and abs(a-b) <= 1e-8+1e-6*max(abs(a),abs(b)),
                'Observed inference no longer matches saved MaAsLin2: '+field)


def hc3_flags(row):
    estimable = row['hc3_estimable'] == 'TRUE'
    positive = estimable and float(row['beta']) > 0 and float(row['hc3_q']) <= .05
    return dict(estimable='TRUE' if estimable else 'FALSE',
                positive_discovery='TRUE' if positive else 'FALSE')


def collect(root, out):
    fresh(out,[root]); verify(root)
    identity=dict(plan_sha256=digest(root/'SHA256SUMS'),image_sha256=json.loads((root/'identity.json').read_text())['image_sha256'])
    observed={(r['context_id'],r['target_label']):r for r in table(root/'observed_targets.tsv')}
    contexts={c['context_id']:c for c in json.loads((root/'contexts.json').read_text())}
    catalog=json.loads((root/'catalog.json').read_text())
    comparisons=[]; seen=set()
    for task in table(root/'tasks.tsv'):
        folder=root/'results'/('batch_%04d'%int(task['index'])); verify(folder)
        require(json.loads((folder/'identity.json').read_text()) == identity, 'Batch identity differs')
        expected=json.loads((root/task['file']).read_text())
        validate_features(folder,expected,catalog)
        rows=table(folder/'targets.tsv'); index={(r['context_id'],r['target_label'],r['arm']):r for r in rows}
        require(len(index) == len(rows) == len(expected)*20, 'Target result coverage differs')
        for c in expected:
            require(c['context_id'] not in seen, 'Context completed twice'); seen.add(c['context_id'])
            for label in sorted({k[1] for k in index if k[0] == c['context_id']}):
                key=(c['context_id'],label); o=observed[key]
                r=index[key+('reference',)]; h=index[key+('observed',)]
                require(o['feature'] == r['feature'] == h['feature'], 'Target feature identity differs')
                require_observed_equivalence(o,h)
                comparisons.append(dict(context_id=c['context_id'],cohort=c['cohort'],profiler=c['profiler'],
                    background=c['background'],n=c['n'],exposure_arm=c['arm'],nominal_total_dose=c['anchor'],
                    allocation_id=c['allocation_id'],target_label=label,category=classify(o,r),
                    hc3_category=classify(hc3_flags(h),hc3_flags(r)),
                    **{'observed_'+k:o[k] for k in ('beta','stderr','raw_p','wrapper_q','status','estimable','positive_discovery')},
                    **{'reference_'+k:r[k] for k in ('beta','stderr','raw_p','wrapper_q','status','estimable','positive_discovery','hc3_p','hc3_q','hc3_stderr','hc3_estimable')},
                    observed_hc3_p=h['hc3_p'],observed_hc3_q=h['hc3_q'],
                    observed_hc3_stderr=h['hc3_stderr'],observed_hc3_estimable=h['hc3_estimable']))
    require(seen == set(contexts) and len(comparisons) == len(contexts)*10, 'Incomplete collection')
    groups=defaultdict(list)
    for r in comparisons:
        groups[(r['cohort'],r['profiler'],r['n'],r['exposure_arm'],r['nominal_total_dose'],r['target_label'])].append(r)
    summary=[]
    for key,rows in sorted(groups.items()):
        for inference,column in (('primary_MaAsLin2','category'),('sensitivity_HC3','hc3_category')):
            for category,count in sorted(Counter(r[column] for r in rows).items()):
                summary.append(dict(zip(('cohort','profiler','n','exposure_arm','nominal_total_dose','target_label'),key),
                    inference=inference,category=category,allocations=len(rows),count=count,conditional_frequency=count/len(rows)))
    out.mkdir(parents=True); write_table(out/'target_comparisons.tsv',comparisons)
    write_table(out/'conditional_categories.tsv',summary)
    (out/'status.json').write_text(json.dumps(dict(status='COMPLETE_PENDING_SCIENTIFIC_REVIEW',
        contexts=len(contexts),target_rows=len(comparisons),production_authorized=False,
        inference='MaAsLin2_1.18.0_primary_HC3_sensitivity_not_randomization_proof'),indent=2)+'\n')
    seal(out,'PASS_REFERENCE_COMPARISON_NOT_PRODUCTION_AUTHORIZATION')


def main():
    p=argparse.ArgumentParser(description=__doc__); s=p.add_subparsers(dest='command',required=True)
    q=s.add_parser('prepare')
    for name in ('plan','inventory','geff-root','observed-report','repo','out'):
        q.add_argument('--'+name,type=Path,required=True)
    q.add_argument('--image',required=True); q.add_argument('--allocations',type=int,default=20)
    q.add_argument('--arms',default='U'); q.add_argument('--batch-size',type=int,default=10)
    q=s.add_parser('run'); q.add_argument('--root',type=Path,required=True); q.add_argument('--index',type=int,required=True)
    q=s.add_parser('collect'); q.add_argument('--root',type=Path,required=True); q.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.command == 'prepare':
        require(a.batch_size > 0, 'Positive batch size required'); prepare(a)
    elif a.command == 'run': run(a.root.resolve(),a.index)
    else: collect(a.root.resolve(),a.out.resolve())


if __name__ == '__main__':
    main()
