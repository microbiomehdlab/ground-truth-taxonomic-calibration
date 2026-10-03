#!/usr/bin/env python3
"""Small direct-MaAsLin2 cost/equivalence trial on immutable DA3 allocations."""
from __future__ import annotations

import argparse
import itertools
import json
import math
import statistics
import subprocess
import time
from collections import Counter
from decimal import Decimal
from pathlib import Path

from audit_bracken_denominators import digest, require, table, write_table
from build_biomarker_abundance_input import parse_profile
from prepare_da_pilot import verify
from plan_da3_canary import finish

COHORTS = ('yachida', 'feng', 'zeller')
TOOLS = ('kraken2_bracken', 'metaphlan4')
KEYS = set(itertools.product(COHORTS, ('Adenoma', 'CRC'), TOOLS,
                             (5, 20), ('0.0001', '0.001')))
CODE = ('scripts/maaslin_timing.py', 'scripts/time_maaslin_context.R',
        'lib/maaslin_context.R', 'lib/maaslin_contract.R', 'lib/unpaired_null.R',
        'scripts/audit_bracken_denominators.py', 'scripts/build_biomarker_abundance_input.py',
        'scripts/prepare_da_pilot.py', 'scripts/plan_da3_canary.py')


def checked(folder, required):
    verify(folder)
    members = {line.split(None, 1)[1] for line in (folder / 'SHA256SUMS').read_text().splitlines()}
    require(set(required) <= members, 'Consumed file absent from checksums')


def fresh(out, sources):
    out = out.resolve()
    require(not out.exists(), 'Fresh output required: ' + str(out))
    for source in sources:
        source = source.resolve()
        require(out != source and out not in source.parents and source not in out.parents,
                'Overlapping input/output')


def prepare(source, out):
    source, out = source.resolve(), out.resolve()
    fresh(out, [source])
    checked(source, ('tasks.tsv', 'catalog.json'))
    require('status\tPASS_DRAFT_BATCH_PLAN_NOT_PRODUCTION_AUTHORIZATION' in
            (source / 'SUCCESS').read_text().splitlines(), 'Wrong DA3 plan status')
    tasks = table(source / 'tasks.tsv')
    members = {line.split(None, 1)[1] for line in (source / 'SHA256SUMS').read_text().splitlines()}
    selected, population, seen = {}, Counter(), set()
    for task in tasks:
        require(task['file'] in members, 'Unsealed batch')
        contexts = json.loads((source / task['file']).read_text())
        require(len(contexts) == int(task['contexts']), 'Batch count mismatch')
        for c in contexts:
            require(c['context_id'] not in seen, 'Duplicate DA3 context')
            seen.add(c['context_id'])
            population[(c['cohort'], c['background'], c['profiler'], int(c['n']))] += 1
            if c['arm'] != 'U':
                continue
            key = (c['cohort'], c['background'], c['profiler'], int(c['n']), str(Decimal(c['anchor'])))
            if key in KEYS and (key not in selected or c['context_id'] < selected[key]['context_id']):
                selected[key] = c
    require(set(selected) == KEYS, 'Missing prespecified timing strata')
    require(len(seen) == 120000, 'Expected the existing complete 120000-context plan')
    require(set(population) == set(itertools.product(COHORTS, ('Adenoma', 'CRC'), TOOLS, (5, 10, 15, 20)))
            and set(population.values()) == {2500}, 'Unexpected full DA3 design counts')
    contexts = [selected[k] for k in sorted(KEYS)]
    catalog = json.loads((source / 'catalog.json').read_text())
    profile_keys = {r['key'] for c in contexts for r in c['observations']}
    subset = dict(profiles={k: catalog['profiles'][k] for k in sorted(profile_keys)},
                  families=catalog['families'], targets=catalog['targets'])
    # Check design and catalog joins before creating any output.
    for c in contexts:
        validate_context(c, subset)
    out.mkdir(parents=True)
    (out / 'contexts.json').write_text(json.dumps(contexts, sort_keys=True) + '\n')
    (out / 'catalog.json').write_text(json.dumps(subset, sort_keys=True) + '\n')
    write_table(out / 'contexts.tsv', [dict(index=i, **{k: v for k, v in c.items()
                                                    if k != 'observations'}) for i, c in enumerate(contexts)])
    write_table(out / 'population.tsv', [dict(cohort=k[0], background=k[1], profiler=k[2],
                n=k[3], contexts=v) for k, v in sorted(population.items())])
    write_table(out / 'input_hashes.tsv', [dict(path=str(source / 'SHA256SUMS'),
                                               sha256=digest(source / 'SHA256SUMS'))])
    finish(out, 'PASS_MAASLIN_TIMING_PLAN')
    print('[PASS] Selected 48 contexts by design only; no results inspected')


def validate_context(c, catalog):
    obs = c['observations']
    require(len(obs) == 2 * c['n'] and len({r['sample_id'] for r in obs}) == len(obs),
            'Biological overlap or wrong sample count')
    require(Counter(r['group'] for r in obs) == {0: c['n'], 1: c['n']}, 'Unbalanced groups')
    features = catalog['families']['|'.join((c['cohort'], c['profiler'], c['background']))]
    require(features and len(features) == len(set(features)), 'Invalid frozen feature family')
    require(set(catalog['targets'][c['profiler']]) <= set(features), 'Targets missing from family')
    for r in obs:
        key = '|'.join((c['cohort'], c['profiler'], r['sample_id'], str(Decimal(r['dose']))))
        require(r['key'] == key, 'Wrong observation profile key')
        p = catalog['profiles'][key]
        require(p['sample_id'] == r['sample_id'] and p['cohort'] == c['cohort'] and
                p['profiler'] == c['profiler'] and p['condition'] == c['background'] and
                p['analysis_population'] == 'community' and
                Decimal(p['nominal_total_dose']) == Decimal(r['dose']), 'Catalog identity mismatch')
        require(Decimal(r['dose']) == (Decimal(c['anchor']) if r['group'] else 0),
                'Uniform timing context dose mismatch')
    return features


def run(plan, results, index, repo):
    checked(plan, ('contexts.json', 'catalog.json', 'contexts.tsv'))
    contexts = json.loads((plan / 'contexts.json').read_text())
    require(len(contexts) == 48 and 0 <= index < len(contexts), 'Invalid timing task')
    context = contexts[index]
    catalog = json.loads((plan / 'catalog.json').read_text())
    features = validate_context(context, catalog)
    out = results / ('task_%03d' % index)
    fresh(out, [plan])
    # A failed task retains inputs/logs but cannot acquire a SUCCESS marker.
    out.mkdir(parents=True)
    start = time.monotonic()
    rows, metadata = [], []
    hashes = {str(repo / 'analysis_v2' / name): digest(repo / 'analysis_v2' / name) for name in CODE}
    for r in context['observations']:
        p = catalog['profiles'][r['key']]
        path = Path(p['source_profile'])
        require(digest(path) == p['sha256'], 'Profile checksum changed: ' + str(path))
        hashes[str(path)] = p['sha256']
        values = parse_profile(path, context['profiler'])
        require(values and all(math.isfinite(v) and 0 <= v <= 1.00001 for v in values.values()),
                'Invalid native abundance')
        rows.append(dict(observation_id=r['sample_id'], **{f: values.get(f, 0) for f in features}))
        metadata.append(dict(observation_id=r['sample_id'], group=r['group']))
    write_table(out / 'abundance.tsv', rows)
    write_table(out / 'metadata.tsv', metadata)
    (out / 'context.json').write_text(json.dumps(context, sort_keys=True) + '\n')
    load_seconds = time.monotonic() - start
    with (out / 'R.stdout').open('w') as stdout, (out / 'R.stderr').open('w') as stderr:
        subprocess.run(['Rscript', str(repo / 'analysis_v2/scripts/time_maaslin_context.R'),
                        str(out.resolve()), str(repo.resolve())], stdout=stdout, stderr=stderr, check=True)
    require(all(digest(Path(p)) == h for p, h in hashes.items()), 'Input changed during fitting')
    timing = table(out / 'timing.tsv')
    require(len(timing) == 1 and timing[0]['maaslin_version'] == '1.18.0', 'Missing pinned backend timing')
    write_table(out / 'summary.tsv', [dict(index=index, context_id=context['context_id'],
                cohort=context['cohort'], background=context['background'], profiler=context['profiler'],
                n=context['n'], anchor=context['anchor'], family_n=len(features),
                input_seconds=load_seconds, task_seconds=time.monotonic() - start,
                output_bytes=sum(p.stat().st_size for p in out.rglob('*') if p.is_file()), **timing[0])])
    write_table(out / 'input_hashes.tsv', [dict(path=p, sha256=h) for p, h in sorted(hashes.items())] +
                [dict(path=str(plan / 'SHA256SUMS'), sha256=digest(plan / 'SHA256SUMS'))])
    finish(out, 'PASS_MAASLIN_TIMING_COMPUTATION')
    print('[PASS] Direct MaAsLin2 timing:', context['context_id'])


def collect(plan, results, out):
    fresh(out, [plan, results])
    failures, summaries, inputs = [], [], []
    try:
        checked(plan, ('contexts.json', 'population.tsv'))
        contexts = json.loads((plan / 'contexts.json').read_text())
        require(len(contexts) == 48, 'Expected 48 timing tasks')
        inputs.append(dict(path=str(plan / 'SHA256SUMS'), sha256=digest(plan / 'SHA256SUMS')))
        for i, c in enumerate(contexts):
            task = results / ('task_%03d' % i)
            try:
                checked(task, ('summary.tsv', 'comparison.tsv', 'input_hashes.tsv'))
                rows = table(task / 'summary.tsv')
                require(len(rows) == 1 and rows[0]['context_id'] == c['context_id'] and
                        int(rows[0]['index']) == i, 'Wrong task summary')
                require(all(rows[0][k] == str(c[k]) for k in ('cohort', 'background', 'profiler', 'n', 'anchor')),
                        'Timing design metadata differs')
                require(rows[0]['maaslin_version'] == '1.18.0', 'Wrong backend version')
                for field in ('task_seconds', 'maaslin_seconds', 'fast_seconds', 'output_bytes'):
                    require(math.isfinite(float(rows[0][field])) and float(rows[0][field]) >= 0,
                            'Invalid timing/resource value')
                summaries.extend(rows)
                inputs.append(dict(path=str(task / 'SHA256SUMS'), sha256=digest(task / 'SHA256SUMS')))
            except Exception as exc:
                failures.append(dict(index=i, error=str(exc)))
    except Exception as exc:
        failures.append(dict(index='plan', error=str(exc)))
    disagreements = sum(int(r['mismatched_features']) for r in summaries)
    out.mkdir(parents=True)
    status = dict(status='INCOMPLETE' if failures else ('REVIEW_DIFFERENCES' if disagreements else
                  'PASS_TIMING_AND_EQUIVALENCE_PENDING_COST_REVIEW'), expected_tasks=48,
                  completed_tasks=len(summaries), mismatched_features=disagreements,
                  failures=failures, production_authorized=False)
    (out / 'status.json').write_text(json.dumps(status, indent=2) + '\n')
    if summaries:
        write_table(out / 'timings.tsv', summaries)
    if len(summaries) == 48 and not failures:
        estimates = []
        population = table(plan / 'population.tsv')
        for tool in TOOLS:
            rows = [r for r in summaries if r['profiler'] == tool]
            for n in (5, 20):
                selected = [r for r in rows if int(r['n']) == n]
                write_rows = dict(profiler=tool, n=n, measured_contexts=len(selected),
                    median_maaslin_seconds=statistics.median(float(r['maaslin_seconds']) for r in selected),
                    min_task_seconds=min(float(r['task_seconds']) for r in selected),
                    mean_task_seconds=statistics.mean(float(r['task_seconds']) for r in selected),
                    max_task_seconds=max(float(r['task_seconds']) for r in selected),
                    mean_output_bytes=statistics.mean(int(r['output_bytes']) for r in selected))
                estimates.append(write_rows)
        write_table(out / 'cost_by_profiler_and_n.tsv', estimates)
        count = sum(int(r['contexts']) for r in population)
        total_seconds = total_bytes = 0.0
        for p in population:
            # Endpoints interpolate n=10/15. Arms/doses and allocations are not fully sampled.
            r = [s for s in summaries if all(s[k] == p[k] for k in ('cohort', 'background', 'profiler'))]
            weight = (int(p['n']) - 5) / 15.0
            for field in ('task_seconds', 'output_bytes'):
                value = sum(w * statistics.mean(float(s[field]) for s in r if int(s['n']) == n)
                            for n, w in ((5, 1 - weight), (20, weight))) * int(p['contexts'])
                if field == 'task_seconds':
                    total_seconds += value
                else:
                    total_bytes += value
        write_table(out / 'rough_projection.tsv', [dict(contexts=count,
                    single_worker_hours=total_seconds / 3600,
                    ideal_hours_at_30_workers=total_seconds / 3600 / 30,
                    model_and_input_output_GiB=total_bytes / 1024**3,
                    interpretation='ROUGH_SINGLE_CONTEXT_JOB_PROJECTION_NOT_SCHEDULER_FORECAST')])
    if inputs:
        write_table(out / 'input_hashes.tsv', inputs)
    # Retain launch provenance in the small report downloaded for review.
    for name in ('source_commit.txt', 'source_files.sha256', 'image.sha256', 'jobs.tsv'):
        source = plan.parent / name
        if source.is_file():
            (out / ('launch_' + name)).write_bytes(source.read_bytes())
    (out / 'README.txt').write_text(
        '48 full-family group-only LM fits, not clinical disease contrasts. Same frozen input, '
        'log2(1+a/1e-8), no extra normalization, and full-family BH as existing DA3. '
        'Native MaAsLin2 q and wrapper full-family q are distinct.\n'
        'MaAsLin2 elapsed includes package loading, model fitting and retained model/table I/O; '
        'fast elapsed averages 50 in-memory calls including transformation. Task elapsed also includes '
        'profile loading/checks and R startup, but excludes Slurm wait, container startup and final sealing.\n'
        'The projection interpolates n=10/15 between n=5/20; only the first allocation and two uniform '
        'doses were sampled. Other arms, cold caches and storage contention can change cost. '
        'Batching may amortize startup. Native models are saved, so inspect disk cost as well as time. '
        'Agreement checks implementation, not statistical calibration or scientific validity. '
        'No production is automatically authorized.\n')
    finish(out, 'PASS_TIMING_REPORT_COLLECTION_NOT_PRODUCTION_AUTHORIZATION')
    print(json.dumps(status, indent=2))
    return 1 if failures else 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p = sub.add_parser('run')
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--results', type=Path, required=True)
    p.add_argument('--index', type=int, required=True)
    p.add_argument('--repo', type=Path, required=True)
    p = sub.add_parser('collect')
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--results', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = vars(parser.parse_args())
    command = args.pop('command')
    result = globals()[command](**args)
    raise SystemExit(result or 0)
