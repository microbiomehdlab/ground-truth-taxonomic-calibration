#!/usr/bin/env python3
"""Original controls versus community-spiked adenomas; frozen DA1 clinical design."""
from __future__ import annotations
import argparse
import json
import math
import shutil
import subprocess
import tempfile
from decimal import Decimal
from pathlib import Path
from audit_bracken_denominators import table, digest, require, write_table
from build_biomarker_abundance_input import parse_profile
from export_da1_clinical import verify

COHORTS = ('yachida', 'feng', 'zeller')
TOOLS = ('kraken2_bracken', 'metaphlan4')
DOSES = ('0.0001', '0.0005', '0.001', '0.005', '0.01', '0.05', '0.1')

def seal(root, status):
    (root/'SUCCESS').write_text('status\t'+status+'\n')
    (root/'SHA256SUMS').write_text(''.join(
        digest(p)+'  '+str(p.relative_to(root))+'\n'
        for p in sorted(root.rglob('*')) if p.is_file() and p.name != 'SHA256SUMS'))

def prepare(root, clinical, inventory):
    root, clinical, inventory = [p.resolve() for p in (root, clinical, inventory)]
    require(not root.joinpath('plan').exists(), 'Plan already exists')
    require(all(root != p and root not in p.parents and p not in root.parents
                for p in (clinical, inventory)), 'Output overlaps inputs')
    covered = verify(clinical, {'SUCCESS', 'status.json', 'context_summary.tsv', 'target_results.tsv'})
    require(json.loads((clinical/'status.json').read_text())['status'] == 'PASS_DA1_CLINICAL_EXPORT',
            'Require completed DA1 clinical export')
    baseline = [r for r in table(clinical/'context_summary.tsv')
                if r['comparison'] == 'Adenoma_vs_Control']
    require(len(baseline) == 6 and {(r['cohort'],r['profiler']) for r in baseline}
            == {(c,t) for c in COHORTS for t in TOOLS}, 'Incomplete clinical baseline grid')
    evidence = {str(clinical/'SHA256SUMS'): digest(clinical/'SHA256SUMS')}
    pending = []
    targets = table(clinical/'target_results.tsv')
    for b in baseline:
        cohort, tool, cid = b['cohort'], b['profiler'], b['context_id']
        inv = inventory/cohort
        verify(inv, {'SUCCESS', 'profile_inventory.tsv'})
        require('status\tPASS_DESIGN_INVENTORY' in (inv/'SUCCESS').read_text().splitlines(),
                'Require sealed design inventory')
        evidence[str(inv/'SHA256SUMS')] = digest(inv/'SHA256SUMS')
        source = clinical/'contexts'/cid/'source'
        require(all(str((source/name).relative_to(clinical)) in covered
                    for name in ('context.tsv','metadata.tsv','abundance.tsv')), 'Unsealed clinical source')
        context = table(source/'context.tsv')
        require(len(context) == 1 and context[0]['paired'] == '0'
                and context[0]['covariates'] == 'age,sex'
                and Decimal(context[0]['nominal_total_dose']) == 0, 'Unexpected clinical model')
        metadata = table(source/'metadata.tsv')
        abundance = table(source/'abundance.tsv')
        family = list(abundance[0])[1:]
        require(len(family) == int(b['family_n']) and len(set(family)) == len(family), 'Family mismatch')
        require(len({m['biological_sample_id'] for m in metadata}) == len(metadata)
                and len(metadata) == len(abundance), 'Duplicate/missing clinical people')
        require({m['group'] for m in metadata} == {'0','1'}
                and all(m['spike_state'] == 'original' and m['sex'] in ('Female','Male')
                        and math.isfinite(float(m['age'])) for m in metadata), 'Invalid clinical metadata')
        require({m['observation_id'] for m in metadata} == {a['observation_id'] for a in abundance},
                'Metadata/matrix mismatch')
        target_map = {r['feature']: r['target_label'] for r in targets if r['context_id'] == cid}
        require(len(target_map) == 10 and set(target_map) <= set(family), 'Incomplete target family')
        profiles = [r for r in table(inv/'profile_inventory.tsv') if r['profiler'] == tool
                    and r['condition'] == 'Adenoma' and r['analysis_population'] == 'community']
        cases = {m['biological_sample_id'] for m in metadata if m['group'] == '1'}
        for dose in DOSES:
            selected = [r for r in profiles if Decimal(r['nominal_total_dose']) == Decimal(dose)]
            # Inventory may contain excluded people; retain exactly DA1's eligible population.
            selected = [r for r in selected if r['sample_id'] in cases]
            require(len(selected) == len(cases) and {r['sample_id'] for r in selected} == cases,
                    'Missing/duplicate eligible adenoma profiles: '+cohort+' '+tool+' '+dose)
            pending.append(dict(cohort=cohort, profiler=tool, baseline_context=cid,
                                dose=dose, profiles=selected, targets=target_map, source=str(source)))
    require(len(pending) == 42, 'Expected 42 fits')
    require(all(digest(Path(p)) == h for p,h in evidence.items()), 'Inputs changed during preparation')
    plan = root/'plan'; plan.mkdir(parents=True)
    tasks = []
    for i, task in enumerate(pending):
        folder = plan/('clinical_spike_%02d'%i); folder.mkdir()
        source = Path(task.pop('source'))
        for name in ('metadata.tsv', 'abundance.tsv'):
            shutil.copyfile(source/name, folder/name)
        task['context_id'] = folder.name
        (folder/'task.json').write_text(json.dumps(task, indent=2)+'\n')
        seal(folder, 'PASS_CLINICAL_SPIKE_PLAN')
        tasks.append(dict(index=i, context_id=folder.name, cohort=task['cohort'],
                          profiler=task['profiler'], nominal_total_dose=task['dose']))
    write_table(plan/'tasks.tsv', tasks)
    write_table(plan/'input_hashes.tsv', [dict(path=p, sha256=h) for p,h in sorted(evidence.items())])
    # Reuse, never refit, zero-dose clinical inference.
    write_table(plan/'baseline_targets.tsv', [r for r in targets if r['context_id'] in
                                             {b['context_id'] for b in baseline}])
    (plan/'baseline_identity.json').write_text(json.dumps(dict(path=str(clinical),
                                                              sha256=evidence[str(clinical/'SHA256SUMS')])))
    seal(plan, 'PASS_CLINICAL_SPIKE_PLAN')
    print('[PASS] Prepared 42 clinical spike fits; six zero-dose references reused', flush=True)

def run(root, index, repo):
    plan = root/'plan'; verify(plan, {'tasks.tsv', 'SUCCESS'})
    tasks = table(plan/'tasks.tsv'); require(0 <= index < len(tasks), 'Invalid task index')
    source = plan/tasks[index]['context_id']; verify(source, {'task.json','metadata.tsv','abundance.tsv'})
    task = json.loads((source/'task.json').read_text())
    destination = root/'results'/task['context_id']
    require(not destination.exists(), 'Result exists; verify it before attempting a rerun')
    attempts = root/'attempts'; attempts.mkdir(exist_ok=True)
    attempt = Path(tempfile.mkdtemp(prefix=task['context_id']+'_', dir=attempts))
    metadata, abundance = table(source/'metadata.tsv'), table(source/'abundance.tsv')
    family = list(abundance[0])[1:]
    by_id = {r['observation_id']: r for r in abundance}
    profiles = {r['sample_id']:r for r in task['profiles']}
    for m in metadata:
        if m['group'] == '1':
            profile = profiles[m['biological_sample_id']]; path = Path(profile['source_profile'])
            require(digest(path) == profile['sha256'], 'Changed spike profile: '+str(path))
            values = parse_profile(path, task['profiler'])
            require(all(math.isfinite(v) for v in values.values()), 'Nonfinite profile abundance')
            by_id[m['observation_id']].update({f:values.get(f,0.0) for f in family})
            m['spike_state'] = 'spiked'
            require(digest(path) == profile['sha256'], 'Profile changed while reading')
    write_table(attempt/'abundance.tsv', abundance)
    write_table(attempt/'metadata.tsv', metadata)
    shutil.copyfile(source/'task.json', attempt/'task.json')
    subprocess.run(['Rscript',str(repo/'analysis_v2/scripts/fit_clinical_spike.R'),str(attempt),str(repo)], check=True)
    write_table(attempt/'profile_evidence.tsv', task['profiles'])
    (attempt/'plan_sha256.txt').write_text(digest(source/'SHA256SUMS')+'\n')
    seal(attempt, 'PASS_CLINICAL_SPIKE_FIT')
    destination.parent.mkdir(exist_ok=True); attempt.rename(destination)
    print('[PASS]', task['context_id'], flush=True)

def collect(root):
    plan=root/'plan'; verify(plan, {'tasks.tsv','baseline_targets.tsv'})
    rows=[]; failures=[]
    for task in table(plan/'tasks.tsv'):
        folder=root/'results'/task['context_id']
        try:
            verify(folder, {'SUCCESS','clinical_results.tsv','plan_sha256.txt'})
            require('status\tPASS_CLINICAL_SPIKE_FIT' in (folder/'SUCCESS').read_text().splitlines(),
                    'Result did not pass clinical fit contract')
            require((folder/'plan_sha256.txt').read_text().strip() ==
                    digest(plan/task['context_id']/'SHA256SUMS'), 'Wrong result plan')
            result=table(folder/'clinical_results.tsv')
            family=list(table(plan/task['context_id']/'abundance.tsv')[0])[1:]
            require(len(result)==len(family) and [r['feature'] for r in result]==family,
                    'Incomplete/reordered full family')
            targets=json.loads((plan/task['context_id']/'task.json').read_text())['targets']
            for r in result:
                if r['feature'] in targets:
                    rows.append(dict(context_id=task['context_id'],cohort=task['cohort'],
                                     profiler=task['profiler'],nominal_total_dose=task['nominal_total_dose'],
                                     nominal_target_dose=str(Decimal(task['nominal_total_dose'])/10),
                                     target_label=targets[r['feature']], **r))
        except Exception as exc: failures.append(dict(context_id=task['context_id'],error=str(exc)))
    require(not failures, 'Incomplete fits: '+json.dumps(failures))
    require(len(rows)==420, 'Incomplete target grid')
    report=root/'REPORT'; require(not report.exists(), 'Report exists'); report.mkdir()
    write_table(report/'spiked_target_results.tsv',rows)
    shutil.copyfile(plan/'baseline_targets.tsv',report/'baseline_target_results.tsv')
    write_table(report/'full_family_inventory.tsv',[
        dict(context_id=t['context_id'],path=str(root/'results'/t['context_id']/'clinical_results.tsv'),
             sha256=digest(root/'results'/t['context_id']/'clinical_results.tsv')) for t in table(plan/'tasks.tsv')])
    (report/'status.json').write_text(json.dumps(dict(status='PASS_CLINICAL_SPIKE_ASSOCIATION',
        fits=42,spiked_target_rows=420,baseline_target_rows=60,scientific_interpretation='PENDING_REVIEW'),indent=2)+'\n')
    seal(report,'PASS_CLINICAL_SPIKE_ASSOCIATION')
    print('[PASS] Clinical spike report:', report, flush=True)

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['prepare','run','collect'])
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--clinical',type=Path); parser.add_argument('--inventory',type=Path)
    parser.add_argument('--repo',type=Path,required=True); parser.add_argument('--index',type=int)
    args=parser.parse_args()
    if args.stage=='prepare': prepare(args.root,args.clinical,args.inventory)
    elif args.stage=='run': run(args.root,args.index,args.repo)
    else: collect(args.root)
