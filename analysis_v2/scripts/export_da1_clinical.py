#!/usr/bin/env python3
"""Verify and export all 12 existing clinical MaAsLin2 contexts, without refitting."""
from __future__ import annotations
import argparse
import csv
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from audit_bracken_denominators import table, digest, require, write_table

COHORTS = ('yachida', 'feng', 'zeller')
PROFILERS = ('kraken2_bracken', 'metaphlan4')
BACKGROUNDS = ('Adenoma', 'CRC')


def verify(root, required):
    covered = {}
    for line in (root/'SHA256SUMS').read_text().splitlines():
        sha, name = line.split(maxsplit=1)
        p = Path(name)
        key = p.as_posix()
        require(not p.is_absolute() and '..' not in p.parts and key not in covered,
                'Unsafe/duplicate checksum member: '+name)
        require(root.resolve() in (root/p).resolve().parents, 'Checksum escapes bundle')
        require(digest(root/p)==sha, 'Checksum mismatch: '+str(root/p))
        covered[key] = sha
    require(set(required) <= set(covered), 'Unsealed required files: '+str(set(required)-set(covered)))
    return covered


def clinical_checkpoint(hashes, cohort):
    """Distinguish clinical person eligibility from the DA inventory's n grid."""
    pinned = {Path(r['path']):r['sha256'] for r in hashes}
    candidates = [p for p in pinned if p.name=='eligibility.tsv' and p.parent.name==cohort]
    clinical = []
    fields = {'cohort','comparison','sample_id','condition','group','age','sex','eligible','exclusion_reason'}
    for path in candidates:
        require(digest(path)==pinned[path], 'Changed eligibility evidence: '+str(path))
        rows = table(path)
        if rows and fields <= set(rows[0]): clinical.append((path,rows))
    require(len(clinical)==1, 'Expected exactly one clinical person-eligibility table for '+cohort+
            '; clinical candidates='+str([str(p) for p,_ in clinical])+'; all eligibility paths='+str([str(p) for p in candidates]))
    path,rows = clinical[0]
    for member in ('SUCCESS','SHA256SUMS'):
        p = path.parent/member
        require(p in pinned and digest(p)==pinned[p], 'Clinical checkpoint identity absent/changed: '+str(p))
    verify(path.parent, {'SUCCESS','eligibility.tsv'})
    require('status\tPASS_DA1_CHECKPOINT' in (path.parent/'SUCCESS').read_text().splitlines(),
            'Selected eligibility is not a passed DA1 checkpoint')
    require(all(r['cohort']==cohort for r in rows), 'Wrong clinical eligibility cohort')
    return path,rows


def build(plan, results, out, repo):
    plan, results, out, repo = [p.resolve() for p in (plan, results, out, repo)]
    require(not out.exists(), 'Choose a fresh report directory')
    for source in (plan, results):
        require(out!=source and out not in source.parents and source not in out.parents, 'Output overlaps inputs')
    plan_members = verify(plan, {'SUCCESS','contexts.tsv','families.tsv','observations.tsv','input_hashes.tsv'})
    contexts = [r for r in table(plan/'contexts.tsv') if r['analysis']=='DA1']
    grid = {(c,p,b) for c in COHORTS for p in PROFILERS for b in BACKGROUNDS}
    require(len(contexts)==12 and {(r['cohort'],r['profiler'],r['background']) for r in contexts}==grid,
            'Require exactly 12 clinical contexts')
    require(len({r['context_id'] for r in contexts})==12, 'Duplicate context IDs')
    families = table(plan/'families.tsv'); observations = table(plan/'observations.tsv')
    hashes = table(plan/'input_hashes.tsv')
    require(len({r['path'] for r in hashes})==len(hashes), 'Duplicate source hashes')
    # Recover the exact sealed checkpoint, without asking the user to rediscover it.
    checkpoints = {}; evidence = {str(plan/'SHA256SUMS'):digest(plan/'SHA256SUMS')}
    for name in ('analysis_v2/scripts/export_da1_clinical.py','analysis_v2/scripts/audit_da1_context.R',
                 'analysis_v2/DA1_CLINICAL_RESULTS.md','analysis_v2/taxon_identity_freeze.sha256'):
        path = repo/name; evidence[str(path)] = digest(path)
    for c in COHORTS:
        path,checkpoints[c] = clinical_checkpoint(hashes,c)
        for p in (path,path.parent/'SUCCESS',path.parent/'SHA256SUMS'): evidence[str(p)] = digest(p)
    alias_path = repo/'examples/spike_taxon_aliases.csv'; panel_path = repo/'spikes/spike_panel.tsv'
    for p in (alias_path,):
        expected = {r['sha256'] for r in hashes if Path(r['path']).name==p.name}
        require(expected=={digest(p)}, 'Target identities differ from pilot: '+p.name)
        evidence[str(p)] = digest(p)
    frozen = {}
    for line in (repo/'analysis_v2/taxon_identity_freeze.sha256').read_text().splitlines():
        sha,name = line.split(maxsplit=1); frozen[name] = sha
    require(frozen.get('spikes/spike_panel.tsv')==digest(panel_path), 'Panel differs from frozen identity policy')
    evidence[str(panel_path)] = digest(panel_path)
    with alias_path.open(newline='') as h:
        aliases = list(csv.DictReader(h))
    targets = {}
    for p in PROFILERS:
        targets[p] = {}
        for r in table(panel_path):
            matches = [a['alias'] for a in aliases if a['tool']==p and a['canonical']==r['taxon_name']]
            require(len(matches)==1, 'Ambiguous target alias')
            require(matches[0] not in targets[p], 'Duplicate target feature')
            targets[p][matches[0]] = r['label']
        require(len(targets[p])==10, 'Require ten targets')
    out.parent.mkdir(parents=True,exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=out.name+'.attempt_',dir=str(out.parent)))
    try:
        full = []; target_rows = []; summary = []; eligibility_rows = []
        for c in COHORTS:
            eligibility_rows.extend(checkpoints[c])
        for context in contexts:
            cid = context['context_id']; cohort = context['cohort']; profiler = context['profiler']
            require(Path(cid).name==cid and cid not in ('.','..',''), 'Unsafe context ID')
            source = results/cid
            needed = {'SUCCESS','provenance.json','context.tsv','metadata.tsv','abundance.tsv',
                      'fit/context_results.tsv','fit/feature_map.tsv','fit/SUCCESS','sessionInfo.txt'}
            members = verify(source, needed)
            require('status\tPASS_ENGINEERING_CONTEXT' in (source/'SUCCESS').read_text().splitlines(), 'Wrong context status')
            require('status\tPASS_CONTEXT' in (source/'fit/SUCCESS').read_text().splitlines(), 'Fit not successful')
            require('Maaslin2_1.18.0' in (source/'sessionInfo.txt').read_text(), 'Missing pinned MaAsLin2 session evidence')
            saved = table(source/'context.tsv')
            require(len(saved)==1 and all(saved[0].get(k)==v for k,v in context.items()) and
                    saved[0].get('inference_method','maaslin2')=='maaslin2', 'Changed context contract')
            provenance = json.loads((source/'provenance.json').read_text())
            require(provenance.get(str(plan/'SHA256SUMS'))==digest(plan/'SHA256SUMS'), 'Different pilot plan provenance')
            family = [r['feature'] for r in families if r['context_id']==cid]
            require(len(family)==len(set(family))==int(context['family_n']), 'Invalid frozen family')
            require({*targets[profiler]} <= set(family), 'Missing forced target')
            obs = [r for r in observations if r['context_id']==cid]
            keys = ('observation_id','biological_sample_id','group','spike_state','age','sex')
            meta = table(source/'metadata.tsv')
            require(meta==[{k:r[k] for k in keys} for r in obs], 'Metadata differs from frozen plan')
            eligible = [r for r in checkpoints[cohort] if r['comparison']==context['background']+'_vs_Control']
            require(eligible and len({r['sample_id'] for r in eligible})==len(eligible), 'Invalid eligibility ledger')
            included = {r['sample_id']:r for r in eligible if r['eligible']=='1'}
            require(len(obs)==len(included) and {r['biological_sample_id'] for r in obs}==set(included), 'Eligible people differ')
            for r in obs:
                e = included[r['biological_sample_id']]
                require(r['group']==e['group'] and float(r['age'])==float(e['age']) and r['sex']==e['sex'], 'Eligibility metadata mismatch')
            raw = table(source/'fit/context_results.tsv')
            require([r['feature'] for r in raw]==family, 'Fit family differs from frozen plan')
            if any(r['estimable']=='TRUE' for r in raw):
                require({'fit/maaslin_native/all_results.tsv','fit/maaslin_native/fits/models.rds'} <= set(members), 'Native fits not sealed')
            dest = stage/'contexts'/cid; dest.mkdir(parents=True)
            # Archive sufficient inputs for reproduction, not thousands of redundant models.
            for name in sorted(needed | {'SHA256SUMS'} | ({'fit/maaslin_native/all_results.tsv'} if
                          'fit/maaslin_native/all_results.tsv' in members else set())):
                target = dest/'source'/name; target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(source/name,target)
            # This is an inventory of the ORIGINAL complete source bundle, not a
            # checksum manifest claiming the reduced copied subset is complete.
            (dest/'source/SHA256SUMS').rename(dest/'source/original_bundle_checksums.txt')
            subprocess.run(['Rscript',str(repo/'analysis_v2/scripts/audit_da1_context.R'),str(dest/'source'),
                            str(dest/'clinical_results.tsv')],check=True)
            rows = table(dest/'clinical_results.tsv')
            for r in rows:
                record = dict(context_id=cid,cohort=cohort,profiler=profiler,
                              comparison=context['background']+'_vs_Control',target_label=targets[profiler].get(r['feature'],''),**r)
                full.append(record)
                if record['target_label']: target_rows.append(record)
            summary.append(dict(context_id=cid,cohort=cohort,profiler=profiler,comparison=context['background']+'_vs_Control',
                control_n=sum(r['group']=='0' for r in obs),disease_n=sum(r['group']=='1' for r in obs),
                excluded_n=len(eligible)-len(obs),family_n=len(rows),estimable_n=sum(r['estimable']=='TRUE' for r in rows),
                significant_two_sided_n=sum(r['significant_two_sided']=='TRUE' for r in rows),
                target_significant_n=sum(r['significant_two_sided']=='TRUE' and r['feature'] in targets[profiler] for r in rows),
                residual_df=len(obs)-4,status='PASS_SAVED_INFERENCE_AUDIT'))
            evidence[str(source/'SHA256SUMS')] = digest(source/'SHA256SUMS')
            # Detect mutation during the export of the required evidence.
            copied = needed | ({'fit/maaslin_native/all_results.tsv'} if
                               'fit/maaslin_native/all_results.tsv' in members else set())
            require(all(digest(source/name)==members[name] and digest(dest/'source'/name)==members[name]
                        for name in copied), 'Source changed during export')
        require(len(target_rows)==120, 'Expected 120 clinical target contrasts')
        write_table(stage/'full_family_results.tsv',full)
        write_table(stage/'target_results.tsv',target_rows)
        write_table(stage/'context_summary.tsv',summary)
        write_table(stage/'eligibility.tsv',eligibility_rows)
        write_table(stage/'source_evidence.tsv',[dict(path=p,sha256=s) for p,s in sorted(evidence.items())])
        for p,s in evidence.items(): require(digest(Path(p))==s, 'Evidence changed during export')
        require(all(digest(plan/name)==sha for name,sha in plan_members.items()), 'Plan changed during export')
        (stage/'status.json').write_text(json.dumps(dict(status='PASS_DA1_CLINICAL_EXPORT',contexts=12,
            target_rows=120,full_family_rows=len(full),new_maaslin_fits=0,scientific_interpretation='PENDING_REVIEW'),indent=2)+'\n')
        (stage/'SUCCESS').write_text('status\tPASS_DA1_CLINICAL_EXPORT\ncontexts\t12\nnew_maaslin_fits\t0\n')
        shutil.copy2(repo/'analysis_v2/DA1_CLINICAL_RESULTS.md',stage/'README.md')
        (stage/'SHA256SUMS').write_text(''.join(digest(p)+'  '+str(p.relative_to(stage))+'\n'
            for p in sorted(stage.rglob('*')) if p.is_file()))
        os.rename(stage,out)
        print('[PASS] DA1 clinical results: 12 verified contexts, 120 target contrasts:',out)
    except Exception as exc:
        (stage/'FAILED.json').write_text(json.dumps({'error':str(exc)}))
        print('[FAIL] Audit diagnostics retained:',stage)
        raise


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan','results','outdir','repo'): parser.add_argument('--'+name,type=Path,required=True)
    a = parser.parse_args(); build(a.plan,a.results,a.outdir,a.repo)
