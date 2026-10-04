#!/usr/bin/env python3
"""Full DA2 grid: descriptive target response plus frozen full-family paired tests."""
from __future__ import annotations
import argparse
import csv
import fcntl
import gzip
import json
import math
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from audit_bracken_denominators import table, digest, require, write_table
from build_biomarker_abundance_input import parse_profile
from export_da1_clinical import verify
from link_da3_recovery import normalize, close, summary_values

COHORTS = ('yachida','feng','zeller')
TOOLS = ('kraken2_bracken','metaphlan4')
CONDITIONS = ('Control','Adenoma','CRC')
DOSES = ('0.0001','0.0005','0.001','0.005','0.01','0.05')
COUNTS = {'yachida':201,'feng':154,'zeller':156}


def seal(root, status):
    (root/'SUCCESS').write_text('status\t'+status+'\n')
    (root/'SHA256SUMS').write_text(''.join(digest(p)+'  '+str(p.relative_to(root))+'\n'
        for p in sorted(root.rglob('*')) if p.is_file() and p.name!='SHA256SUMS'))


def gzwrite(path, rows):
    require(bool(rows),'Empty output table')
    with gzip.open(path,'wt',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)


def gzread(path):
    with gzip.open(path,'rt',newline='') as h: return list(csv.DictReader(h,delimiter='\t'))


def identity(row):
    return (row['profiler'],row['sample_id'],row['analysis_population'],row['profile_id'])


def make_contexts(profiles, families, targets, cohort, fixture=False):
    index={identity(r):r for r in profiles}
    require(len(index)==len(profiles),'Duplicate inventory profile')
    output=[]
    for tool in TOOLS:
        for condition in CONDITIONS:
            base={r['sample_id']:r for r in profiles if r['profiler']==tool and
                  r['condition']==condition and r['analysis_population']=='community' and Decimal(r['nominal_total_dose'])==0}
            require(base,'Missing condition baseline')
            family=[r['feature'] for r in families if r['profiler']==tool and r['condition']==condition]
            require(len(family)==len(set(family)) and set(targets[tool])<=set(family),'Invalid family')
            arms=defaultdict(dict)
            for r in profiles:
                if r['profiler']!=tool or r['condition']!=condition or Decimal(r['nominal_total_dose'])==0: continue
                sid=r['sample_id']; pop=r['analysis_population']
                require(pop in ('community','independent') and sid in base,'Missing baseline/invalid population')
                prefix=sid+'_'; stem=r['profile_id']
                require(stem.startswith(prefix),'Profile ID/sample mismatch')
                label,encoded=stem[len(prefix):].rsplit('_f',1)
                d=Decimal(r['nominal_total_dose'])
                require(Decimal(encoded.replace('p','.'))==d,'Profile ID/dose mismatch')
                require(label=='CRCpanel' if pop=='community' else label in targets[tool].values(),'Invalid spike label')
                key=(pop,label,str(d.normalize()))
                require(sid not in arms[key],'Duplicate dose/person')
                arms[key][sid]=r
            independent_sets=[]
            for (pop,label,d),spikes in sorted(arms.items()):
                if pop=='community': require(set(spikes)==set(base),'Incomplete community people')
                else:
                    independent_sets.append(set(spikes))
                    if not fixture: require(len(spikes)==10,'Require ten independent people per condition')
                output.append(dict(cohort=cohort,profiler=tool,condition=condition,population=pop,
                    spike_label=label,nominal_total_dose=d,family=family,
                    pairs=[dict(sample_id=s,baseline=base[s],spiked=spikes[s]) for s in sorted(spikes)]))
            require(independent_sets and all(s==independent_sets[0] for s in independent_sets),'Independent subset changes across series')
            if not fixture:
                expected={('community','CRCpanel',str(Decimal(d).normalize())) for d in DOSES+('0.1',)}
                expected|={('independent',label,str(Decimal(d).normalize())) for label in targets[tool].values() for d in DOSES}
                require(set(arms)==expected,'Incomplete dose/target grid')
    for i,c in enumerate(output): c['context_id']=cohort+'_da2_%04d'%i
    return output


def prepare(root, inventory, bracken, metaphlan, repo, fixture=False):
    root,inventory,bracken,metaphlan,repo=[p.resolve() for p in (root,inventory,bracken,metaphlan,repo)]
    require(all(root!=p and root not in p.parents and p not in root.parents for p in (inventory,bracken,metaphlan)),
            'Output overlaps source bundles')
    plan=root/'plan'; require(not plan.exists(),'Plan exists')
    evidence={}
    def bundle(p,required):
        verify(p,required); evidence[str(p/'SHA256SUMS')]=digest(p/'SHA256SUMS')
    bundle(inventory,{'SUCCESS'})
    bundle(metaphlan,{'SUCCESS'}|{c+'/endpoints/paired_endpoints.tsv' for c in COHORTS})
    require('status\tPASS_ENDPOINTS' in (metaphlan/'SUCCESS').read_text().splitlines(),'Wrong MetaPhlAn bundle')
    with (repo/'examples/spike_taxon_aliases.csv').open() as h: aliases=list(csv.DictReader(h))
    panel=table(repo/'spikes/spike_panel.tsv')
    require(panel and len({Decimal(r['weight']) for r in panel})==1 and Decimal(panel[0]['weight'])>0,
            'Equal positive community weights required')
    for line in (repo/'analysis_v2/taxon_identity_freeze.sha256').read_text().splitlines():
        sha,name=line.split(maxsplit=1)
        require(digest(repo/name)==sha,'Frozen taxon identity changed')
    targets={t:{} for t in TOOLS}
    for t in TOOLS:
        for p in panel:
            matches=[a['alias'] for a in aliases if a['tool']==t and a['canonical']==p['taxon_name']]
            require(len(matches)==1 and matches[0] not in targets[t],'Invalid target aliases')
            targets[t][matches[0]]=p['label']
        if not fixture: require(len(targets[t])==10,'Ten targets required')
    pending=[]; total_contexts=0; total_points=0
    for cohort in COHORTS:
        inv=inventory/cohort
        bundle(inv,{'SUCCESS','profile_inventory.tsv','feature_families.tsv','input_hashes.tsv'})
        require('status\tPASS_DESIGN_INVENTORY' in (inv/'SUCCESS').read_text().splitlines(),'Wrong inventory status')
        profiles=table(inv/'profile_inventory.tsv'); families=table(inv/'feature_families.tsv')
        require(all(r['cohort']==cohort for r in profiles+families),'Wrong inventory cohort')
        if not fixture: require(len({r['sample_id'] for r in profiles})==COUNTS[cohort],'Wrong cohort size')
        alias_hash={r['sha256'] for r in table(inv/'input_hashes.tsv') if Path(r['path']).name=='spike_taxon_aliases.csv'}
        require(alias_hash=={digest(repo/'examples/spike_taxon_aliases.csv')},'Frozen aliases changed')
        contexts=make_contexts(profiles,families,targets,cohort,fixture)
        expected={}
        for c in contexts:
            labels=targets[c['profiler']].values() if c['population']=='community' else [c['spike_label']]
            for pair in c['pairs']:
                for label in labels:
                    key=(c['profiler'],pair['sample_id'],c['population'],pair['spiked']['profile_id'],label)
                    require(key not in expected,'Repeated target response'); expected[key]=(c,pair)
        b=bracken/cohort; bundle(b,{'summary.json','bracken_recovery_comparison.tsv'})
        status=json.loads((b/'summary.json').read_text())
        require(status['cohort']==cohort and status['status']=='PASS_ENDPOINT_CONSTRUCTION','Wrong Bracken bundle')
        sources=[('kraken2_bracken',table(b/'bracken_recovery_comparison.tsv')),
                 ('metaphlan4',[r for r in table(metaphlan/cohort/'endpoints/paired_endpoints.tsv') if r['profiler']=='metaphlan4'])]
        points={}
        for tool,rows in sources:
            for r in rows:
                key=(tool,r['sample_id'],r['analysis_population'],r['profile_id'],r['target_label'])
                require(key in expected and key not in points and r['cohort']==cohort,'Unexpected/duplicate endpoint')
                c,pair=expected[key]; point=normalize(r,tool); point['population']=c['population']
                require(point['condition']==c['condition'] and Decimal(point['nominal_total_dose'])==Decimal(c['nominal_total_dose']), 'Endpoint condition/dose mismatch')
                require(Path(point['source_profile'])==Path(pair['spiked']['source_profile']) and
                        Path(point['source_baseline'])==Path(pair['baseline']['source_profile']),'Endpoint physical source mismatch')
                point.update(context_id=c['context_id'],target_feature=next(f for f,l in targets[tool].items() if l==r['target_label']))
                if tool=='kraken2_bracken':
                    require(r['target_alias']==point['target_feature'],'Wrong Bracken target alias')
                    old=r['legacy_recovery_ratio']
                else: old=r['response_ratio']
                old=float(old); require(math.isfinite(old),'Nonfinite legacy ratio')
                point['legacy_reference_ratio']=old
                point['transformed_native_difference']=math.log2(1+point['spiked_native']/1e-8)-math.log2(1+point['baseline_native']/1e-8)
                point['nominal_target_dose']=str(Decimal(c['nominal_total_dose'])/(10 if c['population']=='community' else 1))
                points[key]=point
        require(set(points)==set(expected),'Missing target endpoints')
        for tool in TOOLS:
            for condition in CONDITIONS:
                cc=[c for c in contexts if c['profiler']==tool and c['condition']==condition]
                pp=[p for p in points.values() if p['profiler']==tool and p['condition']==condition]
                pending.append((cc,pp))
        total_contexts+=len(contexts); total_points+=len(points)
    if not fixture: require(total_contexts==1206 and total_points==82340,'Incomplete production grid')
    require(all(digest(Path(p))==s for p,s in evidence.items()),'Inputs changed during preparation')
    plan.mkdir()
    tasks=[]
    for i,(contexts,points) in enumerate(pending):
        task=plan/('task_%02d'%i); task.mkdir()
        (task/'contexts.json').write_text(json.dumps(contexts))
        (task/'targets.json').write_text(json.dumps(targets))
        gzwrite(task/'target_responses.tsv.gz',points)
        seal(task,'PASS_DA2_TASK_PLAN')
        tasks.append(dict(index=i,folder=task.name,contexts=len(contexts),target_rows=len(points)))
    write_table(plan/'tasks.tsv',tasks)
    write_table(plan/'input_hashes.tsv',[dict(path=p,sha256=s) for p,s in sorted(evidence.items())])
    seal(plan,'PASS_DA2_PLAN')
    print('[PASS] DA2 planned:',total_contexts,'contexts;',total_points,'target responses',flush=True)


def run(root,index,repo):
    (root/'locks').mkdir(exist_ok=True)
    with (root/'locks'/('%02d.lock'%index)).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        return run_locked(root,index,repo)


def run_locked(root,index,repo):
    task=root/'plan'/('task_%02d'%index); final=root/'results'/task.name
    final.parent.mkdir(exist_ok=True)
    verify(task,{'SUCCESS','contexts.json','targets.json','target_responses.tsv.gz'})
    identity_hash=digest(task/'SHA256SUMS')
    code_names=('scripts/da2_dose_response.py','scripts/run_da2_grid.R','lib/maaslin_contract.R',
                'lib/paired_difference_context.R','scripts/build_biomarker_abundance_input.py')
    code_hashes={n:digest(repo/'analysis_v2'/n) for n in code_names}
    expected_identity=dict(plan_sha256=identity_hash,code_hashes=code_hashes)
    if final.exists():
        verify(final,{'SUCCESS','identity.json','full_family_results.tsv.gz','target_results.tsv','dose_response_summary.tsv'})
        require(json.loads((final/'identity.json').read_text())==expected_identity,'Wrong resumed task inputs/code')
        print('[PASS] Existing task:',index); return
    stage=Path(tempfile.mkdtemp(prefix=task.name+'.attempt_',dir=str(final.parent)))
    try:
        contexts=json.loads((task/'contexts.json').read_text()); points=gzread(task/'target_responses.tsv.gz')
        grouped=defaultdict(list)
        for p in points: grouped[(p['context_id'],p['sample_id'])].append(p)
        sources={}; baselines={}; summary=[]
        with tempfile.TemporaryDirectory(prefix='da2_inputs_') as temporary:
            scratch=Path(temporary)
            for c in contexts:
                cid=c['context_id']; family=c['family']; folder=scratch/cid; folder.mkdir()
                meta=[]
                with (folder/'abundance.tsv').open('w',newline='') as h:
                    w=csv.writer(h,delimiter='\t'); w.writerow(['observation_id']+family)
                    for pair in c['pairs']:
                        values={}
                        for role in ('baseline','spiked'):
                            item=pair[role]; path=Path(item['source_profile']); sha=item['sha256']
                            require(digest(path)==sha,'Changed profile: '+str(path)); sources[str(path)]=sha
                            if role=='baseline' and str(path) in baselines: value=baselines[str(path)]
                            else: value=parse_profile(path,c['profiler'])
                            require(value and all(math.isfinite(v) and 0<=v<=1 for v in value.values()),'Invalid native profile')
                            if role=='baseline': baselines[str(path)]=value
                            values[role]=value
                            oid=pair['sample_id']+'_'+role
                            w.writerow([oid]+[value.get(f,0) for f in family])
                            meta.append(dict(observation_id=oid,biological_sample_id=pair['sample_id'],
                                spike_state='original' if role=='baseline' else 'spiked',group=int(role=='spiked')))
                        for p in grouped[(cid,pair['sample_id'])]:
                            f=p['target_feature']
                            require(close(values['baseline'].get(f,0),p['baseline_native']) and
                                    close(values['spiked'].get(f,0),p['spiked_native']),'Endpoint/native-profile mismatch')
                write_table(folder/'metadata.tsv',meta)
            with (stage/'model.out').open('w') as out,(stage/'model.err').open('w') as err:
                subprocess.run(['Rscript',str(repo/'analysis_v2/scripts/run_da2_grid.R'),str(task),str(scratch),str(stage),str(repo)],
                               stdout=out,stderr=err,check=True)
        cells=defaultdict(list)
        for p in points: cells[(p['context_id'],p['target_label'])].append(p)
        stats={(r['context_id'],r['target_label']):r for r in table(stage/'target_results.tsv')}
        require(set(cells)==set(stats),'Paired target results missing')
        for key,rows in sorted(cells.items()):
            r=rows[0]; fit=stats[key]
            require(len({p['sample_id'] for p in rows})==len(rows)==int(fit['n_pairs']),'Wrong pair count')
            require(close(sum(float(p['transformed_native_difference']) for p in rows)/len(rows),fit['beta']),'Paired target effect mismatch')
            s={k:r[k] for k in ('context_id','cohort','profiler','condition','population','target_label','nominal_total_dose','nominal_target_dose','reference_type')}
            s.update(n_pairs=len(rows),positive_native_changes=sum(float(p['native_increment'])>0 for p in rows),
                     negative_native_changes=sum(float(p['native_increment'])<0 for p in rows),
                     zero_native_changes=sum(float(p['native_increment'])==0 for p in rows),
                     negative_recovery=sum(float(p['recovery_ratio'])<0 for p in rows),
                     baseline_positive=sum(float(p['baseline_native'])>0 for p in rows),
                     spiked_positive=sum(float(p['spiked_native'])>0 for p in rows))
            for field in ('baseline_native','spiked_native','native_increment','transformed_native_difference',
                          'recovery_ratio','legacy_reference_ratio','expected_added_signal','recovered_added_signal',
                          'achieved_total_read_fraction','achieved_target_read_fraction'):
                s.update(summary_values(field,[float(p[field]) for p in rows]))
            s.update({k:fit[k] for k in ('beta','ci_low','ci_high','raw_p','wrapper_q','status','family_n','positive_discovery','significant_two_sided')})
            summary.append(s)
        write_table(stage/'dose_response_summary.tsv',summary)
        shutil.copy2(task/'target_responses.tsv.gz',stage/'target_responses.tsv.gz')
        require(all(digest(Path(p))==s for p,s in sources.items()),'Profiles changed during analysis')
        require(digest(task/'SHA256SUMS')==identity_hash,'Plan changed during analysis')
        require(all(digest(repo/'analysis_v2'/n)==s for n,s in code_hashes.items()),'Code changed during analysis')
        write_table(stage/'profile_hashes.tsv',[dict(path=p,sha256=s) for p,s in sorted(sources.items())])
        (stage/'identity.json').write_text(json.dumps(expected_identity))
        seal(stage,'PASS_DA2_TASK'); os.rename(stage,final)
        print('[PASS] DA2 task:',index,'contexts:',len(contexts),flush=True)
    except Exception as e:
        (stage/'FAILED.json').write_text(json.dumps(dict(error=str(e)))); raise


def collect(root,repo):
    out=root/'REPORT'; require(not out.exists(),'Report exists')
    plan=root/'plan'; verify(plan,{'SUCCESS','tasks.tsv'})
    tasks=table(plan/'tasks.tsv'); failures=[]; summaries=[]; targets=[]; full=[]; context_rows=[]; pairing=[]
    for t in tasks:
        folder=root/'results'/t['folder']
        try:
            verify(folder,{'SUCCESS','identity.json','full_family_results.tsv.gz','target_results.tsv','dose_response_summary.tsv','target_responses.tsv.gz'})
            require(json.loads((folder/'identity.json').read_text())['plan_sha256']==digest(plan/t['folder']/'SHA256SUMS'),'Task belongs to other plan')
            summary=table(folder/'dose_response_summary.tsv'); target=table(folder/'target_results.tsv')
            contexts=json.loads((plan/t['folder']/'contexts.json').read_text())
            context_lookup={c['context_id']:c for c in contexts}
            expected={(c['context_id'],f) for c in contexts for f in c['family']}
            seen=set()
            with gzip.open(folder/'full_family_results.tsv.gz','rt',newline='') as h:
                for r in csv.DictReader(h,delimiter='\t'):
                    key=(r['context_id'],r['feature'])
                    require(key in expected and key not in seen,'Duplicate/unexpected full-family row'); seen.add(key)
                    c=context_lookup[r['context_id']]
                    require(int(r['n_pairs'])==len(c['pairs']) and int(r['family_n'])==len(c['family']), 'Wrong reported family/pair size')
            require(seen==expected,'Missing full-family results')
            ids={c['context_id'] for c in contexts}
            require({r['context_id'] for r in summary}==ids and len(ids)==int(t['contexts']),'Missing task contexts')
            require(len(summary)==len(target)==len({(r['context_id'],r['target_label']) for r in summary}) and
                    {(r['context_id'],r['target_label']) for r in summary}=={(r['context_id'],r['target_label']) for r in target},'Target summary coverage mismatch')
            require(len(gzread(folder/'target_responses.tsv.gz'))==int(t['target_rows']),'Missing person-target rows')
            summaries.extend(summary); targets.extend(target)
            full.append(folder)
            for c in contexts:
                context_rows.append(dict({k:v for k,v in c.items() if k not in ('family','pairs')},
                                         family_n=len(c['family']),n_pairs=len(c['pairs'])))
                for p in c['pairs']:
                    pairing.append(dict(context_id=c['context_id'],sample_id=p['sample_id'],
                        baseline_profile=p['baseline']['source_profile'],baseline_sha256=p['baseline']['sha256'],
                        spiked_profile=p['spiked']['source_profile'],spiked_sha256=p['spiked']['sha256']))
        except Exception as e: failures.append(dict(task=t['folder'],error=str(e)))
    if failures:
        (root/'collection_failures.json').write_text(json.dumps(failures,indent=2)); raise ValueError('Incomplete DA2; see collection_failures.json')
    stage=Path(tempfile.mkdtemp(prefix='report.attempt_',dir=str(root)))
    write_table(stage/'dose_response_summary.tsv',summaries); write_table(stage/'target_results.tsv',targets)
    write_table(stage/'contexts.tsv',context_rows); gzwrite(stage/'pairing_manifest.tsv.gz',pairing)
    shutil.copy2(plan/'input_hashes.tsv',stage/'input_bundle_hashes.tsv')
    for name in ('full_family_results.tsv.gz','target_responses.tsv.gz'):
        with gzip.open(stage/name,'wt',newline='') as h:
            writer=None
            for folder in full:
                with gzip.open(folder/name,'rt',newline='') as src:
                    reader=csv.DictReader(src,delimiter='\t')
                    if writer is None: writer=csv.DictWriter(h,fieldnames=reader.fieldnames,delimiter='\t'); writer.writeheader()
                    require(reader.fieldnames==writer.fieldnames,'Task schema drift'); writer.writerows(reader)
    shutil.copy2(repo/'analysis_v2/DA2_DOSE_RESPONSE.md',stage/'README.md')
    write_table(stage/'task_evidence.tsv',[dict(task=p.name,path=str(p/'SHA256SUMS'),sha256=digest(p/'SHA256SUMS')) for p in full])
    for name in ('source_commit.txt','image.sha256'):
        if (root/name).is_file(): shutil.copy2(root/name,stage/name)
    subprocess.run(['Rscript',str(repo/'analysis_v2/scripts/plot_da2_response.R'),str(stage)],check=True)
    (stage/'status.json').write_text(json.dumps(dict(status='PASS_DA2_DOSE_RESPONSE',tasks=len(tasks),
        contexts=len({r['context_id'] for r in summaries}),summary_cells=len(summaries),
        failures=[],scientific_interpretation='PENDING_REVIEW'),indent=2)+'\n')
    seal(stage,'PASS_DA2_DOSE_RESPONSE'); os.rename(stage,out)
    print((out/'status.json').read_text())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('action',choices=('prepare','run','collect'))
    p.add_argument('--root',type=Path,required=True); p.add_argument('--repo',type=Path,required=True)
    for name in ('inventory','bracken','metaphlan'): p.add_argument('--'+name,type=Path)
    p.add_argument('--index',type=int)
    a=p.parse_args()
    if a.action=='prepare': prepare(a.root,a.inventory,a.bracken,a.metaphlan,a.repo)
    elif a.action=='run': run(a.root,a.index,a.repo)
    else: collect(a.root,a.repo)
