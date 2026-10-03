#!/usr/bin/env python3
"""Join immutable DA3 allocations/results to exact-person quantitative recovery."""
from __future__ import annotations
import argparse, csv, gzip, json, math
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from statistics import mean, median, variance
from audit_bracken_denominators import digest, require, table, write_table
from summarize_bracken_recovery import quantile
from plan_da3_canary import finish

COHORTS = ('yachida', 'feng', 'zeller')

def number(value):
    v=float(value); require(math.isfinite(v), 'Nonfinite numeric evidence'); return v

def dose(value):
    d=Decimal(str(value)); require(d.is_finite() and d>=0, 'Invalid dose'); return str(d.normalize())

def close(a,b):
    return math.isclose(float(a),float(b),rel_tol=1e-7,abs_tol=1e-10)

def manifest(root, required, evidence):
    root=root.resolve(); entries={}
    ledger=root/'SHA256SUMS'
    for line in ledger.read_text().splitlines():
        sha,name=line.split(maxsplit=1); p=(root/name).resolve()
        require(root in p.parents and p not in entries, 'Unsafe/duplicate checksum member')
        require(digest(p)==sha, 'Checksum mismatch: '+str(p));entries[p]=sha
    require(all((root/n).resolve() in entries for n in required), 'Unchecksummed required input')
    evidence.update(entries);evidence[ledger]=digest(ledger)
    return entries

def normalize(row, tool):
    if tool=='kraken2_bracken':
        fields=('legacy_baseline_native','legacy_observed_native','total_fraction_exact',
                'implanted_fraction_target_exact','recovered_signal_all_input',
                'recovery_ratio_all_input','signed_error_all_input')
        nominal=row['nominal_total_dose'];signal=row['implanted_fraction_target_exact']
        reference='all_input_pairs';base_path=row['source_baseline']
    else:
        require(row['reference_type']=='genome_equivalent','Wrong MetaPhlAn recovery scale')
        fields=('baseline_abundance_fraction','observed_abundance_fraction','spike_fraction_total',
                'spike_fraction_target','recovered_spike_signal_profiler_scale',
                'response_ratio_profiler_scale','response_residual_profiler_scale')
        nominal=row['profile_id'].rsplit('_f',1)[1].replace('p','.')
        signal=row['implanted_signal_profiler_scale']
        reference='genome_equivalent';base_path=row['source_baseline_profile']
    b,o,total,target,recovered,ratio,residual=[number(row[f]) for f in fields]
    signal=number(signal)
    require(0<=b<=1.00001 and 0<=o<=1.00001 and 0<target<=total<1 and signal>0, 'Invalid endpoint units')
    require(close(recovered/signal,ratio) and close(recovered-signal,residual), 'Inconsistent recovery arithmetic')
    require(Decimal(dose(nominal))>0,'Endpoint dose must be positive')
    return dict(cohort=row['cohort'],profiler=tool,sample_id=row['sample_id'],
                condition=row['condition'],population='community',target_label=row['target_label'],
                nominal_total_dose=dose(nominal),baseline_native=b,spiked_native=o,native_increment=o-b,
                achieved_total_read_fraction=total,achieved_target_read_fraction=target,
                reference_type=reference,expected_added_signal=signal,recovered_added_signal=recovered,
                recovery_ratio=ratio,response_residual=residual,source_profile=row['source_profile'],
                source_baseline=base_path,profile_id=row['profile_id'],
                source_target_alias=row.get('target_alias',''))

def load_endpoints(bracken, metaphlan, evidence):
    points={}
    manifest(metaphlan,['SUCCESS']+[c+'/endpoints/paired_endpoints.tsv' for c in COHORTS],evidence)
    require('status\tPASS_ENDPOINTS' in (metaphlan/'SUCCESS').read_text().splitlines(),'Wrong MetaPhlAn bundle')
    for c in COHORTS:
        b=bracken/c
        manifest(b,['summary.json','bracken_recovery_comparison.tsv'],evidence)
        s=json.loads((b/'summary.json').read_text())
        require(s['cohort']==c and s['status']=='PASS_ENDPOINT_CONSTRUCTION','Wrong Bracken bundle')
        sources=[('kraken2_bracken',table(b/'bracken_recovery_comparison.tsv')),
                 ('metaphlan4',[r for r in table(metaphlan/c/'endpoints/paired_endpoints.tsv') if r['profiler']=='metaphlan4'])]
        for tool, rows in sources:
            for row in rows:
                require(row['cohort']==c,'Endpoint cohort mismatch')
                if row['analysis_population']!='community':continue
                point=normalize(row,tool)
                key=(c,tool,row['sample_id'],point['nominal_total_dose'],row['target_label'])
                require(key not in points,'Duplicate person/dose/target endpoint')
                points[key]=point
    require(points,'No community recovery endpoints')
    return points

def validate_points(points,catalog):
    baseline={}
    for key,r in points.items():
        c,t,s,d,label=key
        require(label in catalog['targets'][t].values(),'Unrecognized endpoint target')
        if t=='kraken2_bracken':
            require(catalog['targets'][t].get(r['source_target_alias'])==label,'Bracken endpoint alias mismatch')
        # The catalog normalizes numeric doses with Decimal, which preserves scale.
        pk=(c,t,s,Decimal(d));bk=(c,t,s,Decimal(0))
        profiles=catalog['_profiles_by_identity']
        require(pk in profiles and bk in profiles,'Endpoint missing from sealed profile catalog')
        require(Path(r['source_profile'])==Path(profiles[pk]['source_profile']) and
                Path(r['source_baseline'])==Path(profiles[bk]['source_profile']),'Physical profile identity mismatch')
        basekey=(c,t,s,label)
        if basekey in baseline:
            require(close(baseline[basekey]['baseline_native'],r['baseline_native']) and
                    baseline[basekey]['condition']==r['condition'],'Conflicting repeated baseline')
        else:baseline[basekey]=r
    return baseline

def summary_values(prefix, values):
    if not values:
        return {prefix+k:'' for k in ('_median','_q1','_q3','_min','_max')}
    return {prefix+'_median':median(values),prefix+'_q1':quantile(values,.25),
            prefix+'_q3':quantile(values,.75),prefix+'_min':min(values),prefix+'_max':max(values)}

def link_context(context, fits, family_n, catalog, points, baseline):
    c,t=context['cohort'],context['profiler'];obs=context['observations'];n=int(context['n'])
    require(len(obs)==2*n and len({o['sample_id'] for o in obs})==len(obs),'Repeated/missing biological person')
    require(all(o['group'] in (0,1) for o in obs) and sum(o['group'] for o in obs)==n,'Invalid group sizes')
    require(all(Decimal(str(o['dose']))==0 for o in obs if not o['group']),'Spiked artificial control')
    for o in obs:
        require(o['key'] in catalog['profiles'],'Missing selected profile')
        require(tuple(o['key'].split('|')[:3])==(c,t,o['sample_id']) and
                Decimal(o['key'].split('|')[3])==Decimal(str(o['dose'])),'Plan observation identity mismatch')
    targets=catalog['targets'][t]
    require(len(fits)==len(targets) and len({f['target_label'] for f in fits})==len(fits) and
            {f['target_label'] for f in fits}==set(targets.values()),'Incomplete/duplicate target statistics')
    result=[]
    for f in fits:
        label=f['target_label'];require(targets.get(f['feature'])==label,'Target alias mismatch')
        cases=[];controls=[];casebase=[];exposed=[]
        for o in obs:
            sid=o['sample_id'];base=baseline[(c,t,sid,label)]
            require(base['condition']==context['background'],'Wrong clinical background')
            b=base['baseline_native'];d=dose(o['dose'])
            if Decimal(d)>0:
                r=points[(c,t,sid,d,label)];exposed.append(r);v=r['spiked_native']
            else:v=b
            if o['group']:cases.append(v);casebase.append(b)
            else:controls.append(v)
        transformed_cases=[math.log2(1+v/1e-8) for v in cases]
        transformed_controls=[math.log2(1+v/1e-8) for v in controls]
        beta=mean(transformed_cases)-mean(transformed_controls)
        require(close(beta,number(f['beta'])),'Reconstructed native effect differs from saved DA3 result')
        q=number(f['parametric_bh_q']);fw=number(f['max_statistic_fwer_p'])
        require(0<=q<=1 and 0<=fw<=1,'Invalid adjusted result')
        p=None if f['parametric_p'] in ('','NA') else number(f['parametric_p'])
        require(p is None or 0<=p<=1,'Invalid raw p')
        require(p is not None or q==1,'Nonestimable target not retained conservatively')
        se=math.sqrt(((n-1)*variance(transformed_cases)+(n-1)*variance(transformed_controls))/(2*n-2)*(2/n))
        r=dict(context_id=context['context_id'],target_label=label,target_feature=f['feature'],
               cohort=c,profiler=t,background=context['background'],n_per_group=n,
               arm=context['arm'],anchor=context['anchor'],allocation_id=context['allocation_id'],
               family_n=family_n,exposed_cases=len(exposed),unexposed_cases=n-len(exposed),
               achieved_exposure_fraction=len(exposed)/n,requested_exposure=context['requested_exposure'],
               beta=beta,parametric_se=se if p is not None else '',parametric_raw_p='' if p is None else p,
               parametric_bh_q=q,parametric_estimable=int(p is not None),
               positive_parametric_discovery=int(p is not None and beta>0 and q<=.05),
               raw_positive_without_bh=int(p is not None and beta>0 and p<=.05 and q>.05),
               permutation_fwer_p=fw,positive_permutation_discovery=int(beta>0 and fw<=.05),
               reference_type=baseline[(c,t,obs[0]['sample_id'],label)]['reference_type'],
               case_baseline_prevalence=mean(v>0 for v in casebase),
               control_baseline_prevalence=mean(v>0 for v in controls),
               observed_case_prevalence=mean(v>0 for v in cases),
               case_baseline_sd=math.sqrt(variance(casebase)),control_baseline_sd=math.sqrt(variance(controls)),
               exposed_negative_recovery=sum(v['recovered_added_signal']<0 for v in exposed),
               exposed_zero_recovery=sum(v['recovered_added_signal']==0 for v in exposed),
               exposed_positive_native_change=sum(v['native_increment']>0 for v in exposed),
               exposed_negative_native_change=sum(v['native_increment']<0 for v in exposed),
               exposed_zero_native_change=sum(v['native_increment']==0 for v in exposed),
               recovery_scope='EXACT_EXPOSED_PEOPLE_AT_ASSIGNED_DOSES',
               interpretation='DESCRIPTIVE_LINK_NOT_CAUSAL_ATTRIBUTION')
        for prefix,values in [('case_baseline',casebase),('control_baseline',controls),('observed_cases',cases)]:
            r.update(summary_values(prefix,values))
        for field in ('recovery_ratio','native_increment','response_residual','expected_added_signal',
                      'recovered_added_signal','achieved_target_read_fraction','achieved_total_read_fraction'):
            r.update(summary_values('exposed_'+field,[v[field] for v in exposed]))
        result.append(r)
    return result

def build(plan, results, bracken, metaphlan, out, fixture=False):
    plan,results,bracken,metaphlan,out=[p.resolve() for p in (plan,results,bracken,metaphlan,out)]
    require(not out.exists() and all(out!=p and out not in p.parents and p not in out.parents
                                    for p in (plan,results,bracken,metaphlan)),'Fresh nonoverlapping output required')
    evidence={};manifest(plan,['SUCCESS','tasks.tsv','catalog.json'],evidence)
    code_root=Path(__file__).resolve().parent
    for name in ('link_da3_recovery.py','audit_bracken_denominators.py','summarize_bracken_recovery.py','plan_da3_canary.py'):
        p=code_root/name;evidence[p]=digest(p)
    tasks=table(plan/'tasks.tsv')
    require(len({r['index'] for r in tasks})==len(tasks),'Duplicate batch task')
    if not fixture:
        require(len(tasks)==4800 and [int(r['index']) for r in tasks]==list(range(4800)) and
                all(int(r['contexts'])==25 for r in tasks),'Expected complete 4800-batch DA3 plan')
        require('status\tPASS_DRAFT_BATCH_PLAN_NOT_PRODUCTION_AUTHORIZATION' in (plan/'SUCCESS').read_text().splitlines(),'Wrong DA3 plan status')
    catalog=json.loads((plan/'catalog.json').read_text());catalog['_profiles_by_identity']={}
    for k,v in catalog['profiles'].items():
        c,t,s,d=k.split('|');key=(c,t,s,Decimal(d))
        require(key not in catalog['_profiles_by_identity'],'Duplicate catalog profile identity')
        catalog['_profiles_by_identity'][key]=v
    if not fixture:require(all(len(v)==10 and len(set(v.values()))==10 for v in catalog['targets'].values()),'Expected ten unique targets')
    points=load_endpoints(bracken,metaphlan,evidence);baseline=validate_points(points,catalog)
    # Check every batch before producing an apparently complete explanation.
    for task in tasks:
        f=results/('batch_%05d'%int(task['index']))
        manifest(f,['SUCCESS','resume_identity.json','summary.tsv','targets.tsv.gz'],evidence)
        require(json.loads((f/'resume_identity.json').read_text())['plan_sha256']==digest(plan/'SHA256SUMS'),'Results belong to another plan')
        require('status\tPASS_BATCH_COMPUTATION_NOT_PRODUCTION_AUTHORIZATION' in (f/'SUCCESS').read_text().splitlines(),'Wrong batch status')
        require((plan/task['file']).resolve() in evidence,'Unchecksummed allocation file')
    out.mkdir(parents=True)
    with gzip.open(out/'person_target_response.tsv.gz','wt',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(next(iter(points.values()))),delimiter='\t');w.writeheader()
        w.writerows(points[k] for k in sorted(points))
    grouped=defaultdict(lambda:[0,0,0,0]);count=0;seen=set()
    with gzip.open(out/'context_target_explanation.tsv.gz','wt',newline='') as h:
        writer=None
        for task in tasks:
            folder=results/('batch_%05d'%int(task['index']))
            contexts=json.loads((plan/task['file']).read_text())
            require(len(contexts)==int(task['contexts']),'Plan batch count mismatch')
            summaries=table(folder/'summary.tsv')
            require(len(summaries)==len(contexts) and len({s['context_id'] for s in summaries})==len(summaries),'Duplicate/missing batch summaries')
            summaries={s['context_id']:s for s in summaries}
            fits=defaultdict(list)
            with gzip.open(folder/'targets.tsv.gz','rt',newline='') as source:
                for r in csv.DictReader(source,delimiter='\t'):fits[r['context_id']].append(r)
            ids={c['context_id'] for c in contexts}
            require(len(ids)==len(contexts) and set(fits)==ids and set(summaries)==ids and not seen&ids,'Context coverage mismatch')
            seen.update(ids)
            for context in contexts:
                cid=context['context_id'];family=catalog['families']['|'.join((context['cohort'],context['profiler'],context['background']))]
                require(len(family)==len(set(family))==int(summaries[cid]['family_n']),'Feature-family mismatch')
                require(set(catalog['targets'][context['profiler']])<=set(family),'Targets absent from full family')
                rows=link_context(context,fits[cid],len(family),catalog,points,baseline)
                if writer is None:writer=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader()
                writer.writerows(rows);count+=len(rows)
                for row in rows:
                    key=tuple(row[k] for k in ('cohort','profiler','background','target_label','n_per_group','arm','anchor'))
                    v=grouped[key];v[0]+=1;v[1]+=row['positive_parametric_discovery'];v[2]+=row['positive_permutation_discovery'];v[3]+=row['raw_positive_without_bh']
            if int(task['index'])%200==0:print('[INFO] Linked',len(seen),'contexts',flush=True)
    if not fixture:
        require(len(seen)==120000 and count==1200000,'Incomplete final experiment grid')
        require(len(grouped)==12000 and all(v[0]==100 for v in grouped.values()),'Incomplete conditional frequency cells')
    fields=('cohort','profiler','background','target_label','n_per_group','arm','anchor')
    write_table(out/'conditional_significance_frequencies.tsv',[
        dict(zip(fields,k),allocations=v[0],standard_positive_frequency=v[1]/v[0],
             permutation_positive_frequency=v[2]/v[0],raw_positive_without_bh_frequency=v[3]/v[0],
             interpretation='CONDITIONAL_ON_FINITE_COHORT_ALLOCATIONS') for k,v in sorted(grouped.items())])
    require(all(digest(p)==sha for p,sha in evidence.items()),'Inputs changed while linking')
    write_table(out/'input_hashes.tsv',[dict(path=str(p),sha256=sha) for p,sha in sorted(evidence.items())])
    write_table(out/'code_hashes.tsv',[dict(path=str(Path(__file__).resolve()),sha256=digest(Path(__file__)))])
    (out/'status.json').write_text(json.dumps(dict(status='PASS_DESCRIPTIVE_LINK',contexts=len(seen),
        context_target_rows=count,person_target_rows=len(points),production_authorized=False),indent=2)+'\n')
    (out/'README.md').write_text('Exact exposed-person recovery joined to saved DA3 allocations and statistics.\n'
        'Recovery summaries exclude unexposed cases; null arms have blank recovery summaries.\n'
        'Native DA fractions and profiler-specific quantitative recovery are separate scales.\n'
        'Frequencies are conditional on reused cohort members, not population-level probabilities.\n'
        'These associations do not isolate a causal fraction attributable to profiling.\n'
        'Person rows are unique community sample/dose/target responses; context rows identify the immutable allocation.\n')
    finish(out,'PASS_DESCRIPTIVE_LINK')
    print('[PASS] Exact-person recovery/significance link:',count,'target-context rows')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('plan','results','bracken','metaphlan','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();build(a.plan,a.results,a.bracken,a.metaphlan,a.out)
