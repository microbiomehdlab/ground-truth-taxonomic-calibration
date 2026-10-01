#!/usr/bin/env python3
"""Receipt-verified per-cohort input inventory and DA2/DA3 baseline families.

DA1 eligibility/families remain a separate covariate-review gate. No models.
"""
from __future__ import annotations
import argparse
import csv
import math
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from audit_bracken_denominators import table, digest, require, verify_file, write_table
from build_biomarker_abundance_input import parse_profile
from da_allocations import allocations


def build(canonical,state,out,cohort,expected,repetitions=1000):
    canonical,state,out=[p.resolve() for p in (canonical,state,out)]
    require(not out.exists(), 'Output exists')
    require(out!=state and state not in out.parents and out not in state.parents and out!=canonical,
            'Output overlaps inputs')
    seal=state/'production_seal_v2'
    require(not (seal/'AUDIT_IN_PROGRESS').exists(),'Audit in progress')
    covered=set()
    for line in (seal/'production_seal.sha256').read_text().splitlines():
        sha,name=line.split(maxsplit=1)
        name=name.lstrip('*')
        require(Path(name).name==name and name not in covered,'Unsafe seal member')
        require(digest(seal/name)==sha,'Seal checksum mismatch')
        covered.add(name)
    require({'SUCCESS','production_manifest.tsv','sample_flow.tsv'}<=covered,'Missing seal coverage')
    success=dict(line.split('\t',1) for line in (seal/'SUCCESS').read_text().splitlines())
    require(success.get('status')=='PASS' and success.get('cohort')==cohort,'Wrong/failed seal')
    manifest=table(seal/'production_manifest.tsv')
    samples={r['sample_id']:r for r in manifest}
    require(len(samples)==len(manifest)==expected==int(success['samples']),'Sample count mismatch')
    flow=table(seal/'sample_flow.tsv')
    require(len(flow)==expected and {r['sample_id'] for r in flow}==set(samples) and
            all(r['status']=='PASS' for r in flow),'Invalid sample flow')
    rows=table(canonical)
    profiles={}; targets=defaultdict(set); coverage=defaultdict(set); seen=set()
    receipts={}; hashes={canonical:digest(canonical)}
    for name in covered: hashes[seal/name]=digest(seal/name)
    for sid in samples:
        require(sid not in ('','.','..') and Path(sid).name==sid,'Unsafe sample')
        receipt_path=state/'samples'/(sid+'.retained_outputs.tsv')
        hashes[receipt_path]=digest(receipt_path)
        rr={}
        for item in table(receipt_path):
            path=Path(item['path'])
            require(path.is_absolute(),'Nonabsolute receipt path')
            path=path.resolve()
            require(path not in rr,'Duplicate receipt path')
            rr[path]=item
        receipts[sid]=rr
    for row in rows:
        sid=row['sample_id']; profiler=row['profiler']; pop=row['analysis_population']
        require(row['cohort']==cohort and sid in samples and row['include']=='1','Wrong/excluded canonical input')
        require(row['condition']==samples[sid]['condition'],'Clinical condition mismatch')
        require(profiler in ('metaphlan4','kraken2_bracken') and pop in ('community','independent'),'Invalid context')
        path=Path(row['source_profile']).resolve()
        if path not in hashes:
            verify_file(path,receipts[sid]); hashes[path]=digest(path)
        else: require(path in receipts[sid],'Profile absent from sample receipt')
        f=Decimal(row['spike_fraction_total'])
        require(f.is_finite() and 0<=f<1,'Invalid fraction')
        nominal=Decimal('0') if f==0 else Decimal(row['profile_id'].rsplit('_f',1)[1].replace('p','.'))
        key=(sid,pop,profiler,row['profile_id'],row['target_label'])
        require(key not in seen,'Duplicate canonical endpoint'); seen.add(key)
        targets[profiler].add(row['target_taxon'])
        pk=(sid,pop,profiler,row['profile_id'])
        rec=dict(cohort=cohort,sample_id=sid,condition=row['condition'],analysis_population=pop,
                 profiler=profiler,profile_id=row['profile_id'],nominal_total_dose=str(nominal),
                 source_profile=str(path),sha256=hashes[path])
        require(pk not in profiles or profiles[pk]==rec,'Conflicting physical profile')
        profiles[pk]=rec
        coverage[(sid,pop,profiler,row['target_label'])].add(nominal)
    require({r['sample_id'] for r in rows}==set(samples),'Canonical sample coverage mismatch')
    for sid in samples:
        for profiler in ('metaphlan4','kraken2_bracken'):
            for pop in ('community','independent'):
                labels={k[3] for k in coverage if k[:3]==(sid,pop,profiler)}
                needed=10 if pop=='community' or samples[sid]['independent_subset']=='1' else 0
                require(len(labels)==needed,'Missing/unexpected target series: '+sid)
    for (sid,pop,profiler,target),doses in coverage.items():
        required={Decimal(x) for x in (('0','.0001','.0005','.001','.005','.01','.05','.1') if pop=='community'
                                      else ('0','.0001','.0005','.001','.005','.01','.05'))}
        require(doses==required,'Incomplete dose series')
    families=[]; eligibility=[]; ledger=[]
    # Use frozen exact target aliases, not name normalization or inferred synonyms.
    # Repository location is based on this script, not a cluster work layout.
    aliases=Path(__file__).resolve().parents[2]/'examples/spike_taxon_aliases.csv'
    hashes[aliases]=digest(aliases)
    with aliases.open(newline='') as h: alias_rows=list(csv.DictReader(h))
    for condition in ('Control','Adenoma','CRC'):
        pool=sorted(s for s in samples if samples[s]['condition']==condition)
        require(pool,'Empty clinical pool')
        subset=sorted(s for s in pool if samples[s]['independent_subset']=='1')
        require(len(subset)==10,'Independent subset must be ten per condition')
        for profiler in ('metaphlan4','kraken2_bracken'):
            base=[r for r in profiles.values() if r['condition']==condition and r['profiler']==profiler and
                  r['analysis_population']=='community' and Decimal(r['nominal_total_dose'])==0]
            require(len(base)==len(pool),'Missing baseline')
            counts=defaultdict(int)
            for r in base:
                values=parse_profile(Path(r['source_profile']),profiler)
                require(values and all(math.isfinite(v) for v in values.values()),'Invalid/empty species profile')
                for feature,v in values.items(): counts[feature]+=int(v>0)
            target_names={r['alias'] for r in alias_rows if r['tool']==profiler and r['canonical'] in targets[profiler]}
            require(len(target_names)==10,'Target alias coverage mismatch')
            for feature in sorted(set(counts)|target_names):
                if counts[feature]/len(pool)>=.1 or feature in target_names:
                    families.append(dict(cohort=cohort,condition=condition,profiler=profiler,feature=feature,
                        baseline_positive=counts[feature],baseline_n=len(pool),target=int(feature in target_names)))
        for n in (5,10,15,20):
            eligibility.append(dict(cohort=cohort,condition=condition,pool_n=len(pool),n_per_group=n,
                                    feasible=int(2*n<=len(pool)),independent_n=len(subset)))
        if condition=='Control': continue
        for exhaustive,p in ((False,pool),(True,subset)):
            for draw in allocations(p,cohort,condition,repetitions=repetitions,exhaustive=exhaustive):
                for group in ('cases','controls'):
                    for position,sid in enumerate(draw[group],1):
                        ledger.append(dict(cohort=cohort,condition=condition,pool='independent_subset' if exhaustive else 'full',
                            allocation_id=draw['allocation_id'],n=draw['n'],group=group,sample_id=sid,
                            position=position,pool_hash=draw['pool_hash'],seed_namespace=draw['seed_namespace']))
    require(all(digest(p)==sha for p,sha in hashes.items()),'Inputs changed')
    out.mkdir(parents=True)
    for name,data in [('profile_inventory.tsv',list(profiles.values())),('feature_families.tsv',families),
                      ('eligibility.tsv',eligibility),('allocations.tsv',ledger),
                      ('input_hashes.tsv',[dict(path=str(p),sha256=sha) for p,sha in sorted(hashes.items())])]:
        write_table(out/name,data)
    (out/'SUCCESS').write_text('status\tPASS_DESIGN_INVENTORY\nDA1_families\tPENDING_COVARIATE_REVIEW\n')
    (out/'SHA256SUMS').write_text(''.join(digest(p)+'  '+p.name+'\n' for p in sorted(out.iterdir())))
    print('[PASS] Design inventory; no models:',cohort,len(profiles))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ('canonical','state-dir','outdir'): p.add_argument('--'+flag,type=Path,required=True)
    p.add_argument('--cohort',choices=['yachida','feng','zeller'],required=True)
    p.add_argument('--expected-samples',type=int,required=True)
    p.add_argument('--repetitions',type=int,default=1000)
    a=p.parse_args()
    try: build(a.canonical,a.state_dir,a.outdir,a.cohort,a.expected_samples,a.repetitions)
    except (ValueError,OSError,KeyError,IndexError) as e: p.exit(1,'[ERROR] '+str(e)+'\n')
