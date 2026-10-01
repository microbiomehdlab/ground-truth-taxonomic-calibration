#!/usr/bin/env python3
"""Prepare a fixed small engineering pilot, not definitive DA results."""
from __future__ import annotations
import argparse
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from audit_bracken_denominators import table, digest, require, write_table


def verify(folder):
    covered=set()
    for line in (folder/'SHA256SUMS').read_text().splitlines():
        sha,name=line.split(maxsplit=1)
        path=Path(name)
        require(not path.is_absolute() and '..' not in path.parts and name not in covered,'Unsafe checksum member')
        require(digest(folder/path)==sha,'Checksum mismatch: '+str(folder/path))
        covered.add(name)
    require('SUCCESS' in covered,'SUCCESS not checksummed')


def build(inventory_root, da1_root, out, cohorts=('yachida','feng','zeller')):
    inventory_root,da1_root,out=[p.resolve() for p in (inventory_root,da1_root,out)]
    require(not out.exists(),'Output exists')
    for root in (inventory_root,da1_root):
        require(out!=root and out not in root.parents and root not in out.parents,'Overlapping output')
        verify(root)
    contexts=[]; observations=[]; family_rows=[]; achieved=[]; matches=[]; hashes={}
    profilers=('kraken2_bracken','metaphlan4')
    caches={}
    for cohort in cohorts:
        inv=inventory_root/cohort; da1=da1_root/cohort
        for folder,status in ((inv,'PASS_DESIGN_INVENTORY'),(da1,'PASS_DA1_CHECKPOINT')):
            verify(folder)
            require('status\t'+status in (folder/'SUCCESS').read_text().splitlines(),'Wrong checkpoint status')
            for p in folder.iterdir():
                if p.is_file(): hashes[p]=digest(p)
            for r in table(folder/'input_hashes.tsv'):
                p=Path(r['path']).resolve()
                require(p not in hashes or hashes[p]==r['sha256'],'Conflicting source hashes')
                hashes[p]=r['sha256']
        profiles=table(inv/'profile_inventory.tsv')
        index={}
        for r in profiles:
            require(r['cohort']==cohort,'Cohort mismatch')
            key=(r['sample_id'],r['profiler'],r['analysis_population'],Decimal(r['nominal_total_dose']),r['profile_id'])
            require(key not in index,'Duplicate profile'); index[key]=r
        canon=[p for p in hashes if p.name=='canonical_input.tsv' and p in
               {Path(r['path']).resolve() for r in table(inv/'input_hashes.tsv')}]
        require(len(canon)==1,'Ambiguous canonical input')
        endpoints=table(canon[0]); labels={r['profile_id']:r['target_label'] for r in endpoints if r['analysis_population']=='independent'}
        for r in endpoints:
            if Decimal(r['spike_fraction_total'])==0: continue
            achieved.append(dict(cohort=cohort,sample_id=r['sample_id'],profiler=r['profiler'],
                population=r['analysis_population'],profile_id=r['profile_id'],target_label=r['target_label'],
                achieved_total_fraction=r['spike_fraction_total'],achieved_target_fraction=r['spike_fraction_target']))
        # Pair by nominal target dose only; save actual achieved values separately.
        # Actual totals have rounding: nominal dose is recovered from profile identity.
        matched=defaultdict(dict)
        for r in endpoints:
            if Decimal(r['spike_fraction_total'])==0: continue
            nominal=Decimal(r['profile_id'].rsplit('_f',1)[1].replace('p','.'))
            target_nominal=nominal/10 if r['analysis_population']=='community' else nominal
            matched[(r['sample_id'],r['profiler'],r['target_label'],target_nominal)][r['analysis_population']]=r
        for key,arms in sorted(matched.items()):
            if set(arms)!= {'community','independent'}: continue
            a,b=arms['independent'],arms['community']
            matches.append(dict(cohort=cohort,sample_id=key[0],profiler=key[1],target_label=key[2],
                nominal_target_fraction=str(key[3]),independent_profile=a['profile_id'],community_profile=b['profile_id'],
                independent_achieved=a['spike_fraction_target'],community_achieved=b['spike_fraction_target'],
                achieved_difference=str(Decimal(b['spike_fraction_target'])-Decimal(a['spike_fraction_target'])),
                match_policy='NOMINAL_ONLY_NOT_EXACT_ACHIEVED'))
        caches[cohort]=(profiles,labels,table(inv/'feature_families.tsv'),table(da1/'feature_families.tsv'),
                        table(da1/'eligibility.tsv'),table(inv/'allocations.tsv'))

    def add(cohort,profiler,analysis,background,pop,dose,assignment,paired=False,label='',allocation=''):
        profiles,labels,fam,clinical,eligible,_=caches[cohort]
        context_id='pilot_%04d'%len(contexts)
        family=[r['feature'] for r in (clinical if analysis=='DA1' else fam)
                if r['profiler']==profiler and r.get('comparison',r.get('condition'))==
                (background+'_vs_Control' if analysis=='DA1' else background)]
        require(family and len(family)==len(set(family)),'Invalid family')
        family_rows.extend(dict(context_id=context_id,feature=f) for f in family)
        for sid,group,spiked,age,sex in assignment:
            use_dose=Decimal(dose) if spiked else Decimal(0)
            chosen=[r for r in profiles if r['sample_id']==sid and r['profiler']==profiler and
                    r['analysis_population']==pop and Decimal(r['nominal_total_dose'])==use_dose and
                    (not spiked or pop=='community' or labels.get(r['profile_id'])==label)]
            require(len(chosen)==1,'Missing/ambiguous selected profile: '+str((sid,pop,use_dose,label)))
            r=chosen[0]
            observations.append(dict(context_id=context_id,observation_id=sid+('_spiked' if spiked else '_original'),
                biological_sample_id=sid,group=group,spike_state='spiked' if spiked else 'original',age=age,sex=sex,
                source_profile=r['source_profile'],sha256=r['sha256']))
        contexts.append(dict(context_id=context_id,cohort=cohort,profiler=profiler,analysis=analysis,
            background=background,population=pop,nominal_total_dose=dose,target_label=label,
            allocation_id=allocation,paired=int(paired),covariates='age,sex' if analysis=='DA1' else '',
            family_n=len(family),observations_n=len(assignment),stage='ENGINEERING_PILOT_NOT_DEFINITIVE'))

    for cohort in cohorts:
        profiles,labels,fam,clinical,eligible,alloc=caches[cohort]
        for profiler in profilers:
            for background in ('Adenoma','CRC'):
                rows=[r for r in eligible if r['comparison']==background+'_vs_Control' and r['eligible']=='1']
                add(cohort,profiler,'DA1',background,'community','0',
                    [(r['sample_id'],int(r['group']),False,r['age'],r['sex']) for r in rows])
            # Full community paired and ten-person Fnuc paired smoke contexts.
            for pop,dose,label in (('community','.001',''),('independent','.0001','Fnuc')):
                ids=sorted({r['sample_id'] for r in profiles if r['condition']=='Adenoma' and
                            r['profiler']==profiler and r['analysis_population']==pop})
                add(cohort,profiler,'DA2','Adenoma',pop,dose,
                    [(sid,g,bool(g),'','') for sid in ids for g in (0,1)],True,label)
    # Twenty allocation contexts balanced cyclically over twelve strata. Same
    # allocation within a cohort/background across profilers; no outcome selection.
    strata=[(c,p,b) for c in cohorts for b in ('Adenoma','CRC') for p in profilers]
    for i in range(20):
        cohort,profiler,background=strata[i%len(strata)]
        allocation='full_%04d'%(i//len(strata))
        rows=[r for r in caches[cohort][5] if r['pool']=='full' and r['condition']==background
              and r['allocation_id']==allocation and r['n']=='5']
        require(len(rows)==10 and len({r['sample_id'] for r in rows})==10,'Invalid allocation')
        for dose in ('0','.0001','.005'):
            add(cohort,profiler,'NULL' if dose=='0' else 'DA3',background,'community',dose,
                [(r['sample_id'],int(r['group']=='cases'),dose!='0' and r['group']=='cases','','') for r in rows],
                allocation=allocation)
    require(len(contexts)==8*len(cohorts)+60,'Unexpected pilot context count')
    require(all(digest(p)==sha for p,sha in hashes.items()),'Changed source before pilot preparation')
    require(all(out!=p and out not in p.parents for p in hashes),'Output overlaps source')
    out.mkdir(parents=True)
    for name,rows in [('contexts.tsv',contexts),('observations.tsv',observations),('families.tsv',family_rows),
                      ('achieved_doses.tsv',achieved),('nominal_matches.tsv',matches),
                      ('input_hashes.tsv',[dict(path=str(p),sha256=s) for p,s in sorted(hashes.items())])]:
        require(rows,'Empty output: '+name); write_table(out/name,rows)
    (out/'SUCCESS').write_text('status\tPASS_PILOT_PREPARATION\ncontexts\t%d\nmodels_run\t0\n'%len(contexts))
    (out/'SHA256SUMS').write_text(''.join(digest(p)+'  '+p.name+'\n' for p in sorted(out.iterdir())))
    print('[PASS] Prepared engineering contexts; no models:',len(contexts),out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('inventory-root','da1-root','outdir'): p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    try: build(a.inventory_root,a.da1_root,a.outdir)
    except (ValueError,OSError,KeyError) as e: p.exit(1,'[ERROR] '+str(e)+'\n')
