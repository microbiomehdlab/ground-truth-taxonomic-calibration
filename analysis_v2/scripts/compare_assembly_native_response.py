#!/usr/bin/env python3
"""Matched assembly-choice response, without legacy recovery or new DA fits."""
from __future__ import annotations
import argparse
import json
import math
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from statistics import median
from audit_bracken_denominators import table, require, digest, write_table
from build_biomarker_abundance_input import parse_profile
from reference_response_pilot import fresh, seal


def pair_rows(rows, designs):
    baselines={};cells=defaultdict(dict)
    for r in rows:
        require(r['include']=='1' and r['cohort']=='yachida' and
                r['analysis_population']=='independent', 'Wrong/excluded assembly panel')
        require(r['target_label'] in ('Pana','Pint') and r['assembly_arm'] in ('original','clean'), 'Wrong assembly target/arm')
        value=float(r['abundance_fraction'])
        require(math.isfinite(value) and 0<=value<=1.00001,'Invalid native abundance')
        key=(r['sample_id'],r['target_label'],r['profiler'])
        if Decimal(r['spike_fraction_target'])==0:
            # Historical canonical builder repeats the same baseline per arm.
            old=baselines.get(key)
            require(old is None or (old['source_profile']==r['source_profile'] and
                    old['abundance_fraction']==r['abundance_fraction']),'Conflicting baseline')
            baselines[key]=r;continue
        d=[d for d in designs[r['source_design']] if d['sample_id']==r['sample_id'] and
           int(d['N_inserted'])==int(r['implanted_read_pairs_target']) and
           abs(float(d['f_hat'])-float(r['spike_fraction_target']))<1e-12]
        require(len(d)==1,'Ambiguous matched assembly design')
        d=d[0]; cell=key+(str(Decimal(d['fraction'])),)
        require(r['assembly_arm'] not in cells[cell], 'Duplicate assembly cell')
        cells[cell][r['assembly_arm']]=(r,d)
    result=[]
    for key,arms in sorted(cells.items()):
        require(set(arms)=={'original','clean'} and key[:3] in baselines,'Missing paired arm/baseline')
        (o,od),(c,cd)=arms['original'],arms['clean']
        require(o['condition']==c['condition'], 'Condition differs')
        for field in ('R','R1','R2','N_inserted','seed'):
            require(field in od and field in cd and od[field]==cd[field], 'Unmatched design: '+field)
        require(abs(float(o['spike_fraction_target'])-float(c['spike_fraction_target']))<1e-12,'Achieved dose differs')
        b=float(baselines[key[:3]]['abundance_fraction']);ov=float(o['abundance_fraction']);cv=float(c['abundance_fraction'])
        result.append(dict(sample_id=key[0],target_label=key[1],profiler=key[2],nominal_target_fraction=key[3],
            condition=o['condition'],achieved_target_fraction=o['spike_fraction_target'],baseline_native=b,
            original_native=ov,replacement_native=cv,original_native_change=ov-b,replacement_native_change=cv-b,
            replacement_minus_original_native=cv-ov,original_nonzero=int(ov>0),replacement_nonzero=int(cv>0),
            original_positive_change=int(ov>b),replacement_positive_change=int(cv>b),
            comparison='MATCHED_NATIVE_RESPONSE_NOT_GENOME_CORRECTED_RECOVERY_OR_CONTAMINATION_CAUSAL_EFFECT'))
    require(result,'No assembly pairs')
    return result


def build(canonical,aliases,out):
    fresh(out,[canonical.parent,aliases]); rows=table(canonical)
    require((canonical.parent/'validation/SUCCESS').is_file(), 'Canonical validation absent')
    hashes={str(canonical):digest(canonical),str(aliases):digest(aliases)}
    with aliases.open(newline='') as h:
        import csv
        alias={(r['canonical'],r['tool']):r['alias'] for r in csv.DictReader(h)}
    designs={};profiles={}
    for r in rows:
        path=Path(r['source_profile'])
        if str(path) not in profiles:
            hashes[str(path)]=digest(path);profiles[str(path)]=parse_profile(path,r['profiler'])
        feature=alias[(r['target_taxon'],r['profiler'])]
        require(abs(profiles[str(path)].get(feature,0)-float(r['abundance_fraction']))<1e-12,'Canonical native value differs')
        if r['source_design']!='BASELINE' and r['source_design'] not in designs:
            path=Path(r['source_design']);hashes[str(path)]=digest(path);designs[str(path)]=table(path)
    pairs=pair_rows(rows,designs)
    require(len(pairs)==720 and len({r['sample_id'] for r in pairs})==30,'Expected 720 matched pairs in 30 people')
    groups=defaultdict(list)
    for r in pairs:
        for scope in ('all_conditions',r['condition']):
            groups[(r['target_label'],r['profiler'],r['nominal_target_fraction'],scope)].append(r)
    summary=[]
    for key,values in sorted(groups.items()):
        n=len(values)
        require(n==(30 if key[3]=='all_conditions' else 10),'Incomplete assembly dose/condition cell')
        summary.append(dict(zip(('target_label','profiler','nominal_target_fraction','condition_scope'),key),
            people=n,median_replacement_minus_original_native=median(r['replacement_minus_original_native'] for r in values),
            min_replacement_minus_original_native=min(r['replacement_minus_original_native'] for r in values),
            max_replacement_minus_original_native=max(r['replacement_minus_original_native'] for r in values),
            original_positive_change_frequency=sum(r['original_positive_change'] for r in values)/n,
            replacement_positive_change_frequency=sum(r['replacement_positive_change'] for r in values)/n))
    require(all(digest(Path(p))==h for p,h in hashes.items()),'Assembly sources changed')
    out.mkdir(parents=True);write_table(out/'person_assembly_pairs.tsv',pairs)
    write_table(out/'assembly_response_summary.tsv',summary)
    write_table(out/'input_hashes.tsv',[dict(path=p,sha256=h) for p,h in sorted(hashes.items())])
    (out/'status.json').write_text(json.dumps(dict(status='PASS_MATCHED_NATIVE_RESPONSE_PENDING_SOURCE_SEAL_AND_SCIENTIFIC_REVIEW',
        pairs=len(pairs),people=30,new_model_fits=0,production_authorized=False,
        scope='Native measurement comparison only; no genome-corrected recovery ratios'),indent=2)+'\n')
    seal(out,'PASS_MATCHED_NATIVE_ASSEMBLY_RESPONSE_PENDING_REVIEW')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--canonical',type=Path,required=True);p.add_argument('--aliases',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args();build(a.canonical.resolve(),a.aliases.resolve(),a.out.resolve())
