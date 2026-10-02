#!/usr/bin/env python3
"""Export native target baselines from a checksummed DA inventory; no models."""
from __future__ import annotations
import argparse
import csv
import math
from decimal import Decimal
from pathlib import Path
from statistics import median
from audit_bracken_denominators import digest,require,table,write_table
from build_biomarker_abundance_input import parse_profile
from prepare_da_pilot import verify

EXPECTED={'yachida':201,'feng':154,'zeller':156}
TOOLS=('kraken2_bracken','metaphlan4')
CONDITIONS=('Control','Adenoma','CRC')

def quantile(values,p):
    v=sorted(values);position=(len(v)-1)*p;low=int(position)
    return v[low]+(v[min(low+1,len(v)-1)]-v[low])*(position-low)

def build(inventory,out,expected=None):
    expected=EXPECTED if expected is None else expected
    inventory,out=inventory.resolve(),out.resolve()
    require(not out.exists(),'Use a fresh output directory')
    require(out!=inventory and out not in inventory.parents and inventory not in out.parents,'Output overlaps inventory')
    repo=Path(__file__).resolve().parents[2]
    panel=table(repo/'spikes/spike_panel.tsv');labels={r['taxon_name']:r['label'] for r in panel}
    with (repo/'examples/spike_taxon_aliases.csv').open(newline='') as h:aliases=list(csv.DictReader(h))
    targets={}
    for r in aliases:
        if r['canonical'] not in labels:continue
        key=(r['tool'],labels[r['canonical']])
        require(key not in targets,'Duplicate target alias')
        targets[key]=r['alias']
    require(len(targets)==20 and len(set(labels.values()))==10,'Ten exact targets per profiler required')
    hashes={p:digest(p) for p in (repo/'spikes/spike_panel.tsv',repo/'examples/spike_taxon_aliases.csv',Path(__file__).resolve())}
    long=[];summary=[];groups={}
    for cohort,n in expected.items():
        folder=inventory/cohort;verify(folder)
        require('status\tPASS_DESIGN_INVENTORY' in (folder/'SUCCESS').read_text().splitlines(),'Wrong inventory status')
        for p in folder.iterdir():
            if p.is_file():hashes[p]=digest(p)
        profiles=[r for r in table(folder/'profile_inventory.tsv') if r['analysis_population']=='community' and Decimal(r['nominal_total_dose'])==0]
        identities=set();conditions={};source_paths=set()
        for r in profiles:
            require(r['cohort']==cohort and r['profiler'] in TOOLS and r['condition'] in CONDITIONS,'Wrong baseline context')
            sid=r['sample_id'];tool=r['profiler'];key=(sid,tool)
            require(key not in identities,'Duplicate baseline identity');identities.add(key)
            require(sid not in conditions or conditions[sid]==r['condition'],'Conflicting sample condition');conditions[sid]=r['condition']
            path=Path(r['source_profile']).resolve()
            require(path not in source_paths,'Reused physical baseline');source_paths.add(path)
            require(path.is_file() and digest(path)==r['sha256'],'Baseline profile checksum mismatch: '+str(path))
            require(out!=path and out not in path.parents,'Output contains source profile')
            hashes[path]=r['sha256'];values=parse_profile(path,tool)
            require(values and all(math.isfinite(v) and 0<=v<=1.00001 for v in values.values()),'Invalid native species profile')
            for label in sorted(labels.values()):
                feature=targets[(tool,label)];value=values.get(feature,0.0)
                long.append(dict(cohort=cohort,condition=r['condition'],sample_id=sid,profiler=tool,target_label=label,target_feature=feature,native_fraction=value,reported_positive=int(value>0),feature_row_present=int(feature in values)))
                groups.setdefault((cohort,tool,r['condition'],label),[]).append(value)
        require(len(conditions)==n and identities=={(sid,tool) for sid in conditions for tool in TOOLS},'Incomplete sample/profiler baseline coverage')
        require(set(conditions.values())==set(CONDITIONS),'Missing clinical condition')
    for (cohort,tool,condition,label),values in sorted(groups.items()):
        positive=[v for v in values if v>0]
        summary.append(dict(cohort=cohort,profiler=tool,condition=condition,target_label=label,n=len(values),positive_n=len(positive),reported_prevalence=len(positive)/len(values),mean_native_fraction=sum(values)/len(values),median_native_fraction=median(values),q25_native_fraction=quantile(values,.25),q75_native_fraction=quantile(values,.75),q95_native_fraction=quantile(values,.95),positive_only_median=median(positive) if positive else ''))
    lookup={(r['cohort'],r['profiler'],r['condition'],r['target_label']):r for r in summary}
    contrasts=[]
    for cohort in expected:
        for tool in TOOLS:
            for label in sorted(labels.values()):
                control=lookup[(cohort,tool,'Control',label)]
                for condition in ('Adenoma','CRC'):
                    case=lookup[(cohort,tool,condition,label)]
                    contrasts.append(dict(cohort=cohort,profiler=tool,target_label=label,contrast=condition+'-Control',case_n=case['n'],control_n=control['n'],reported_prevalence_difference=case['reported_prevalence']-control['reported_prevalence'],mean_native_fraction_difference=case['mean_native_fraction']-control['mean_native_fraction'],median_native_fraction_difference=case['median_native_fraction']-control['median_native_fraction'],interpretation='DESCRIPTIVE_UNADJUSTED_NOT_DA_INFERENCE'))
    require(all(digest(p)==sha for p,sha in hashes.items()),'Input changed during export')
    out.mkdir(parents=True)
    for name,rows in [('baseline_targets_long.tsv',long),('baseline_summary.tsv',summary),('clinical_descriptive_contrasts.tsv',contrasts),('input_hashes.tsv',[dict(path=str(p),sha256=h) for p,h in sorted(hashes.items())])]:write_table(out/name,rows)
    (out/'README.md').write_text('# Clinical target reference\n\nNative baseline abundance only, from checksum-verified profiles identified by the existing receipt-derived DA inventory. No models or new upstream processing.\n\nReported positive abundance is not proven biological presence. Missing target rows are reported as zero with feature_row_present=0. Quantiles include zeros; positive-only medians are separate. Clinical contrasts are descriptive, unadjusted and not differential-abundance results. Native fractions are profiler-specific, not interchangeable with inserted read fractions or quantitative recovery truth. No biological arms or production methods are selected by this export.\n')
    (out/'SUCCESS').write_text('status\tPASS_DESCRIPTIVE_CLINICAL_REFERENCE\nsamples\t%d\ntarget_rows\t%d\nproduction_authorized\t0\n'%(sum(expected.values()),len(long)))
    (out/'SHA256SUMS').write_text(''.join(digest(p)+'  '+p.name+'\n' for p in sorted(out.iterdir()) if p.is_file()))
    print('[PASS] Clinical reference:',sum(expected.values()),'people;',len(long),'target rows;',out)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inventory-root',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args();build(a.inventory_root,a.outdir)
