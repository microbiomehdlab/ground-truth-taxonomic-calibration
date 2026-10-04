#!/usr/bin/env python3
"""Tiny all-cohort fixtures exercise real paired R inference, not cluster data."""
import csv
import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path
REPO=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'analysis_v2/scripts'))
import da2_dose_response as m


def fixture(root):
    inv=root/'inventory'; br=root/'bracken'; mp=root/'metaphlan'; run=root/'run'; run.mkdir()
    with (REPO/'examples/spike_taxon_aliases.csv').open() as h: aliases=list(csv.DictReader(h))
    panel=m.table(REPO/'spikes/spike_panel.tsv')
    targets={tool:{p['label']:next(a['alias'] for a in aliases if a['canonical']==p['taxon_name'] and a['tool']==tool)
                   for p in panel} for tool in m.TOOLS}
    for cohort in m.COHORTS:
        profiles=[]; families=[]; brrows=[]; mprows=[]
        for tool in m.TOOLS:
            for condition in m.CONDITIONS:
                families.extend(dict(cohort=cohort,profiler=tool,condition=condition,feature=f) for f in targets[tool].values())
                for i in range(6):
                    sid=condition+str(i); base=.02*(i+1)
                    folder=root/'profiles'/cohort/tool/sid; folder.mkdir(parents=True)
                    def profile(pop,label,d):
                        pid=sid if d==0 else sid+'_'+label+'_f'+str(d).replace('.','p')
                        path=folder/(pid+'.tsv')
                        values={f:base for f in targets[tool].values()}
                        labels=list(targets[tool]) if pop=='community' else [label]
                        if d:
                            signal=d/(10 if pop=='community' else 1)
                            for l in labels: values[targets[tool][l]]=base*(1-d)+signal*(-.5 if i==0 else i/3)
                        if tool=='kraken2_bracken':
                            m.write_table(path,[dict(name=f,taxonomy_lvl='S',fraction_total_reads=v) for f,v in values.items()])
                        else:
                            path.write_text(''.join('s__'+f.replace(' ','_')+'\t1\t'+str(v*100)+'\n' for f,v in values.items()))
                        profiles.append(dict(cohort=cohort,sample_id=sid,condition=condition,analysis_population=pop,
                            profiler=tool,profile_id=pid,nominal_total_dose=str(d),source_profile=str(path),sha256=m.digest(path)))
                        if not d: return
                        for l in labels:
                            signal=d/(10 if pop=='community' else 1); rec=signal*(-.5 if i==0 else i/3)
                            common=dict(cohort=cohort,sample_id=sid,condition=condition,analysis_population=pop,
                                profile_id=pid,target_label=l,source_profile=str(path))
                            basepath=str(folder/(sid+'.tsv')); o=values[targets[tool][l]]
                            if tool=='kraken2_bracken':
                                brrows.append(dict(common,source_baseline=basepath,target_alias=targets[tool][l],nominal_total_dose=d,
                                    legacy_baseline_native=base,legacy_observed_native=o,total_fraction_exact=d,
                                    implanted_fraction_target_exact=signal,recovered_signal_all_input=rec,recovery_ratio_all_input=rec/signal,
                                    signed_error_all_input=rec-signal,legacy_recovery_ratio=rec/signal))
                            else:
                                mprows.append(dict(common,profiler=tool,source_baseline_profile=basepath,reference_type='genome_equivalent',
                                    baseline_abundance_fraction=base,observed_abundance_fraction=o,spike_fraction_total=d,
                                    spike_fraction_target=signal,implanted_signal_profiler_scale=signal,
                                    recovered_spike_signal_profiler_scale=rec,response_ratio_profiler_scale=rec/signal,
                                    response_residual_profiler_scale=rec-signal,response_ratio=rec/signal))
                    profile('community','CRCpanel',0)
                    for d in (.001,.005):
                        profile('community','CRCpanel',d)
                        if i<3: profile('independent','Fnuc',d)
        folder=inv/cohort; folder.mkdir(parents=True)
        m.write_table(folder/'profile_inventory.tsv',profiles); m.write_table(folder/'feature_families.tsv',families)
        p=REPO/'examples/spike_taxon_aliases.csv'
        m.write_table(folder/'input_hashes.tsv',[dict(path=str(p),sha256=m.digest(p))])
        m.seal(folder,'PASS_DESIGN_INVENTORY')
        folder=br/cohort; folder.mkdir(parents=True)
        m.write_table(folder/'bracken_recovery_comparison.tsv',brrows)
        (folder/'summary.json').write_text(json.dumps(dict(cohort=cohort,status='PASS_ENDPOINT_CONSTRUCTION')))
        m.seal(folder,'PASS_ENDPOINT_CONSTRUCTION')
        folder=mp/cohort/'endpoints'; folder.mkdir(parents=True)
        m.write_table(folder/'paired_endpoints.tsv',mprows)
    m.seal(inv,'PASS'); m.seal(mp,'PASS_ENDPOINTS')
    return run,inv,br,mp


class Tests(unittest.TestCase):
    def test_all_stages_real_r_and_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,inv,br,mp=fixture(Path(tmp))
            m.prepare(root,inv,br,mp,REPO,fixture=True)
            self.assertEqual(len(m.table(root/'plan/tasks.tsv')),18)
            for i in range(18): m.run(root,i,REPO)
            m.run(root,0,REPO)
            m.collect(root,REPO)
            report=root/'REPORT'; m.verify(report,{'SUCCESS','status.json','target_results.tsv'})
            status=json.loads((report/'status.json').read_text())
            self.assertEqual(status['contexts'],72); self.assertEqual(status['summary_cells'],396)
            self.assertEqual(len(list(report.glob('*.pdf'))),6)
            rows=m.table(report/'dose_response_summary.tsv')
            self.assertTrue(any(int(r['negative_recovery'])>0 for r in rows))
            self.assertEqual({r['n_pairs'] for r in rows if r['population']=='independent'},{'3'})
            self.assertEqual({r['family_n'] for r in rows},{'10'})
            bad=root/'results/task_00/target_results.tsv'; bad.write_text('corrupt')
            with self.assertRaises(ValueError): m.run(root,0,REPO)

    def test_missing_or_duplicate_endpoint(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,inv,br,mp=fixture(Path(tmp))
            p=br/'feng/bracken_recovery_comparison.tsv'; original=m.table(p)
            for rows in (original[:-1],original+[original[0]]):
                m.write_table(p,rows); m.seal(p.parent,'PASS_ENDPOINT_CONSTRUCTION')
                with self.assertRaises(ValueError): m.prepare(root,inv,br,mp,REPO,fixture=True)
                self.assertFalse((root/'plan').exists())

    def test_changed_source_and_missing_result(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,inv,br,mp=fixture(Path(tmp))
            m.prepare(root,inv,br,mp,REPO,fixture=True)
            context=json.loads((root/'plan/task_00/contexts.json').read_text())[0]
            path=Path(context['pairs'][0]['spiked']['source_profile']); path.write_text('changed')
            with self.assertRaises(ValueError): m.run(root,0,REPO)
            self.assertFalse((root/'results/task_00').exists())
            with self.assertRaises(ValueError): m.collect(root,REPO)
            self.assertFalse((root/'REPORT').exists())

    def test_production_counts_are_not_optional(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,inv,br,mp=fixture(Path(tmp))
            with self.assertRaises(ValueError): m.prepare(root,inv,br,mp,REPO)

    def test_complete_grid_and_truncated_dose_series(self):
        targets={t:{'feature'+str(i):'label'+str(i) for i in range(10)} for t in m.TOOLS}
        profiles=[]; families=[]
        for tool in m.TOOLS:
            for condition in m.CONDITIONS:
                families.extend(dict(profiler=tool,condition=condition,feature=f) for f in targets[tool])
                for i in range(10):
                    sid=condition+str(i)
                    common=dict(profiler=tool,condition=condition,sample_id=sid)
                    profiles.append(dict(common,analysis_population='community',nominal_total_dose='0',profile_id=sid))
                    for pop in ('community','independent'):
                        for label in (['CRCpanel'] if pop=='community' else targets[tool].values()):
                            for d in (m.DOSES+('0.1',) if pop=='community' else m.DOSES):
                                profiles.append(dict(common,analysis_population=pop,nominal_total_dose=d,
                                                     profile_id=sid+'_'+label+'_f'+d.replace('.','p')))
        contexts=m.make_contexts(profiles,families,targets,'feng')
        self.assertEqual(len(contexts),402)
        with self.assertRaises(ValueError): m.make_contexts(profiles[:-1],families,targets,'feng')
        truncated=[p for p in profiles if p['nominal_total_dose']!='0.1']
        with self.assertRaises(ValueError): m.make_contexts(truncated,families,targets,'feng')
        duplicated=profiles+[copy.deepcopy(profiles[-1])]
        with self.assertRaises(ValueError): m.make_contexts(duplicated,families,targets,'feng')


if __name__=='__main__': unittest.main()
