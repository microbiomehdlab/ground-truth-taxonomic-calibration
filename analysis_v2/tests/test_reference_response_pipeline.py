"""Whole prepare/run/resume/collect fixture; R fitting is independently tested."""
from __future__ import annotations
import gzip
import json
import sys
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import reference_response_pilot as module
from audit_bracken_denominators import digest,write_table,table


class PipelineTests(unittest.TestCase):
    def test_complete_pipeline_and_tampering(self):
        repo=Path(__file__).resolve().parents[2]
        with tempfile.TemporaryDirectory() as temporary:
            base=Path(temporary);plan=base/'plan';plan.mkdir(); inventory=base/'inventory';inventory.mkdir()
            geff=base/'geff';geff.mkdir();report=base/'observed';report.mkdir()
            image=base/'image.sif';image.write_text('fixture image')
            labels={f'T{i}':f'L{i}' for i in range(10)};features=list(labels)+['bystander']
            sizes=geff/'targets';sizes.mkdir()
            write_table(sizes/'target_genome_sizes.tsv',[dict(target_label=l,genome_size_bp=2e6) for l in labels.values()])
            catalog=dict(profiles={},families={},targets={t:labels for t in module.TOOLS})
            contexts=[]
            for cohort in module.COHORTS:
                inv=inventory/cohort;inv.mkdir();source=base/(cohort+'_profiles');source.mkdir()
                canonical=[];gfolder=geff/cohort/'geff_primary_cov95';gfolder.mkdir(parents=True)
                write_table(gfolder/'effective_genome_size.tsv',[dict(sample_id=str(i),effective_genome_size_bp=2e6,mapping_coverage=1) for i in range(40)])
                for tool in module.TOOLS:
                    catalog['families']['|'.join((cohort,tool,'Adenoma'))]=features
                    for sid in map(str,range(40)):
                        for dose in ('0',)+module.DOSES:
                            path=source/(tool+'_'+sid+'_'+dose+'.tsv')
                            if tool==module.TOOLS[0]:
                                write_table(path,[dict(name=f,taxonomy_lvl='S',fraction_total_reads=.01 if f!='bystander' else .9,new_est_reads=10 if f!='bystander' else 900) for f in features])
                            else:
                                path.write_text(''.join('s__'+f+'\t0\t'+('1' if f!='bystander' else '90')+'\n' for f in features))
                            key='|'.join((cohort,tool,sid,dose))
                            catalog['profiles'][key]=dict(cohort=cohort,profiler=tool,sample_id=sid,condition='Adenoma',
                                analysis_population='community',nominal_total_dose=dose,source_profile=str(path),sha256=digest(path))
                            if dose=='0':continue
                            added=int(float(dose)*1000000);F=added/(1000000+added)
                            design=source/('design_'+sid+'_'+dose+'.tsv')
                            write_table(design,[dict(sample_id=sid,fraction=dose,R=1000000,N_total=added)])
                            for label in labels.values():
                                canonical.append(dict(analysis_population='community',spike_fraction_total=F,
                                    spike_fraction_target=(added//10)/(1000000+added),implanted_read_pairs_target=added//10,
                                    profiler=tool,sample_id=sid,source_profile=str(path),source_design=str(design),target_label=label))
                    for n in (10,20):
                        for dose in module.DOSES:
                            for allocation in range(20):
                                contexts.append(dict(context_id='da3_%07d'%len(contexts),cohort=cohort,profiler=tool,
                                    background='Adenoma',n=n,arm='U',anchor=dose,allocation_id='full_%04d'%allocation,
                                    observations=[dict(sample_id=str(i),group=int(i<n),dose=dose if i<n else '0',
                                        key='|'.join((cohort,tool,str(i),dose if i<n else '0'))) for i in range(2*n)]))
                canon=inv/'canonical_input.tsv';write_table(canon,canonical)
                write_table(inv/'input_hashes.tsv',[dict(path=str(canon),sha256=digest(canon))])
                module.seal(inv,'PASS_DESIGN_INVENTORY')
            module.seal(geff,'PASS')
            (plan/'catalog.json').write_text(json.dumps(catalog));(plan/'batch.json').write_text(json.dumps(contexts))
            write_table(plan/'tasks.tsv',[dict(index=0,file='batch.json',contexts=len(contexts))])
            module.seal(plan,'PASS_DRAFT_BATCH_PLAN_NOT_PRODUCTION_AUTHORIZATION')
            summaries=[];observed=[]
            for c in contexts:
                summaries.append(dict(context_id=c['context_id'],**{k:c[k] for k in ('cohort','profiler','background','n','arm','anchor','allocation_id')},family_n=11))
                for f,l in labels.items():
                    observed.append(dict(context_id=c['context_id'],target_label=l,feature=f,
                        beta='NA',stderr='NA',raw_p='NA',wrapper_q='1',estimable='FALSE',status='NON_ESTIMABLE_CONSTANT',positive_discovery='FALSE'))
            write_table(report/'context_summary.tsv',summaries)
            with gzip.open(report/'targets.tsv.gz','wt',newline='') as h:
                import csv
                w=csv.DictWriter(h,fieldnames=list(observed[0]),delimiter='\t');w.writeheader();w.writerows(observed)
            (report/'status.json').write_text(json.dumps(dict(completed_contexts=120000,mismatched_features=0,failures=[])))
            (report/'identity.json').write_text(json.dumps(dict(plan_sha256=digest(plan/'SHA256SUMS'),image_sha256=digest(image),
                code={n:digest(repo/'analysis_v2'/n) for n in ('lib/maaslin_contract.R','lib/maaslin_context.R','scripts/build_biomarker_abundance_input.py')})))
            module.seal(report,'PASS_DIRECT_MAASLIN')
            root=base/'extension'
            a=Namespace(out=root,repo=repo,plan=plan,inventory=inventory,geff_root=geff,observed_report=report,
                        image=str(image),allocations=20,arms='U',batch_size=720)
            module.prepare(a)
            self.assertEqual(json.loads((root/'identity.json').read_text())['contexts'],720)
            reference=json.loads((root/'reference_profiles.json').read_text())
            self.assertAlmostEqual(reference['yachida|kraken2_bracken|0|0.001']['T0'],.055)
            def fake_r(command,check):
                inp,out=Path(command[2]),Path(command[3]);out.mkdir()
                rows=[];targets=[]
                for c in table(inp/'contexts.tsv'):
                    for arm in ('reference','observed'):
                        for f in features:
                            r=dict(context_id=c['context_id'],arm=arm,feature=f,beta='NA',stderr='NA',raw_p='NA',
                                wrapper_q='1',estimable='FALSE',status='NON_ESTIMABLE_CONSTANT',positive_discovery='FALSE',
                                hc3_p='NA',hc3_q='1',hc3_stderr='NA',hc3_estimable='FALSE')
                            rows.append(r)
                            if f in labels:targets.append(dict(r,target_label=labels[f]))
                with gzip.open(out/'features.tsv.gz','wt',newline='') as h:
                    w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
                write_table(out/'targets.tsv',targets)
            with patch.object(module.subprocess,'run',side_effect=fake_r):module.run(root,0)
            # Completed batch resume does not fit again.
            with patch.object(module.subprocess,'run',side_effect=AssertionError('Unexpected refit')):module.run(root,0)
            module.collect(root,root.parent/'report')
            status=json.loads((root.parent/'report/status.json').read_text())
            self.assertEqual(status['target_rows'],7200)
            self.assertFalse(status['production_authorized'])
            with (root/'results/batch_0000/targets.tsv').open('a') as h:h.write('tampered\n')
            with self.assertRaises(ValueError):module.collect(root,root.parent/'bad-report')


if __name__=='__main__':unittest.main()
