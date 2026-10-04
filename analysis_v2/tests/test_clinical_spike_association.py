#!/usr/bin/env python3
from __future__ import annotations
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
REPO=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'analysis_v2/scripts'))
import clinical_spike_association as c

class ClinicalSpikeTest(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name); self.clinical=self.root/'clinical'; self.clinical.mkdir()
        self.inventory=self.root/'inventory'; self.inventory.mkdir()
        summary=[]; targets=[]
        for cohort in c.COHORTS:
            inv=self.inventory/cohort; inv.mkdir(); profiles=[]
            for tool in c.TOOLS:
                cid=cohort+'_'+tool
                source=self.clinical/'contexts'/cid/'source'; source.mkdir(parents=True)
                metadata=[]; abundance=[]
                for i in range(12):
                    sid='person_%02d'%i; oid=sid+'_original'
                    metadata.append(dict(observation_id=oid,biological_sample_id=sid,
                                         group=str(i>=6 and 1 or 0),spike_state='original',
                                         age=30+i,sex='Female' if i%2 else 'Male'))
                    abundance.append(dict(observation_id=oid,**{'taxon%d'%j:(i+j+1)*1e-6 for j in range(10)}))
                    if i>=6:
                        for dose in c.DOSES:
                            path=self.root/(cohort+'_'+tool+'_'+sid+'_'+dose+'.tsv')
                            if tool=='kraken2_bracken':
                                c.write_table(path,[dict(name='taxon%d'%j,taxonomy_lvl='S',fraction_total_reads=.001+j*1e-5) for j in range(10)])
                            else:
                                path.write_text(''.join('k__Bacteria|s__taxon%d\t1\t%s\n'%(j,.1+j*.001) for j in range(10)))
                            profiles.append(dict(cohort=cohort,profiler=tool,condition='Adenoma',
                                analysis_population='community',sample_id=sid,nominal_total_dose=dose,
                                source_profile=str(path),sha256=c.digest(path)))
                c.write_table(source/'metadata.tsv',metadata); c.write_table(source/'abundance.tsv',abundance)
                c.write_table(source/'context.tsv',[dict(paired=0,covariates='age,sex',nominal_total_dose=0)])
                summary.append(dict(context_id=cid,cohort=cohort,profiler=tool,comparison='Adenoma_vs_Control',family_n=10))
                targets.extend(dict(context_id=cid,feature='taxon%d'%j,target_label='label%d'%j) for j in range(10))
            c.write_table(inv/'profile_inventory.tsv',profiles); c.seal(inv,'PASS_DESIGN_INVENTORY')
        c.write_table(self.clinical/'context_summary.tsv',summary)
        c.write_table(self.clinical/'target_results.tsv',targets)
        (self.clinical/'status.json').write_text(json.dumps(dict(status='PASS_DA1_CLINICAL_EXPORT')))
        c.seal(self.clinical,'PASS_DA1_CLINICAL_EXPORT')
        self.out=self.root/'output'

    def test_grid_and_unchanged_controls(self):
        c.prepare(self.out,self.clinical,self.inventory)
        self.assertEqual(len(c.table(self.out/'plan/tasks.tsv')),42)
        task=self.out/'plan/clinical_spike_00'
        original=c.table(task/'abundance.tsv')
        def fake_fit(args,check):
            folder=Path(args[-2])
            matrix=c.table(folder/'abundance.tsv'); meta=c.table(folder/'metadata.tsv')
            self.assertEqual(matrix[:6],original[:6])
            self.assertNotEqual(matrix[6:],original[6:])
            self.assertTrue(all(r['spike_state']=='original' for r in meta[:6]))
            self.assertTrue(all(r['spike_state']=='spiked' for r in meta[6:]))
            c.write_table(folder/'clinical_results.tsv',[dict(feature='taxon%d'%j,beta=1) for j in range(10)])
        with patch.object(c.subprocess,'run',side_effect=fake_fit): c.run(self.out,0,REPO)
        c.verify(self.out/'results/clinical_spike_00',{'clinical_results.tsv','profile_evidence.tsv'})
        with self.assertRaises(ValueError): c.collect(self.out)

    def test_missing_person_refused(self):
        inv=self.inventory/'feng'
        rows=c.table(inv/'profile_inventory.tsv'); c.write_table(inv/'profile_inventory.tsv',rows[1:]); c.seal(inv,'PASS_DESIGN_INVENTORY')
        with self.assertRaisesRegex(ValueError,'Missing/duplicate'): c.prepare(self.out,self.clinical,self.inventory)
        self.assertFalse((self.out/'plan').exists())

    def test_complete_collection(self):
        c.prepare(self.out,self.clinical,self.inventory)
        for t in c.table(self.out/'plan/tasks.tsv'):
            folder=self.out/'results'/t['context_id']; folder.mkdir(parents=True)
            c.write_table(folder/'clinical_results.tsv',[
                dict(feature='taxon%d'%j,beta=1,stderr=.2,raw_p=.001,wrapper_q=.01,
                     significant_two_sided='TRUE') for j in range(10)])
            (folder/'plan_sha256.txt').write_text(c.digest(self.out/'plan'/t['context_id']/'SHA256SUMS'))
            c.seal(folder,'PASS_CLINICAL_SPIKE_FIT')
        c.collect(self.out)
        c.verify(self.out/'REPORT',{'status.json','spiked_target_results.tsv','baseline_target_results.tsv'})
        self.assertEqual(len(c.table(self.out/'REPORT/spiked_target_results.tsv')),420)
        self.assertEqual(len(c.table(self.out/'REPORT/baseline_target_results.tsv')),60)

    def test_corrupt_profile_refused(self):
        c.prepare(self.out,self.clinical,self.inventory)
        task=json.loads((self.out/'plan/clinical_spike_00/task.json').read_text())
        Path(task['profiles'][0]['source_profile']).write_text('corrupted')
        with self.assertRaisesRegex(ValueError,'Changed spike profile'): c.run(self.out,0,REPO)
        self.assertFalse((self.out/'results/clinical_spike_00').exists())

    def test_metadata_corruption_refused(self):
        path=self.clinical/'contexts/yachida_kraken2_bracken/source/metadata.tsv'
        path.write_text(path.read_text()+'bad\n')
        with self.assertRaisesRegex(ValueError,'Checksum mismatch'): c.prepare(self.out,self.clinical,self.inventory)

if __name__=='__main__': unittest.main()
