#!/usr/bin/env python3
"""Synthetic saved-fit fixtures. No real clinical data or cluster access."""
import csv
import json
import os
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'analysis_v2/scripts'))
from export_da1_clinical import build, verify, COHORTS, PROFILERS, BACKGROUNDS
from audit_bracken_denominators import write_table, table, digest


def seal(folder):
    (folder/'SHA256SUMS').write_text(''.join(digest(p)+'  '+str(p.relative_to(folder))+'\n'
        for p in sorted(folder.rglob('*')) if p.is_file() and p.name!='SHA256SUMS'))


class ExportTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(); cls.root = Path(cls.tmp.name)
        cls.plan = cls.root/'plan'; cls.results = cls.root/'results'
        cls.plan.mkdir(); cls.results.mkdir()
        contexts=[]; observations=[]; families=[]; hashes=[]
        with (REPO/'examples/spike_taxon_aliases.csv').open() as h: aliases=list(csv.DictReader(h))
        rng=random.Random(12)
        for cohort in COHORTS:
            eligibility=[]
            for bg in BACKGROUNDS:
                eligibility.extend(dict(cohort=cohort,comparison=bg+'_vs_Control',sample_id=bg+str(i),
                    condition='Control' if i<10 else bg,group=str(int(i>=10)),age=str(30+(i*7)%43),
                    sex='Male' if i%2 else 'Female',eligible='1',exclusion_reason='',bmi='') for i in range(20))
                eligibility.append(dict(cohort=cohort,comparison=bg+'_vs_Control',sample_id=bg+'excluded',
                    condition=bg,group='1',age='',sex='Male',eligible='0',exclusion_reason='missing_or_invalid_age',bmi=''))
            checkpoint=cls.root/'checkpoint'/cohort; checkpoint.mkdir(parents=True)
            write_table(checkpoint/'eligibility.tsv',eligibility)
            hashes.append(dict(path=str(checkpoint/'eligibility.tsv'),sha256=digest(checkpoint/'eligibility.tsv')))
            for profiler in PROFILERS:
                features=[a['alias'] for a in aliases if a['tool']==profiler]
                for bg in BACKGROUNDS:
                    cid='pilot_%04d'%len(contexts)
                    c=dict(context_id=cid,cohort=cohort,profiler=profiler,analysis='DA1',background=bg,
                        population='community',nominal_total_dose='0',paired='0',covariates='age,sex',
                        family_n=str(len(features)),observations_n='20')
                    contexts.append(c); dest=cls.results/cid; dest.mkdir()
                    write_table(dest/'context.tsv',[dict(c,inference_method='maaslin2')])
                    obs=[dict(context_id=cid,observation_id=e['sample_id']+'_original',biological_sample_id=e['sample_id'],
                        group=e['group'],spike_state='original',age=e['age'],sex=e['sex']) for e in eligibility
                        if e['comparison']==bg+'_vs_Control' and e['eligible']=='1']
                    observations.extend(obs)
                    write_table(dest/'metadata.tsv',[{k:v for k,v in r.items() if k!='context_id'} for r in obs])
                    write_table(dest/'abundance.tsv',[dict(observation_id=r['observation_id'],**{
                        f:1e-8*(2**(8+int(r['group'])*.3+rng.random())-1) for f in features}) for r in obs])
                    families.extend(dict(context_id=cid,feature=f) for f in features)
                    (dest/'SUCCESS').write_text('status\tPASS_ENGINEERING_CONTEXT\n')
        alias=REPO/'examples/spike_taxon_aliases.csv'
        hashes.append(dict(path=str(alias),sha256=digest(alias)))
        for name,rows in [('contexts',contexts),('observations',observations),('families',families),('input_hashes',hashes)]:
            write_table(cls.plan/(name+'.tsv'),rows)
        (cls.plan/'SUCCESS').write_text('status\tPASS\n'); seal(cls.plan)
        subprocess.run(['Rscript','analysis_v2/tests/test_da1_clinical_export.R','--fixtures',str(cls.results)],
                       cwd=REPO,check=True,env=dict(os.environ,TZ='UTC'),stdout=subprocess.DEVNULL)
        for dest in cls.results.iterdir():
            (dest/'provenance.json').write_text(json.dumps({str(cls.plan/'SHA256SUMS'):digest(cls.plan/'SHA256SUMS')}))
            seal(dest)

    @classmethod
    def tearDownClass(cls): cls.tmp.cleanup()

    def test_full_export_and_safe_rerun(self):
        out=self.root/'complete'
        build(self.plan,self.results,out,REPO)
        verify(out,{'SUCCESS','target_results.tsv','full_family_results.tsv','eligibility.tsv'})
        self.assertEqual(len(table(out/'target_results.tsv')),120)
        self.assertEqual({r['excluded_n'] for r in table(out/'context_summary.tsv')},{'1'})
        self.assertEqual({r['residual_df'] for r in table(out/'target_results.tsv')},{'16'})
        self.assertTrue((out/'contexts/pilot_0000/source/original_bundle_checksums.txt').is_file())
        self.assertFalse((out/'contexts/pilot_0000/source/SHA256SUMS').exists())
        with self.assertRaises(Exception): build(self.plan,self.results,out,REPO)

    def test_uncovered_input_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/'SUCCESS').write_text('ok'); seal(root)
            (root/'metadata.tsv').write_text('unsealed')
            with self.assertRaises(Exception): verify(root,{'SUCCESS','metadata.tsv'})

    def test_changed_input_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/'SUCCESS').write_text('ok'); seal(root)
            (root/'SUCCESS').write_text('changed')
            with self.assertRaises(Exception): verify(root,{'SUCCESS'})

    def test_missing_context_rejected(self):
        source=self.results/'pilot_0000'; hidden=self.root/'hidden'
        source.rename(hidden)
        try:
            out=self.root/'incomplete'
            with self.assertRaises(Exception): build(self.plan,self.results,out,REPO)
            self.assertFalse(out.exists())
        finally: hidden.rename(source)

    def test_overlap_rejected(self):
        with self.assertRaises(Exception): build(self.plan,self.results,self.results/'new',REPO)
        self.assertFalse((self.results/'new').exists())

    def test_unsafe_and_duplicate_checksums(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/'SUCCESS').write_text('ok'); seal(root)
            manifest=(root/'SHA256SUMS').read_text()
            for invalid in (manifest+manifest, digest(root/'SUCCESS')+'  ../SUCCESS\n'):
                (root/'SHA256SUMS').write_text(invalid)
                with self.assertRaises(Exception): verify(root,{'SUCCESS'})

    def test_sealed_but_wrong_metadata_rejected(self):
        source=self.results/'pilot_0000'
        original=(source/'metadata.tsv').read_text()
        try:
            rows=table(source/'metadata.tsv'); rows[0]['age']='100'
            write_table(source/'metadata.tsv',rows); seal(source)
            out=self.root/'wrong_metadata'
            with self.assertRaisesRegex(Exception,'Metadata differs'):
                build(self.plan,self.results,out,REPO)
            self.assertFalse(out.exists())
        finally:
            (source/'metadata.tsv').write_text(original); seal(source)

    def test_sealed_but_wrong_family_rejected(self):
        source=self.results/'pilot_0000'
        original=(source/'fit/context_results.tsv').read_text()
        try:
            rows=table(source/'fit/context_results.tsv'); rows[0]['feature']='wrong_taxon'
            write_table(source/'fit/context_results.tsv',rows); seal(source)
            out=self.root/'wrong_family'
            with self.assertRaisesRegex(Exception,'Fit family differs'):
                build(self.plan,self.results,out,REPO)
            self.assertFalse(out.exists())
        finally:
            (source/'fit/context_results.tsv').write_text(original); seal(source)


if __name__=='__main__': unittest.main()
