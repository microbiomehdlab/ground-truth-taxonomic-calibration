#!/usr/bin/env python3
from __future__ import annotations
import contextlib
import csv
import gzip
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

REPO=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'analysis_v2/scripts'))
import maaslin_production as p


def fixture(root):
    plan=root/'plan'; (plan/'batches').mkdir(parents=True)
    features=['taxon%d'%i for i in range(11)]
    targets={f:'T%d'%i for i,f in enumerate(features[:10])}
    catalog=dict(profiles={},families={'feng|kraken2_bracken|Adenoma':features},targets={'kraken2_bracken':targets})
    contexts=[]
    for idx in range(2):
        c=dict(context_id='da3_%07d'%idx,cohort='feng',profiler='kraken2_bracken',background='Adenoma',
               n=5,arm='P50' if idx==0 else 'N',anchor='0.001' if idx==0 else '',allocation_id='a000',observations=[])
        for i in range(10):
            sid='person%d'%i; group=int(i>=5); dose='0.001' if idx==0 and i>=8 else '0'
            key='|'.join(('feng','kraken2_bracken',sid,dose))
            c['observations'].append(dict(sample_id=sid,group=group,dose=dose,key=key))
            path=root/('%s_%s.tsv'%(sid,dose))
            path.write_text('name\ttaxonomy_lvl\tfraction_total_reads\n'+''.join(f+'\tS\t0.001\n' for f in features))
            catalog['profiles'][key]=dict(cohort='feng',profiler='kraken2_bracken',sample_id=sid,condition='Adenoma',
                nominal_total_dose=dose,analysis_population='community',source_profile=str(path),sha256=p.digest(path))
        contexts.append(c)
    (plan/'batches/batch_00000.json').write_text(json.dumps(contexts))
    (plan/'catalog.json').write_text(json.dumps(catalog))
    p.write_table(plan/'tasks.tsv',[dict(index=0,file='batches/batch_00000.json',contexts=2)])
    p.finish(plan,'PASS_DRAFT_BATCH_PLAN_NOT_PRODUCTION_AUTHORIZATION')
    state=dict(plan_sha256=p.digest(plan/'SHA256SUMS'),contexts=2,batches=1,
               code={name:p.digest(REPO/'analysis_v2'/name) for name in p.CODE})
    (root/'identity.json').write_text(json.dumps(state))
    (root/'identity.sha256').write_text(p.digest(root/'identity.json')+'  identity.json\n')
    return contexts,catalog


def fake_fit(command,**kwargs):
    inp,out=Path(command[2]),Path(command[3])
    rows=p.table(inp/'contexts.tsv'); summaries=[]; targets=[]; features=[]
    for row in rows:
        row.update(family_n=11,estimable=11,target_rows=10,mismatched_features=0,backend_version='1.18.0')
        summaries.append(row)
        names=list(p.table(inp/row['context_id']/'abundance.tsv')[0])[1:]
        for j,name in enumerate(names):
            values=dict(context_id=row['context_id'],feature=name,model_feature='feature_%06d'%(j+1),
                beta=1,stderr=1,raw_p=.5,native_q=.5,estimable='TRUE',status='ESTIMABLE',p_for_BH=.5,
                wrapper_q=.5,positive_discovery='FALSE')
            features.append(values)
        for r in p.table(inp/row['context_id']/'targets.tsv'):
            targets.append(dict(context_id=row['context_id'],**r))
    p.write_table(out/'summary.tsv',summaries)
    with gzip.open(out/'targets.tsv.gz','wt',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(targets[0]),delimiter='\t'); w.writeheader(); w.writerows(targets)
    with gzip.open(out/'features.tsv.gz','wt',newline='') as h:
        w=csv.DictWriter(h,fieldnames=p.FEATURE_FIELDS,delimiter='\t'); w.writeheader(); w.writerows(features)
    return subprocess.CompletedProcess(command,0)


class ProductionTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(); self.root=Path(self.temp.name)
        self.contexts,self.catalog=fixture(self.root)

    def tearDown(self): self.temp.cleanup()

    def test_run_resume_collection_and_immutability(self):
        before=p.digest(self.root/'plan/SHA256SUMS')
        with patch.object(p.subprocess,'run',side_effect=fake_fit) as backend,contextlib.redirect_stdout(io.StringIO()):
            p.run(self.root,0,REPO)
            p.run(self.root,0,REPO)
            self.assertEqual(backend.call_count,1)
            self.assertEqual(p.collect(self.root,REPO,self.root/'REPORT'),0)
        self.assertEqual(p.digest(self.root/'plan/SHA256SUMS'),before)
        p.verify(self.root/'results/batch_00000'); p.verify(self.root/'REPORT')
        status=json.loads((self.root/'REPORT/status.json').read_text())
        self.assertEqual(status['completed_contexts'],2); self.assertFalse(status['failures'])
        self.assertFalse(list((self.root/'attempts').iterdir()))
        path=next(iter(self.catalog['profiles'].values()))['source_profile']
        Path(path).write_text('changed')
        with self.assertRaisesRegex(ValueError,'Completed input changed'):
            p.run(self.root,0,REPO)

    def test_failure_preserves_scratch_and_missing_batch_report(self):
        with patch.object(p.subprocess,'run',side_effect=subprocess.CalledProcessError(1,'Rscript')):
            with self.assertRaises(subprocess.CalledProcessError): p.run(self.root,0,REPO)
        self.assertFalse((self.root/'results/batch_00000').exists())
        self.assertEqual(len(list((self.root/'attempts').glob('*/failed_inputs'))),1)
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(p.collect(self.root,REPO,self.root/'REPORT'),1)
        status=json.loads((self.root/'REPORT/status.json').read_text())
        self.assertEqual(status['completed_contexts'],0); self.assertEqual(len(status['failures']),1)

    def test_corrupt_final_never_overwritten(self):
        with patch.object(p.subprocess,'run',side_effect=fake_fit),contextlib.redirect_stdout(io.StringIO()):
            p.run(self.root,0,REPO)
        output=self.root/'results/batch_00000/summary.tsv'; output.write_text('corrupt')
        with self.assertRaises(ValueError): p.run(self.root,0,REPO)
        self.assertEqual(output.read_text(),'corrupt')

    def test_wrong_family_people_or_identity_rejected(self):
        self.contexts[0]['observations'][0]['sample_id']='person1'
        with self.assertRaisesRegex(ValueError,'disjoint'): p.context_features(self.contexts[0],self.catalog)
        self.catalog['targets']['kraken2_bracken'].pop('taxon0')
        with self.assertRaisesRegex(ValueError,'ten targets'): p.context_features(self.contexts[1],self.catalog)
        (self.root/'identity.json').write_text('{}')
        with self.assertRaisesRegex(ValueError,'identity changed'): p.run(self.root,0,REPO)

    def test_incomplete_source_plan_rejected_before_copy(self):
        target=self.root/'new'; target.mkdir()
        with self.assertRaisesRegex(ValueError,'4800 batches'):
            p.prepare(self.root/'plan',target,REPO,self.root/'image.sif')
        self.assertFalse((target/'plan').exists())

    def test_missing_full_family_rows_cannot_be_sealed(self):
        def broken(command,**kwargs):
            fake_fit(command,**kwargs)
            path=Path(command[3])/'features.tsv.gz'
            with gzip.open(path,'wt') as h: h.write('\t'.join(p.FEATURE_FIELDS)+'\n')
            return subprocess.CompletedProcess(command,0)
        with patch.object(p.subprocess,'run',side_effect=broken):
            with self.assertRaisesRegex(ValueError,'full-family result'): p.run(self.root,0,REPO)
        self.assertFalse((self.root/'results/batch_00000').exists())

    def test_launcher_dependencies_total_throttle_and_active_guard(self):
        # Local fake scheduler only: never calls Slurm or connects to a cluster.
        bins=self.root/'bin'; bins.mkdir()
        mock_git=bins/'git'
        mock_git.write_text('#!/bin/sh\nif [ "$1" = archive ]; then\n'
            'tar -cf - -C "$PROJECT" analysis_v2/run_maaslin_production.sbatch analysis_v2/submit_maaslin_production.sh\n'
            'else printf "fixture-commit\\n"; fi\n')
        mock_git.chmod(0o755)
        mock_sbatch=bins/'sbatch'
        mock_sbatch.write_text('#!'+sys.executable+'\nimport json,os,sys\nfrom pathlib import Path\n'
            'p=Path(os.environ["MOCK_CALLS"])\nrows=p.read_text().splitlines() if p.exists() else []\n'
            'with p.open("a") as h: h.write(json.dumps(sys.argv[1:])+"\\n")\nprint(1001+len(rows))\n')
        mock_sbatch.chmod(0o755)
        mock_squeue=bins/'squeue'; mock_squeue.write_text('#!/bin/sh\nprintf "1002_0\\n"\n'); mock_squeue.chmod(0o755)
        image=self.root/'image.sif'; image.write_text('fixture image')
        runroot=self.root/'launch'; calls=self.root/'calls.jsonl'
        env=dict(os.environ,PROJECT=str(REPO),ANALYSIS_SIF=str(image),DA3_BATCH_PLAN=str(self.root/'plan'),
            MAASLIN_RUN_ROOT=str(runroot),MAASLIN_CONCURRENCY='31',MOCK_CALLS=str(calls),
            PATH=str(bins)+os.pathsep+os.environ['PATH'])
        script=REPO/'analysis_v2/submit_maaslin_production.sh'
        result=subprocess.run(['bash',str(script)],env=env,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        submitted=[json.loads(line) for line in calls.read_text().splitlines()]
        self.assertEqual(len(submitted),8)
        arrays=submitted[1:7]
        throttles=[int(next(a for a in row if a.startswith('--array=')).split('%')[1]) for row in arrays]
        self.assertEqual(sum(throttles),31)
        self.assertTrue(all('--dependency=afterok:1001' in row for row in arrays))
        self.assertIn('--dependency=afterany:1002:1003:1004:1005:1006:1007',submitted[-1])
        (runroot/'identity.json').write_text('{}')
        resumed=subprocess.run(['bash',str(script),'--resume'],env=env,capture_output=True,text=True)
        self.assertNotEqual(resumed.returncode,0)
        self.assertIn('still active',resumed.stderr)
        self.assertEqual(len(calls.read_text().splitlines()),8)


if __name__=='__main__': unittest.main()
