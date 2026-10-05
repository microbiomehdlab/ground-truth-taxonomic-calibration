from __future__ import annotations
import itertools
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from reference_response_pilot import COHORTS,TOOLS,DOSES,selected_contexts,seal
from audit_bracken_denominators import write_table

class ExtensionTests(unittest.TestCase):
    def test_full_grid_contains_pilot_without_changing_allocations(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);contexts=[]
            for cohort,tool,n,dose,allocation in itertools.product(COHORTS,TOOLS,(10,20),DOSES,range(100)):
                contexts.append(dict(context_id='c%05d'%len(contexts),cohort=cohort,profiler=tool,
                    background='Adenoma',n=n,arm='U',anchor=dose,allocation_id='full_%04d'%allocation,
                    observations=[dict(sample_id=str(i),group=int(i<n),dose=dose if i<n else '0') for i in range(2*n)]))
            (root/'batch.json').write_text(json.dumps(contexts))
            write_table(root/'tasks.tsv',[dict(index=0,file='batch.json')])
            seal(root,'PASS_DRAFT_BATCH_PLAN_NOT_PRODUCTION_AUTHORIZATION')
            full=selected_contexts(root,100,('U',));pilot=selected_contexts(root,20,('U',))
            self.assertEqual(len(full),3600);self.assertEqual(len(pilot),720)
            full_by_id={r['context_id']:r for r in full}
            self.assertTrue(all(full_by_id[r['context_id']]==r for r in pilot))
    def test_shell_submissions_and_invalid_allocation_guard(self):
        repo=Path(__file__).resolve().parents[2]
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp);bin=base/'bin';bin.mkdir();calls=base/'calls'
            sbatch=bin/'sbatch'
            sbatch.write_text('#!/usr/bin/env bash\nprintf "%s\\n" "$*" >> "$CALLS"\necho 12345\n')
            sbatch.chmod(0o755)
            for name in ('plan/SHA256SUMS','inventory/yachida/SHA256SUMS','geff/SUCCESS','observed/status.json','image.sif'):
                p=base/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('fixture')
            env=dict(os.environ,PROJECT=str(repo),ANALYSIS_SIF=str(base/'image.sif'),
                DA3_BATCH_PLAN=str(base/'plan'),DA_INVENTORY_ROOT=str(base/'inventory'),GEFF_ROOT=str(base/'geff'),
                DA3_OBSERVED_REPORT=str(base/'observed'),CALLS=str(calls),PATH=str(bin)+os.pathsep+os.environ['PATH'],REFERENCE_CONCURRENCY='72')
            for allocations,last in (('20','71'),('100','359')):
                env.update(REFERENCE_ALLOCATIONS=allocations,REFERENCE_ROOT=str(base/('run'+allocations)))
                result=subprocess.run(['bash',str(repo/'analysis_v2/submit_reference_extension.sh')],env=env,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
                commands=calls.read_text().splitlines();self.assertEqual(len(commands),3)
                self.assertIn('--array=0-'+last+'%72',commands[1]);self.assertIn('--dependency=afterok:12345',commands[1])
                self.assertIn('--dependency=afterok:12345',commands[2]);calls.unlink()
            env.update(REFERENCE_ALLOCATIONS='30',REFERENCE_ROOT=str(base/'bad'))
            result=subprocess.run(['bash',str(repo/'analysis_v2/submit_reference_extension.sh')],env=env,capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0);self.assertFalse(calls.exists())

if __name__=='__main__':unittest.main()
