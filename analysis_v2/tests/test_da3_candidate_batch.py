#!/usr/bin/env python3
from __future__ import annotations
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

REPO=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'analysis_v2/scripts'))
from run_da3_candidate import run
from collect_da3_candidate import collect
from prepare_da_pilot import verify
from audit_bracken_denominators import digest,write_table

def seal(folder):
    (folder/'SUCCESS').write_text('status\tPASS\n')
    (folder/'SHA256SUMS').write_text(''.join(digest(p)+'  '+p.name+'\n'
        for p in sorted(folder.iterdir()) if p.is_file() and p.name!='SHA256SUMS'))

class CandidateTests(unittest.TestCase):
    def test_real_positive_runner_and_failure_aware_collector(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);plan=root/'plan';plan.mkdir();pilot=root/'pilot';source=pilot/'pilot_0000';source.mkdir(parents=True)
            contexts=[dict(context_id='pilot_%04d'%i,analysis='DA3') for i in range(40)]
            write_table(plan/'contexts.tsv',contexts);seal(plan)
            write_table(source/'context.tsv',[dict(context_id='pilot_0000',analysis='DA3',paired='0',cohort='fixture')])
            write_table(source/'metadata.tsv',[dict(observation_id='s%d'%i,biological_sample_id='s%d'%i,group=int(i>=5)) for i in range(10)])
            write_table(source/'abundance.tsv',[dict(observation_id='s%d'%i,feature=0 if i<5 else (i-4)*.00001) for i in range(10)])
            seal(source);output=root/'output'
            run(REPO,output,root/'validation',pilot,plan,12)
            verify(output/'candidate_012')
            with self.assertRaises(ValueError):run(REPO,output,root/'validation',pilot,plan,12)
            null_source=root/'validation/task_012';null_source.mkdir(parents=True)
            for member in ('context.tsv','abundance.tsv'):
                (null_source/member).write_bytes((source/member).read_bytes())
            allocations=[dict(pool=pool,allocation_id=str(j),n='5',sample_id='s%d'%i,
                group='cases' if i>=5 else 'controls')
                for pool,count in (('full',1000),('independent_subset',252))
                for j in range(count) for i in range(10)]
            write_table(null_source/'allocations.tsv',allocations);seal(null_source)
            run(REPO,output,root/'validation',pilot,plan,0)
            verify(output/'candidate_000')
            with patch('builtins.print'):collect(output)
            import json
            status=json.loads((output/'REPORT/status.json').read_text())
            self.assertEqual(status['completed_tasks'],2)
            self.assertEqual(len(status['failures']),50)
            self.assertFalse(status['production_authorized'])

if __name__=='__main__':unittest.main()
