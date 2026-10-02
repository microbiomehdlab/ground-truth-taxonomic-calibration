#!/usr/bin/env python3
from __future__ import annotations
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

REPO=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'analysis_v2/scripts'))
from run_da_robustness import run
from collect_da_robustness import collect
from audit_bracken_denominators import table,write_table,digest
from prepare_da_pilot import verify

def seal(folder):
    (folder/'SUCCESS').write_text('status\tPASS\n')
    (folder/'SHA256SUMS').write_text(''.join(digest(p)+'  '+p.name+'\n'
        for p in sorted(folder.iterdir()) if p.is_file() and p.name!='SHA256SUMS'))

class Tests(unittest.TestCase):
    def test_sharding_and_failure_collection(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);validation=root/'validation';folder=validation/'task_012';folder.mkdir(parents=True)
            write_table(folder/'context.tsv',[dict(cohort='yachida',profiler='kraken2_bracken',background='Adenoma')])
            (folder/'abundance.tsv').write_text('observation_id\tf\ns0\t0\n')
            rows=[dict(pool='full',n='10',allocation_id='full_%04d'%j,sample_id='s%d'%i,
                group='cases' if i>=10 else 'controls') for j in range(1000) for i in range(20)]
            write_table(folder/'allocations.tsv',rows);seal(folder)
            output=root/'output'
            def backend(cmd,**kwargs):
                task=Path(cmd[2]);alloc=table(task/'allocations.tsv');keys=sorted({r['allocation_id'] for r in alloc})
                self.assertEqual(len(keys),100)
                write_table(task/'draws.tsv',[dict(allocation_id=k,parametric_bh=0,mc_bh=0,mc_max_fwer=0) for k in keys])
            with patch('run_da_robustness.subprocess.run',side_effect=backend):
                run(REPO,output,validation,0);run(REPO,output,validation,1)
            verify(output/'robust_000');verify(output/'robust_001')
            first=table(output/'robust_000/allocations.tsv');second=table(output/'robust_001/allocations.tsv')
            self.assertFalse({r['allocation_id'] for r in first}&{r['allocation_id'] for r in second})
            self.assertTrue(all(0<int(r['seed'])<2147483647 for r in first))
            with patch('builtins.print'):collect(output)
            import json
            status=json.loads((output/'REPORT/status.json').read_text())
            self.assertEqual(status['completed_tasks'],2)
            self.assertEqual(status['status'],'INCOMPLETE')
            self.assertFalse(status['production_authorized'])
            summaries=table(output/'REPORT/comparison_summary.tsv')
            self.assertTrue(all(r['repetitions']=='200' for r in summaries))

if __name__=='__main__':unittest.main()
