#!/usr/bin/env python3
from __future__ import annotations
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

REPO=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'analysis_v2/scripts'))
from audit_bracken_denominators import digest, write_table
from run_final_validation_task import run
from collect_final_validation import collect
from prepare_da_pilot import verify

def seal(folder):
    (folder/'SUCCESS').write_text('status\tPASS\n')
    (folder/'SHA256SUMS').write_text(''.join(digest(p)+'  '+str(p.relative_to(folder))+'\n'
        for p in sorted(folder.rglob('*')) if p.is_file() and p.name!='SHA256SUMS'))

class BatchTests(unittest.TestCase):
    def test_modes_and_failure_collection(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); paired=root/'paired'; inventory=root/'inventory'; output=root/'output'
            source=paired/'pilot_0002'; source.mkdir(parents=True)
            write_table(source/'context.tsv',[dict(context_id='pilot_0002',analysis='DA2',
                inference_method='paired_differences',cohort='yachida',profiler='kraken2_bracken')])
            (source/'abundance.tsv').write_text('observation_id\tf\na\t0\n')
            (source/'metadata.tsv').write_text('observation_id\tgroup\na\t0\n'); seal(source)
            folder=inventory/'yachida'; folder.mkdir(parents=True)
            profiles=[]
            for i in range(40):
                profile=root/('native_%d.tsv'%i)
                profile.write_text('name\ttaxonomy_lvl\tfraction_total_reads\nf\tS\t0.1\n')
                profiles.append(dict(sample_id='s%02d'%i,condition='Adenoma',profiler='kraken2_bracken',
                    analysis_population='community',nominal_total_dose='0',source_profile=str(profile),sha256=digest(profile)))
            write_table(folder/'profile_inventory.tsv',profiles)
            write_table(folder/'feature_families.tsv',[dict(condition='Adenoma',profiler='kraken2_bracken',feature='f')])
            allocations=[]
            for pool,repetitions,sizes in [('full',1000,(5,10,15,20)),('independent_subset',252,(5,))]:
                for repetition in range(repetitions):
                    for n in sizes:
                        for i in range(2*n):
                            allocations.append(dict(condition='Adenoma',pool=pool,n=str(n),
                                allocation_id=str(repetition),group='cases' if i<n else 'controls',sample_id='s%02d'%i))
            write_table(folder/'allocations.tsv',allocations); seal(folder)
            def fake_backend(cmd,**kwargs):
                task=Path(cmd[2]); write_table(task/'summary.tsv',[dict(features=1)])
                if cmd[3]=='sensitivity':
                    write_table(task/'sensitivity_results.tsv',[dict(feature='f',pseudocount=str(p),
                        wrapper_q='.1',positive_discovery='FALSE',status='ESTIMABLE') for p in (1e-9,1e-8,1e-7)])
            with patch('run_final_validation_task.subprocess.run',side_effect=fake_backend):
                for index in (0,12,24):
                    run(REPO,output,paired,inventory,index)
                    verify(output/('task_%03d'%index))
                with self.assertRaises(ValueError): run(REPO,output,paired,inventory,0)
                with self.assertRaises(ValueError): run(REPO,output,paired,inventory,40)
            # Collector must classify failure, not silently discard missing tasks.
            with patch('builtins.print'): collect(output)
            import json
            status=json.loads((output/'REPORT/status.json').read_text())
            self.assertEqual(status['status'],'INCOMPLETE')
            self.assertFalse(status['production_authorized'])
            self.assertEqual(len(status['failures']),37)
            self.assertEqual(status['completed_tasks'],3)

    def test_complete_report_still_requires_review(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for i in range(40):
                folder=root/('task_%03d'%i); folder.mkdir()
                write_table(folder/'task.tsv',[dict(index=i,mode='synthetic')])
                write_table(folder/'summary.tsv',[dict(any_discovery_rate='.1' if i==0 else '0',
                    review_flag='TRUE' if i==0 else 'FALSE')])
                seal(folder)
            with patch('builtins.print'): collect(root)
            import json
            status=json.loads((root/'REPORT/status.json').read_text())
            self.assertEqual(status['status'],'REVIEW_FLAGS')
            self.assertEqual(status['completed_tasks'],40)
            self.assertEqual(status['review_flags'],1)
            self.assertFalse(status['production_authorized'])

if __name__=='__main__': unittest.main()
