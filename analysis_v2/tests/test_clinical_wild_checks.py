#!/usr/bin/env python3
import sys, tempfile, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from clinical_wild_checks import plan, collect, run
from plan_da3_canary import finish
from audit_bracken_denominators import write_table, table

class Checks(unittest.TestCase):
    def test_grid_and_missing_results(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); pilot=root/'pilot';pilot.mkdir();i=0
            for c in ('yachida','feng','zeller'):
                for t in ('kraken2_bracken','metaphlan4'):
                    for b in ('Adenoma','CRC'):
                        source=pilot/('pilot_%04d'%i);source.mkdir();i+=1
                        write_table(source/'context.tsv',[dict(analysis='DA1',cohort=c,profiler=t,background=b)])
                        write_table(source/'metadata.tsv',[
                            dict(observation_id='s%d'%k,biological_sample_id='p%d'%k,group=int(k>=12),age=30+(k*7)%19,sex='Male' if k%2 else 'Female') for k in range(24)])
                        write_table(source/'abundance.tsv',[
                            dict(observation_id='s%d'%k,**{'f%02d'%j:((k*13+j*7)%31+1)/10000 for j in range(20)}) for k in range(24)])
                        finish(source,'fixture')
            plan(pilot,root/'plan');rows=table(root/'plan/tasks.tsv')
            self.assertEqual(len(rows),156)
            self.assertEqual(sum(r['scenario']=='actual' for r in rows),12)
            self.assertEqual(sum(int(r['simulations']) for r in rows),7200)
            # Small end-to-end fixtures exercise real R execution, copying and sealing.
            for r in rows:r['draws']=39;r['simulations']=0 if r['scenario']=='actual' else 2
            write_table(root/'plan/tasks.tsv',rows);finish(root/'plan','PASS_WILD_ENGINEERING_PLAN')
            repo=Path(__file__).resolve().parents[2]
            for index in (0,1):
                run(root/'plan',root/'results',index,repo)
                self.assertTrue((root/'results'/('task_%03d'%index)/'SUCCESS').is_file())
            with self.assertRaises(SystemExit):collect(root/'plan',root/'absent',root/'report')
            import json
            status=json.loads((root/'report/status.json').read_text())
            self.assertEqual(len(status['failures']),156)
            self.assertFalse(status['production_authorized'])
            with self.assertRaises(ValueError):plan(pilot,pilot/'overlap')
            for row in rows:
                folder=root/'complete'/('task_%03d'%int(row['index']));folder.mkdir(parents=True)
                if row['scenario']=='actual':write_table(folder/'comparison.tsv',[dict(feature='f00',wild_p=.5)])
                else:
                    write_table(folder/'draws.tsv',[
                        dict(simulation=int(row['chunk'])*50+k,method=m,discoveries=0,pointwise_rejections=1)
                        for k in range(1,51) for m in ('ordinary','HC3','wild')])
                finish(folder,'PASS_WILD_ENGINEERING_TASK')
            collect(root/'plan',root/'complete',root/'complete_report')
            summary=table(root/'complete_report/summary.tsv')
            self.assertEqual(len(summary),108)
            self.assertTrue(all(r['simulations']=='200' for r in summary))
            self.assertFalse(json.loads((root/'complete_report/status.json').read_text())['production_authorized'])

if __name__=='__main__':unittest.main()
