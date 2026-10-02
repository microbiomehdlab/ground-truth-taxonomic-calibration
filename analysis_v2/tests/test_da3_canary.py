import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from plan_da3_canary import assign,DOSES,assemble,finish,write_table,table,digest,plan,collect

class CanaryTests(unittest.TestCase):
    def test_complete_matrix_and_failure_aware_collector(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);inventory=root/'inventory'
            for cohort in ('yachida','feng','zeller'):
                folder=inventory/cohort;folder.mkdir(parents=True)
                profiles=[];families=[]
                for background in ('Adenoma','CRC'):
                    for tool in ('kraken2_bracken','metaphlan4'):
                        families.append(dict(cohort=cohort,condition=background,profiler=tool,feature='Target'))
                        for i in range(40):
                            for dose in ('0',)+DOSES:
                                profiles.append(dict(cohort=cohort,condition=background,sample_id=background+str(i),profiler=tool,analysis_population='community',nominal_total_dose=dose,source_profile='/fixture/native',sha256='fixture'))
                write_table(folder/'profile_inventory.tsv',profiles);write_table(folder/'feature_families.tsv',families)
                finish(folder,'PASS_DESIGN_INVENTORY')
            out=root/'plan';plan(inventory,out)
            contexts=table(out/'contexts.tsv');obs=table(out/'observations.tsv')
            self.assertEqual(len(contexts),384)
            self.assertEqual({r['arm'] for r in contexts},{'U','P25','P50','P75','V','PV','N'})
            self.assertEqual({r['n'] for r in contexts},{'5','10','15','20'})
            for r in contexts:
                selected=[p for p in obs if p['context_id']==r['context_id']]
                self.assertEqual(len(selected),2*int(r['n']))
                self.assertEqual(len({p['biological_sample_id'] for p in selected}),len(selected))
                self.assertEqual(sum(p['exposed']=='1' for p in selected),int(r['exposed_n']))
            collect(out,root/'missing_results',root/'report')
            import json
            status=json.loads((root/'report/status.json').read_text())
            self.assertEqual(len(status['failures']),384)
            self.assertFalse(status['production_authorized'])

    def test_assignments_are_stable_and_shared(self):
        cases=['s'+str(i) for i in range(20)]
        for arm,p in [('P25','.25'),('P50','.5'),('P75','.75'),('V','1'),('PV','.5'),('N','0')]:
            a=assign(cases,'seed',arm,p,'.001')
            self.assertEqual(a,assign(list(reversed(cases)),'seed',arm,p,'.001'))
            self.assertEqual(sum(v!='0' for v in a.values()),int(float(p)*20))
            self.assertTrue(set(a.values())<=set(DOSES)|{'0','.001'})
        self.assertEqual(sum(v!='0' for v in assign(cases[:5],'seed','P50','.5','.001').values()),3)
        a=assign(cases,'seed','P25','.25','.001');b=assign(cases,'seed','P75','.75','.005')
        self.assertTrue({s for s,v in a.items() if v!='0'}<={s for s,v in b.items() if v!='0'})

    def test_real_profile_assembly_and_tamper(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);plan=root/'plan';plan.mkdir();profile=root/'native.tsv'
            profile.write_text('name\ttaxonomy_lvl\tfraction_total_reads\nTarget\tS\t0.01\n')
            write_table(plan/'contexts.tsv',[dict(context_id='c',n=5,arm='U',exposed_n=5,profiler='kraken2_bracken')])
            write_table(plan/'families.tsv',[dict(context_id='c',feature='Target')])
            write_table(plan/'observations.tsv',[dict(context_id='c',observation_id='s'+str(i),biological_sample_id='s'+str(i),group=int(i<5),exposed=int(i<5),source_profile=str(profile),sha256=digest(profile)) for i in range(10)])
            finish(plan,'FIXTURE');assemble(plan,root/'out',0,fit=True)
            self.assertEqual(len(table(root/'out'/'abundance.tsv')),10)
            self.assertEqual(table(root/'out'/'candidate_results.tsv')[0]['exact_p'],'1')
            profile.write_text(profile.read_text()+'\n')
            with self.assertRaises(ValueError):assemble(plan,root/'bad',0)
            self.assertFalse((root/'bad').exists())

if __name__=='__main__':unittest.main()
