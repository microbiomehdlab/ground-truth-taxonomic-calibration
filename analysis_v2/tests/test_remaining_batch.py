import sys,tempfile,unittest,json,subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from plan_da3_batches import run,finish,digest,write_table,table
from remaining_validation import collect
REPO=Path(__file__).resolve().parents[2]

class RemainingTests(unittest.TestCase):
 def test_compact_batch_actual_exact_backend(self):
  with tempfile.TemporaryDirectory() as d:
   root=Path(d);plan=root/'plan';plan.mkdir();(plan/'batches').mkdir();profile=root/'native.tsv'
   profile.write_text('name\ttaxonomy_lvl\tfraction_total_reads\nTarget\tS\t0.01\n')
   profiles={str(i):dict(source_profile=str(profile),sha256=digest(profile)) for i in range(10)}
   (plan/'catalog.json').write_text(json.dumps(dict(profiles=profiles,families={'feng|kraken2_bracken|Adenoma':['Target']},targets={'kraken2_bracken':{'Target':'Fnuc'}})))
   context=dict(context_id='test',cohort='feng',profiler='kraken2_bracken',background='Adenoma',n=5,arm='N',anchor='0',allocation_id='a',permutation_seed=1,observations=[dict(sample_id='s'+str(i),key=str(i),group=int(i<5)) for i in range(10)])
   (plan/'batches/batch_00000.json').write_text(json.dumps([context]));write_table(plan/'tasks.tsv',[dict(index=0,file='batches/batch_00000.json',contexts=1)]);finish(plan,'TEST')
   run(plan,root/'results',0,REPO)
   self.assertEqual(table(root/'results/batch_00000/summary.tsv')[0]['positive_fwer_discoveries'],'0')
   import gzip,csv
   with gzip.open(root/'results/batch_00000/targets.tsv.gz','rt') as h:self.assertEqual(next(csv.DictReader(h,delimiter='\t'))['exact_p'],'1')
   with self.assertRaises(ValueError):run(plan,root/'results',0,REPO)

 def test_partial_null_runs_and_collector_keeps_missing(self):
  with tempfile.TemporaryDirectory() as d:
   root=Path(d);task=root/'task';task.mkdir()
   write_table(task/'settings.tsv',[dict(n=5,family_n=20,scenario='correlated')])
   subprocess.run(['Rscript',str(REPO/'analysis_v2/scripts/run_remaining_validation.R'),str(task),'partial_null',str(REPO)],check=True,stdout=subprocess.DEVNULL)
   self.assertEqual(len(table(task/'draws.tsv')),200)
   p=root/'plan';p.mkdir();write_table(p/'tasks.tsv',[dict(index=0,mode='partial_null')]);finish(p,'TEST')
   collect(p,root/'missing',root/'report')
   status=json.loads((root/'report/status.json').read_text());self.assertEqual(len(status['failures']),1);self.assertFalse(status['production_authorized'])

 def test_actual_design_synthetic_covariate_null(self):
  with tempfile.TemporaryDirectory() as d:
   task=Path(d)
   write_table(task/'settings.tsv',[dict(scenario='heteroskedastic')])
   write_table(task/'metadata.tsv',[dict(observation_id='s'+str(i),group=int(i>=10),age=30+i,sex='Male' if i%2 else 'Female') for i in range(20)])
   write_table(task/'abundance.tsv',[dict(observation_id='s'+str(i),target1=.001,target2=.002) for i in range(20)])
   subprocess.run(['Rscript',str(REPO/'analysis_v2/scripts/run_remaining_validation.R'),str(task),'clinical_null',str(REPO)],check=True,stdout=subprocess.DEVNULL)
   self.assertEqual(len(table(task/'draws.tsv')),200)

if __name__=='__main__':unittest.main()
