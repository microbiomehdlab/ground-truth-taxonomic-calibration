from __future__ import annotations
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import reference_scaling_sensitivity as m
from reference_response_pilot import CODE,seal,expected_profile
from audit_bracken_denominators import digest,write_table,table

class ScalingTests(unittest.TestCase):
    def test_sealed_preparation_collection_and_snapshot_worker(self):
        with tempfile.TemporaryDirectory() as tmp:
            parent=Path(tmp);base=parent/'pilot';base.mkdir();image=parent/'image.sif';image.write_text('fixture')
            targets={f'T{i}':f'L{i}' for i in range(10)}
            catalog=dict(profiles={},targets={'metaphlan4':targets})
            baseline={f:.01 for f in targets};baseline['background']=.9
            native={};audits={};contexts=[];saved=[];comparisons=[]
            for cohort in ('yachida','feng','zeller'):
                for n in (10,20):
                    for dose in ('.0001','.001','.01'):
                        for allocation in range(20):
                            cid='c%04d'%len(contexts);obs=[]
                            for person in range(2*n):
                                group=int(person<n);actual=dose if group else '0';sid=str(person)
                                key='|'.join((cohort,'metaphlan4',sid,actual));bkey='|'.join((cohort,'metaphlan4',sid,'0'))
                                for k in (key,bkey):
                                    catalog['profiles'][k]=dict(source_profile='/saved/'+k)
                                    native['/saved/'+k]=baseline
                                obs.append(dict(key=key,sample_id=sid,group=group,dose=actual))
                                if actual!='0':
                                    counts={label:int(float(dose)*1e6)//10 for label in targets.values()};total=sum(counts.values())
                                    audits[key]=dict(key=key,cohort=cohort,sample_id=sid,profiler='metaphlan4',original_pairs=1000000,
                                        inserted_pairs=total,achieved_total_fraction=total/(1e6+total),effective_genome_size_bp=2000000,
                                        inserted_counts_json=json.dumps(counts),target_genome_sizes_json=json.dumps({l:2000000 for l in counts}))
                            contexts.append(dict(context_id=cid,cohort=cohort,profiler='metaphlan4',n=n,anchor=dose,arm='U',observations=obs))
                            for feature,label in targets.items():
                                saved.append(dict(context_id=cid,target_label=label,feature=feature,beta=1,stderr=.1,raw_p=.01,
                                    wrapper_q=.02,status='ESTIMABLE',estimable='TRUE',positive_discovery='TRUE'))
                                comparisons.append(dict(context_id=cid,target_label=label,cohort=cohort,profiler='metaphlan4',n=n,
                                    nominal_total_dose=dose,observed_beta='1',observed_estimable='TRUE',category='BOTH_POSITIVE'))
            all_contexts=contexts+[dict(c,context_id=c['context_id']+'kb',profiler='kraken2_bracken') for c in contexts]
            all_comparisons=comparisons+[dict(r,context_id=r['context_id']+'kb',profiler='kraken2_bracken') for r in comparisons]
            for name,value in (('contexts.json',all_contexts),('catalog.json',catalog),('native_profiles.json',native),
                ('identity.json',dict(image_path=str(image),image_sha256=digest(image)))):
                (base/name).write_text(json.dumps(value))
            write_table(base/'reference_audit.tsv',list(audits.values()));write_table(base/'observed_targets.tsv',saved)
            repo=Path(__file__).resolve().parents[1]
            for name in CODE:
                dest=base/'source/analysis_v2'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(repo/name,dest)
            seal(base,'PASS')
            original_hash=digest(base/'SHA256SUMS');original_report=base/'REPORT';original_report.mkdir()
            (original_report/'status.json').write_text(json.dumps(dict(contexts=720,target_rows=7200,status='COMPLETE_PENDING_SCIENTIFIC_REVIEW')))
            (original_report/'collector_provenance.json').write_text(json.dumps(dict(plan_sha256=original_hash)))
            write_table(original_report/'target_comparisons.tsv',all_comparisons)
            seal(original_report,'PASS')
            out=parent/'stress';m.prepare(base,out)
            self.assertEqual(digest(base/'SHA256SUMS'),original_hash)
            self.assertEqual(len(table(out/'scale_050/tasks.tsv')),36)
            self.assertEqual(table(out/'scale_050/observed_targets.tsv'),table(base/'observed_targets.tsv'))
            first=next(iter(audits));a=audits[first]
            lower=json.loads((out/'scale_050/reference_profiles.json').read_text())[first]
            upper=json.loads((out/'scale_200/reference_profiles.json').read_text())[first]
            nominal,_=expected_profile(baseline,'metaphlan4',json.loads(a['inserted_counts_json']),
                a['achieved_total_fraction'],targets,geff=2e6,genome_sizes=json.loads(a['target_genome_sizes_json']))
            self.assertLess(lower['T0'],nominal['T0']);self.assertGreater(upper['T0'],nominal['T0'])
            bkey=contexts[0]['observations'][-1]['key']
            self.assertEqual(json.loads((out/'scale_050/reference_profiles.json').read_text())[bkey],baseline)
            with patch.object(m.subprocess,'run') as call:
                m.run(out,36)
                self.assertIn('scale_200/source/analysis_v2/scripts/reference_response_pilot.py',call.call_args[0][0][2])
                self.assertEqual(call.call_args[0][0][-1],'0')
            # Exact report schema: no feature column is invented in comparisons.
            for scenario in m.VARIANTS:
                folder=out/scenario;report=folder/'REPORT';report.mkdir()
                write_table(report/'target_comparisons.tsv',comparisons)
                (report/'collector_provenance.json').write_text(json.dumps(dict(plan_sha256=digest(folder/'SHA256SUMS'))))
                seal(report,'PASS')
            m.report(out)
            self.assertEqual(len(table(out/'REPORT/target_comparisons.tsv')),10800)
            self.assertEqual(json.loads((out/'REPORT/status.json').read_text())['observed_refits'],0)
            with self.assertRaises(ValueError):m.prepare(base,out)
            with self.assertRaises(ValueError):m.run(out,72)
            # Changing an original input is detected before another plan exists.
            (base/'contexts.json').write_text('[]')
            with self.assertRaises(ValueError):m.prepare(base,parent/'bad')

if __name__=='__main__':unittest.main()
