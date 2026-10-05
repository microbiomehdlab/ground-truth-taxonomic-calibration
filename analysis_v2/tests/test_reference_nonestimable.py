from __future__ import annotations
import sys
import json
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_reference_nonestimable import describe,build
from reference_response_pilot import seal
from audit_bracken_denominators import digest,write_table,table

class InputTests(unittest.TestCase):
    def test_zero_and_nonzero_constants(self):
        groups=[0,0,1,1]
        zero=describe([0]*4,groups)
        self.assertEqual(zero['constant_kind'],'ALL_ZERO')
        self.assertEqual(zero['native_zero_count'],4)
        other=describe([.01]*4,groups)
        self.assertEqual(other['arithmetic_input_class'],'NON_ESTIMABLE_CONSTANT')
        self.assertEqual(other['constant_kind'],'NONZERO_CONSTANT_OR_NUMERICALLY_CONSTANT')
    def test_perfect_fit_is_not_constant_or_missed(self):
        x=describe([0,0,.001,.001],[0,0,1,1])
        self.assertEqual(x['arithmetic_input_class'],'NON_ESTIMABLE_PERFECT_FIXED_FIT')
        self.assertEqual(x['group0_zero_count'],2)
        self.assertEqual(x['group1_zero_count'],0)
        self.assertGreater(x['group1_transformed_mean'],x['group0_transformed_mean'])
    def test_variable_and_invalid(self):
        self.assertEqual(describe([0,.001,.002,.003],[0,0,1,1])['arithmetic_input_class'],'ESTIMABLE')
        for values,groups in (([0,0,0,0],[0]*4),([0,-1,0,0],[0,0,1,1]),
                              ([0,float('nan'),0,0],[0,0,1,1])):
            with self.assertRaises(ValueError):describe(values,groups)
    def test_sealed_audit_and_input_preservation(self):
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp);root=base/'run';root.mkdir();report=root/'REPORT';out=base/'audit'
            contexts=[];profiles={};native={};rows=[]
            for i,tool in enumerate(('metaphlan4','kraken2_bracken')):
                cid='c'+str(i);observations=[]
                for person in range(4):
                    key=cid+'_'+str(person);path='/source/'+key
                    profiles[key]={'source_profile':path};group=person//2
                    native[path]={'feature':.001 if i==1 and group==1 else 0}
                    observations.append(dict(sample_id=str(person),key=key,group=group,dose='.0001' if group else '0'))
                contexts.append(dict(context_id=cid,cohort='feng',profiler=tool,n=2,observations=observations))
                rows.append(dict(context_id=cid,cohort='feng',profiler=tool,n=2,nominal_total_dose='.0001',
                    allocation_id='full_0000',target_label='Fnuc',feature='feature',observed_estimable='FALSE',
                    observed_status='NON_ESTIMABLE_CONSTANT' if i==0 else 'NON_ESTIMABLE_PERFECT_FIXED_FIT'))
            (root/'contexts.json').write_text(json.dumps(contexts))
            (root/'catalog.json').write_text(json.dumps(dict(profiles=profiles,
                targets={t:{'feature':'Fnuc'} for t in ('metaphlan4','kraken2_bracken')})))
            (root/'native_profiles.json').write_text(json.dumps(native));seal(root,'PASS')
            before=digest(root/'SHA256SUMS');report.mkdir()
            (report/'collector_provenance.json').write_text(json.dumps(dict(plan_sha256=before)))
            (report/'status.json').write_text(json.dumps(dict(status='COMPLETE_PENDING_SCIENTIFIC_REVIEW',contexts=2,target_rows=2)))
            write_table(report/'target_comparisons.tsv',rows);seal(report,'PASS')
            build(root,report,out)
            self.assertEqual(digest(root/'SHA256SUMS'),before)
            self.assertEqual(len(table(out/'input_vectors.tsv')),8)
            self.assertEqual(json.loads((out/'status.json').read_text())['new_fits'],0)
            self.assertEqual(table(out/'input_summary.tsv')[0]['constant_kind'],'ALL_ZERO')
            with self.assertRaises(ValueError):build(root,report,out)
            native['/source/c0_0']['feature']=.002
            (root/'native_profiles.json').write_text(json.dumps(native))
            with self.assertRaises(ValueError):build(root,report,base/'bad')

if __name__=='__main__':unittest.main()
