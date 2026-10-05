import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from summarize_reference_response import summarize,quantile

class SummaryTests(unittest.TestCase):
    def test_ne_denominator_and_matched_effects(self):
        base=dict(cohort='feng',profiler='metaphlan4',n='10',nominal_total_dose='.0001',target_label='Fnuc',
                  background='Adenoma',exposure_arm='U',observed_estimable='TRUE',reference_estimable='TRUE',
                  observed_beta='1',reference_beta='2',observed_stderr='.5',reference_stderr='.4',
                  observed_positive_discovery='FALSE',reference_positive_discovery='TRUE')
        rows=[dict(base,context_id='a',category='REFERENCE_ONLY',hc3_category='NEITHER_POSITIVE'),
              dict(base,context_id='b',category='OBSERVED_NON_ESTIMABLE',hc3_category='OBSERVED_NON_ESTIMABLE',observed_estimable='FALSE',observed_beta='NA')]
        f,e,s=summarize(rows)
        self.assertEqual(next(r['frequency'] for r in s if r['arm']=='reference'),1)
        self.assertEqual(sum(r['frequency'] for r in f if r['inference']=='primary_MaAsLin2'),1)
        self.assertEqual(next(r['frequency'] for r in f if r['category']=='REFERENCE_ONLY'),.5)
        self.assertTrue(all(r['matched_estimable']==1 for r in e))
        self.assertEqual(next(r['median'] for r in e if r['quantity']=='beta_difference'),-1)
        with self.assertRaises(ValueError):summarize(rows+rows)
        self.assertEqual(quantile([], .5),'NA')
        self.assertEqual(quantile([0,4],.25),1)

if __name__=='__main__':unittest.main()
