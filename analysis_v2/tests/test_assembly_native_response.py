from __future__ import annotations
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from compare_assembly_native_response import pair_rows

class AssemblyTests(unittest.TestCase):
    def fixture(self):
        base=dict(include='1',cohort='yachida',analysis_population='independent',sample_id='s',
                  target_label='Pana',profiler='metaphlan4',condition='Adenoma',source_profile='baseline',
                  abundance_fraction='.1',spike_fraction_target='0',assembly_arm='original')
        rows=[base]
        design=dict(sample_id='s',fraction='.01',f_hat='.01',N_inserted='10',R='990',R1='990',R2='990',seed='7')
        designs={'original':[design.copy()],'clean':[design.copy()]}
        for arm,value in (('original','.11'),('clean','.12')):
            rows.append(dict(base,assembly_arm=arm,abundance_fraction=value,spike_fraction_target='.01',
                             implanted_read_pairs_target='10',source_design=arm))
        return rows,designs
    def test_matched_delta(self):
        rows,designs=self.fixture();r=pair_rows(rows,designs)[0]
        self.assertAlmostEqual(r['replacement_minus_original_native'],.01)
        self.assertAlmostEqual(r['original_native_change'],.01)
    def test_seed_and_dose_mismatch(self):
        for field in ('seed','R','N_inserted'):
            rows,designs=self.fixture();designs['clean'][0][field]='99'
            with self.assertRaises(ValueError):pair_rows(rows,designs)
    def test_missing_and_duplicate(self):
        rows,designs=self.fixture()
        with self.assertRaises(ValueError):pair_rows(rows[:-1],designs)
        with self.assertRaises(ValueError):pair_rows(rows+[rows[-1]],designs)

if __name__=='__main__':unittest.main()
