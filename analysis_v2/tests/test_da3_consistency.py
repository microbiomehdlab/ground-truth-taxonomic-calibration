from __future__ import annotations
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from summarize_da3_consistency import category

class AgreementTests(unittest.TestCase):
    def test_all_categories(self):
        r=lambda e,d:dict(estimable=e,positive_discovery=d)
        yes,no,ne=r('TRUE','TRUE'),r('TRUE','FALSE'),r('FALSE','FALSE')
        self.assertEqual(category(yes,yes),'BOTH')
        self.assertEqual(category(yes,no),'BRACKEN_ONLY')
        self.assertEqual(category(no,yes),'METAPHLAN_ONLY')
        self.assertEqual(category(no,no),'NEITHER')
        self.assertEqual(category(ne,yes),'AT_LEAST_ONE_NON_ESTIMABLE')
        with self.assertRaises(ValueError):category(r('FALSE','TRUE'),yes)

if __name__=='__main__':unittest.main()
