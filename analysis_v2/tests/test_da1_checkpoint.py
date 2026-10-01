#!/usr/bin/env python3
from __future__ import annotations
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import build_da1_checkpoint as d
import test_da_design_inventory as fixture


class Tests(unittest.TestCase):
    def test_eligibility_and_rank(self):
        self.assertEqual(d.eligibility(dict(age='50',sex='Female')), [])
        for age in ('','NaN','inf','-1','121'):
            self.assertIn('missing_or_invalid_age',d.eligibility(dict(age=age,sex='Male')))
        self.assertIn('missing_or_invalid_sex',d.eligibility(dict(age='50',sex='unknown')))
        rows=[dict(age=40+i,group=i//4,sex='Male' if i%2 else 'Female') for i in range(8)]
        self.assertEqual(d.design_rank(rows),4)
        for r in rows: r['sex']='Male' if r['group'] else 'Female'
        self.assertEqual(d.design_rank(rows),3)

    def test_full_checkpoint_and_corruption(self):
        original_write=fixture.m.write_table
        original_build=fixture.m.build
        def write(path,rows):
            if path.name=='production_manifest.tsv':
                rows=[dict(r,age=str(40+i%10),sex='Male' if i%2 else 'Female') for i,r in enumerate(rows)]
            original_write(path,rows)
        def build(*args,**kwargs):
            result=original_build(*args,**kwargs)
            out=args[2]
            d.build(out,out.parent/'da1','feng')
            self.assertEqual(len(d.table(out.parent/'da1'/'feature_families.tsv')),40)
            self.assertTrue(all(r['design_rank']=='4' and r['excluded']=='0'
                                for r in d.table(out.parent/'da1'/'model_design.tsv')))
            self.assertEqual(len(d.table(out.parent/'da1'/'eligibility.tsv')),40)
            with self.assertRaises(ValueError): d.build(out,out.parent/'da1','feng')
            source=Path(d.table(out/'profile_inventory.tsv')[0]['source_profile'])
            source.write_text('corrupt\n')
            with self.assertRaises(ValueError): d.build(out,out.parent/'refused','feng')
            self.assertFalse((out.parent/'refused').exists())
            return result
        with patch.object(fixture.m,'write_table',write),patch.object(fixture.m,'build',build):
            fixture.Tests('test_receipt_verified_inventory_all_doses_families_and_allocations').test_receipt_verified_inventory_all_doses_families_and_allocations()


if __name__=='__main__': unittest.main()
