#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from da_allocations import allocations


class Tests(unittest.TestCase):
    def test_nested_disjoint_deterministic_and_independent_namespaces(self):
        pool=['S'+str(i) for i in range(47)]
        a=list(allocations(pool,'feng','Adenoma',repetitions=3))
        self.assertEqual(a,list(allocations(list(reversed(pool)),'feng','Adenoma',repetitions=3)))
        self.assertEqual(len(a),12)
        for row in a:
            self.assertEqual(len(set(row['cases']+row['controls'])),2*row['n'])
        for start in range(0,12,4):
            for smaller,larger in zip(a[start:start+3],a[start+1:start+4]):
                self.assertEqual(smaller['cases'],larger['cases'][:smaller['n']])
                self.assertEqual(smaller['controls'],larger['controls'][:smaller['n']])
        other=list(allocations(pool,'zeller','Adenoma',repetitions=3))
        self.assertNotEqual(a[0]['seed_namespace'],other[0]['seed_namespace'])

    def test_exact_subset_enumeration(self):
        a=list(allocations([str(i) for i in range(10)],'feng','CRC',exhaustive=True))
        self.assertEqual(len(a),252)
        self.assertEqual(len({r['cases'] for r in a}),252)
        for r in a: self.assertEqual(len(set(r['cases']+r['controls'])),10)

    def test_feasibility_and_duplicates(self):
        self.assertEqual({r['n'] for r in allocations([str(i) for i in range(15)],'feng','CRC',repetitions=1)},{5})
        with self.assertRaises(ValueError): list(allocations(['A','A'],'feng','CRC'))
        with self.assertRaises(ValueError): list(allocations(['A'],'feng','CRC',exhaustive=True))


if __name__=='__main__': unittest.main()
