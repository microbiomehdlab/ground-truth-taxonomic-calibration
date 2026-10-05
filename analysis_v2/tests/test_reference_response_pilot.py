#!/usr/bin/env python3
"""Reference scale, complete-family dilution, frozen selections and error gates."""
from __future__ import annotations
import json
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from reference_response_pilot import expected_profile, selected_contexts, classify, seal, fresh, COHORTS, TOOLS, DOSES, require_observed_equivalence, hc3_flags, reconcile_fractions


class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.targets={'T'+str(i):'L'+str(i) for i in range(10)}
        self.inserted={label:10 for label in self.targets.values()}
        self.baseline={'T0':.1,'bystander':.9}

    def test_bracken_native_not_all_input(self):
        p,scale=expected_profile(self.baseline,TOOLS[0],self.inserted,.01,self.targets,species_total=900)
        self.assertAlmostEqual(scale,.9)
        self.assertAlmostEqual(p['T0'],.1)
        self.assertAlmostEqual(p['T1'],.01)
        self.assertAlmostEqual(p['bystander'],.81)
        self.assertAlmostEqual(sum(p.values()),1)
        # Changing read fraction alone cannot change the native count reference.
        q,_=expected_profile(self.baseline,TOOLS[0],self.inserted,.1,self.targets,species_total=900)
        self.assertEqual(p,q)

    def test_genome_equivalent_and_full_composition(self):
        sizes={label:2e6 for label in self.inserted}
        p,scale=expected_profile(self.baseline,TOOLS[1],self.inserted,.1,self.targets,geff=2e6,genome_sizes=sizes)
        self.assertAlmostEqual(scale,.9); self.assertAlmostEqual(p['bystander'],.81)
        self.assertAlmostEqual(p['T0'],.1); self.assertAlmostEqual(sum(p.values()),1)
        sizes['L0']=1e6
        q,_=expected_profile(self.baseline,TOOLS[1],self.inserted,.1,self.targets,geff=2e6,genome_sizes=sizes)
        self.assertGreater(q['T0'],p['T0'])
        self.assertLess(q['bystander'],p['bystander'])

    def test_invalid_references(self):
        for counts in ({'L0':10},{**self.inserted,'L0':0},{**self.inserted,'L0':-1}):
            with self.assertRaises(ValueError):
                expected_profile(self.baseline,TOOLS[0],counts,.1,self.targets,species_total=900)
        for invalid in (0,-1,float('nan'),1):
            with self.assertRaises(ValueError):
                expected_profile(self.baseline,TOOLS[0],self.inserted,invalid,self.targets,species_total=900)
        with self.assertRaises(ValueError):
            expected_profile(self.baseline,TOOLS[1],self.inserted,.1,self.targets,geff=0,genome_sizes={})

    def test_nonestimability_not_nonsignificance(self):
        r=lambda e,d:dict(estimable=e,positive_discovery=d)
        self.assertEqual(classify(r('TRUE','FALSE'),r('TRUE','TRUE')),'REFERENCE_ONLY')
        self.assertEqual(classify(r('FALSE','FALSE'),r('TRUE','TRUE')),'OBSERVED_NON_ESTIMABLE')
        self.assertEqual(classify(r('TRUE','FALSE'),r('FALSE','FALSE')),'REFERENCE_NON_ESTIMABLE')
        self.assertEqual(classify(r('FALSE','FALSE'),r('FALSE','FALSE')),'BOTH_NON_ESTIMABLE')

    def test_complete_deterministic_720_context_selection(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/'plan';root.mkdir()
            contexts=[]
            for cohort in COHORTS:
                for tool in TOOLS:
                    for n in (10,20):
                        for dose in DOSES:
                            for allocation in range(20):
                                contexts.append(dict(context_id='da3_%07d'%len(contexts),cohort=cohort,
                                    profiler=tool,n=n,background='Adenoma',arm='U',anchor=dose,
                                    allocation_id='full_%04d'%allocation,
                                    observations=[dict(sample_id=str(i),group=int(i<n),dose=dose if i<n else '0') for i in range(2*n)]))
            (root/'batch.json').write_text(json.dumps(contexts))
            (root/'tasks.tsv').write_text('index\tfile\tcontexts\n0\tbatch.json\t720\n')
            seal(root,'PASS_DRAFT_BATCH_PLAN_NOT_PRODUCTION_AUTHORIZATION')
            self.assertEqual(len(selected_contexts(root)),720)
            with self.assertRaises(ValueError):selected_contexts(root,100)
            contexts.pop();(root/'batch.json').write_text(json.dumps(contexts))
            # Hash gate, then coverage gate even with an updated legitimate seal.
            with self.assertRaises(ValueError):selected_contexts(root)
            seal(root,'PASS_DRAFT_BATCH_PLAN_NOT_PRODUCTION_AUTHORIZATION')
            with self.assertRaises(ValueError):selected_contexts(root)

    def test_overlap_gate(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            with self.assertRaises(ValueError):fresh(root/'output',[root])

    def test_observed_equivalence_and_hc3(self):
        r=dict(estimable='TRUE',beta='1',stderr='.1',raw_p='.01',wrapper_q='.02')
        require_observed_equivalence(r,r)
        with self.assertRaises(ValueError):require_observed_equivalence(r,dict(r,beta='2'))
        self.assertEqual(hc3_flags(dict(hc3_estimable='FALSE',beta='NA',hc3_q='1'))['estimable'],'FALSE')
        self.assertEqual(hc3_flags(dict(hc3_estimable='TRUE',beta='1',hc3_q='.01'))['positive_discovery'],'TRUE')

    def historical_rows(self):
        exact=100/(1000003+100)
        recorded=format(exact,'.8f')
        design=dict(R='1000003',N_total='100',f_hat=recorded)
        rows=[dict(target_label=label,implanted_read_pairs_target='10',spike_fraction_total=recorded,
                   spike_fraction_target=repr(10/1000103)) for label in self.inserted]
        return rows,design,exact

    def test_historical_eight_decimal_total_is_accepted_without_changing_reference(self):
        rows,design,exact=self.historical_rows()
        self.assertGreater(abs(float(design['f_hat'])-exact),1e-10)
        original,added,inserted,observed=reconcile_fractions(rows,design)
        self.assertEqual((original,added),(1000003,100))
        self.assertEqual(observed,exact)
        self.assertNotEqual(observed,float(design['f_hat']))
        self.assertEqual(inserted,self.inserted)

    def test_rounding_does_not_allow_bad_design_or_canonical_totals(self):
        rows,design,_=self.historical_rows()
        with self.assertRaises(ValueError):reconcile_fractions(rows,dict(design,f_hat='.00011'))
        rows[0]['spike_fraction_total']=str(float(design['f_hat'])+1e-9)
        with self.assertRaises(ValueError):reconcile_fractions(rows,design)

    def test_target_fractions_remain_exact_and_counts_reconcile(self):
        rows,design,_=self.historical_rows()
        rows[0]['spike_fraction_target']=format(float(rows[0]['spike_fraction_target']),'.8f')
        with self.assertRaises(ValueError):reconcile_fractions(rows,design)
        rows,design,_=self.historical_rows()
        rows[0]['implanted_read_pairs_target']='11'
        with self.assertRaises(ValueError):reconcile_fractions(rows,design)


if __name__=='__main__':unittest.main()
