#!/usr/bin/env python3
import json
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import summarize_bracken_recovery as m


def fixture(root,demo=False):
    for cohort in m.COHORTS:
        folder=root/cohort
        folder.mkdir(parents=True)
        rows=[]
        targets=['Bfrag','Csym','Dpne','Fnuc','Hhat','Pmic','Pana','Psto','Porp','Pint'] if demo else ['Fnuc']
        conditions=['Control','Adenoma','CRC'] if demo else ['Adenoma']
        doses=['.0001','.001','.01'] if demo else ['.001']
        populations=['community','independent'] if demo else ['community']
        for condition in conditions:
            for population in populations:
                for target_index,target in enumerate(targets):
                    for dose in doses:
                        for i in range(5 if demo else 2):
                            old = (target_index-2)*.3+i*.3+1 if demo else [0,4][i]
                            new = old*.75 if demo else 1
                            rows.append(dict(cohort=cohort,sample_id=condition+str(i),condition=condition,
                              analysis_population=population,target_label=target,nominal_total_dose=dose,
                              profile_id=condition+str(i)+population+target+dose,
                              implanted_fraction_target_exact=float(dose)/10,total_fraction_exact=float(dose),
                              detected_native_nonzero=str(int(i>0)),estimated_count_nonzero=str(int(i>0)),
                              legacy_recovery_ratio=old,recovery_ratio_all_input=new))
        m.write_table(folder/'bracken_recovery_comparison.tsv',rows)
        for name in ['pair_reconciliation.tsv','input_hashes.tsv']:
            (folder/name).write_text('placeholder\nfixture\n')
        summary=dict(cohort=cohort,status='PASS_ENDPOINT_CONSTRUCTION',target_endpoints=len(rows),
                     samples=len({r['sample_id'] for r in rows}),
                     populations=dict(m.Counter(r['analysis_population'] for r in rows)))
        (folder/'summary.json').write_text(json.dumps(summary))
        (folder/'SHA256SUMS').write_text(''.join(m.digest(p)+'  '+p.name+'\n' for p in sorted(folder.iterdir())))


class Tests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)/'input'
        self.out=Path(self.temp.name)/'summary'
        fixture(self.root)

    def test_paired_difference_and_all_zero_flags(self):
        m.build(self.root,self.out,True)
        rows=m.table(self.out/'paired_method_differences.tsv')
        self.assertEqual(len(rows),3)
        for r in rows:
            self.assertEqual(float(r['paired_error_delta_median']),-2)
            self.assertEqual(int(r['samples_error_lower']),2)
            self.assertEqual(int(r['native_zero_samples']),1)
        summary=m.table(self.out/'recovery_summary.tsv')
        self.assertEqual(float(summary[0]['recovery_ratio_q1']),1)
        self.assertEqual(float(summary[0]['recovery_ratio_q3']),3)

    def test_negative_values_retained(self):
        row=m.table(self.root/'feng/bracken_recovery_comparison.tsv')[0]
        row['recovery_ratio_all_input']='-2'
        ss,_=m.summarize_group(('feng','community','Adenoma','Fnuc','.001'),[row])
        self.assertEqual(ss[1]['recovery_ratio_median'],-2)
        self.assertEqual(ss[1]['negative_recovery_samples'],1)

    def test_hash_mismatch(self):
        with (self.root/'feng/summary.json').open('a') as h: h.write(' ')
        with self.assertRaises(ValueError): m.build(self.root,self.out,True)
        self.assertFalse(self.out.exists())

    def test_real_counts_required(self):
        with self.assertRaises(ValueError): m.build(self.root,self.out)

    def test_no_overwrite(self):
        self.out.mkdir()
        with self.assertRaises(ValueError): m.build(self.root,self.out,True)

    def test_no_overlap(self):
        with self.assertRaises(ValueError): m.build(self.root,self.root/'summary',True)


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='--demo':
        root=Path(sys.argv[2])
        fixture(root/'input',True)
        m.build(root/'input',root/'tables',True)
    else:
        unittest.main()
