#!/usr/bin/env python3
from __future__ import annotations
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import summarize_metaphlan_recovery as m


def fixture(root):
    for cohort in m.COHORTS:
        folder=root/cohort/'endpoints'
        folder.mkdir(parents=True)
        rows=[]
        for i in range(2):
            rows.append(dict(profiler='metaphlan4',cohort=cohort,sample_id='S'+str(i),
                condition='Adenoma',analysis_population='community',target_label='Fnuc',
                profile_id='S'+str(i)+'_CRCpanel_f0p001',reference_type='genome_equivalent',
                response_ratio=[0,4][i],response_ratio_profiler_scale=1,
                spike_fraction_target='.0001',spike_fraction_total='.001',
                observed_detected_native_nonzero=str(i)))
        m.write_table(folder/'paired_endpoints.tsv',rows)
    (root/'SUCCESS').write_text('status\tPASS_ENDPOINTS\n')
    seal(root)


def seal(root):
    (root/'SHA256SUMS').write_text(''.join(m.digest(p)+'  ./'+str(p.relative_to(root))+'\n'
        for p in sorted(root.rglob('*')) if p.is_file() and p.name!='SHA256SUMS'))


class Tests(unittest.TestCase):
    def test_summary_pairing_labels_and_no_bracken_count_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,out=Path(tmp)/'input',Path(tmp)/'out'
            fixture(root)
            m.build(root,out,True)
            rows=m.table(out/'paired_method_differences.tsv')
            self.assertEqual(float(rows[0]['paired_error_delta_median']),-2)
            self.assertNotIn('estimated_count_zero_samples',rows[0])
            summary=m.table(out/'recovery_summary.tsv')
            self.assertEqual(summary[1]['method'],m.LABELS[1])
            self.assertEqual(float(summary[0]['recovery_ratio_q1']),1)

    def test_corruption_counts_and_overwrite_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,out=Path(tmp)/'input',Path(tmp)/'out'
            fixture(root)
            with self.assertRaises(ValueError): m.build(root,out)
            self.assertFalse(out.exists())
            with (root/'feng/endpoints/paired_endpoints.tsv').open('a') as h: h.write('corrupt\n')
            with self.assertRaises(ValueError): m.build(root,out,True)
            self.assertFalse(out.exists())
            with self.assertRaises(ValueError): m.build(root,root/'out',True)

    def test_negative_retained_and_nonfinite_refused(self):
        row=dict(reference_type='genome_equivalent',profile_id='S_CRCpanel_f0p001',
                 response_ratio='-2',response_ratio_profiler_scale='-1',
                 observed_detected_native_nonzero='0',spike_fraction_target='.0001',spike_fraction_total='.001')
        self.assertEqual(m.convert(row)['recovery_ratio_all_input'],-1)
        row['response_ratio']='nan'
        with self.assertRaises(ValueError): m.convert(row)


if __name__=='__main__': unittest.main()
