#!/usr/bin/env python3
import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location('audit', Path(__file__).resolve().parents[1]/'scripts/audit_bracken_denominators.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class AuditTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.state, self.results, self.out = [self.root/x for x in ('state', 'results', 'audit')]
        self.seal = self.state/'production_seal_v2'
        self.seal.mkdir(parents=True)
        (self.state/'samples').mkdir()
        self.receipt = []
        self.paths = []
        for category, total, unclassified in [('baseline', 100, 20), ('community', 110, 20)]:
            parent = self.results/'study'/'sample'/'profiles'/category/'profile'
            parent.mkdir(parents=True)
            report = parent/'profile.kraken2.report'
            report.write_text('20\t%d\t%d\tU\t0\tunclassified\n80\t%d\t0\tR\t1\troot\n' % (unclassified, unclassified, total-unclassified))
            bracken = parent/'profile.bracken.S.tsv'
            bracken.write_text('taxonomy_id\ttaxonomy_lvl\tnew_est_reads\tfraction_total_reads\n2\tS\t70\t1.00000\n')
            self.paths.append(bracken)
            for path in (report, bracken):
                self.receipt.append(dict(path=str(path), sha256=m.digest(path), bytes=path.stat().st_size))
        m.write_table(self.state/'samples/sample.retained_outputs.tsv', self.receipt)
        m.write_table(self.seal/'production_manifest.tsv', [dict(sample_id='sample', study='study', condition='CRC')])
        m.write_table(self.seal/'sample_flow.tsv', [dict(sample_id='sample', status='PASS', expected_baseline_profiles=1, expected_community_profiles=1, expected_independent_profiles=0)])
        (self.seal/'SUCCESS').write_text('cohort\tfeng\nsamples\t1\nstatus\tPASS\n')
        self.reseal()

    def reseal(self):
        (self.seal/'production_seal.sha256').write_text(''.join(m.digest(self.seal/name)+'  '+name+'\n' for name in ('production_manifest.tsv', 'sample_flow.tsv', 'SUCCESS')))

    def run_audit(self):
        m.audit('feng', self.state, self.results, self.out)

    def test_complete_and_no_source_mutation(self):
        before = {p: m.digest(p) for p in self.root.rglob('*') if p.is_file()}
        self.run_audit()
        rows = m.table(self.out/'profile_denominators.tsv')
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]['total_observations'], '100')
        self.assertEqual(rows[0]['classified_not_in_species'], '10')
        self.assertEqual(rows[0]['represented_fraction'], '0.7')
        self.assertEqual(rows[1]['delta_total_observations'], '10')
        self.assertEqual(rows[1]['spike_pair_reconciliation'], 'NOT_YET_VERIFIED')
        for p, h in before.items():
            self.assertEqual(m.digest(p), h)

    def test_corrupt_profile_rejected_before_output(self):
        self.paths[0].write_text('bad')
        with self.assertRaises(ValueError): self.run_audit()
        self.assertFalse(self.out.exists())

    def test_missing_report(self):
        self.paths[0].with_name('profile.kraken2.report').unlink()
        with self.assertRaises(ValueError): self.run_audit()

    def test_wrong_cohort(self):
        with self.assertRaises(ValueError): m.audit('zeller', self.state, self.results, self.out)

    def test_output_overlap(self):
        with self.assertRaises(ValueError): m.audit('feng', self.state, self.results, self.results/'audit')

    def test_no_overwrite(self):
        self.out.mkdir()
        with self.assertRaises(ValueError): self.run_audit()

    def test_expected_count_mismatch(self):
        path = self.seal/'sample_flow.tsv'
        rows = m.table(path)
        rows[0]['expected_independent_profiles'] = 60
        m.write_table(path, rows)
        self.reseal()
        with self.assertRaises(ValueError): self.run_audit()

    def test_seal_hash_mismatch(self):
        with (self.seal/'SUCCESS').open('a') as h: h.write('extra\tx\n')
        with self.assertRaises(ValueError): self.run_audit()

    def test_fraction_is_not_all_input(self):
        row = m.parse_counts(self.paths[0].with_name('profile.kraken2.report'), self.paths[0])
        self.assertEqual(row['native_fraction_sum'], 1)
        self.assertEqual(row['represented_fraction'], .7)

    def test_duplicate_taxid(self):
        with self.paths[0].open('a') as h: h.write('2\tS\t1\t0.1\n')
        with self.assertRaises(ValueError): m.parse_counts(self.paths[0].with_name('profile.kraken2.report'), self.paths[0])


if __name__ == '__main__':
    unittest.main()
