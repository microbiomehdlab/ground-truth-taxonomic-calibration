#!/usr/bin/env python3
import csv
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

SCRIPTS = Path(__file__).resolve().parents[1]/'scripts'
sys.path.insert(0, str(SCRIPTS))
import derive_bracken_all_input_recovery as m


class RecoveryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        root = Path(self.temp.name)
        self.project, self.audit, self.out = root/'project', root/'audit', root/'output'
        self.project.mkdir()
        self.audit.mkdir()
        (self.project/'spikes').mkdir()
        (self.project/'examples').mkdir()
        (self.project/'analysis_v2/scripts').mkdir(parents=True)
        shutil.copyfile(SCRIPTS/'build_crc_cohort_canonical_input.py', self.project/'analysis_v2/scripts/build_crc_cohort_canonical_input.py')
        self.panel = [dict(label='T'+str(i), taxon_name='Species '+str(i), weight='1') for i in range(10)]
        m.write_table(self.project/'spikes/spike_panel.tsv', self.panel)
        with (self.project/'examples/spike_taxon_aliases.csv').open('w', newline='') as h:
            w = csv.DictWriter(h, fieldnames=['canonical','alias','tool'])
            w.writeheader()
            for r in self.panel:
                w.writerow(dict(canonical=r['taxon_name'],alias=r['taxon_name'],tool='kraken2_bracken'))
        frozen = ['spikes/spike_panel.tsv', 'examples/spike_taxon_aliases.csv']
        (self.project/'analysis_v2/taxon_identity_freeze.sha256').write_text(''.join(m.digest(self.project/n)+'  '+n+'\n' for n in frozen))
        self.design = root/'design/CRCpanel.tsv'
        self.design.parent.mkdir()
        m.write_table(self.design, [dict(sample_id='S',community='CRCpanel',fraction='.01',R1=10000,R2=10000,R=10000,N_total=101,f_hat=format(101/10101,'.8f'))])
        allocation = m.allocate(101, self.panel)
        records, sources, canonical = [], [], []
        self.profile_paths = []
        for category, total in [('baseline',10000),('community',10101)]:
            pid = 'S' if category == 'baseline' else 'S_CRCpanel_f0p01'
            parent = root/'profiles'/category
            parent.mkdir(parents=True)
            path = parent/(pid+'.bracken.S.tsv')
            report = parent/(pid+'.kraken2.report')
            self.profile_paths.append(path)
            report.write_text('20\t2000\t2000\tU\t0\tunclassified\n80\t%d\t0\tR\t1\troot\n' % (total-2000))
            features = []
            for i, target in enumerate(self.panel):
                n = 100+(allocation[target['label']] if category == 'community' else 0)
                features.append(dict(name=target['taxon_name'],taxonomy_id=i+2,taxonomy_lvl='S',new_est_reads=n,fraction_total_reads=n/2000))
            m.write_table(path, features)
            records.append(dict(cohort='feng',sample_id='S',condition='CRC',category=category,profile=str(path),**m.parse_counts(report,path)))
            for source in (path,report): sources.append(dict(path=str(source),sha256=m.digest(source),role='profile'))
            for target, feature in zip(self.panel,features):
                baseline = category == 'baseline'
                canonical.append(dict(profiler='kraken2_bracken',include='1',cohort='feng',assembly_arm='original',
                    source_profile=str(path),analysis_population='community',target_label=target['label'],
                    sample_id='S',condition='CRC',target_taxon=target['taxon_name'],source_design=str(self.design),
                    implanted_read_pairs_target=0 if baseline else allocation[target['label']],
                    spike_fraction_total=0 if baseline else format(101/10101,'.8f'),
                    spike_fraction_target=0 if baseline else allocation[target['label']]/10101,
                    abundance_fraction=feature['fraction_total_reads']))
        (root/'canonical').mkdir()
        self.canonical = root/'canonical/canonical.tsv'
        m.write_table(self.canonical,canonical)
        m.write_table(self.audit/'profile_denominators.tsv',records)
        m.write_table(self.audit/'input_hashes.tsv',sources)
        m.write_table(self.audit/'design_evidence_inventory.tsv',[dict(path=str(self.design),sha256=m.digest(self.design))])
        (self.audit/'summary.json').write_text(json.dumps(dict(cohort='feng',status='PASS_DENOMINATOR_AUDIT_ONLY',profiles=2,samples=1)))
        self.reseal()
        self.args = SimpleNamespace(audit=self.audit,canonical=self.canonical,project=self.project,
                                    outdir=self.out,cohort='feng',community_allocation='frozen-panel-reconstruction')

    def reseal(self):
        (self.audit/'SHA256SUMS').write_text(''.join(m.digest(p)+'  '+p.name+'\n' for p in sorted(self.audit.iterdir()) if p.name != 'SHA256SUMS'))

    def change_canonical(self, change):
        rows=m.table(self.canonical)
        change(rows)
        m.write_table(self.canonical,rows)

    def test_integration_exact_recovery_and_legacy(self):
        before={p:m.digest(p) for p in self.project.parent.rglob('*') if p.is_file()}
        m.run(self.args)
        rows=m.table(self.out/'bracken_recovery_comparison.tsv')
        self.assertEqual(len(rows),10)
        self.assertEqual(sum(int(r['target_pairs']) for r in rows),101)
        self.assertNotEqual(rows[-1]['target_pairs'], rows[0]['target_pairs'])
        for r in rows:
            self.assertAlmostEqual(float(r['recovery_ratio_all_input']),1)
            self.assertGreater(float(r['legacy_recovery_ratio']),1)
        for p,h in before.items(): self.assertEqual(m.digest(p),h)

    def test_missing_canonical_target(self):
        self.change_canonical(lambda rows: rows.pop())
        with self.assertRaises(ValueError): m.run(self.args)
        self.assertFalse(self.out.exists())

    def test_wrong_target_count(self):
        self.change_canonical(lambda rows: rows[-1].update(implanted_read_pairs_target='999'))
        with self.assertRaises(ValueError): m.run(self.args)

    def test_tampered_design(self):
        with self.design.open('a') as h: h.write('\n')
        with self.assertRaises(ValueError): m.run(self.args)

    def test_tampered_profile(self):
        with self.profile_paths[0].open('a') as h: h.write('\n')
        with self.assertRaises(ValueError): m.run(self.args)

    def test_changed_panel(self):
        with (self.project/'spikes/spike_panel.tsv').open('a') as h: h.write('\n')
        with self.assertRaises(ValueError): m.run(self.args)

    def test_output_exists(self):
        self.out.mkdir()
        with self.assertRaises(ValueError): m.run(self.args)

    def test_output_overlaps_inputs(self):
        self.args.outdir=self.audit/'new'
        with self.assertRaises(ValueError): m.run(self.args)

    def test_wrong_cohort(self):
        self.args.cohort='zeller'
        with self.assertRaises(ValueError): m.run(self.args)

    def test_rounding_detection_difference(self):
        r=m.endpoints(0,1,1000000,1000010,10,0,0,.00001,.00001)
        self.assertEqual(r['estimated_count_nonzero'],1)
        self.assertEqual(r['detected_native_nonzero'],0)

    def test_negative_response_retained(self):
        r=m.endpoints(100,90,1000,1100,100,.1,.09,.1,.1)
        self.assertLess(r['recovery_ratio_all_input'],0)

    def test_zero_target_refused(self):
        with self.assertRaises(ValueError): m.endpoints(0,0,100,110,0,0,0,.1,.1)

    def test_individual_exact_fraction_vs_rounded_legacy(self):
        r=m.endpoints(0,7,100,107,7,0,.1,.06542056,.06542056)
        self.assertAlmostEqual(r['recovery_ratio_all_input'],1)
        self.assertAlmostEqual(r['legacy_recovery_ratio'],.1/.06542056)


    def test_individual_integration(self):
        root=self.project.parent
        design=root/'design/T0.tsv'
        m.write_table(design,[dict(sample_id='S',label='T0',fraction='.001',R1=10000,R2=10000,R=10000,N_inserted=10,f_hat=format(10/10010,'.8f'))])
        path=root/'profiles/independent/S_T0_f0p001.bracken.S.tsv'
        path.parent.mkdir()
        m.write_table(path,[dict(name='Species 0',taxonomy_id=2,taxonomy_lvl='S',new_est_reads=110,fraction_total_reads=.055)])
        report=path.with_name('S_T0_f0p001.kraken2.report')
        report.write_text('20\t2000\t2000\tU\t0\tunclassified\n80\t8010\t0\tR\t1\troot\n')
        rows=m.table(self.audit/'profile_denominators.tsv')
        rows.append(dict(cohort='feng',sample_id='S',condition='CRC',category='independent',profile=str(path),**m.parse_counts(report,path)))
        m.write_table(self.audit/'profile_denominators.tsv',rows)
        sources=m.table(self.audit/'input_hashes.tsv')
        sources.extend(dict(path=str(p),sha256=m.digest(p),role='profile') for p in (path,report))
        m.write_table(self.audit/'input_hashes.tsv',sources)
        designs=m.table(self.audit/'design_evidence_inventory.tsv')
        designs.append(dict(path=str(design),sha256=m.digest(design)))
        m.write_table(self.audit/'design_evidence_inventory.tsv',designs)
        rows=m.table(self.canonical)
        base=dict(rows[0],analysis_population='independent')
        spike=dict(base,source_profile=str(path),source_design=str(design),
                   implanted_read_pairs_target='10',spike_fraction_total=format(10/10010,'.8f'),
                   spike_fraction_target=format(10/10010,'.8f'),abundance_fraction='.055')
        rows.extend([base,spike])
        m.write_table(self.canonical,rows)
        summary=json.loads((self.audit/'summary.json').read_text())
        summary['profiles']=3
        (self.audit/'summary.json').write_text(json.dumps(summary))
        self.reseal()
        m.run(self.args)
        output=m.table(self.out/'bracken_recovery_comparison.tsv')
        individual=[r for r in output if r['analysis_population']=='independent']
        self.assertEqual(len(individual),1)
        self.assertAlmostEqual(float(individual[0]['recovery_ratio_all_input']),1)

    def test_pair_mismatch_even_with_valid_design_hash(self):
        rows=m.table(self.design)
        rows[0]['R']='9999'
        m.write_table(self.design,rows)
        m.write_table(self.audit/'design_evidence_inventory.tsv',[dict(path=str(self.design),sha256=m.digest(self.design))])
        self.reseal()
        with self.assertRaises(ValueError): m.run(self.args)
        self.assertFalse(self.out.exists())


if __name__ == '__main__': unittest.main()
