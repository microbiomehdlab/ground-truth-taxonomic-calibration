#!/usr/bin/env python3
"""Timing pilot orchestration fixtures; no cluster or real biological files."""
from __future__ import annotations

import contextlib
import io
import itertools
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'analysis_v2/scripts'))
import maaslin_timing as m


def source_fixture(root):
    root.mkdir()
    (root / 'batches').mkdir()
    catalog = dict(profiles={}, families={}, targets={t: {'target': 'T'} for t in m.TOOLS})
    tasks = []
    for j, (cohort, background, tool, n) in enumerate(itertools.product(
            m.COHORTS, ('Adenoma', 'CRC'), m.TOOLS, (5, 10, 15, 20))):
        catalog['families']['|'.join((cohort, tool, background))] = ['target', 'background']
        contexts = []
        for k in range(2500):
            c = dict(context_id='da3_%07d' % (2500*j+k), cohort=cohort,
                     background=background, profiler=tool, n=n, arm='N', anchor='', allocation_id='a0000')
            if k < 2 and n in (5, 20):
                c.update(arm='U', anchor=('0.0001', '0.001')[k], observations=[])
                for i in range(2*n):
                    sid = background + '_s%d' % i
                    group = int(i >= n)
                    dose = c['anchor'] if group else '0'
                    key = '|'.join((cohort, tool, sid, dose))
                    c['observations'].append(dict(sample_id=sid, group=group, dose=dose, key=key))
                    catalog['profiles'][key] = dict(sample_id=sid, cohort=cohort, condition=background,
                        profiler=tool, analysis_population='community', nominal_total_dose=dose,
                        source_profile='/fixture-only/' + key, sha256='unavailable-fixture-profile')
            contexts.append(c)
        name = 'batches/batch_%05d.json' % j
        (root / name).write_text(json.dumps(contexts))
        tasks.append(dict(index=j, file=name, contexts=len(contexts)))
    m.write_table(root / 'tasks.tsv', tasks)
    (root / 'catalog.json').write_text(json.dumps(catalog))
    m.finish(root, 'PASS_DRAFT_BATCH_PLAN_NOT_PRODUCTION_AUTHORIZATION')


class TimingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        cls.source, cls.plan = cls.root / 'source', cls.root / 'plan'
        source_fixture(cls.source)
        with contextlib.redirect_stdout(io.StringIO()):
            m.prepare(cls.source, cls.plan)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_selection_and_immutability(self):
        contexts = json.loads((self.plan / 'contexts.json').read_text())
        self.assertEqual(len(contexts), 48)
        self.assertEqual({(c['cohort'], c['background'], c['profiler'], c['n'], c['anchor'])
                          for c in contexts}, m.KEYS)
        self.assertTrue(all(c['arm'] == 'U' for c in contexts))
        m.verify(self.source)
        m.verify(self.plan)
        with self.assertRaisesRegex(ValueError, 'Fresh'):
            m.prepare(self.source, self.plan)
        with self.assertRaisesRegex(ValueError, 'Overlap'):
            m.fresh(self.source / 'new', [self.source])

    def test_unsealed_or_corrupt_inputs_refused(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d)
            (folder / 'consumed.tsv').write_text('content')
            (folder / 'SUCCESS').write_text('pass')
            (folder / 'SHA256SUMS').write_text(m.digest(folder / 'SUCCESS')+'  SUCCESS\n')
            with self.assertRaisesRegex(ValueError, 'absent'):
                m.checked(folder, ['consumed.tsv'])
            (folder / 'SUCCESS').write_text('changed')
            with self.assertRaisesRegex(ValueError, 'Checksum'):
                m.checked(folder, [])

    def test_sample_overlap_and_profile_identity(self):
        c = json.loads((self.plan / 'contexts.json').read_text())[0]
        catalog = json.loads((self.plan / 'catalog.json').read_text())
        c['observations'][1]['sample_id'] = c['observations'][0]['sample_id']
        with self.assertRaisesRegex(ValueError, 'Biological overlap'):
            m.validate_context(c, catalog)
        c = json.loads((self.plan / 'contexts.json').read_text())[0]
        catalog['profiles'][c['observations'][0]['key']]['condition'] = 'Control'
        with self.assertRaisesRegex(ValueError, 'identity mismatch'):
            m.validate_context(c, catalog)

    def test_incomplete_collection_and_projection(self):
        with tempfile.TemporaryDirectory() as d, contextlib.redirect_stdout(io.StringIO()):
            root = Path(d)
            self.assertEqual(m.collect(self.plan, root/'results', root/'incomplete'), 1)
            status = json.loads((root/'incomplete/status.json').read_text())
            self.assertEqual(len(status['failures']), 48)
            self.assertFalse((root/'incomplete/rough_projection.tsv').exists())
            contexts = json.loads((self.plan/'contexts.json').read_text())
            for i, c in enumerate(contexts):
                folder = root/'results'/('task_%03d' % i)
                folder.mkdir(parents=True)
                row = {k: c[k] for k in ('context_id', 'cohort', 'background', 'profiler', 'n', 'anchor')}
                row.update(index=i, task_seconds=c['n'], maaslin_seconds=1, fast_seconds=.01,
                           output_bytes=1024, maaslin_version='1.18.0', mismatched_features=int(i==0))
                m.write_table(folder/'summary.tsv', [row])
                (folder/'comparison.tsv').write_text('fixture\n')
                (folder/'input_hashes.tsv').write_text('fixture\n')
                m.finish(folder, 'PASS_MAASLIN_TIMING_COMPUTATION')
            self.assertEqual(m.collect(self.plan, root/'results', root/'complete'), 0)
            status = json.loads((root/'complete/status.json').read_text())
            self.assertEqual(status['status'], 'REVIEW_DIFFERENCES')
            self.assertEqual(status['completed_tasks'], 48)
            self.assertFalse(status['production_authorized'])
            projection = m.table(root/'complete/rough_projection.tsv')[0]
            self.assertEqual(int(projection['contexts']), 120000)
            self.assertAlmostEqual(float(projection['single_worker_hours']), 120000*12.5/3600)
            m.verify(root/'complete')

    def test_worker_failed_backend_retains_inputs_without_success(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            # Use a private modified plan and real tiny Bracken profiles.
            plan = root/'plan'; plan.mkdir()
            contexts = json.loads((self.plan/'contexts.json').read_text())
            catalog = json.loads((self.plan/'catalog.json').read_text())
            i = next(i for i,c in enumerate(contexts) if c['profiler']=='kraken2_bracken')
            for j, obs in enumerate(contexts[i]['observations']):
                path = root/('profile_%02d.tsv' % j)
                path.write_text('name\ttaxonomy_lvl\tfraction_total_reads\ntarget\tS\t0.001\nbackground\tS\t0.05\n')
                catalog['profiles'][obs['key']].update(source_profile=str(path),sha256=m.digest(path))
            (plan/'contexts.json').write_text(json.dumps(contexts))
            (plan/'catalog.json').write_text(json.dumps(catalog))
            (plan/'contexts.tsv').write_bytes((self.plan/'contexts.tsv').read_bytes())
            m.finish(plan,'PASS_MAASLIN_TIMING_PLAN')
            with patch.object(m.subprocess,'run',side_effect=subprocess.CalledProcessError(1,'Rscript')):
                with self.assertRaises(subprocess.CalledProcessError):
                    m.run(plan,root/'results',i,REPO)
            folder = root/'results'/('task_%03d' % i)
            self.assertTrue((folder/'abundance.tsv').is_file())
            self.assertFalse((folder/'SUCCESS').exists())
            self.assertEqual(len(m.table(folder/'abundance.tsv')), 2*contexts[i]['n'])

            def backend_fixture(command, **kwargs):
                target = Path(command[2])
                m.write_table(target/'timing.tsv', [dict(maaslin_version='1.18.0',
                    maaslin_seconds=1, fast_seconds=.01, mismatched_features=0)])
                (target/'comparison.tsv').write_text('fixture\n')
                return subprocess.CompletedProcess(command, 0)

            with patch.object(m.subprocess,'run',side_effect=backend_fixture), contextlib.redirect_stdout(io.StringIO()):
                m.run(plan,root/'successful_results',i,REPO)
            success = root/'successful_results'/('task_%03d' % i)
            m.checked(success, ('summary.tsv', 'comparison.tsv', 'input_hashes.tsv'))
            self.assertEqual(m.table(success/'summary.tsv')[0]['family_n'], '2')


if __name__ == '__main__':
    unittest.main()
