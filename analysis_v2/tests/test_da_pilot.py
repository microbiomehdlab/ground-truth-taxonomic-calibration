#!/usr/bin/env python3
from __future__ import annotations
import shutil
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import prepare_da_pilot as p
import run_da_pilot_task as w
import build_da1_checkpoint as d
import test_da_design_inventory as fixture


class Tests(unittest.TestCase):
    def test_plan_memberships_doses_and_resumable_worker(self):
        original_write=fixture.m.write_table; original_build=fixture.m.build
        repo=Path(__file__).resolve().parents[2]
        panel=fixture.m.table(repo/'spikes/spike_panel.tsv')
        label={r['taxon_name']:r['label'] for r in panel}
        def write(path,rows):
            if path.name=='production_manifest.tsv':
                rows=[dict(r,age=str(40+i%10),sex='Male' if i%2 else 'Female') for i,r in enumerate(rows)]
            if path.name=='canonical.tsv':
                rows=[dict(r,target_label=label[r['target_taxon']],
                           spike_fraction_target=str(float(r['spike_fraction_total'])/(10 if r['analysis_population']=='community' else 1))) for r in rows]
            original_write(path,rows)
        def build(*args,**kwargs):
            if args[0].name=='canonical.tsv':
                canonical=args[0].with_name('canonical_input.tsv')
                if args[0].exists(): args[0].rename(canonical)
                args=(canonical,)+args[1:]
            result=original_build(*args,**kwargs)
            out=args[2]; root=out.parent
            inv=root/'inventories'; da1=root/'clinical'
            shutil.copytree(out,inv/'feng')
            d.build(out,da1/'feng','feng')
            for folder in (inv,da1):
                (folder/'SUCCESS').write_text('status\tPASS\n')
                members=sorted(x for x in folder.rglob('*') if x.is_file())
                (folder/'SHA256SUMS').write_text(''.join(p.digest(x)+'  '+str(x.relative_to(folder))+'\n' for x in members))
            plan=root/'plan'; p.build(inv,da1,plan,cohorts=('feng',))
            contexts=p.table(plan/'contexts.tsv'); obs=p.table(plan/'observations.tsv')
            self.assertEqual(len(contexts),68)
            self.assertEqual(sum(c['analysis']=='NULL' for c in contexts),20)
            self.assertTrue(p.table(plan/'nominal_matches.tsv'))
            for c in contexts:
                rows=[r for r in obs if r['context_id']==c['context_id']]
                if c['analysis'] in ('NULL','DA3'):
                    self.assertEqual(len({r['biological_sample_id'] for r in rows}),10)
                    self.assertEqual(sum(r['group']=='1' for r in rows),5)
                if c['analysis']=='NULL':
                    self.assertTrue(all(r['spike_state']=='original' for r in rows))
            def fake_fit(command,**kwargs):
                task=Path(command[2]); (task/'fit').mkdir()
                (task/'fit/SUCCESS').write_text('status\tPASS_CONTEXT\n')
            results=root/'results'
            with patch.object(w.subprocess,'run',side_effect=fake_fit) as backend:
                w.run(plan,results,0,repo)
                w.run(plan,results,0,repo)
                self.assertEqual(backend.call_count,1)
            w.verify(results/'pilot_0000')
            with self.assertRaises(ValueError): w.run(plan,results,999,repo)
            # A failed fit retains diagnostics and does not publish a success.
            with patch.object(w.subprocess,'run',side_effect=RuntimeError('mock failure')):
                with self.assertRaises(RuntimeError): w.run(plan,results,1,repo)
            self.assertFalse((results/'pilot_0001').exists())
            self.assertTrue(list((results/'attempts').glob('pilot_0001_*/FAILED.json')))
            with self.assertRaises(ValueError): p.build(inv,da1,plan,cohorts=('feng',))
            return result
        with patch.object(fixture.m,'write_table',write),patch.object(fixture.m,'build',build):
            with patch.object(fixture.Tests,'repetitions',6,create=True):
                fixture.Tests('test_receipt_verified_inventory_all_doses_families_and_allocations').test_receipt_verified_inventory_all_doses_families_and_allocations()


if __name__=='__main__': unittest.main()
