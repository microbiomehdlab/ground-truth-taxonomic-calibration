#!/usr/bin/env python3
from __future__ import annotations
import csv
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import build_da_design_inventory as m


class Tests(unittest.TestCase):
    def test_receipt_verified_inventory_all_doses_families_and_allocations(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); state=root/'state'; seal=state/'production_seal_v2'
            seal.mkdir(parents=True)
            aliases=Path(__file__).resolve().parents[2]/'examples/spike_taxon_aliases.csv'
            with aliases.open() as h: alias=list(csv.DictReader(h))
            names=sorted({r['canonical'] for r in alias})
            maps={(r['canonical'],r['tool']):r['alias'] for r in alias}
            manifest=[]; canonical=[]
            for condition in ('Control','Adenoma','CRC'):
                for i in range(10):
                    sid=condition+str(i)
                    manifest.append(dict(sample_id=sid,condition=condition,study='Study',independent_subset='1'))
                    receipt=[]
                    for profiler in ('metaphlan4','kraken2_bracken'):
                        for pop in ('community','independent'):
                            doses=('0','.0001','.0005','.001','.005','.01','.05','.1') if pop=='community' else ('0','.0001','.0005','.001','.005','.01','.05')
                            for dose in doses:
                                # Shared files per sample/profiler sufficient for inventory fixture;
                                # physical identity uniqueness is checked in canonical/profile IDs.
                                path=root/'profiles'/sid/(profiler+'.tsv')
                                if not path.exists():
                                    path.parent.mkdir(parents=True,exist_ok=True)
                                    if profiler=='kraken2_bracken':
                                        m.write_table(path,[dict(name=maps[n,profiler],taxonomy_lvl='S',fraction_total_reads='.001') for n in names])
                                    else:
                                        path.write_text(''.join('k__Bacteria|s__'+maps[n,profiler].replace(' ','_')+'\t2\t0.1\n' for n in names))
                                    receipt.append(dict(path=str(path),sha256=m.digest(path),bytes=path.stat().st_size))
                                for j,n in enumerate(names):
                                    tag='f'+dose.replace('.','p') if dose!='0' else ''
                                    pid=sid if dose=='0' else sid+'_'+('CRCpanel' if pop=='community' else 'target'+str(j))+'_'+tag
                                    canonical.append(dict(cohort='feng',sample_id=sid,condition=condition,profiler=profiler,
                                        analysis_population=pop,include='1',source_profile=str(path),spike_fraction_total=dose,
                                        profile_id=pid,target_label='target'+str(j),target_taxon=n))
                    m.write_table(state/'samples'/(sid+'.retained_outputs.tsv'),receipt) if (state/'samples').exists() else None
                    if not (state/'samples').exists():
                        (state/'samples').mkdir()
                        m.write_table(state/'samples'/(sid+'.retained_outputs.tsv'),receipt)
            m.write_table(seal/'production_manifest.tsv',manifest)
            m.write_table(seal/'sample_flow.tsv',[dict(sample_id=r['sample_id'],status='PASS') for r in manifest])
            (seal/'SUCCESS').write_text('cohort\tfeng\nstatus\tPASS\nsamples\t30\n')
            (seal/'production_seal.sha256').write_text(''.join(m.digest(p)+'  '+p.name+'\n' for p in sorted(seal.iterdir())))
            c=root/'canonical.tsv'; m.write_table(c,canonical)
            out=root/'out'
            m.build(c,state,out,'feng',30,repetitions=2)
            self.assertTrue((out/'SUCCESS').exists())
            families=m.table(out/'feature_families.tsv')
            self.assertEqual(len(families),60)
            self.assertTrue(all(r['baseline_n']=='10' for r in families))
            allocations=m.table(out/'allocations.tsv')
            self.assertEqual(len(allocations),2*(252+2)*10)
            # Corrupt source must fail before output is created.
            source=Path(canonical[0]['source_profile']); source.write_text('corrupt\n')
            with self.assertRaises(ValueError): m.build(c,state,root/'bad','feng',30,repetitions=2)
            self.assertFalse((root/'bad').exists())


if __name__=='__main__': unittest.main()
