import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from export_clinical_target_reference import build,digest,table,write_table,TOOLS,CONDITIONS

class ClinicalReferenceTests(unittest.TestCase):
    def fixture(self,root):
        inventory=root/'inventory';folder=inventory/'feng';folder.mkdir(parents=True)
        rows=[]
        for i,condition in enumerate(CONDITIONS):
            for tool in TOOLS:
                path=root/('%s_%s.tsv'%(condition,tool))
                value=(i+1)*.001
                if tool=='kraken2_bracken':
                    path.write_text('name\ttaxonomy_lvl\tfraction_total_reads\nFusobacterium nucleatum\tS\t'+str(value)+'\n')
                else:
                    path.write_text('# test\nk__Bacteria|s__Fusobacterium_nucleatum_subsp._nucleatum\t1\t'+str(value*100)+'\n')
                rows.append(dict(cohort='feng',sample_id='s'+str(i),condition=condition,profiler=tool,analysis_population='community',nominal_total_dose='0',source_profile=str(path),sha256=digest(path)))
        write_table(folder/'profile_inventory.tsv',rows)
        (folder/'SUCCESS').write_text('status\tPASS_DESIGN_INVENTORY\n')
        self.seal(folder)
        return inventory,folder,rows

    def seal(self,folder):
        (folder/'SHA256SUMS').write_text(''.join(digest(folder/name)+'  '+name+'\n' for name in ('SUCCESS','profile_inventory.tsv')))

    def test_complete_export_and_missing_rows(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);inventory,folder,rows=self.fixture(root)
            out=root/'out';build(inventory,out,{'feng':3})
            self.assertEqual(len(table(out/'baseline_targets_long.tsv')),60)
            self.assertEqual(len(table(out/'baseline_summary.tsv')),60)
            self.assertEqual(len(table(out/'clinical_descriptive_contrasts.tsv')),40)
            absent=[r for r in table(out/'baseline_targets_long.tsv') if r['target_label']=='Bfrag']
            self.assertTrue(all(r['native_fraction']=='0.0' and r['feature_row_present']=='0' for r in absent))
            with self.assertRaises(ValueError):build(inventory,out,{'feng':3})

    def test_tampered_profile_and_duplicate_rejected_without_output(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);inventory,folder,rows=self.fixture(root)
            path=Path(rows[0]['source_profile']);original=path.read_text();path.write_text(original+'\n')
            out=root/'out'
            with self.assertRaises(ValueError):build(inventory,out,{'feng':3})
            self.assertFalse(out.exists());path.write_text(original)
            write_table(folder/'profile_inventory.tsv',rows+[rows[0]]);self.seal(folder)
            with self.assertRaises(ValueError):build(inventory,out,{'feng':3})
            self.assertFalse(out.exists())

if __name__=='__main__':unittest.main()
