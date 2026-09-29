#!/usr/bin/env python3
import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from check_spike_panel_identity import check


class SpikePanelIdentityTests(unittest.TestCase):
    def test_absolute_and_relative_fasta_are_equivalent(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "references/genomes").mkdir(parents=True)
            frozen = root / "frozen.tsv"
            runtime = root / "runtime.tsv"
            frozen.write_text("label\ttaxon_name\tassembly\tfasta\tweight\turl\nFnuc\tFusobacterium nucleatum\tA1\treferences/genomes/Fnuc.fa\t1\t\n")
            runtime.write_text(f"label\ttaxon_name\tassembly\tfasta\tweight\turl\nFnuc\tFusobacterium nucleatum\tA1\t{root}/references/genomes/Fnuc.fa\t1\t\n")
            check(runtime, frozen, root)
            runtime.write_text(runtime.read_text().replace("A1", "A2"))
            with self.assertRaisesRegex(ValueError, "assembly"):
                check(runtime, frozen, root)
            runtime.write_text(runtime.read_text().replace("A2", "A1").replace("Fnuc.fa", "Other.fa"))
            with self.assertRaisesRegex(ValueError, "FASTA"):
                check(runtime, frozen, root)


if __name__ == "__main__":
    unittest.main()
