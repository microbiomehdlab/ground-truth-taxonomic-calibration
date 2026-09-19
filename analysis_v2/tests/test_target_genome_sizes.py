#!/usr/bin/env python3
"""Fixture tests for target genome-size generation from the frozen spike panel."""

from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis_v2/scripts/build_target_genome_sizes.py"
LABELS = ["Bfrag", "Csym", "Dpne", "Fnuc", "Hhat", "Pmic", "Pana", "Psto", "Porp", "Pint"]


def panel(path: Path, entries):
    lines = ["label\ttaxon_name\tassembly\tfasta\tweight\turl"]
    for label, taxon, assembly, fasta in entries:
        lines.append(f"{label}\t{taxon}\t{assembly}\t{fasta}\t1\t")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def fasta(path: Path, length: int, line_width: int = 70):
    seq = "ACGT" * (length // 4) + "A" * (length % 4)
    body = "\n".join(seq[i:i + line_width] for i in range(0, len(seq), line_width))
    path.write_text(f">{path.stem} test\n{body}\n", encoding="utf-8")


def run(spike_panel, outdir, root=None, extra=()):
    cmd = ["python3", str(SCRIPT), "--spike-panel", str(spike_panel),
           "--outdir", str(outdir), *map(str, extra)]
    if root:
        cmd += ["--fasta-root", str(root)]
    return subprocess.run(cmd, text=True, capture_output=True)


def read(path):
    rows = [l.split("\t") for l in Path(path).read_text().splitlines()]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


class TargetGenomeSizeTest(unittest.TestCase):

    def _full_panel(self, root, lengths=None):
        entries = []
        for i, label in enumerate(LABELS):
            rel = f"references/genomes/{label}.fa"
            target = root / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            fasta(target, (lengths or {}).get(label, 1_000_000 + i * 1000))
            entries.append((label, f"Species {label}", f"GCF_00{i}.1", rel))
        p = root / "spike_panel.tsv"
        panel(p, entries)
        return p

    def test_measures_lengths_from_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            lengths = {"Dpne": 1_247_407, "Hhat": 5_697_783}
            p = self._full_panel(root, lengths)
            out = root / "out"
            done = run(p, out, root)
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            rows = {r["target_label"]: r for r in read(out / "target_genome_sizes.tsv")}
            self.assertEqual(len(rows), 10)
            # Measured, not remembered.
            self.assertEqual(int(rows["Dpne"]["genome_size_bp"]), 1_247_407)
            self.assertEqual(int(rows["Hhat"]["genome_size_bp"]), 5_697_783)
            self.assertTrue(all(len(r["fasta_sha256"]) == 64 for r in rows.values()))
            self.assertTrue(all(int(r["genome_size_bp"]) > 0 for r in rows.values()))
            self.assertEqual([r["target_label"] for r in
                              read(out / "target_genome_sizes.tsv")], sorted(LABELS))
            self.assertTrue((out / "SUCCESS").is_file())
            prov = (out / "target_genome_sizes_provenance.tsv").read_text()
            self.assertIn("spike_panel_sha256", prov)
            self.assertIn("targets\t10", prov)

    def test_missing_fasta_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            p = self._full_panel(root)
            (root / "references/genomes/Fnuc.fa").unlink()
            done = run(p, root / "out", root)
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("not found", done.stdout + done.stderr)
            self.assertIn("never substituted", done.stdout + done.stderr)

    def test_empty_fasta_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            p = self._full_panel(root)
            (root / "references/genomes/Fnuc.fa").write_text("", encoding="utf-8")
            done = run(p, root / "out", root)
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("empty", done.stdout + done.stderr)

    def test_duplicate_label_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            p = self._full_panel(root)
            with p.open("a", encoding="utf-8") as handle:
                handle.write("Fnuc\tSpecies Fnuc\tGCF_dup.1\treferences/genomes/Fnuc.fa\t1\t\n")
            done = run(p, root / "out", root)
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("duplicate target label", done.stdout + done.stderr)

    def test_unexpected_label_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            p = self._full_panel(root)
            text = p.read_text().replace("Porp\t", "Wrong\t")
            p.write_text(text, encoding="utf-8")
            done = run(p, root / "out", root)
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("do not match the expected set", done.stdout + done.stderr)

    def test_implausible_length_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            p = self._full_panel(root)
            fasta(root / "references/genomes/Fnuc.fa", 500)
            done = run(p, root / "out", root)
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("outside", done.stdout + done.stderr)

    def test_headerless_file_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            p = self._full_panel(root)
            (root / "references/genomes/Fnuc.fa").write_text(
                "ACGT" * 300_000, encoding="utf-8")
            done = run(p, root / "out", root)
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("no FASTA header", done.stdout + done.stderr)


if __name__ == "__main__":
    unittest.main()
