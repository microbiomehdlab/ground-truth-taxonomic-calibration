#!/usr/bin/env python3
"""Safety and provenance gates for the provisional paper figure book."""
import csv
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis_v2/scripts/assemble_provisional_figure_book.py"
INVENTORY = ROOT / "analysis_v2/PROVISIONAL_FIGURE_INVENTORY.tsv"
FIELDS = ("figure", "panel", "title", "state", "asset", "seal", "note")


def inventory(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class FigureBook(unittest.TestCase):
    def invoke(self, root, manifest):
        return subprocess.run(["python3", str(SCRIPT), "--inventory", str(manifest),
                               "--project-root", str(root), "--outdir", str(root / "book")],
                              capture_output=True, text=True)

    def test_project_inventory_builds_with_explicit_placeholders(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            result = self.invoke(root, INVENTORY)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with (root / "book/FIGURE_MANIFEST.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 19)
            self.assertTrue(all(row["preview_status"] in
                                {"PLACEHOLDER", "MISSING PREVIEW"} for row in rows))
            self.assertIn("DEVELOPMENT ONLY", (root / "book/index.html").read_text())

    def test_available_preview_is_copied_and_hashed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "plot.png").write_bytes(b"png fixture")
            (root / "seal").write_text("PASS\n")
            manifest = root / "inventory.tsv"
            inventory(manifest, [dict(figure="Figure 2", panel="B", title="Recovery",
                                      state="checked", asset="plot.png", seal="seal", note="fixture")])
            result = self.invoke(root, manifest)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with (root / "book/FIGURE_MANIFEST.tsv").open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["preview_status"], "CHECKED")
            self.assertEqual(len(row["copied_sha256"]), 64)
            self.assertEqual((root / "book" / row["copied_asset"]).read_bytes(), b"png fixture")

    def test_checked_missing_seal_fails_closed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "plot.png").write_bytes(b"png fixture")
            manifest = root / "inventory.tsv"
            inventory(manifest, [dict(figure="Figure 2", panel="B", title="Recovery",
                                      state="checked", asset="plot.png", seal="missing", note="fixture")])
            result = self.invoke(root, manifest)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((root / "book").exists())


if __name__ == "__main__":
    unittest.main()
