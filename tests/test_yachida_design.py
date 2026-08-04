#!/usr/bin/env python3
from __future__ import annotations

import csv
import pathlib
import subprocess
import tempfile
import unittest


ROOT = pathlib.Path(__file__).parents[1]


class YachidaDesignTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = pathlib.Path(self.temp.name)
        self.fields = ["sample_id", "Target_Condition", "age", "sex"]
        self.rows = [
            {"sample_id": "A1", "Target_Condition": "Adenoma", "age": "45", "sex": "Female"},
            {"sample_id": "A2", "Target_Condition": "Adenoma", "age": "65", "sex": "Male"},
            {"sample_id": "C1", "Target_Condition": "Control", "age": "40", "sex": "Female"},
            {"sample_id": "C2", "Target_Condition": "Control", "age": "42", "sex": "Female"},
            {"sample_id": "C3", "Target_Condition": "Control", "age": "61", "sex": "Male"},
            {"sample_id": "C4", "Target_Condition": "Control", "age": "68", "sex": "Male"},
            {"sample_id": "R1", "Target_Condition": "CRC", "age": "47", "sex": "Female"},
            {"sample_id": "R2", "Target_Condition": "CRC", "age": "49", "sex": "Female"},
            {"sample_id": "R3", "Target_Condition": "CRC", "age": "62", "sex": "Male"},
            {"sample_id": "R4", "Target_Condition": "CRC", "age": "69", "sex": "Male"},
        ]

    def tearDown(self):
        self.temp.cleanup()

    def write_manifest(self, path: pathlib.Path, rows):
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=self.fields, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)

    def run_pilot(self, source: pathlib.Path, output: pathlib.Path):
        subprocess.run([
            "python3", str(ROOT / "datasets/yachida/build_pilot_design.py"),
            "--manifest", str(source), "--output", str(output), "--per-condition", "2",
            "--selection-seed", "test-seed",
        ], check=True, capture_output=True, text=True)

    def read_rows(self, path: pathlib.Path):
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_selection_and_batches_are_order_invariant(self):
        forward = self.root / "forward.tsv"
        reverse = self.root / "reverse.tsv"
        self.write_manifest(forward, self.rows)
        self.write_manifest(reverse, reversed(self.rows))
        pilot_a = self.root / "pilot_a.tsv"
        pilot_b = self.root / "pilot_b.tsv"
        self.run_pilot(forward, pilot_a)
        self.run_pilot(reverse, pilot_b)
        ids_a = {row["sample_id"] for row in self.read_rows(pilot_a)}
        ids_b = {row["sample_id"] for row in self.read_rows(pilot_b)}
        self.assertEqual(ids_a, ids_b)
        self.assertEqual(len(ids_a), 6)

        mappings = []
        for label, pilot in (("a", pilot_a), ("b", pilot_b)):
            output = self.root / f"batched_{label}.tsv"
            batch_dir = self.root / f"batches_{label}"
            subprocess.run([
                "python3", str(ROOT / "scripts/assign_processing_batches.py"),
                "--manifest", str(pilot), "--output", str(output),
                "--batch-dir", str(batch_dir), "--max-batch-size", "2", "--batch-seed", "test-batch",
            ], check=True, capture_output=True, text=True)
            rows = self.read_rows(output)
            self.assertTrue(all(int(row["batch_size"]) <= 2 for row in rows))
            mappings.append({row["sample_id"]: (row["processing_order"], row["batch_id"]) for row in rows})
        self.assertEqual(mappings[0], mappings[1])


if __name__ == "__main__":
    unittest.main()
