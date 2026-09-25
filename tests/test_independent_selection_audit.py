import csv
import hashlib
import pathlib
import subprocess
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]


class IndependentSelectionAuditTests(unittest.TestCase):
    def write_manifest(self, path: pathlib.Path, reverse: bool = False) -> None:
        fields = ["Name", "Study condition", "Age", "Sex", "BMI"]
        rows = []
        for condition in ("Control", "Adenoma", "CRC"):
            for index in range(12):
                rows.append([
                    f"{condition}_{index:02d}", condition, str(45 + index),
                    "Female" if index % 2 else "Male", str(20 + index / 10),
                ])
        if reverse:
            rows.reverse()
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(fields)
            writer.writerows(rows)

    def select(self, source: pathlib.Path, output: pathlib.Path) -> None:
        subprocess.run([
            "python3", str(ROOT / "scripts/select_samples_deterministically.py"),
            "--manifest", str(source), "--output", str(output),
            "--per-condition", "10", "--selection-seed", "test-seed",
            "--id-column", "Name", "--condition-column", "Study condition",
        ], check=True)

    def test_selection_is_row_order_invariant_and_auditable(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            first_source = root / "first.tsv"
            second_source = root / "second.tsv"
            first_output = root / "first-selection.tsv"
            second_output = root / "second-selection.tsv"
            self.write_manifest(first_source)
            self.write_manifest(second_source, reverse=True)
            self.select(first_source, first_output)
            self.select(second_source, second_output)
            self.assertEqual(first_output.read_bytes(), second_output.read_bytes())

            audit = root / "audit"
            subprocess.run([
                "python3", str(ROOT / "datasets/audit_independent_selection.py"),
                "--eligible-manifest", str(first_source),
                "--selection", str(first_output), "--output-dir", str(audit),
            ], check=True)
            self.assertTrue((audit / "selection_balance.tsv").is_file())
            methods = (audit / "selection_methods.txt").read_text(encoding="utf-8")
            self.assertIn("no profiler, recovery, spike, or biomarker outcome", methods)
            expected = hashlib.sha256((audit / "selection_balance.tsv").read_bytes()).hexdigest()
            self.assertTrue((audit / "selection_balance.tsv.sha256").read_text().startswith(expected))


class SelectionProvenanceCollisionTests(unittest.TestCase):
    """The selector must never silently append a duplicate selection column.

    Appending a second `selection_rank` beside an inherited one is what produced
    the ambiguous historical Yachida independent manifest.
    """

    SELECTION_FIELDS = ("selection_rank", "selection_hash", "selection_seed")

    def base_manifest(self, path: pathlib.Path, with_selection: bool = False):
        fields = ["Name", "Study condition", "Age"]
        if with_selection:
            fields = fields + list(self.SELECTION_FIELDS)
        rows = []
        for condition in ("Control", "Adenoma", "CRC"):
            for index in range(12):
                row = [f"{condition}_{index:02d}", condition, str(45 + index)]
                if with_selection:
                    row += [str(index + 1), f"hash_{condition}_{index}", "pilot-seed"]
                rows.append(row)
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(fields)
            writer.writerows(rows)

    def select(self, source, output, extra=()):
        return subprocess.run([
            "python3", str(ROOT / "scripts/select_samples_deterministically.py"),
            "--manifest", str(source), "--output", str(output),
            "--per-condition", "10", "--selection-seed", "independent-seed",
            "--id-column", "Name", "--condition-column", "Study condition",
            *extra,
        ], text=True, capture_output=True)

    def header_of(self, path):
        with path.open(encoding="utf-8") as handle:
            return handle.readline().rstrip("\n").split("\t")

    def test_non_colliding_output_is_unchanged(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            source, output = root / "in.tsv", root / "out.tsv"
            self.base_manifest(source)
            done = self.select(source, output)
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            self.assertEqual(self.header_of(output),
                             ["Name", "Study condition", "Age",
                              "selection_rank", "selection_hash",
                              "selection_seed"])

    def test_collision_fails_without_explicit_prefixes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            source, output = root / "in.tsv", root / "out.tsv"
            self.base_manifest(source, with_selection=True)
            done = self.select(source, output)
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("already carries", done.stdout + done.stderr)
            self.assertIn("--existing-selection-prefix", done.stdout + done.stderr)
            self.assertFalse(output.exists())

    def test_prefixed_output_has_unique_headers(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            source, output = root / "in.tsv", root / "out.tsv"
            self.base_manifest(source, with_selection=True)
            done = self.select(source, output, extra=[
                "--existing-selection-prefix", "pilot",
                "--new-selection-prefix", "independent"])
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            header = self.header_of(output)
            self.assertEqual(len(header), len(set(header)), header)
            for field in self.SELECTION_FIELDS:
                self.assertIn("pilot_%s" % field, header)
                self.assertIn("independent_%s" % field, header)
                self.assertNotIn(field, header)
            with output.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 30)
            # Inherited provenance is preserved verbatim, new provenance added.
            self.assertTrue(all(row["pilot_selection_seed"] == "pilot-seed"
                                for row in rows))
            self.assertTrue(all(row["independent_selection_seed"]
                                == "independent-seed" for row in rows))
            self.assertTrue(all(row["pilot_selection_hash"].startswith("hash_")
                                for row in rows))

    def test_row_selection_is_unchanged_by_the_prefixes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            plain, annotated = root / "plain.tsv", root / "annotated.tsv"
            self.base_manifest(plain)
            self.base_manifest(annotated, with_selection=True)
            first, second = root / "first.tsv", root / "second.tsv"
            self.assertEqual(self.select(plain, first).returncode, 0)
            self.assertEqual(self.select(annotated, second, extra=[
                "--existing-selection-prefix", "pilot",
                "--new-selection-prefix", "independent"]).returncode, 0)
            with first.open(newline="", encoding="utf-8") as handle:
                chosen = [row["Name"] for row in csv.DictReader(handle, delimiter="\t")]
            with second.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["Name"] for row in rows], chosen)
            self.assertEqual([row["independent_selection_rank"] for row in rows],
                             [str(rank) for _ in range(3) for rank in range(1, 11)])

    def test_partial_prefixing_that_stays_ambiguous_is_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            source, output = root / "in.tsv", root / "out.tsv"
            self.base_manifest(source, with_selection=True)
            done = self.select(source, output, extra=[
                "--new-selection-prefix", ""])
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("already carries", done.stdout + done.stderr)

    def test_duplicate_input_header_is_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            source, output = root / "in.tsv", root / "out.tsv"
            self.base_manifest(source)
            text = source.read_text(encoding="utf-8").split("\n")
            text[0] = text[0] + "\tAge"
            source.write_text("\n".join(
                [text[0]] + [row + "\t1" if row else row for row in text[1:]]),
                encoding="utf-8")
            done = self.select(source, output)
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("duplicate column(s): Age", done.stdout + done.stderr)


if __name__ == "__main__":
    unittest.main()
