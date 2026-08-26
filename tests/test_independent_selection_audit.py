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


if __name__ == "__main__":
    unittest.main()
