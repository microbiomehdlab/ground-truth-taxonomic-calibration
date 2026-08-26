import csv
import pathlib
import subprocess
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]


class EligibleManifestTests(unittest.TestCase):
    def test_ineligible_rows_are_retained_in_exclusion_ledger(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            source = root / "metadata.tsv"
            output = root / "eligible.tsv"
            fields = ["Name", "Sample Source ID", "Organism", "Body Site", "Study name", "Study condition"]
            rows = [
                ["good", "source-good", "Homo sapiens", "feces", "Study", "Control"],
                ["wrong-body", "source-body", "Homo sapiens", "blood", "Study", "CRC"],
                ["wrong-condition", "source-condition", "Homo sapiens", "feces", "Study", "Other"],
            ]
            with source.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(fields)
                writer.writerows(rows)
            subprocess.run([
                "python3", str(ROOT / "datasets/build_crc_eligible_manifest.py"),
                "--metadata", str(source), "--output", str(output),
                "--study-name", "Study",
            ], check=True)
            with output.open(newline="", encoding="utf-8") as handle:
                eligible = list(csv.DictReader(handle, delimiter="\t"))
            with (root / "eligible.excluded.tsv").open(newline="", encoding="utf-8") as handle:
                excluded = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["Name"] for row in eligible], ["good"])
            reasons = {row["Name"]: row["metadata_exclusion_reason"] for row in excluded}
            self.assertEqual(reasons["wrong-body"], "not_fecal")
            self.assertEqual(reasons["wrong-condition"], "unsupported_condition")
            self.assertTrue((root / "eligible.provenance.tsv").is_file())
            self.assertTrue((root / "eligible.tsv.sha256").is_file())


if __name__ == "__main__":
    unittest.main()
