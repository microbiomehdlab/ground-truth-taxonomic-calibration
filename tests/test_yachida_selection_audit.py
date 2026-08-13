import csv
import pathlib
import subprocess
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]


class SelectionAuditTests(unittest.TestCase):
    def test_selection_audit_namespaces_both_selection_stages(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            pilot = root / "pilot.tsv"
            nested = root / "nested.tsv"
            fields = ["sample_id", "Target_Condition", "age", "sex", "matching_age_bin", "matching_stratum",
                      "selection_rank", "selection_hash", "selection_seed", "matching_reference"]
            rows = []
            for condition in ("Control", "Adenoma", "CRC"):
                rows.append([condition, condition, "60", "Male", "60-69", "Male|60-69", "1", "pilot", "pilot-seed", "Adenoma"])
            with pilot.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(fields); writer.writerows(rows)
            with nested.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["sample_id", "Target_Condition", "selection_rank", "selection_hash", "selection_seed"])
                for condition in ("Control", "Adenoma", "CRC"):
                    writer.writerow([condition, condition, "1", "nested", "nested-seed"])
            output = root / "audit"
            subprocess.run(["python3", str(ROOT / "datasets/yachida/audit_selection_balance.py"),
                            "--pilot-manifest", str(pilot), "--independent-manifest", str(nested),
                            "--output-dir", str(output)], check=True)
            header = (output / "independent_selection_provenance.tsv").read_text().splitlines()[0]
            self.assertIn("pilot_selection_seed", header)
            self.assertIn("independent_selection_seed", header)


if __name__ == "__main__":
    unittest.main()
