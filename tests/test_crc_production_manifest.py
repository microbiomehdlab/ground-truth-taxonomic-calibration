import csv
import pathlib
import subprocess
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]


class ProductionManifestTests(unittest.TestCase):
    def test_multiple_runs_are_preserved_in_stable_order(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            eligible = root / "eligible.tsv"
            inventory = root / "wgets.sh"
            ena = root / "ena.tsv"
            output = root / "production.tsv"
            with eligible.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["Name", "Study condition", "Study name"])
                writer.writerow(["sample", "Control", "Study"])
            selection = root / "selection.tsv"
            with selection.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["Name", "Study condition", "selection_rank",
                                 "selection_hash", "selection_seed"])
                writer.writerow(["sample", "Control", "1", "frozen-hash", "frozen-seed"])
            inventory.write_text(
                "# Sample: sample\n"
                "wget ftp://ftp.sra.ebi.ac.uk/x/ERR2_1.fastq.gz\n"
                "wget ftp://ftp.sra.ebi.ac.uk/x/ERR2_2.fastq.gz\n"
                "wget ftp://ftp.sra.ebi.ac.uk/x/ERR1_1.fastq.gz\n"
                "wget ftp://ftp.sra.ebi.ac.uk/x/ERR1_2.fastq.gz\n",
                encoding="utf-8",
            )
            md5_a = "a" * 32
            md5_b = "b" * 32
            with ena.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["run_accession", "sample_accession", "secondary_sample_accession",
                                 "fastq_ftp", "fastq_md5", "fastq_bytes"])
                for run in ("ERR1", "ERR2"):
                    writer.writerow([run, "SAMEA", "ERS", f"ftp.sra.ebi.ac.uk/x/{run}_1.fastq.gz;ftp.sra.ebi.ac.uk/x/{run}_2.fastq.gz",
                                     f"{md5_a};{md5_b}", "10;20"])
            subprocess.run([
                "python3", str(ROOT / "datasets/build_crc_production_manifest.py"),
                "--eligible-manifest", str(eligible), "--url-inventory", str(inventory),
                "--ena-report", str(ena), "--output", str(output),
                "--independent-selection", str(selection),
                "--ena-study-accession", "PRJTEST",
                "--ena-report-retrieved-date", "2026-08-26",
            ], check=True)
            with output.open(newline="", encoding="utf-8") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["run_count"], "2")
            self.assertEqual(row["run_accessions"], "ERR1;ERR2")
            self.assertEqual(row["total_download_bytes"], "60")
            self.assertEqual(row["run_combination_rule"],
                             "concatenate mates separately in ascending run_accession order")
            independent = output.with_name("production.independent.tsv")
            with independent.open(newline="", encoding="utf-8") as handle:
                independent_row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(independent_row["sample_id"], "sample")
            self.assertEqual(independent_row["independent_subset"], "1")
            self.assertTrue(independent.with_suffix(".tsv.sha256").is_file())


if __name__ == "__main__":
    unittest.main()
