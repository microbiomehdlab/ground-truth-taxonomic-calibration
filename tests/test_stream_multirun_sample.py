import csv
import gzip
import hashlib
import os
import pathlib
import subprocess
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]


def md5(path: pathlib.Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


class MultiRunStreamingTests(unittest.TestCase):
    def test_verified_runs_are_assembled_once_in_accession_order(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            source = root / "source"
            scratch = root / "scratch"
            state = root / "state"
            retained = root / "retained"
            source.mkdir()
            retained.mkdir()
            files = {}
            for run, marker in (("ERR2", b"second"), ("ERR1", b"first")):
                for mate in (1, 2):
                    path = source / f"{run}_{mate}.fastq.gz"
                    with gzip.open(path, "wb") as handle:
                        handle.write(marker + f"-mate{mate}\n".encode())
                    files[(run, mate)] = path

            manifest = root / "manifest.tsv"
            fields = ["sample_id", "condition", "study", "run_count", "run_accessions",
                      "fastq1_urls", "fastq2_urls", "fastq1_md5s", "fastq2_md5s",
                      "fastq1_bytes", "fastq2_bytes"]
            ordered = ("ERR1", "ERR2")
            row = {
                "sample_id": "sample", "condition": "Control", "study": "Study", "run_count": "2",
                "run_accessions": ";".join(ordered),
            }
            for mate in (1, 2):
                row[f"fastq{mate}_urls"] = ";".join(files[(run, mate)].as_uri() for run in ordered)
                row[f"fastq{mate}_md5s"] = ";".join(md5(files[(run, mate)]) for run in ordered)
                row[f"fastq{mate}_bytes"] = ";".join(str(files[(run, mate)].stat().st_size) for run in ordered)
            with manifest.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader(); writer.writerow(row)

            runner = root / "runner.py"
            runner.write_text(
                "#!/usr/bin/env python3\n"
                "import gzip, hashlib, os, pathlib\n"
                "out=pathlib.Path(os.environ['RETAINED'])/'assembled.txt'\n"
                "with gzip.open(os.environ['RAW_R1'],'rb') as h: out.write_bytes(h.read())\n"
                "d=hashlib.sha256(out.read_bytes()).hexdigest()\n"
                "r=pathlib.Path(os.environ['RECEIPT']); r.write_text('path\\tsha256\\tbytes\\n'+f'{out}\\t{d}\\t{out.stat().st_size}\\n')\n",
                encoding="utf-8",
            )
            runner.chmod(0o755)
            environment = {**os.environ, "RETAINED": str(retained), "CRC_DOWNLOAD_ATTEMPTS": "1"}
            subprocess.run([
                "python3", str(ROOT / "datasets/stream_multirun_sample.py"),
                "--manifest", str(manifest), "--sample-id", "sample",
                "--scratch-root", str(scratch), "--state-dir", str(state),
                "--runner", str(runner),
            ], check=True, env=environment)
            self.assertEqual((retained / "assembled.txt").read_bytes(), b"first-mate1\nsecond-mate1\n")
            provenance = (state / "sample.input_provenance.tsv").read_text(encoding="utf-8")
            self.assertLess(provenance.index("ERR1"), provenance.index("ERR2"))
            self.assertTrue((state / "sample.input_provenance.tsv.sha256").is_file())
            self.assertTrue((state / "sample.verified").is_file())


if __name__ == "__main__":
    unittest.main()
