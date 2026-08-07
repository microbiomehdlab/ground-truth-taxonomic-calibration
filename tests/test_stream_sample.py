#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import pathlib
import subprocess
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "datasets/yachida/stream_sample.py"


class StreamSampleResumeTests(unittest.TestCase):
    def make_verified_case(self, root: pathlib.Path, *, sentinel: bool = True) -> list[str]:
        sample = "SAMPLE1"
        scratch = root / "scratch"
        work = scratch / sample
        state = root / "state"
        retained = root / "persistent" / "result.tsv"
        work.mkdir(parents=True)
        state.mkdir()
        retained.parent.mkdir()
        retained.write_text("result\n", encoding="utf-8")
        if sentinel:
            (work / ".ground_truth_stream_sample").write_text("managed\n", encoding="utf-8")
        (work / "large.fastq.gz").write_text("disposable\n", encoding="utf-8")
        manifest = root / "manifest.tsv"
        manifest.write_text(
            "sample_id\tfastq1_url\tfastq1_md5\tfastq1_bytes\tfastq2_url\tfastq2_md5\tfastq2_bytes\tTarget_Condition\tStudy\n"
            f"{sample}\tunused\tunused\t1\tunused\tunused\t1\tControl\tYachidaS_2019\n",
            encoding="utf-8",
        )
        runner = root / "runner.sh"
        runner.write_text("#!/usr/bin/env bash\nexit 99\n", encoding="utf-8")
        runner.chmod(0o755)
        digest = hashlib.sha256(retained.read_bytes()).hexdigest()
        (state / f"{sample}.retained_outputs.tsv").write_text(
            f"path\tsha256\tbytes\n{retained}\t{digest}\t{retained.stat().st_size}\n",
            encoding="utf-8",
        )
        (state / f"{sample}.verified").write_text("outputs\t1\n", encoding="utf-8")
        return [
            "python3", str(SCRIPT), "--manifest", str(manifest), "--sample-id", sample,
            "--scratch-root", str(scratch), "--state-dir", str(state), "--runner", str(runner),
            "--delete-inputs-after-verification",
        ]

    def test_verified_resume_only_deletes_sentinel_protected_scratch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            process = subprocess.run(self.make_verified_case(root), text=True, capture_output=True)
            self.assertEqual(process.returncode, 0, process.stderr)
            self.assertFalse((root / "scratch/SAMPLE1").exists())
            self.assertTrue((root / "persistent/result.tsv").is_file())
            self.assertIn("already verified", process.stdout)
            self.assertIn("deleted previously verified", process.stdout)

    def test_verified_resume_refuses_cleanup_without_sentinel(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            process = subprocess.run(self.make_verified_case(root, sentinel=False), text=True, capture_output=True)
            self.assertNotEqual(process.returncode, 0)
            self.assertTrue((root / "scratch/SAMPLE1").is_dir())
            self.assertIn("cleanup sentinel missing", process.stderr)


if __name__ == "__main__":
    unittest.main()
