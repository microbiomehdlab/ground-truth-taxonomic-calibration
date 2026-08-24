#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import http.client
import importlib.util
import io
import pathlib
import subprocess
import tempfile
import unittest
from unittest import mock


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "datasets/yachida/stream_sample.py"
SPEC = importlib.util.spec_from_file_location("yachida_stream_sample", SCRIPT)
STREAM_SAMPLE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(STREAM_SAMPLE)


class FakeResponse(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class StreamSampleResumeTests(unittest.TestCase):
    def test_download_retries_disconnect_and_checksum_failure(self) -> None:
        payload = b"verified payload"
        expected_md5 = hashlib.md5(payload).hexdigest()
        with tempfile.TemporaryDirectory() as directory, mock.patch.dict(
            STREAM_SAMPLE.os.environ,
            {
                "YACHIDA_DOWNLOAD_ATTEMPTS": "4",
                "YACHIDA_DOWNLOAD_RETRY_SECONDS": "1",
                "YACHIDA_DOWNLOAD_TIMEOUT_SECONDS": "9",
            },
        ), mock.patch.object(STREAM_SAMPLE.time, "sleep") as sleep, mock.patch.object(
            STREAM_SAMPLE.urllib.request,
            "urlopen",
            side_effect=[
                http.client.RemoteDisconnected("transient"),
                FakeResponse(b"corrupt"),
                FakeResponse(payload),
            ],
        ) as urlopen:
            destination = pathlib.Path(directory) / "reads.fastq.gz"
            STREAM_SAMPLE.download("https://example.invalid/reads", destination, expected_md5, len(payload))
            self.assertEqual(destination.read_bytes(), payload)
            self.assertFalse(destination.with_suffix(".gz.partial").exists())
            self.assertEqual(urlopen.call_count, 3)
            self.assertEqual(sleep.call_args_list, [mock.call(1), mock.call(2)])
            self.assertTrue(all(call.kwargs["timeout"] == 9 for call in urlopen.call_args_list))

    def test_download_fails_closed_after_bounded_attempts(self) -> None:
        with tempfile.TemporaryDirectory() as directory, mock.patch.dict(
            STREAM_SAMPLE.os.environ,
            {
                "YACHIDA_DOWNLOAD_ATTEMPTS": "2",
                "YACHIDA_DOWNLOAD_RETRY_SECONDS": "1",
            },
        ), mock.patch.object(STREAM_SAMPLE.time, "sleep"), mock.patch.object(
            STREAM_SAMPLE.urllib.request,
            "urlopen",
            side_effect=[FakeResponse(b"bad"), FakeResponse(b"still bad")],
        ):
            destination = pathlib.Path(directory) / "reads.fastq.gz"
            with self.assertRaisesRegex(SystemExit, "download failed after 2 attempts"):
                STREAM_SAMPLE.download("https://example.invalid/reads", destination, "0" * 32, 10)
            self.assertFalse(destination.exists())
            self.assertFalse(destination.with_suffix(".gz.partial").exists())

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
