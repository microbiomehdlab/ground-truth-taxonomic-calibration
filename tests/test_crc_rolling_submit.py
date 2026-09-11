import os
import pathlib
import subprocess
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]


class RollingSubmitTests(unittest.TestCase):
    def test_builds_bounded_download_generations_and_correlated_lanes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            manifest = root / "production.tsv"
            manifest.write_text(
                "sample_id\n" + "".join(f"sample-{index}\n" for index in range(1, 101)),
                encoding="utf-8",
            )
            crc_env = root / "crc.env"
            crc_env.write_text("CRC_SCRATCH_ROOT=/tmp/scratch\n", encoding="utf-8")
            binary = root / "bin"
            binary.mkdir()
            calls = root / "calls.txt"
            counter = root / "counter.txt"
            fake_sbatch = binary / "sbatch"
            fake_sbatch.write_text(
                "#!/usr/bin/env bash\n"
                "set -euo pipefail\n"
                f"printf '%s\\n' \"$*\" >> {calls!s}\n"
                f"n=$(cat {counter!s} 2>/dev/null || echo 7000)\n"
                "n=$((n + 1))\n"
                f"printf '%s\\n' \"$n\" > {counter!s}\n"
                "printf '%s\\n' \"$n\"\n",
                encoding="utf-8",
            )
            fake_sbatch.chmod(0o755)
            environment = {
                **os.environ,
                "PATH": f"{binary}:{os.environ['PATH']}",
                "PROJECT": str(ROOT),
                "CRC_ENV": str(crc_env),
                "CRC_ROLLING_LOG_ROOT": str(root / "logs"),
            }
            completed = subprocess.run(
                [
                    "bash", str(ROOT / "datasets/submit_crc_rolling.sh"),
                    "--manifest", str(manifest), "--lanes", "42",
                    "--download-concurrent", "8", "--exclude", "compute-1,compute-10",
                    "--delete-verified-inputs",
                ],
                check=True, capture_output=True, text=True, env=environment,
            )
            self.assertIn("At most 8 downloads and 42 compute lanes", completed.stdout)
            records = calls.read_text(encoding="utf-8").splitlines()
            self.assertEqual(len(records), 6)
            self.assertIn("--array=1-42%8", records[0])
            self.assertIn("--dependency=aftercorr:7001", records[1])
            self.assertIn("--dependency=afterany:7001,aftercorr:7002", records[2])
            self.assertIn("--dependency=aftercorr:7003", records[3])
            self.assertIn("--array=1-16%8", records[4])
            self.assertIn("--exclude=compute-1,compute-10", records[5])


if __name__ == "__main__":
    unittest.main()
