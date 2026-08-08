#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import pathlib
import subprocess
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "spikes/scripts/spikein/stage_finalized_pools.sh"


class StagePoolsTests(unittest.TestCase):
    def test_stage_and_verify_finalized_pool_fixture(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            source = root / "source"
            destination = root / "destination"
            source.mkdir()
            pool_names = ["Fnuc.pool_1.fq", "Fnuc.pool_2.fq"]
            for mate, name in enumerate(pool_names, 1):
                (source / name).write_text(f"@read/{mate}\nACGT\n+\nIIII\n", encoding="utf-8")
            with (source / "pool_files.sha256").open("w", encoding="utf-8") as handle:
                for name in pool_names:
                    digest = hashlib.sha256((source / name).read_bytes()).hexdigest()
                    handle.write(f"{digest}  {name}\n")
            pair_index = source / "pool_pair_counts.tsv"
            pair_index.write_text(
                "label\tmate1_pairs\tmate2_pairs\tpool_pairs\nFnuc\t1\t1\t1\n",
                encoding="utf-8",
            )
            digest = hashlib.sha256(pair_index.read_bytes()).hexdigest()
            # Older indexes recorded an absolute source path. Staging must
            # verify the digest itself and rewrite a destination-local record.
            (source / "pool_pair_counts.tsv.sha256").write_text(
                f"{digest}  {pair_index}\n", encoding="utf-8"
            )

            process = subprocess.run(
                [str(SCRIPT), "--source", str(source), "--destination", str(destination)],
                text=True,
                capture_output=True,
            )
            self.assertEqual(process.returncode, 0, process.stderr)
            self.assertTrue((destination / "STAGED_SUCCESS").is_file())
            self.assertTrue(
                (destination / "pool_pair_counts.tsv.sha256").read_text(encoding="utf-8").endswith(
                    "  pool_pair_counts.tsv\n"
                )
            )
            for name in pool_names:
                self.assertEqual((destination / name).read_bytes(), (source / name).read_bytes())


if __name__ == "__main__":
    unittest.main()
