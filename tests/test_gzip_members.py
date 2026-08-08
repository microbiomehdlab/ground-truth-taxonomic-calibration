#!/usr/bin/env python3
from __future__ import annotations

import gzip
import hashlib
import pathlib
import tempfile
import unittest


class ConcatenatedGzipTests(unittest.TestCase):
    def test_concatenated_members_equal_recompressed_stream(self) -> None:
        background = b"@background/1\nACGT\n+\nIIII\n"
        spike = b"@spike/1\nTGCA\n+\nIIII\n"
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            background_gz = root / "background.fq.gz"
            spike_gz = root / "spike.fq.gz"
            recompressed = root / "recompressed.fq.gz"
            concatenated = root / "concatenated.fq.gz"
            with gzip.open(background_gz, "wb") as handle:
                handle.write(background)
            with gzip.open(spike_gz, "wb") as handle:
                handle.write(spike)
            with gzip.open(recompressed, "wb") as handle:
                handle.write(background + spike)
            concatenated.write_bytes(background_gz.read_bytes() + spike_gz.read_bytes())

            with gzip.open(recompressed, "rb") as handle:
                left = handle.read()
            with gzip.open(concatenated, "rb") as handle:
                right = handle.read()
            self.assertEqual(left, background + spike)
            self.assertEqual(right, left)
            self.assertEqual(hashlib.sha256(right).hexdigest(), hashlib.sha256(left).hexdigest())


if __name__ == "__main__":
    unittest.main()
