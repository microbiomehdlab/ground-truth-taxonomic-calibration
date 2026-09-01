#!/usr/bin/env python3
"""Select synchronized pairs at deterministic midpoints spanning a full FASTQ pool."""
from __future__ import annotations

import argparse
import gzip
import io
from pathlib import Path


def record(handle, path: Path, index: int):
    lines = [handle.readline() for _ in range(4)]
    if not lines[0]:
        return None
    if any(line == "" for line in lines[1:]):
        raise SystemExit(f"[ERROR] Truncated FASTQ record {index + 1} in {path}")
    if not lines[0].startswith("@") or not lines[2].startswith("+"):
        raise SystemExit(f"[ERROR] Invalid FASTQ record {index + 1} in {path}")
    if len(lines[1].rstrip("\r\n")) != len(lines[3].rstrip("\r\n")):
        raise SystemExit(f"[ERROR] Sequence/quality mismatch at record {index + 1} in {path}")
    return lines


def read_id(header: str) -> str:
    value = header[1:].split(maxsplit=1)[0]
    for suffix in ("/1", "/2", ".1", ".2", "_1", "_2"):
        if value.endswith(suffix):
            return value[: -len(suffix)]
    return value


def gzip_text(path: Path):
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    return raw, compressed, io.TextIOWrapper(compressed, encoding="utf-8", newline="")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--output-r1", required=True, type=Path)
    parser.add_argument("--output-r2", required=True, type=Path)
    parser.add_argument("--total-pairs", required=True, type=int)
    parser.add_argument("--sample-pairs", required=True, type=int)
    args = parser.parse_args()
    if not (1 <= args.sample_pairs <= args.total_pairs):
        parser.error("require 1 <= --sample-pairs <= --total-pairs")

    args.output_r1.parent.mkdir(parents=True, exist_ok=True)
    args.output_r2.parent.mkdir(parents=True, exist_ok=True)
    raw1, gz1, out1 = gzip_text(args.output_r1)
    raw2, gz2, out2 = gzip_text(args.output_r2)
    selected = 0
    next_index = ((2 * selected + 1) * args.total_pairs) // (2 * args.sample_pairs)
    try:
        with args.r1.open(encoding="utf-8") as left, args.r2.open(encoding="utf-8") as right:
            for index in range(args.total_pairs):
                rec1 = record(left, args.r1, index)
                rec2 = record(right, args.r2, index)
                if rec1 is None or rec2 is None:
                    raise SystemExit(f"[ERROR] Pool ended before declared pair count at pair {index + 1}")
                if read_id(rec1[0]) != read_id(rec2[0]):
                    raise SystemExit(f"[ERROR] Mate identifiers differ at pair {index + 1}")
                if index == next_index:
                    out1.writelines(rec1); out2.writelines(rec2)
                    selected += 1
                    if selected < args.sample_pairs:
                        next_index = ((2 * selected + 1) * args.total_pairs) // (2 * args.sample_pairs)
            if record(left, args.r1, args.total_pairs) is not None or record(right, args.r2, args.total_pairs) is not None:
                raise SystemExit("[ERROR] Pool contains more pairs than declared")
    finally:
        out1.close(); out2.close()
        gz1.close(); gz2.close()
        raw1.close(); raw2.close()
    if selected != args.sample_pairs:
        raise SystemExit(f"[ERROR] Selected {selected}, expected {args.sample_pairs}")
    print(f"[PASS] Systematically selected {selected} of {args.total_pairs} synchronized pairs")


if __name__ == "__main__":
    main()
