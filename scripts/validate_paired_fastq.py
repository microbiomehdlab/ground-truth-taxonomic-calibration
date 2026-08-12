#!/usr/bin/env python3
"""Fail unless two FASTQ files contain valid, synchronized read pairs."""
from __future__ import annotations

import argparse
import gzip
import pathlib
from itertools import zip_longest


def open_text(path: pathlib.Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="strict") if path.suffix == ".gz" else path.open(encoding="utf-8")


def records(path: pathlib.Path):
    with open_text(path) as handle:
        while True:
            lines = [handle.readline() for _ in range(4)]
            if not lines[0]:
                if any(lines[1:]):
                    raise ValueError(f"truncated FASTQ record in {path}")
                return
            if any(line == "" for line in lines[1:]):
                raise ValueError(f"truncated FASTQ record in {path}")
            header, sequence, plus, quality = (line.rstrip("\r\n") for line in lines)
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"invalid FASTQ structure in {path}: {header!r}")
            if len(sequence) != len(quality):
                raise ValueError(f"sequence/quality length mismatch in {path}: {header!r}")
            yield header, sequence, quality


def normalized_id(header: str) -> str:
    token = header[1:].split(maxsplit=1)[0]
    for suffix in ("/1", "/2", ".1", ".2", "_1", "_2"):
        if token.endswith(suffix):
            return token[: -len(suffix)]
    return token


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", required=True, type=pathlib.Path)
    parser.add_argument("--r2", required=True, type=pathlib.Path)
    parser.add_argument("--output", required=True, type=pathlib.Path)
    parser.add_argument("--minimum-pairs", type=int, default=1)
    args = parser.parse_args()
    if args.minimum_pairs < 1:
        parser.error("--minimum-pairs must be positive")

    count = bases1 = bases2 = gc1 = gc2 = quality_sum1 = quality_sum2 = 0
    try:
        for index, pair in enumerate(zip_longest(records(args.r1), records(args.r2)), start=1):
            left, right = pair
            if left is None or right is None:
                raise ValueError(f"mate record counts differ at pair {index}")
            if normalized_id(left[0]) != normalized_id(right[0]):
                raise ValueError(f"mate identifiers differ at pair {index}: {left[0]!r} versus {right[0]!r}")
            count += 1
            bases1 += len(left[1]); bases2 += len(right[1])
            gc1 += sum(base in "GCgc" for base in left[1]); gc2 += sum(base in "GCgc" for base in right[1])
            quality_sum1 += sum(ord(char) - 33 for char in left[2])
            quality_sum2 += sum(ord(char) - 33 for char in right[2])
    except (OSError, UnicodeError, ValueError) as error:
        raise SystemExit(f"[ERROR] Paired FASTQ validation failed: {error}") from error

    if count < args.minimum_pairs:
        raise SystemExit(f"[ERROR] Only {count} synchronized pairs; minimum is {args.minimum_pairs}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "metric\tvalue\n"
        f"synchronized_pairs\t{count}\n"
        f"r1_mean_length\t{bases1 / count:.6f}\n"
        f"r2_mean_length\t{bases2 / count:.6f}\n"
        f"r1_gc_pct\t{100 * gc1 / bases1:.6f}\n"
        f"r2_gc_pct\t{100 * gc2 / bases2:.6f}\n"
        f"r1_mean_phred\t{quality_sum1 / bases1:.6f}\n"
        f"r2_mean_phred\t{quality_sum2 / bases2:.6f}\n"
        "identifier_synchronization\tPASS\n"
        "fastq_structure\tPASS\n",
        encoding="utf-8",
    )
    print(f"[PASS] Validated {count} synchronized read pairs")


if __name__ == "__main__":
    main()
