#!/usr/bin/env python3
"""Check a runtime spike panel against the frozen panel, allowing path spelling only."""
import argparse
import csv
import sys
from pathlib import Path


def read_panel(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        if fields is None or len(fields) != len(set(fields)) or "fasta" not in fields:
            raise ValueError(f"Invalid spike panel header: {path}")
        rows = list(reader)
    if not rows or any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"Invalid spike panel rows: {path}")
    return fields, rows


def check(runtime, frozen, root):
    runtime_fields, runtime_rows = read_panel(runtime)
    frozen_fields, frozen_rows = read_panel(frozen)
    if runtime_fields != frozen_fields or len(runtime_rows) != len(frozen_rows):
        raise ValueError("Spike panel schema or row count differs from frozen identity policy")
    for index, (actual, expected) in enumerate(zip(runtime_rows, frozen_rows), 1):
        for field in runtime_fields:
            if field == "fasta":
                actual_path = Path(actual[field])
                expected_path = Path(expected[field])
                actual_path = (actual_path if actual_path.is_absolute() else root / actual_path).resolve()
                expected_path = (expected_path if expected_path.is_absolute() else root / expected_path).resolve()
                if actual_path != expected_path:
                    raise ValueError(f"Spike panel row {index} FASTA differs from frozen identity policy")
            elif actual[field] != expected[field]:
                raise ValueError(f"Spike panel row {index} field {field} differs from frozen identity policy")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime", type=Path, required=True)
    parser.add_argument("--frozen", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    try:
        check(args.runtime, args.frozen, args.root.resolve())
    except (OSError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    print("[PASS] Runtime spike panel matches frozen identity policy")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
