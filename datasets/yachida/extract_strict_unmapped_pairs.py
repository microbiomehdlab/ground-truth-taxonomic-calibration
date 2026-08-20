#!/usr/bin/env python3
"""Write paired FASTQ only when both primary SAM records are unmapped."""

from __future__ import annotations

import argparse
import gzip
import sys


def emit(records, r1, r2):
    primary = [row for row in records if not (row[1] & (0x100 | 0x800))]
    mates = {}
    for row in primary:
        if row[1] & 0x40:
            mates[1] = row
        elif row[1] & 0x80:
            mates[2] = row
    if set(mates) != {1, 2}:
        raise SystemExit(f"[ERROR] Expected one primary record per mate for {records[0][0]}")
    if not all((mates[mate][1] & 0x4) and (mates[mate][1] & 0x8) for mate in (1, 2)):
        return 0
    for mate, output in ((1, r1), (2, r2)):
        qname, _flag, sequence, quality = mates[mate]
        if sequence == "*" or quality == "*":
            raise SystemExit(f"[ERROR] Missing sequence/quality for unmapped record {qname}/{mate}")
        output.write(f"@{qname}/{mate}\n{sequence}\n+\n{quality}\n")
    return 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    args = parser.parse_args()
    current_name = None
    records = []
    kept = 0
    observed = 0
    with gzip.open(args.r1, "wt", encoding="ascii", newline="\n") as r1, gzip.open(
        args.r2, "wt", encoding="ascii", newline="\n"
    ) as r2:
        for line in sys.stdin:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise SystemExit("[ERROR] Malformed SAM record")
            qname = fields[0]
            row = (qname, int(fields[1]), fields[9], fields[10])
            if current_name is not None and qname != current_name:
                kept += emit(records, r1, r2)
                observed += 1
                records = []
            current_name = qname
            records.append(row)
        if records:
            kept += emit(records, r1, r2)
            observed += 1
    print(f"[PASS] Strict extraction retained {kept} of {observed} input pairs", file=sys.stderr)


if __name__ == "__main__":
    main()
