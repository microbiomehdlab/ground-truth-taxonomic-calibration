#!/usr/bin/env python3
"""Allocate an integer community total deterministically and exactly."""
from __future__ import annotations

import argparse
from decimal import Decimal


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--total", required=True, type=int)
    parser.add_argument("--member", action="append", required=True, help="LABEL:WEIGHT; repeat for each member")
    args = parser.parse_args()
    if args.total < 0:
        parser.error("--total cannot be negative")
    members = []
    for index, item in enumerate(args.member):
        try:
            label, raw_weight = item.rsplit(":", 1)
            weight = Decimal(raw_weight)
        except Exception as error:
            parser.error(f"invalid --member {item!r}: {error}")
        if not label or weight <= 0:
            parser.error(f"invalid --member {item!r}")
        members.append((index, label, weight))
    total_weight = sum((member[2] for member in members), Decimal(0))
    shares = [(index, label, Decimal(args.total) * weight / total_weight) for index, label, weight in members]
    # Preserve the published workflow's allocation whenever it is valid:
    # round all but the final member and assign the exact remainder to it.
    allocations = {index: int(share.to_integral_value()) for index, _label, share in shares[:-1]}
    allocations[shares[-1][0]] = args.total - sum(allocations.values())
    # Very small or highly uneven totals can make that final remainder negative.
    # Such cases previously failed downstream; use largest remainder only there.
    if allocations[shares[-1][0]] < 0:
        allocations = {index: int(share) for index, _label, share in shares}
        remaining = args.total - sum(allocations.values())
        order = sorted(shares, key=lambda row: (-(row[2] - int(row[2])), row[0]))
        for index, _label, _share in order[:remaining]:
            allocations[index] += 1
    print("label\tn")
    for index, label, _weight in members:
        print(f"{label}\t{allocations[index]}")
    if sum(allocations.values()) != args.total:
        raise SystemExit("[ERROR] internal allocation error")


if __name__ == "__main__":
    main()
