#!/usr/bin/env python3
"""Compute how many reads to add for an exact final spike fraction.

Primary formula (exogenous spike; spike species absent from background):
    N = f/(1-f) * R

If spike species is already present in background with fraction f0 (0..1), and you
want the *final* fraction to be f_target, then:
    N = (f_target - f0)/(1 - f_target) * R

This second form is useful for "endogenous boosts".

The script prints:
  - N (rounded)
  - achieved fraction after rounding

Examples:
  calc_spike_reads.py --R 500000 --f 0.001
  calc_spike_reads.py --R 500000 --f-target 0.060 --f0 0.050
"""

import argparse


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=int, required=True, help="background read count (pairs for PE, reads for SE)")

    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--f", type=float, help="target FINAL spike fraction for exogenous spike (0..1)")
    g.add_argument("--f-target", type=float, help="target FINAL fraction for endogenous boost (0..1)")

    ap.add_argument("--f0", type=float, default=0.0, help="baseline fraction of the boosted taxon (0..1); only used with --f-target")
    ap.add_argument("--round", choices=["round", "floor", "ceil"], default="round")

    return ap.parse_args()


def round_int(x: float, mode: str) -> int:
    import math

    if mode == "round":
        return int(round(x))
    if mode == "floor":
        return int(math.floor(x))
    if mode == "ceil":
        return int(math.ceil(x))
    raise ValueError(mode)


def main() -> None:
    a = parse_args()

    R = a.R
    if R <= 0:
        raise SystemExit("R must be > 0")

    if a.f is not None:
        f = a.f
        if not (0.0 < f < 1.0):
            raise SystemExit("--f must be in (0,1)")
        N_exact = (f / (1.0 - f)) * R
        N = round_int(N_exact, a.round)
        f_hat = N / (R + N)
        print(f"N\t{N}")
        print(f"f_hat\t{f_hat:.8f}")
        return

    # endogenous target
    f_target = a.f_target
    f0 = a.f0
    if not (0.0 <= f0 < 1.0):
        raise SystemExit("--f0 must be in [0,1)")
    if not (0.0 < f_target < 1.0):
        raise SystemExit("--f-target must be in (0,1)")
    if f_target <= f0:
        raise SystemExit("--f-target must be > --f0")

    N_exact = ((f_target - f0) / (1.0 - f_target)) * R
    N = round_int(N_exact, a.round)

    # achieved final fraction after rounding
    # background already contains f0 of R reads from the boosted taxon
    f_hat = (f0 * R + N) / (R + N)

    print(f"N\t{N}")
    print(f"f_hat\t{f_hat:.8f}")


if __name__ == "__main__":
    main()
