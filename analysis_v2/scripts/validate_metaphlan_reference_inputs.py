#!/usr/bin/env python3
"""Read-only reference validation, using the endpoint stage's exact loaders."""
from __future__ import annotations

import argparse
import math
from pathlib import Path

from derive_paired_endpoints import (
    read_table, load_target_genome_sizes, load_effective_genome_sizes,
    effective_size,
)


def validate(canonical: Path, targets: Path, effective: Path) -> None:
    rows = read_table(canonical, ["profiler", "cohort", "sample_id",
                                  "analysis_population", "target_label"])
    sizes = load_target_genome_sizes(targets)
    geff = load_effective_genome_sizes(effective)
    for value in list(sizes.values()) + list(geff.values()):
        if not math.isfinite(value):
            raise ValueError("genome sizes must be finite")
    for row in rows:
        if row["profiler"] != "metaphlan4":
            continue
        if row["target_label"] not in sizes:
            raise ValueError("no genome size for target: " + row["target_label"])
        if effective_size(geff, row) is None:
            raise ValueError("no effective community genome size for scope: " +
                             "/".join(row[k] for k in
                                      ("cohort", "sample_id", "analysis_population")))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical", type=Path, required=True)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--effective", type=Path, required=True)
    args = parser.parse_args()
    try:
        validate(args.canonical, args.targets, args.effective)
    except (ValueError, OSError) as error:
        raise SystemExit("[ERROR] " + str(error))
    print("[PASS] Reference values and target/population coverage validated")


if __name__ == "__main__":
    main()
