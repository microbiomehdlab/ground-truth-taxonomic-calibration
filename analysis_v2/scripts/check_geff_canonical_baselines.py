#!/usr/bin/env python3
"""Require canonical MetaPhlAn baselines to match the receipt-verified G_eff sources."""
from __future__ import annotations
import argparse
from pathlib import Path
from audit_bracken_denominators import table, digest, require


def check(canonical, selection, cohort, expected):
    manifest = table(selection / "baseline_manifest.tsv")
    sources = {r["sample_id"]: Path(r["source_profile"]).resolve() for r in manifest}
    require(len(sources) == len(manifest) == expected, "Selection count mismatch")
    require(all(r["cohort"] == cohort for r in manifest), "Wrong selection cohort")
    for item in table(selection / "input_hashes.tsv"):
        require(digest(Path(item["path"])) == item["sha256"], "G_eff source changed: " + item["path"])
    seen = set()
    for row in table(canonical):
        require(row["cohort"] == cohort, "Wrong canonical cohort")
        if row["profiler"] != "metaphlan4" or float(row["spike_fraction_total"]) != 0:
            continue
        require(row["include"] == "1", "Excluded canonical baseline")
        sid = row["sample_id"]
        require(sid in sources and Path(row["source_profile"]).resolve() == sources[sid],
                "Canonical/G_eff baseline mismatch: " + sid)
        seen.add(sid)
    require(seen == set(sources), "Canonical baseline coverage mismatch")
    print("[PASS] Canonical and receipt-verified G_eff baselines agree:", cohort, len(seen))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical", type=Path, required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--expected-samples", type=int, required=True)
    args = parser.parse_args()
    try:
        check(args.canonical, args.selection, args.cohort, args.expected_samples)
    except (OSError, ValueError) as error:
        raise SystemExit("[ERROR] " + str(error))
