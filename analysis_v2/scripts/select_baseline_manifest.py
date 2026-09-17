#!/usr/bin/env python3
"""Select unspiked, sample-wide MetaPhlAn baselines from a canonical v2 table.

Header-aware and fail-closed, replacing positional `awk` field indexing, which
breaks silently if the canonical column order changes.

Sample-wide, not per-population. `build_crc_cohort_canonical_input.py` writes
the *same physical* zero-dose MetaPhlAn baseline for both analytical
populations: community and independent rows share `profile_id == sample_id`,
`baseline_profile_id == sample_id`, and one `source_profile`, differing only in
`analysis_population`. G_eff is therefore a property of that single physical
baseline, so rows are collapsed on `(cohort, sample_id)` and the selected
manifest leaves `analysis_population` **empty** on purpose: that blank value is
what makes `derive_paired_endpoints.py` fall back to the generic per-sample
lookup and serve both populations from one value.

Selection keeps rows with the requested cohort, `profiler == metaphlan4`,
`spike_fraction_total == 0`, `spike_fraction_target == 0`,
`profile_id == baseline_profile_id`, and `include == 1`.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

REQUIRED = ["cohort", "study", "sample_id", "condition", "analysis_population",
            "target_label", "profiler", "profile_id", "baseline_profile_id",
            "spike_fraction_total", "spike_fraction_target", "source_profile",
            "include", "exclusion_reason"]
OUTPUT = ["cohort", "sample_id", "analysis_population", "profiler", "profile_id",
          "baseline_profile_id", "spike_fraction_total", "source_profile",
          "include", "exclusion_reason"]
# Fields describing the one physical baseline. analysis_population is
# deliberately absent: it is expected to differ between populations.
PHYSICAL_FIELDS = ["cohort", "sample_id", "profiler", "profile_id",
                   "baseline_profile_id", "source_profile",
                   "spike_fraction_total", "spike_fraction_target", "include"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--profiler", default="metaphlan4")
    parser.add_argument("--cohort", help="Optional cohort filter.")
    parser.add_argument("--expected-profiles", type=int,
                        help="Withhold SUCCESS unless exactly N sample-wide "
                             "baselines are selected.")
    args = parser.parse_args()

    with args.canonical.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"[FAIL] {args.canonical}: no rows")
    missing = [f for f in REQUIRED if f not in rows[0]]
    if missing:
        raise SystemExit(f"[FAIL] {args.canonical}: missing column(s) {', '.join(missing)}")

    # Collision-safe: keyed by cohort and sample, never by bare profile_id,
    # because profile_id is the sample id and repeats across cohorts.
    selected: dict[tuple[str, str], dict[str, str]] = {}
    seen_populations: dict[tuple[str, str], set[str]] = {}
    physical_seen: dict[tuple[str, str], dict[str, str]] = {}
    exclusion_ledger: list[dict[str, str]] = []
    counts = {"total": len(rows), "profiler_mismatch": 0, "positive_dose": 0,
              "cohort_mismatch": 0}

    for index, row in enumerate(rows, start=2):
        # include=0 must always carry a reason, even when the row is dropped.
        if row["include"].strip() == "0":
            if not row["exclusion_reason"].strip():
                raise SystemExit(
                    f"[FAIL] line {index}: include=0 requires a non-empty exclusion_reason")
            exclusion_ledger.append({
                "line": str(index), "cohort": row["cohort"], "sample_id": row["sample_id"],
                "profile_id": row["profile_id"], "analysis_population": row["analysis_population"],
                "profiler": row["profiler"], "reason": row["exclusion_reason"]})
            continue
        if row["include"].strip() != "1":
            raise SystemExit(f"[FAIL] line {index}: include must be 0 or 1")
        if args.cohort and row["cohort"] != args.cohort:
            counts["cohort_mismatch"] += 1
            continue
        if row["profiler"] != args.profiler:
            counts["profiler_mismatch"] += 1
            continue
        try:
            total = float(row["spike_fraction_total"])
            target = float(row["spike_fraction_target"])
        except ValueError:
            raise SystemExit(f"[FAIL] line {index}: spike fraction is not numeric") from None
        if total != 0 or target != 0 or row["profile_id"] != row["baseline_profile_id"]:
            counts["positive_dose"] += 1
            continue

        key = (row["cohort"], row["sample_id"])
        seen_populations.setdefault(key, set()).add(row["analysis_population"])
        physical = {field: row[field] for field in PHYSICAL_FIELDS}
        if key not in selected:
            record = {field: row[field] for field in OUTPUT}
            record["analysis_population"] = ""      # sample-wide, not per-population
            selected[key] = record
            physical_seen[key] = physical
            continue
        # Repetition across targets and across populations is expected; every
        # physical field must still agree or the manifest is inconsistent.
        for field, value in physical.items():
            if physical_seen[key][field] != value:
                raise SystemExit(
                    f"[FAIL] line {index}: repeated baseline rows for "
                    f"{key[0]}/{key[1]} disagree on {field}: "
                    f"{physical_seen[key][field]!r} vs {value!r}")

    if not selected:
        raise SystemExit(
            f"[FAIL] no {args.profiler} zero-dose included rows selected from "
            f"{args.canonical}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    manifest = args.outdir / "baseline_manifest.tsv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for key in sorted(selected):
            writer.writerow(selected[key])

    ledger = args.outdir / "excluded_canonical_rows.tsv"
    with ledger.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["line", "cohort", "sample_id", "profile_id",
                                "analysis_population", "profiler", "reason"],
            delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(exclusion_ledger)

    both = sum(1 for pops in seen_populations.values() if len(pops) > 1)
    ok = args.expected_profiles is None or len(selected) == args.expected_profiles
    audit = args.outdir / "baseline_manifest_audit.tsv"
    audit.write_text("\n".join([
        "metric\tvalue",
        f"canonical_rows\t{counts['total']}",
        f"selected_baseline_samples\t{len(selected)}",
        f"samples_in_both_populations\t{both}",
        f"dropped_other_profiler\t{counts['profiler_mismatch']}",
        f"dropped_positive_dose\t{counts['positive_dose']}",
        f"dropped_other_cohort\t{counts['cohort_mismatch']}",
        f"dropped_include_0\t{len(exclusion_ledger)}",
        f"profiler\t{args.profiler}",
        f"expected_profiles\t{args.expected_profiles if args.expected_profiles is not None else 'not_enforced'}",
        f"status\t{'PASS' if ok else 'FAIL_EXPECTED_PROFILE_COUNT'}",
    ]) + "\n", encoding="utf-8")

    (args.outdir / "baseline_manifest.sha256").write_text(
        "".join(f"{sha256(p)}  {p.resolve()}\n"
                for p in (args.canonical, manifest, ledger, audit)), encoding="utf-8")

    if not ok:
        raise SystemExit(
            f"[FAIL] expected {args.expected_profiles} sample-wide baselines, "
            f"selected {len(selected)}")
    (args.outdir / "SUCCESS").write_text(
        f"selected_baseline_samples\t{len(selected)}\nstatus\tPASS\n", encoding="utf-8")
    print(f"[PASS] selected {len(selected)} sample-wide {args.profiler} baseline(s)")
    print(f"[INFO] {manifest}")


if __name__ == "__main__":
    main()
