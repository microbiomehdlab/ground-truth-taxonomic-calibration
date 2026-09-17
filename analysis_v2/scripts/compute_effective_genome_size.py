#!/usr/bin/env python3
"""Derive sample-specific effective community genome size G_eff,i.

For each UNSPIKED baseline MetaPhlAn profile i:

    G_eff,i = sum_j(c_ij * G_j) / sum_j(c_ij)

over eligible mapped features j, where c_ij is the baseline relative abundance
and G_j the database representative genome length.

Rank policy
-----------
The production MetaPhlAn vJan25 taxonomy is SGB-based: database lineages
terminate in `t__SGB...` (see scripts/build_reference_representation_table.py,
metaphlan_sgb_map). Raw `.metaphlan.tsv` profiles carry the full hierarchy, so
the SGB rows are present (workflows/metaphlan4/postprocess_local.sh exists
precisely to strip them for the species-only derivative).

Default `--profile-rank sgb` therefore weights **terminal** rows at or below
species rank. Terminal means the row has no descendant row in the same profile,
so parent species/genus/phylum rows never contribute and hierarchical mass is
counted exactly once. A terminal `s__` row with no SGB child stays eligible and
is reported as unmapped rather than dropped, so coverage reflects it honestly.

`--profile-rank species` weights every `s__` row instead, matching the canonical
input's species-level estimand, and requires a species-rank mapping table.

Either way the script fails closed if the mapping table carries no identifier of
the rank the chosen mode needs.

Scope
-----
One physical baseline per biological sample yields **one** sample-wide
`G_eff,i`. The canonical builder writes the same zero-dose MetaPhlAn profile for
the community and independent populations, so rows are collapsed on
`(cohort, sample_id)` and `analysis_population` is emitted **blank** on purpose:
that is what makes `derive_paired_endpoints.py` fall back to its generic
per-sample lookup and serve both populations from the same value.

Independence from spike outcomes
--------------------------------
Only unspiked baselines may inform G_eff,i. The manifest must declare
`spike_fraction_total == 0` and `profile_id == baseline_profile_id` for every
row, and every included row must be `profiler == metaphlan4`. There is no
override flag, and nothing here reads observed post-spike abundance.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path

MANIFEST_REQUIRED = ["cohort", "sample_id", "analysis_population", "profiler",
                     "profile_id", "baseline_profile_id", "spike_fraction_total",
                     "source_profile", "include", "exclusion_reason"]
# One physical baseline per biological sample. build_crc_cohort_canonical_input
# writes the same zero-dose MetaPhlAn profile for both analytical populations
# (same profile_id, baseline_profile_id and source_profile; only
# analysis_population differs), so identity is (cohort, sample_id) and
# analysis_population is deliberately NOT an invariant.
IDENTITY_FIELDS = ["cohort", "sample_id"]
INVARIANT_FIELDS = ["cohort", "sample_id", "profiler", "profile_id",
                    "baseline_profile_id", "source_profile",
                    "spike_fraction_total", "include"]
SAMPLE_FIELDS = ["cohort", "sample_id", "profile_id", "analysis_population",
                 "effective_genome_size_bp", "mapping_coverage",
                 "eligible_abundance", "mapped_abundance", "unmapped_abundance",
                 "excluded_abundance", "n_mapped_features", "n_unmapped_features",
                 "profile_rank", "abundance_unit", "source_profile"]
REQUIRED_PROFILER = "metaphlan4"
MIN_GENOME_BP = 100_000
MAX_GENOME_BP = 50_000_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def render(value: float) -> str:
    return format(value, ".17g")


def leaf_rank(clade: str) -> str:
    leaf = clade.split("|")[-1]
    return leaf.split("__", 1)[0] if "__" in leaf else ""


def load_genome_sizes(path: Path, match_on: str, need_rank: str) -> dict[str, int]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"[FAIL] {path}: no genome-size rows")
    column = "clade_lineage" if match_on == "lineage" else "clade_leaf"
    for required in (column, "genome_size_bp"):
        if required not in rows[0]:
            raise SystemExit(f"[FAIL] {path}: missing column {required}")

    sizes: dict[str, int] = {}
    rank_counts: dict[str, int] = {}
    for index, row in enumerate(rows, start=2):
        key = row[column].strip()
        if not key:
            raise SystemExit(f"[FAIL] {path} line {index}: empty {column}")
        try:
            size = int(float(row["genome_size_bp"]))
        except ValueError:
            raise SystemExit(f"[FAIL] {path} line {index}: genome_size_bp not numeric") from None
        if not MIN_GENOME_BP <= size <= MAX_GENOME_BP:
            raise SystemExit(
                f"[FAIL] {path} line {index}: genome_size_bp {size} outside "
                f"[{MIN_GENOME_BP}, {MAX_GENOME_BP}]")
        if key in sizes and sizes[key] != size:
            raise SystemExit(
                f"[FAIL] {path}: {column} {key!r} has conflicting sizes "
                f"{sizes[key]} and {size}")
        sizes[key] = size
        rank = (row.get("rank") or leaf_rank(key)).strip()
        rank_counts[rank] = rank_counts.get(rank, 0) + 1

    # Rank compatibility: an SGB-only mapping cannot serve a species-rank
    # estimator, and vice versa. Failing here beats silently mapping nothing.
    if rank_counts.get(need_rank, 0) == 0:
        observed = ", ".join(f"{r or '?'}={n}" for r, n in sorted(rank_counts.items()))
        raise SystemExit(
            f"[FAIL] {path}: no identifiers of rank {need_rank!r}, required by the "
            f"selected --profile-rank; observed ranks: {observed}")
    return sizes


def read_manifest(path: Path, policy: str) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"[FAIL] {path}: no manifest rows")
    missing = [f for f in MANIFEST_REQUIRED if f not in rows[0]]
    if missing:
        raise SystemExit(f"[FAIL] {path}: missing column(s) {', '.join(missing)}")

    filtered: list[dict[str, str]] = []
    excluded_rows: list[dict[str, str]] = []
    for index, row in enumerate(rows, start=2):
        state = row["include"].strip()
        if state not in {"0", "1"}:
            raise SystemExit(f"[FAIL] {path} line {index}: include must be 0 or 1")
        if state == "0":
            if not row["exclusion_reason"].strip():
                raise SystemExit(
                    f"[FAIL] {path} line {index}: include=0 requires exclusion_reason")
            if policy == "require_included_only":
                raise SystemExit(
                    f"[FAIL] {path} line {index}: excluded row present but policy is "
                    "require_included_only; pre-filter the manifest or pass "
                    "--inclusion-policy filter_excluded")
            excluded_rows.append({"sample_id": row["sample_id"],
                                  "profile_id": row["profile_id"],
                                  "reason": "canonical_include_0",
                                  "detail": row["exclusion_reason"]})
            continue
        filtered.append(row)

    if not filtered:
        raise SystemExit(f"[FAIL] {path}: no included rows")

    # MetaPhlAn only. A Bracken table is a different format entirely and must
    # never be parsed as a MetaPhlAn profile.
    profilers = sorted({r["profiler"] for r in filtered})
    if profilers != [REQUIRED_PROFILER]:
        raise SystemExit(
            f"[FAIL] {path}: G_eff requires {REQUIRED_PROFILER} rows only, found "
            f"profiler(s): {', '.join(profilers)}. Select MetaPhlAn baseline rows "
            "with scripts/select_baseline_manifest.py")

    offenders = []
    for index, row in enumerate(filtered, start=2):
        try:
            dose = float(row["spike_fraction_total"])
        except ValueError:
            raise SystemExit(f"[FAIL] {path}: spike_fraction_total not numeric") from None
        if dose != 0 or row["profile_id"] != row["baseline_profile_id"]:
            offenders.append(row["profile_id"])
    if offenders:
        raise SystemExit(
            f"[FAIL] {path}: post-spike profile(s) present; only unspiked baselines "
            f"may inform G_eff: {', '.join(sorted(set(offenders))[:5])}")

    # Baselines legitimately repeat once per implanted target and once per
    # analytical population. Collapse on (cohort, sample_id) only after proving
    # every physical field agrees; a conflict means the manifest is wrong.
    grouped: dict[tuple[str, ...], list[dict[str, str]]] = {}
    for row in filtered:
        grouped.setdefault(tuple(row[f] for f in IDENTITY_FIELDS), []).append(row)
    unique: list[dict[str, str]] = []
    for identity in sorted(grouped):
        members = grouped[identity]
        first = members[0]
        for other in members[1:]:
            for field in INVARIANT_FIELDS:
                if other[field] != first[field]:
                    raise SystemExit(
                        f"[FAIL] {path}: repeated baseline rows for "
                        f"{'/'.join(identity)} disagree on {field}: "
                        f"{first[field]!r} vs {other[field]!r}")
        record = dict(first)
        record["analysis_population"] = ""   # sample-wide; blank drives the
        unique.append(record)                # generic downstream lookup
    return unique, excluded_rows


def parse_profile(path: Path) -> list[tuple[str, str, float]]:
    records = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            raise SystemExit(f"[FAIL] {path}: malformed row {line[:80]!r}")
        clade = parts[0].strip()
        try:
            value = float(parts[2])
        except ValueError:
            raise SystemExit(f"[FAIL] {path}: non-numeric abundance for {clade}") from None
        if not math.isfinite(value) or value < 0:
            raise SystemExit(f"[FAIL] {path}: invalid abundance for {clade}")
        records.append((clade, clade.split("|")[-1], value))
    if not records:
        raise SystemExit(f"[FAIL] {path}: no data rows")
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--genome-sizes", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--profile-root", type=Path)
    parser.add_argument("--profile-rank", choices=["sgb", "species"], default="sgb",
                        help="sgb (default): terminal rows at or below species rank, "
                             "matching the SGB-based vJan25 taxonomy. species: all s__ rows.")
    parser.add_argument("--match-on", choices=["lineage", "leaf"], default="lineage")
    parser.add_argument("--abundance-unit", choices=["auto", "percent", "fraction"],
                        default="auto")
    parser.add_argument("--inclusion-policy",
                        choices=["require_included_only", "filter_excluded"],
                        default="require_included_only")
    parser.add_argument("--min-coverage", type=float, default=0.95)
    parser.add_argument("--max-excluded-samples", type=int, default=0)
    args = parser.parse_args()

    if not 0 < args.min_coverage <= 1:
        raise SystemExit("[FAIL] --min-coverage must be in (0, 1]")
    args.outdir.mkdir(parents=True, exist_ok=True)

    need_rank = "t" if args.profile_rank == "sgb" else "s"
    sizes = load_genome_sizes(args.genome_sizes, args.match_on, need_rank)
    manifest, exclusions = read_manifest(args.manifest, args.inclusion_policy)

    samples: list[dict[str, str]] = []
    unmapped_tally: dict[str, list[float]] = {}
    ledger: dict[str, dict[str, str]] = {}

    for row in manifest:
        source = Path(row["source_profile"])
        if args.profile_root and not source.is_absolute():
            source = args.profile_root / source
        if not source.is_file():
            exclusions.append({"sample_id": row["sample_id"], "profile_id": row["profile_id"],
                               "reason": "profile_not_found", "detail": str(source)})
            continue

        records = parse_profile(source)
        clades = [c for c, _, _ in records]
        terminal = {c: not any(o.startswith(c + "|") for o in clades) for c in clades}

        if args.profile_rank == "sgb":
            # Terminal rows at or below species rank form a non-overlapping
            # partition, so species mass is never added to its own SGB mass.
            eligible = [(c, leaf, v) for c, leaf, v in records
                        if c.upper() != "UNCLASSIFIED" and terminal[c]
                        and (leaf.startswith("t__") or leaf.startswith("s__"))]
        else:
            eligible = [(c, leaf, v) for c, leaf, v in records
                        if c.upper() != "UNCLASSIFIED" and leaf.startswith("s__")]

        excluded_mass = sum(
            v for c, leaf, v in records
            if c.upper() == "UNCLASSIFIED"
            or (terminal[c] and not leaf.startswith("s__") and not leaf.startswith("t__")))

        eligible_total = sum(v for _, _, v in eligible)
        if eligible_total <= 0:
            exclusions.append({"sample_id": row["sample_id"], "profile_id": row["profile_id"],
                               "reason": "zero_eligible_abundance", "detail": str(source)})
            continue

        unit = args.abundance_unit
        if unit == "auto":
            unit = "percent" if eligible_total > 10 else "fraction"

        mapped_mass = weighted = unmapped_mass = 0.0
        n_mapped = n_unmapped = 0
        for clade, leaf, value in eligible:
            key = clade if args.match_on == "lineage" else leaf
            size = sizes.get(key)          # exact match only; never imputed
            if size is None:
                n_unmapped += 1
                unmapped_mass += value
                unmapped_tally.setdefault(key, []).append(value)
                continue
            n_mapped += 1
            mapped_mass += value
            weighted += value * size

        coverage = mapped_mass / eligible_total
        if mapped_mass <= 0:
            exclusions.append({"sample_id": row["sample_id"], "profile_id": row["profile_id"],
                               "reason": "no_mapped_features", "detail": str(source)})
            continue
        if coverage < args.min_coverage:
            exclusions.append({"sample_id": row["sample_id"], "profile_id": row["profile_id"],
                               "reason": "coverage_below_threshold",
                               "detail": f"{coverage:.6f} < {args.min_coverage}"})
            continue

        resolved = str(source.resolve())
        owner = (row["cohort"], row["sample_id"])
        if resolved in ledger:
            # Hash each physical file once, but never let two different
            # biological samples silently claim the same baseline file.
            previous = (ledger[resolved]["cohort"], ledger[resolved]["sample_id"])
            if previous != owner:
                raise SystemExit(
                    f"[FAIL] source profile {resolved} is claimed by both "
                    f"{previous[0]}/{previous[1]} and {owner[0]}/{owner[1]}")
        else:
            ledger[resolved] = {
                "cohort": row["cohort"], "sample_id": row["sample_id"],
                "profile_id": row["profile_id"], "source_profile": resolved,
                "bytes": str(source.stat().st_size), "sha256": sha256(source)}

        samples.append({
            "cohort": row["cohort"], "sample_id": row["sample_id"],
            "profile_id": row["profile_id"],
            "analysis_population": "",          # sample-wide, serves both populations
            "effective_genome_size_bp": render(weighted / mapped_mass),
            "mapping_coverage": render(coverage),
            "eligible_abundance": render(eligible_total),
            "mapped_abundance": render(mapped_mass),
            "unmapped_abundance": render(unmapped_mass),
            "excluded_abundance": render(excluded_mass),
            "n_mapped_features": str(n_mapped), "n_unmapped_features": str(n_unmapped),
            "profile_rank": args.profile_rank, "abundance_unit": unit,
            "source_profile": resolved,
        })

    samples.sort(key=lambda r: (r["cohort"], r["sample_id"]))
    exclusions.sort(key=lambda r: (r["sample_id"], r["profile_id"], r["reason"]))

    def dump(name, fields, rows):
        path = args.outdir / name
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        return path

    table = dump("effective_genome_size.tsv", SAMPLE_FIELDS, samples)
    excl = dump("excluded_profiles.tsv", ["sample_id", "profile_id", "reason", "detail"],
                exclusions)
    prov = dump("source_profile_checksums.tsv",
                ["cohort", "sample_id", "profile_id", "source_profile", "bytes", "sha256"],
                [ledger[k] for k in sorted(ledger)])

    unmapped = args.outdir / "unmapped_features.tsv"
    with unmapped.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["feature", "profiles_affected", "total_abundance", "max_abundance"])
        for key in sorted(unmapped_tally):
            values = unmapped_tally[key]
            writer.writerow([key, len(values), render(sum(values)), render(max(values))])

    coverages = sorted(float(s["mapping_coverage"]) for s in samples)
    geffs = sorted(float(s["effective_genome_size_bp"]) for s in samples)

    def quantile(values, q):
        if not values:
            return "NA"
        return render(values[min(len(values) - 1, int(q * (len(values) - 1) + 0.5))])

    ok = len(exclusions) <= args.max_excluded_samples
    audit = args.outdir / "effective_genome_size_audit.tsv"
    audit.write_text("\n".join([
        "metric\tvalue",
        f"manifest_profiles\t{len(manifest)}",
        f"profiles_with_geff\t{len(samples)}",
        f"profiles_excluded\t{len(exclusions)}",
        f"genome_size_entries\t{len(sizes)}",
        f"profile_rank\t{args.profile_rank}",
        f"match_on\t{args.match_on}",
        f"inclusion_policy\t{args.inclusion_policy}",
        f"min_coverage_threshold\t{args.min_coverage}",
        f"coverage_min\t{quantile(coverages, 0.0)}",
        f"coverage_median\t{quantile(coverages, 0.5)}",
        f"coverage_max\t{quantile(coverages, 1.0)}",
        f"geff_min_bp\t{quantile(geffs, 0.0)}",
        f"geff_median_bp\t{quantile(geffs, 0.5)}",
        f"geff_max_bp\t{quantile(geffs, 1.0)}",
        f"distinct_unmapped_features\t{len(unmapped_tally)}",
        f"status\t{'PASS' if ok else 'FAIL_EXCLUSIONS'}",
    ]) + "\n", encoding="utf-8")

    (args.outdir / "effective_genome_size.sha256").write_text(
        "".join(f"{sha256(p)}  {p.resolve()}\n"
                for p in [args.manifest, args.genome_sizes, table, excl, prov,
                          unmapped, audit]),
        encoding="utf-8")

    if not ok:
        raise SystemExit(
            f"[FAIL] {len(exclusions)} profile(s) excluded, limit is "
            f"{args.max_excluded_samples}; see {excl}")
    (args.outdir / "SUCCESS").write_text(
        f"profiles_with_geff\t{len(samples)}\nstatus\tPASS\n", encoding="utf-8")
    print(f"[PASS] G_eff for {len(samples)} baseline profile(s) at rank {args.profile_rank}")
    print(f"[INFO] {table}")


if __name__ == "__main__":
    main()
