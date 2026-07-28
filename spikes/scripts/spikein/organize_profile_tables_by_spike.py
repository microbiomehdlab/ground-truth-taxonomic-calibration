#!/usr/bin/env python3
"""
Organize per-sample profiler profile CSVs into per-spike/per-tool/per-fraction matrices.

Raw profiling layout expected:
  <inroot>/<sample_id>/<tool>/<sample_id>.open_world.split.csv
  <inroot>/<sample_id>/<tool>/<sample_id>.closed_world.csv

Organized output:
  <outroot>/<spike_label>/<tool>/spike_f0p0001.open_world.csv
  <outroot>/<spike_label>/<tool>/spike_f0p0001.closed_world.csv

Definitions:
  open_world:
    Profiler-reported scale. Species-level features may be selected, but
    unclassified / above-species mass is not redistributed across retained taxa.
    Recommended primary scale for quantitative spike-in recovery when spike
    fractions are defined relative to all final read-pairs.

  closed_world:
    Retained feature table after renormalization to the retained feature space.
    Useful as a biomarker-style or sensitivity table, but it changes the
    denominator and should not be compared to an all-read spike fraction unless
    the expected value is denominator-adjusted.
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


DEFAULT_PATTERNS = {
    # Matches both <sample>.open_world.csv and <sample>.open_world.split.csv.
    "open_world": "*.open_world*.csv",
    "closed_world": "*.closed_world.csv",
}


def parse_args():
    ap = argparse.ArgumentParser(
        description=(
            "Group per-sample open_world/closed_world CSVs into "
            "per-spike/per-tool/per-fraction tables."
        )
    )
    ap.add_argument(
        "--inroot",
        required=True,
        help=(
            "Root of raw profiling results, e.g. "
            "results/<sample_id>/<tool>/<sample_id>.open_world.split.csv"
        ),
    )
    ap.add_argument(
        "--outroot",
        required=True,
        help=(
            "Output root, e.g. "
            "profile_results/<label>/<tool>/spike_f0p0001.open_world.csv"
        ),
    )
    ap.add_argument(
        "--kind",
        choices=["open_world", "closed_world", "both"],
        default="both",
        help=(
            "Which table kind to organize. Default: both. "
            "Use open_world for primary quantitative spike-in recovery."
        ),
    )
    ap.add_argument(
        "--open-pattern",
        default=DEFAULT_PATTERNS["open_world"],
        help="Filename pattern for per-sample open-world files (default: *.open_world*.csv)",
    )
    ap.add_argument(
        "--closed-pattern",
        default=DEFAULT_PATTERNS["closed_world"],
        help="Filename pattern for per-sample closed-world files (default: *.closed_world.csv)",
    )
    ap.add_argument(
        "--sample-column",
        default="sample",
        help="Name of the sample ID column (default: sample).",
    )
    ap.add_argument(
        "--drop-columns",
        nargs="*",
        default=[],
        help=(
            "Optional columns to remove from output matrices, e.g. "
            "--drop-columns AboveSpecies Unclassified. Default: keep all columns. "
            "Dropping columns does not renormalize the remaining values."
        ),
    )
    return ap.parse_args()


def kinds_to_run(kind_arg: str):
    return ["open_world", "closed_world"] if kind_arg == "both" else [kind_arg]


def pattern_for_kind(args, kind: str):
    if kind == "open_world":
        return args.open_pattern
    if kind == "closed_world":
        return args.closed_pattern
    raise ValueError(f"Unknown kind: {kind}")


def parse_sample_id(sample_id: str):
    """
    Expected sample_id format:
      <original_sample>_<spike_label>_<fraction>

    Example:
      ERR478961_Bfrag_f0p0001
    """
    parts = sample_id.rsplit("_", 2)
    if len(parts) != 3:
        raise ValueError(
            f"Could not parse sample_id '{sample_id}'. "
            f"Expected format <original>_<label>_<fraction>."
        )
    original_id, spike_label, fraction = parts
    return original_id, spike_label, fraction


def read_profile_csv(path: Path, sample_column: str, drop_columns):
    drop_columns = set(drop_columns or [])

    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            raise ValueError(f"No header found in {path}")
        if sample_column not in reader.fieldnames:
            raise ValueError(f"'{sample_column}' column missing in {path}")

        rows = []
        for row in reader:
            sample = (row.get(sample_column) or "").strip()
            if not sample:
                continue
            taxa = {}
            for k, v in row.items():
                if k == sample_column or k in drop_columns:
                    continue
                taxa[k] = v if v not in (None, "") else "0"
            rows.append((sample, taxa))
        return rows


def organize_kind(inroot: Path, outroot: Path, kind: str, pattern: str, sample_column: str, drop_columns):
    files = sorted(inroot.rglob(pattern))
    if not files:
        print(f"[warn] No {kind} files found under {inroot} matching {pattern}", file=sys.stderr)
        return 0

    grouped = defaultdict(lambda: {"taxa": set(), "rows": []})

    for fpath in files:
        rel = fpath.relative_to(inroot)
        parts = rel.parts

        # Expected raw layout: <sample_id>/<tool>/<file>
        if len(parts) < 3:
            print(f"[warn] Skipping unexpected path: {rel}", file=sys.stderr)
            continue

        sample_dir = parts[0]
        tool = parts[1]

        try:
            _original_id, spike_label, fraction = parse_sample_id(sample_dir)
        except Exception as e:
            print(f"[warn] Skipping {rel}: {e}", file=sys.stderr)
            continue

        try:
            rows = read_profile_csv(fpath, sample_column, drop_columns)
        except Exception as e:
            print(f"[warn] Skipping {rel}: {e}", file=sys.stderr)
            continue

        key = (spike_label, tool, fraction)

        for sample, taxa in rows:
            grouped[key]["taxa"].update(taxa.keys())
            grouped[key]["rows"].append((sample, taxa, str(fpath)))

    if not grouped:
        print(f"[warn] No valid {kind} CSVs were parsed.", file=sys.stderr)
        return 0

    written = 0
    for (spike_label, tool, fraction), payload in sorted(grouped.items()):
        taxa_cols = sorted(payload["taxa"])
        rows = payload["rows"]

        seen = {}
        deduped = []
        for sample, taxa, source_file in rows:
            if sample in seen:
                raise ValueError(
                    f"Duplicate sample '{sample}' found for "
                    f"label={spike_label}, tool={tool}, fraction={fraction}, kind={kind}\n"
                    f"First source: {seen[sample]}\n"
                    f"Second source: {source_file}"
                )
            seen[sample] = source_file
            deduped.append((sample, taxa))

        deduped.sort(key=lambda x: x[0])

        outdir = outroot / spike_label / tool
        outdir.mkdir(parents=True, exist_ok=True)
        outfile = outdir / f"spike_{fraction}.{kind}.csv"

        with outfile.open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow([sample_column] + taxa_cols)
            for sample, taxa in deduped:
                writer.writerow([sample] + [taxa.get(t, "0") for t in taxa_cols])

        print(f"[ok] Wrote {outfile} ({len(deduped)} rows, {len(taxa_cols)} features)")
        written += 1

    return written


def main():
    args = parse_args()
    inroot = Path(args.inroot).resolve()
    outroot = Path(args.outroot).resolve()
    outroot.mkdir(parents=True, exist_ok=True)

    if not inroot.exists():
        print(f"[error] Input root does not exist: {inroot}", file=sys.stderr)
        sys.exit(1)

    total_written = 0
    for kind in kinds_to_run(args.kind):
        total_written += organize_kind(
            inroot=inroot,
            outroot=outroot,
            kind=kind,
            pattern=pattern_for_kind(args, kind),
            sample_column=args.sample_column,
            drop_columns=args.drop_columns,
        )

    if total_written == 0:
        print("[error] No organized profile tables were created.", file=sys.stderr)
        sys.exit(1)

    print("[done] Organized requested profile tables by spike label, tool, and fraction.")


if __name__ == "__main__":
    main()
