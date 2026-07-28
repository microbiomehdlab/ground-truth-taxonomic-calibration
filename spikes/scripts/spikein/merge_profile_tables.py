#!/usr/bin/env python3
"""
Merge fraction-level profile CSV files across spike labels, grouped by tool.

Expected organized layout:
  <root>/<spike_label>/<tool>/spike_f0p0001.open_world.csv
  <root>/<spike_label>/<tool>/spike_f0p0001.closed_world.csv

Outputs:
  <outdir>/<tool>.open_world.merged.csv
  <outdir>/<tool>.closed_world.merged.csv

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


KIND_PATTERNS = {
    "open_world": "spike_*.open_world.csv",
    "closed_world": "spike_*.closed_world.csv",
}


def parse_args():
    ap = argparse.ArgumentParser(
        description=(
            "Merge fraction-level open_world/closed_world CSV files across "
            "spike labels, grouped by tool."
        )
    )
    ap.add_argument(
        "--root",
        required=True,
        help=(
            "Root directory containing organized tables as "
            "<spike_label>/<tool>/spike_f*.open_world.csv and/or "
            "<spike_label>/<tool>/spike_f*.closed_world.csv"
        ),
    )
    ap.add_argument(
        "--outdir",
        required=True,
        help="Directory where merged per-tool CSVs will be written.",
    )
    ap.add_argument(
        "--kind",
        choices=["open_world", "closed_world", "both"],
        default="both",
        help=(
            "Which table kind to merge. Default: both. "
            "Use open_world for primary quantitative spike-in recovery."
        ),
    )
    ap.add_argument(
        "--open-pattern",
        default=KIND_PATTERNS["open_world"],
        help="Pattern for organized open-world files (default: spike_*.open_world.csv)",
    )
    ap.add_argument(
        "--closed-pattern",
        default=KIND_PATTERNS["closed_world"],
        help="Pattern for organized closed-world files (default: spike_*.closed_world.csv)",
    )
    ap.add_argument(
        "--sample-column",
        default="sample",
        help="Name of the sample ID column (default: sample).",
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


def infer_spike_label_and_tool(root: Path, fpath: Path):
    rel = fpath.relative_to(root)
    parts = rel.parts
    if len(parts) < 3:
        raise ValueError(f"Expected <spike_label>/<tool>/<file>, got: {rel}")
    return parts[0], parts[1]


def read_profile_csv(fpath: Path, sample_column: str):
    with fpath.open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            raise ValueError(f"No header found in {fpath}")
        if sample_column not in reader.fieldnames:
            raise ValueError(f"'{sample_column}' column not found in {fpath}")

        rows = []
        for row in reader:
            sample = (row.get(sample_column) or "").strip()
            if not sample:
                continue

            taxa = {}
            for k, v in row.items():
                if k == sample_column:
                    continue
                taxa[k] = v if v not in (None, "") else "0"
            rows.append((sample, taxa))
        return rows


def merge_kind(root: Path, outdir: Path, kind: str, pattern: str, sample_column: str):
    files = sorted(root.rglob(pattern))
    if not files:
        print(f"[warn] No {kind} files found under {root} matching {pattern}", file=sys.stderr)
        return 0

    grouped = defaultdict(lambda: {"taxa": set(), "rows": []})

    for fpath in files:
        try:
            spike_label, tool = infer_spike_label_and_tool(root, fpath)
            rows = read_profile_csv(fpath, sample_column)
        except Exception as e:
            print(f"[warn] Skipping {fpath}: {e}", file=sys.stderr)
            continue

        for sample, taxa in rows:
            grouped[tool]["taxa"].update(taxa.keys())
            grouped[tool]["rows"].append(
                {
                    "sample": sample,
                    "spike_label": spike_label,
                    "taxa": taxa,
                    "source_file": str(fpath),
                }
            )

    if not grouped:
        print(f"[warn] No valid {kind} files could be parsed.", file=sys.stderr)
        return 0

    written = 0
    for tool, payload in sorted(grouped.items()):
        taxa_cols = sorted(payload["taxa"])
        rows = payload["rows"]

        seen = {}
        deduped = []
        for r in rows:
            sid = r["sample"]
            if sid in seen:
                raise ValueError(
                    f"Duplicate sample ID detected for tool='{tool}', kind='{kind}': {sid}\n"
                    f"First source: {seen[sid]}\n"
                    f"Second source: {r['source_file']}\n"
                    "Make sure sample names are globally unique."
                )
            seen[sid] = r["source_file"]
            deduped.append(r)

        deduped.sort(key=lambda x: x["sample"])

        out_csv = outdir / f"{tool}.{kind}.merged.csv"
        with out_csv.open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow([sample_column] + taxa_cols)
            for r in deduped:
                writer.writerow([r["sample"]] + [r["taxa"].get(t, "0") for t in taxa_cols])

        print(f"[ok] Wrote {out_csv} ({len(deduped)} rows, {len(taxa_cols)} features)")
        written += 1

    return written


def main():
    args = parse_args()
    root = Path(args.root).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if not root.exists():
        print(f"[error] Root does not exist: {root}", file=sys.stderr)
        sys.exit(1)

    total_written = 0
    for kind in kinds_to_run(args.kind):
        total_written += merge_kind(
            root=root,
            outdir=outdir,
            kind=kind,
            pattern=pattern_for_kind(args, kind),
            sample_column=args.sample_column,
        )

    if total_written == 0:
        print("[error] No merged profile tables were created.", file=sys.stderr)
        sys.exit(1)

    print("[done] All requested merged profile tables created.")


if __name__ == "__main__":
    main()
