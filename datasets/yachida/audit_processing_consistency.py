#!/usr/bin/env python3
"""Report which frozen Yachida samples have final results and paired-read QC."""
from __future__ import annotations

import argparse
import csv
import hashlib
import pathlib


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=pathlib.Path)
    parser.add_argument("--results-root", required=True, type=pathlib.Path)
    parser.add_argument("--qc-root", required=True, type=pathlib.Path)
    parser.add_argument("--study", default="YachidaS_2019")
    parser.add_argument("--output", required=True, type=pathlib.Path)
    args = parser.parse_args()
    with args.manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit("[ERROR] Empty manifest")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output_rows = []
    for row in rows:
        sample = row["sample_id"]
        sample_results = args.results_root / args.study / sample
        sample_qc = args.qc_root / args.study / sample
        complete = (sample_results / "SUCCESS").is_file()
        paired = sample_qc / "paired_fastq_integrity.tsv"
        provenance = sample_qc / "metashotgunprep_provenance.tsv"
        output_rows.append({
            "sample_id": sample,
            "condition": row["Target_Condition"],
            "analysis_complete": "yes" if complete else "no",
            "paired_fastq_integrity": "yes" if paired.is_file() and paired.stat().st_size else "no",
            "preprocessing_provenance": "yes" if provenance.is_file() and provenance.stat().st_size else "no",
            "action": "none" if complete and paired.is_file() and paired.stat().st_size else
                      ("reprocess_with_final_pipeline" if complete else "process_with_final_pipeline"),
        })
    fields = list(output_rows[0])
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(output_rows)
    digest = hashlib.sha256(args.output.read_bytes()).hexdigest()
    args.output.with_suffix(args.output.suffix + ".sha256").write_text(
        f"{digest}  {args.output.name}\n", encoding="utf-8"
    )
    counts: dict[str, int] = {}
    for row in output_rows:
        counts[row["action"]] = counts.get(row["action"], 0) + 1
    print(f"[OK] Processing-consistency report: {args.output}")
    for action in sorted(counts):
        print(f"[INFO] {action}: {counts[action]}")


if __name__ == "__main__":
    main()
