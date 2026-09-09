#!/usr/bin/env python3
"""Adapt historical native Feng/Zeller profiles for v2 development only."""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import subprocess
import sys
from pathlib import Path

from build_assembly_sensitivity_input import HEADER, bracken_abundance, metaphlan_abundance, read_rows
from build_crc_cohort_canonical_input import make_row

PROFILERS = ("kraken2_bracken", "metaphlan4")
SUFFIX = {"kraken2_bracken": ".bracken.S.tsv", "metaphlan4": ".metaphlan.tsv"}
PROFILE_DIRS = {
    "kraken2_bracken": ("kraken2_bracken", "kraken_bracken"),
    "metaphlan4": ("metaphlan4",),
}
FROZEN_DOSES = {"0p0001", "0p0005", "0p001", "0p005", "0p01", "0p05", "0p1"}
NAME = re.compile(r"^(ERR\d+)(?:_([A-Za-z0-9]+)_f(0p[0-9]+))?$")


def args_parser() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results-root", required=True, type=Path)
    p.add_argument("--feng-manifest", required=True, type=Path)
    p.add_argument("--zeller-manifest", required=True, type=Path)
    p.add_argument("--spike-panel", required=True, type=Path)
    p.add_argument("--aliases", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument(
        "--profiler-coverage", choices=("paired", "available"), default="paired",
        help="Use only common Kraken/MetaPhlAn profile pairs (default), or all available profiles.",
    )
    return p.parse_args()


def accession_map(paths: list[tuple[str, Path]]) -> dict[str, dict[str, str]]:
    result = {}
    for cohort, path in paths:
        for row in read_rows(path):
            sample_accessions = sorted(
                accession.strip() for accession in row["run_accessions"].split(";")
                if accession.strip()
            )
            if not sample_accessions:
                raise ValueError(f"sample has no run accession: {row['sample_id']}")
            for accession in sample_accessions:
                accession = accession.strip()
                if accession in result and result[accession]["sample_id"] != row["sample_id"]:
                    raise ValueError(f"run accession maps to multiple samples: {accession}")
                result[accession] = {
                    "cohort": cohort, "sample_id": row["sample_id"],
                    "condition": row["condition"], "study": row["study"],
                }
    return result


def one_profile(root: Path, accession: str, profiler: str) -> Path | None:
    found = sorted(
        path
        for directory in PROFILE_DIRS[profiler]
        for path in (root / directory).glob(accession + "*" + SUFFIX[profiler])
    )
    return found[0] if len(found) == 1 else None


def abundance(path: Path, profiler: str, alias: str) -> float:
    return (bracken_abundance if profiler == "kraken2_bracken" else metaphlan_abundance)(path, alias)


def select_representative_runs(results_root: Path, accessions: dict[str, dict[str, str]],
                               targets: dict[str, dict[str, str]]) -> dict[tuple[str, str], str]:
    """Choose one historical run per biological sample by paired-profile coverage."""
    scores: dict[str, int] = {accession: 0 for accession in accessions}
    for directory in results_root.iterdir():
        if not directory.is_dir():
            continue
        match = NAME.fullmatch(directory.name)
        if not match or match.group(2) is None:
            continue
        run, label, dose_tag = match.groups()
        if (run not in accessions or dose_tag not in FROZEN_DOSES or
                (label != "CRCpanel" and label not in targets)):
            continue
        baseline = results_root / run
        if all(one_profile(baseline, run, profiler) is not None and
               one_profile(directory, directory.name, profiler) is not None
               for profiler in PROFILERS):
            scores[run] += 1

    grouped: dict[tuple[str, str], list[str]] = {}
    for accession, meta in accessions.items():
        grouped.setdefault((meta["cohort"], meta["sample_id"]), []).append(accession)
    return {
        sample: sorted(runs, key=lambda run: (-scores[run], run))[0]
        for sample, runs in grouped.items()
    }


def main() -> None:
    a = args_parser()
    try:
        panel = read_rows(a.spike_panel)
        targets = {row["label"]: row for row in panel}
        aliases = read_rows(a.aliases, delimiter=",")
        alias = {(row["canonical"], row["tool"]): row["alias"] for row in aliases}
        accessions = accession_map([("feng", a.feng_manifest), ("zeller", a.zeller_manifest)])
        representatives = select_representative_runs(a.results_root, accessions, targets)
        positives: list[tuple[dict[str, str], str, str, float, str, Path, Path]] = []
        exclusions: list[list[str]] = []

        for directory in sorted(a.results_root.iterdir()):
            if not directory.is_dir():
                continue
            match = NAME.fullmatch(directory.name)
            if not match or match.group(2) is None:
                continue
            run, label, dose_tag = match.groups()
            if label != "CRCpanel" and label not in targets:
                exclusions.append([directory.name, run, label, dose_tag, "unsupported_target"]); continue
            if dose_tag not in FROZEN_DOSES:
                exclusions.append([directory.name, run, label, dose_tag, "outside_frozen_dose_grid"]); continue
            if run not in accessions:
                exclusions.append([directory.name, run, label, dose_tag, "run_not_in_frozen_manifests"]); continue
            sample_key = (accessions[run]["cohort"], accessions[run]["sample_id"])
            if run != representatives[sample_key]:
                exclusions.append([directory.name, run, label, dose_tag,
                                   "non_primary_run_for_multirun_sample"]); continue
            dose = float(dose_tag.replace("p", "."))
            baseline_dir = a.results_root / run
            paths = {
                profiler: (
                    one_profile(baseline_dir, run, profiler),
                    one_profile(directory, directory.name, profiler),
                )
                for profiler in PROFILERS
            }
            if a.profiler_coverage == "paired" and any(
                base is None or spike is None for base, spike in paths.values()
            ):
                exclusions.append([
                    directory.name, run, label, dose_tag,
                    "incomplete_common_profiler_pair",
                ])
                continue
            for profiler in PROFILERS:
                base, spike = paths[profiler]
                if base is None or spike is None:
                    exclusions.append([directory.name, run, label, dose_tag, f"missing_{profiler}_pair"]); continue
                positives.append((accessions[run], run, label, dose, profiler, base, spike))

        rows = []
        seen_baselines = set()
        for meta, run, label, total_dose, profiler, base, spike in positives:
            population = "community" if label == "CRCpanel" else "independent"
            # Community profiles are expanded once per panel member.
            expanded = panel if population == "community" else [targets[label]]
            for member in expanded:
                member_alias = alias[(member["taxon_name"], profiler)]
                target_dose = total_dose / len(panel) if population == "community" else total_dose
                key = (meta["cohort"], meta["sample_id"], population, member["label"], profiler)
                if key not in seen_baselines:
                    rows.append(make_row(meta, population, member, profiler, run, run, 0, 0, 0,
                                         abundance(base, profiler, member_alias), base, "BASELINE", meta["cohort"]))
                    seen_baselines.add(key)
                # Historical exact pair allocations are unavailable. A positive
                # placeholder satisfies the structural fixture only; dose models
                # use the explicitly nominal fractions, never this count.
                rows.append(make_row(meta, population, member, profiler, spike.parent.parent.name, run,
                                     total_dose, target_dose, 1,
                                     abundance(spike, profiler, member_alias), spike,
                                     "LEGACY_NOMINAL_FRACTION_ONLY", meta["cohort"]))

        if not rows:
            raise ValueError("no complete legacy baseline-spike pairs found")
        a.outdir.mkdir(parents=True, exist_ok=True)
        canonical = a.outdir / "canonical_input.tsv"
        with canonical.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=HEADER, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
        with (a.outdir / "exclusion_ledger.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["directory", "run_accession", "target_label", "dose_tag", "reason"])
            writer.writerows(exclusions)
        (a.outdir / "DEVELOPMENT_ONLY.txt").write_text(
            "status=DEVELOPMENT_ONLY\nexact_implanted_pair_counts_available=0\n"
            f"fractions=historical_nominal\nprofiler_coverage={a.profiler_coverage}\n"
            "prohibited_for_manuscript_results=1\n", encoding="utf-8")
        validator = Path(__file__).with_name("validate_canonical_input.py")
        subprocess.run([sys.executable, str(validator), "--input", str(canonical),
                        "--outdir", str(a.outdir / "validation")], check=True)
        files = [canonical, a.outdir / "exclusion_ledger.tsv", a.outdir / "DEVELOPMENT_ONLY.txt"]
        with (a.outdir / "development_input.sha256").open("w", encoding="utf-8") as handle:
            for path in files:
                handle.write(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n")
        print(f"[PASS] Legacy development input built: {len(rows)} canonical rows")
        print(f"[INFO] Excluded/incomplete directories: {len(exclusions)}")
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        raise SystemExit(f"[ERROR] {error}") from error


if __name__ == "__main__":
    main()
