#!/usr/bin/env python3
"""Fail-closed readiness gate for the definitive Yachida downstream run."""
from __future__ import annotations

import argparse
import csv
import hashlib
from collections import Counter
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty {label}: {path}")


def verify_seal(seal: Path) -> None:
    checksum = seal / "production_seal.sha256"
    require_file(checksum, "production-seal checksum")
    for line in checksum.read_text(encoding="utf-8").splitlines():
        expected, name = line.split(maxsplit=1)
        path = seal / name.strip()
        require_file(path, "sealed artifact")
        if digest(path) != expected:
            raise ValueError(f"production-seal checksum mismatch: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--state-dir", required=True, type=Path)
    parser.add_argument("--canonical", required=True, type=Path)
    parser.add_argument("--canonical-success", required=True, type=Path)
    parser.add_argument("--analysis-sif", required=True, type=Path)
    parser.add_argument("--assembly-sensitivity-success", required=True, type=Path)
    parser.add_argument("--expected-samples", type=int, default=201)
    parser.add_argument("--expected-independent", type=int, default=30)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    try:
        for path, label in ((args.manifest, "manifest"), (args.canonical, "canonical input"),
                            (args.canonical_success, "canonical validation marker"),
                            (args.analysis_sif, "analysis image"),
                            (args.assembly_sensitivity_success, "assembly-sensitivity seal")):
            require_file(path, label)
        manifest = rows(args.manifest)
        if len(manifest) != args.expected_samples:
            raise ValueError(f"expected {args.expected_samples} manifest samples; found {len(manifest)}")
        required = {"sample_id", "Target_Condition"}
        if not manifest or not required.issubset(manifest[0]):
            raise ValueError("manifest lacks sample_id or Target_Condition")
        identifiers = [row["sample_id"] for row in manifest]
        if len(set(identifiers)) != len(identifiers):
            raise ValueError("manifest contains duplicate sample IDs")
        if args.expected_samples == 201 and Counter(row["Target_Condition"] for row in manifest) != {
                "Control": 67, "Adenoma": 67, "CRC": 67}:
            raise ValueError("manifest is not the frozen balanced 67/67/67 design")
        state = args.state_dir.resolve()
        for sample in identifiers:
            require_file(state / "samples" / f"{sample}.verified", f"verified marker for {sample}")
            require_file(state / "samples" / f"{sample}.retained_outputs.tsv", f"receipt for {sample}")
        seal = state / "production_seal"
        require_file(seal / "SUCCESS", "production seal")
        verify_seal(seal)

        canonical = rows(args.canonical)
        needed = {"cohort", "sample_id", "analysis_population", "profiler", "include"}
        if not canonical or not needed.issubset(canonical[0]):
            raise ValueError("canonical input lacks readiness columns")
        included = [row for row in canonical if row["include"] == "1"]
        if {row["cohort"] for row in included} != {"yachida"}:
            raise ValueError("canonical input is not Yachida-only")
        if {row["profiler"] for row in included} != {"kraken2_bracken", "metaphlan4"}:
            raise ValueError("canonical input does not contain both frozen profilers")
        populations = {row["analysis_population"] for row in included}
        if populations != {"community", "independent"}:
            raise ValueError(f"canonical populations are {sorted(populations)}, expected community and independent")
        community = {row["sample_id"] for row in included if row["analysis_population"] == "community"}
        independent = {row["sample_id"] for row in included if row["analysis_population"] == "independent"}
        if community != set(identifiers):
            raise ValueError("community canonical population does not cover the frozen manifest exactly")
        if len(independent) != args.expected_independent or not independent.issubset(community):
            raise ValueError("independent canonical population is not the frozen nested subset")

        args.outdir.mkdir(parents=True, exist_ok=True)
        report = args.outdir / "readiness.tsv"
        report.write_text(
            "metric\tvalue\n"
            f"manifest_samples\t{len(manifest)}\n"
            f"verified_samples\t{len(identifiers)}\n"
            f"canonical_rows\t{len(canonical)}\n"
            f"community_samples\t{len(community)}\n"
            f"independent_samples\t{len(independent)}\n"
            "profilers\tkraken2_bracken;metaphlan4\n"
            "status\tPASS\n", encoding="utf-8")
        (args.outdir / "readiness.sha256").write_text(
            f"{digest(args.manifest)}  {args.manifest.resolve()}\n"
            f"{digest(args.canonical)}  {args.canonical.resolve()}\n"
            f"{digest(args.analysis_sif)}  {args.analysis_sif.resolve()}\n"
            f"{digest(report)}  {report.resolve()}\n", encoding="utf-8")
        (args.outdir / "SUCCESS").write_text("analysis\tyachida_definitive_readiness\nstatus\tPASS\n", encoding="utf-8")
        print(f"[PASS] Definitive Yachida readiness: {len(manifest)} samples; both populations")
    except (ValueError, OSError, KeyError) as error:
        raise SystemExit(f"[ERROR] {error}") from error


if __name__ == "__main__":
    main()
