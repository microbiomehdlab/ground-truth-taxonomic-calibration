#!/usr/bin/env python3
"""Audit and seal a completed Feng or Zeller upstream cohort."""
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
from collections import Counter
from pathlib import Path


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def valid(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def verify_sidecar(path: Path) -> None:
    sidecar = path.with_suffix(path.suffix + ".sha256")
    if not valid(sidecar):
        raise ValueError(f"missing manifest checksum: {sidecar}")
    fields = sidecar.read_text(encoding="utf-8").split()
    if len(fields) != 2 or fields[0] != digest(path) or fields[1] != path.name:
        raise ValueError(f"manifest checksum mismatch: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", required=True, choices=("feng", "zeller"))
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--independent-manifest", required=True, type=Path)
    parser.add_argument("--state-dir", required=True, type=Path)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--expected-samples", required=True, type=int)
    parser.add_argument("--expected-independent", type=int, default=30)
    args = parser.parse_args()
    try:
        verify_sidecar(args.manifest); verify_sidecar(args.independent_manifest)
        manifest = read(args.manifest); independent_rows = read(args.independent_manifest)
        required = {"sample_id", "condition", "study", "age", "sex", "bmi", "independent_subset"}
        if not manifest or not required.issubset(manifest[0]):
            raise ValueError("production manifest lacks required metadata columns")
        if len(manifest) != args.expected_samples:
            raise ValueError(f"expected {args.expected_samples} samples; found {len(manifest)}")
        identifiers = [row["sample_id"] for row in manifest]
        independent = {row["sample_id"] for row in independent_rows}
        if len(set(identifiers)) != len(identifiers) or len(independent) != args.expected_independent or not independent.issubset(identifiers):
            raise ValueError("production/independent sample identity contract failed")
        if args.expected_independent == 30 and Counter(row["condition"] for row in independent_rows) != {"Control": 10, "Adenoma": 10, "CRC": 10}:
            raise ValueError("independent subset is not balanced 10/10/10")
        flow = []
        for row in manifest:
            sample = row["sample_id"]; sample_root = args.results_root / row["study"] / sample
            expected_profiles = 68 if sample in independent else 8
            observed_profiles = sum(1 for path in (sample_root / "profiles").rglob("SUCCESS") if valid(path))
            checks = {
                "verified_marker": valid(args.state_dir / "samples" / f"{sample}.verified"),
                "retained_output_receipt": valid(args.state_dir / "samples" / f"{sample}.retained_outputs.tsv"),
                "input_provenance": valid(args.state_dir / "samples" / f"{sample}.input_provenance.tsv"),
                "sample_success": valid(sample_root / "SUCCESS"),
                "profile_count": observed_profiles == expected_profiles,
            }
            flow.append({"sample_id": sample, "study": row["study"], "condition": row["condition"],
                         "independent_subset": int(sample in independent), "expected_profiles": expected_profiles,
                         "observed_profiles": observed_profiles, **{name: int(ok) for name, ok in checks.items()},
                         "status": "PASS" if all(checks.values()) else "FAIL",
                         "failure_reasons": ";".join(name for name, ok in checks.items() if not ok)})
        args.outdir.mkdir(parents=True, exist_ok=True)
        flow_path = args.outdir / "sample_flow.tsv"
        with flow_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(flow[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(flow)
        covariates = args.outdir / "covariate_audit.tsv"
        with covariates.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["field", "missing", "present", "distinct_nonmissing"])
            for field in ("condition", "age", "sex", "bmi"):
                values = [row[field].strip() for row in manifest]
                writer.writerow([field, sum(not value for value in values), sum(bool(value) for value in values),
                                 len({value for value in values if value})])
        for source, name in ((args.manifest, "production_manifest.tsv"),
                             (args.independent_manifest, "production_manifest.independent.tsv")):
            shutil.copyfile(source, args.outdir / name)
        if any(row["status"] != "PASS" for row in flow):
            raise ValueError(f"{sum(row['status'] != 'PASS' for row in flow)} samples are incomplete; see {flow_path}")
        success = args.outdir / "SUCCESS"
        success.write_text(f"cohort\t{args.cohort}\nsamples\t{len(manifest)}\nindependent_samples\t{len(independent)}\nstatus\tPASS\n", encoding="utf-8")
        artifacts = [flow_path, covariates, args.outdir / "production_manifest.tsv",
                     args.outdir / "production_manifest.independent.tsv", success]
        (args.outdir / "production_seal.sha256").write_text(
            "".join(f"{digest(path)}  {path.name}\n" for path in artifacts), encoding="utf-8")
        print(f"[PASS] Sealed {args.cohort} upstream: {len(manifest)} samples")
    except (ValueError, OSError, KeyError) as error:
        raise SystemExit(f"[ERROR] {error}") from error


if __name__ == "__main__":
    main()
