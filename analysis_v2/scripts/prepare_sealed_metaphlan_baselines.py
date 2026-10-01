#!/usr/bin/env python3
"""Select raw baseline profiles exclusively from verified v2 seals and receipts."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
from audit_bracken_denominators import require, digest, table, verify_file
from compute_effective_genome_size import MANIFEST_REQUIRED


def prepare(cohort, state, results, out, expected):
    state, results, out = [p.resolve() for p in (state, results, out)]
    require(not out.exists(), "Output already exists")
    for source in (state, results):
        require(out != source and source not in out.parents and out not in source.parents,
                "Output overlaps source")
    seal = state / "production_seal_v2"
    require(not (seal / "AUDIT_IN_PROGRESS").exists(), "Audit in progress")
    checked = set()
    inputs = []
    for line in (seal / "production_seal.sha256").read_text().splitlines():
        checksum, name = line.split(maxsplit=1)
        name = name.lstrip("*")
        require(Path(name).name == name and name not in checked, "Unsafe/duplicate member")
        require(digest(seal / name) == checksum, "Seal checksum mismatch: " + name)
        checked.add(name)
        inputs.append(dict(path=str(seal / name), sha256=checksum))
    require({"SUCCESS", "production_manifest.tsv", "sample_flow.tsv"} <= checked,
            "Required seal member not checksummed")
    success = dict(line.split("\t", 1) for line in (seal / "SUCCESS").read_text().splitlines())
    require(success.get("status") == "PASS" and success.get("cohort") == cohort, "Wrong/failed seal")
    samples = table(seal / "production_manifest.tsv")
    ids = [s["sample_id"] for s in samples]
    require(len(ids) == len(set(ids)) == expected == int(success["samples"]), "Sample count mismatch")
    flows = table(seal / "sample_flow.tsv")
    require(len(flows) == len(ids) and {s["sample_id"] for s in flows} == set(ids), "Flow mismatch")
    require(all(s["status"] == "PASS" for s in flows), "Failed sample flow")
    rows = []
    for sample in samples:
        sid, study = sample["sample_id"], sample["study"]
        require(all(x not in ("", ".", "..") and Path(x).name == x for x in (sid, study)), "Unsafe identity")
        receipt_path = state / "samples" / (sid + ".retained_outputs.tsv")
        inputs.append(dict(path=str(receipt_path), sha256=digest(receipt_path)))
        receipt = {}
        for item in table(receipt_path):
            path = Path(item["path"])
            require(path.is_absolute(), "Nonabsolute receipt path")
            path = path.resolve()
            require(path not in receipt, "Duplicate receipt path")
            receipt[path] = item
        baseline = results / study / sid / "profiles" / "baseline"
        candidates = list(baseline.rglob(sid + ".metaphlan.tsv"))
        require(len(candidates) == 1, "Missing/ambiguous raw baseline: " + sid)
        path = candidates[0].resolve()
        require(results in path.parents and baseline.resolve() in path.parents, "Baseline escaped results")
        verify_file(path, receipt)
        inputs.append(dict(path=str(path), sha256=digest(path)))
        rows.append(dict(cohort=cohort, sample_id=sid, analysis_population="",
                         profiler="metaphlan4", profile_id=sid, baseline_profile_id=sid,
                         spike_fraction_total="0", source_profile=str(path), include="1",
                         exclusion_reason=""))
    require(all(digest(Path(r["path"])) == r["sha256"] for r in inputs), "Input changed during validation")
    out.mkdir(parents=True)
    for name, fields, data in (("baseline_manifest.tsv", MANIFEST_REQUIRED, rows),
                               ("input_hashes.tsv", ["path", "sha256"], inputs)):
        with (out / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(data)
    (out / "SHA256SUMS").write_text("".join(digest(out / n) + "  " + n + "\n"
                                            for n in ("baseline_manifest.tsv", "input_hashes.tsv")))
    print("[PASS] Receipt-verified baselines:", cohort, len(rows))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", choices=["yachida", "feng", "zeller"], required=True)
    for flag in ("state-dir", "results-root", "outdir"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--expected-samples", type=int, required=True)
    args = parser.parse_args()
    try:
        prepare(args.cohort, args.state_dir, args.results_root, args.outdir, args.expected_samples)
    except (ValueError, OSError) as error:
        raise SystemExit("[ERROR] " + str(error))
