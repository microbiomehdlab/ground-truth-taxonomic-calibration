#!/usr/bin/env python3
"""Build canonical v2 rows for the additive Pana/Pint assembly sensitivity."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


SCHEMA = "paired-dose-response-v2.1"
HEADER = [
    "schema_version", "cohort", "study", "sample_id", "condition",
    "analysis_population", "target_label", "target_taxon", "assembly_arm",
    "profiler", "profile_id", "baseline_profile_id", "spike_fraction_total",
    "spike_fraction_target", "implanted_read_pairs_target", "native_abundance",
    "native_unit", "abundance_fraction", "detected_native_nonzero",
    "source_profile", "source_design", "include", "exclusion_reason",
]


def read_rows(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def require_columns(rows: list[dict[str, str]], columns: set[str], label: str) -> None:
    present = set(rows[0]) if rows else set()
    missing = columns - present
    if missing:
        raise ValueError("{} missing columns: {}".format(label, ", ".join(sorted(missing))))


def unique_file(root: Path, filename: str) -> Path:
    matches = [path for path in root.rglob(filename) if path.is_file() and path.stat().st_size]
    if len(matches) != 1:
        raise ValueError("expected one {} under {}; found {}".format(filename, root, len(matches)))
    return matches[0]


def profile_id_from_design(row: dict[str, str], sample: str, arm: str) -> str:
    output = Path(row["output_r1"]).name
    for suffix in ("_1.fq.gz", "_1.fastq.gz", "_1.fq", "_1.fastq"):
        if output.endswith(suffix):
            return output[:-len(suffix)]
    # Frozen naming fallback, checked against the actual profile directory.
    fraction = float(row["fraction"])
    tag = ("{:.6f}".format(fraction)).rstrip("0").rstrip(".").replace(".", "p")
    return "{}_{}_f{}".format(sample, arm, tag)


def bracken_abundance(path: Path, alias: str) -> float:
    rows = read_rows(path)
    require_columns(rows, {"name", "taxonomy_lvl", "fraction_total_reads"}, str(path))
    matches = [row for row in rows if row["taxonomy_lvl"].strip().upper() == "S"
               and row["name"].strip() == alias]
    if len(matches) > 1:
        raise ValueError("duplicate Bracken target {} in {}".format(alias, path))
    return float(matches[0]["fraction_total_reads"]) if matches else 0.0


def metaphlan_abundance(path: Path, alias: str) -> float:
    wanted = alias.replace(" ", "_")
    matches = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError("malformed MetaPhlAn row in {}".format(path))
            terminal = fields[0].split("|")[-1]
            if terminal.startswith("s__") and terminal[3:] == wanted:
                matches.append(float(fields[2]))
    if len(matches) > 1:
        raise ValueError("duplicate MetaPhlAn target {} in {}".format(alias, path))
    return matches[0] if matches else 0.0


def locate_profile(root: Path, profile_id: str, profiler: str) -> Path:
    suffix = ".bracken.S.tsv" if profiler == "kraken2_bracken" else ".metaphlan.tsv"
    return unique_file(root, profile_id + suffix)


def canonical_row(meta: dict[str, str], arm: dict[str, str], profiler: str,
                  profile_id: str, baseline_id: str, total_fraction: float,
                  target_fraction: float, pairs: int, abundance: float,
                  profile: Path, design: str, assembly_arm: str) -> dict[str, str]:
    unit = "fraction_total_reads" if profiler == "kraken2_bracken" else "relative_abundance_pct"
    fraction = abundance if profiler == "kraken2_bracken" else abundance / 100.0
    return {
        "schema_version": SCHEMA, "cohort": "yachida", "study": meta["Study"],
        "sample_id": meta["sample_id"], "condition": meta["Target_Condition"],
        "analysis_population": "independent", "target_label": arm["target_label"],
        "target_taxon": arm["target_taxon"], "assembly_arm": assembly_arm,
        "profiler": profiler, "profile_id": profile_id,
        "baseline_profile_id": baseline_id,
        "spike_fraction_total": format(total_fraction, ".17g"),
        "spike_fraction_target": format(target_fraction, ".17g"),
        "implanted_read_pairs_target": str(pairs),
        "native_abundance": format(abundance, ".17g"), "native_unit": unit,
        "abundance_fraction": format(fraction, ".17g"),
        "detected_native_nonzero": "1" if fraction > 0 else "0",
        "source_profile": str(profile.resolve()), "source_design": design,
        "include": "1", "exclusion_reason": "",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--arms", type=Path, required=True)
    parser.add_argument("--spike-panel", type=Path, required=True)
    parser.add_argument("--aliases", type=Path, required=True)
    parser.add_argument("--sensitivity-root", type=Path, required=True)
    parser.add_argument(
        "--original-root", type=Path,
        help="Optional sealed original-assembly run root; adds matched original arms.",
    )
    parser.add_argument("--baseline-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        manifest = read_rows(args.manifest)
        arms = read_rows(args.arms)
        panel = read_rows(args.spike_panel)
        aliases = read_rows(args.aliases, delimiter=",")
        require_columns(manifest, {"sample_id", "Study", "Target_Condition"}, "manifest")
        require_columns(arms, {"arm_label", "target_label", "assembly_accession"}, "arms")
        require_columns(panel, {"label", "taxon_name"}, "spike panel")
        require_columns(aliases, {"canonical", "alias", "tool"}, "aliases")
        taxa = {row["label"]: row["taxon_name"] for row in panel}
        alias_map = {(row["canonical"], row["tool"]): row["alias"] for row in aliases}
        for arm in arms:
            if arm["target_label"] not in taxa:
                raise ValueError("arm target absent from spike panel: " + arm["target_label"])
            arm["target_taxon"] = taxa[arm["target_label"]]

        output_rows = []
        for meta in manifest:
            sample = meta["sample_id"]
            study = meta["Study"]
            sample_sensitivity = args.sensitivity_root / "results" / study / sample
            sample_baseline = args.baseline_root / "results" / study / sample
            for arm in arms:
                sources = [("clean", sample_sensitivity, arm["arm_label"])]
                if args.original_root:
                    original_sample = args.original_root / "results" / study / sample
                    sources.insert(0, ("original", original_sample, arm["target_label"]))
                for assembly_arm, source_sample, source_label in sources:
                    design_path = source_sample / "spike_design" / "independent" / (source_label + ".tsv")
                    design_rows = read_rows(design_path)
                    require_columns(design_rows, {"sample_id", "label", "fraction", "N_inserted",
                                                  "f_hat", "output_r1"}, str(design_path))
                    design_rows = [row for row in design_rows if row["sample_id"] == sample
                                   and row["label"] == source_label]
                    if len(design_rows) != 6:
                        raise ValueError("expected six design rows in {}; found {}".format(
                            design_path, len(design_rows)))
                    for profiler in ("kraken2_bracken", "metaphlan4"):
                        alias = alias_map.get((arm["target_taxon"], profiler))
                        if not alias:
                            raise ValueError("missing alias for {} {}".format(arm["target_taxon"], profiler))
                        baseline_profile = locate_profile(sample_baseline / "profiles" / "baseline", sample, profiler)
                        abundance_reader = bracken_abundance if profiler == "kraken2_bracken" else metaphlan_abundance
                        baseline_abundance = abundance_reader(baseline_profile, alias)
                        output_rows.append(canonical_row(
                            meta, arm, profiler, sample, sample, 0.0, 0.0, 0,
                            baseline_abundance, baseline_profile, "BASELINE", assembly_arm))
                        for design in sorted(design_rows, key=lambda row: float(row["fraction"])):
                            profile_id = profile_id_from_design(design, sample, source_label)
                            profile_root = source_sample / "profiles" / "independent" / source_label
                            profile = locate_profile(profile_root, profile_id, profiler)
                            abundance = abundance_reader(profile, alias)
                            achieved = float(design["f_hat"])
                            output_rows.append(canonical_row(
                                meta, arm, profiler, profile_id, sample, achieved, achieved,
                                int(design["N_inserted"]), abundance, profile,
                                str(design_path.resolve()), assembly_arm))

        args.outdir.mkdir(parents=True, exist_ok=True)
        output = args.outdir / "canonical_input.tsv"
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=HEADER, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(output_rows)
        validator = Path(__file__).with_name("validate_canonical_input.py")
        subprocess.run([sys.executable, str(validator), "--input", str(output),
                        "--outdir", str(args.outdir / "validation")], check=True)
        print("[PASS] Assembly-sensitivity canonical input built: {} rows".format(len(output_rows)))
        print("[INFO] {}".format(output))
    except (ValueError, FileNotFoundError, KeyError, subprocess.CalledProcessError) as error:
        raise SystemExit("[ERROR] {}".format(error)) from error


if __name__ == "__main__":
    main()
