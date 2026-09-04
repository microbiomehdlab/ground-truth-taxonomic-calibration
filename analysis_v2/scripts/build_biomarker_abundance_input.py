#!/usr/bin/env python3
"""Expand canonical profile paths into a species-abundance input for paired DA."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


CONTEXT = ["cohort", "study", "analysis_population", "sample_id", "condition",
           "target_label", "assembly_arm", "profiler", "profile_id",
           "baseline_profile_id", "spike_fraction_target", "dose_level", "source_profile"]
MANIFEST = CONTEXT + ["target_taxon", "target_feature", "include", "exclusion_reason"]
ABUNDANCE = ["profiler", "source_profile", "feature", "abundance_fraction"]


def rows(path: Path, delimiter: str = "\t"):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError("empty table: {}".format(path))
        return list(reader), set(reader.fieldnames)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_profile(path: Path, profiler: str):
    found = {}
    if profiler == "kraken2_bracken":
        data, fields = rows(path)
        needed = {"name", "taxonomy_lvl", "fraction_total_reads"}
        if not needed <= fields:
            raise ValueError("{} missing Bracken columns".format(path))
        for row in data:
            if row["taxonomy_lvl"].strip().upper() != "S":
                continue
            feature = row["name"].strip()
            value = float(row["fraction_total_reads"])
            if feature in found:
                raise ValueError("duplicate species {} in {}".format(feature, path))
            found[feature] = value
    elif profiler == "metaphlan4":
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    raise ValueError("malformed MetaPhlAn row in {}".format(path))
                terminal = fields[0].split("|")[-1]
                if not terminal.startswith("s__"):
                    continue
                feature = terminal[3:].replace("_", " ")
                value = float(fields[2]) / 100.0
                if feature in found:
                    raise ValueError("duplicate species {} in {}".format(feature, path))
                found[feature] = value
    else:
        raise ValueError("unsupported profiler: {}".format(profiler))
    if any(value < 0 or value > 1.00001 for value in found.values()):
        raise ValueError("abundance outside fraction scale in {}".format(path))
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical", type=Path, required=True)
    parser.add_argument("--aliases", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    try:
        canonical, canonical_fields = rows(args.canonical)
        aliases, alias_fields = rows(args.aliases, delimiter=",")
        required = set(CONTEXT + ["target_taxon", "include", "exclusion_reason"])
        required.remove("dose_level")
        if not required <= canonical_fields:
            raise ValueError("canonical table missing: {}".format(
                ", ".join(sorted(required - canonical_fields))))
        if not {"canonical", "alias", "tool"} <= alias_fields:
            raise ValueError("alias table lacks canonical, alias, or tool")
        alias_map = {(row["canonical"], row["tool"]): row["alias"] for row in aliases}
        manifest = []
        identities = set()
        profiles = {}
        dose_groups = {}
        rank_key = ["cohort", "study", "analysis_population", "sample_id", "condition",
                    "target_label", "assembly_arm", "profiler"]
        for row in canonical:
            dose_groups.setdefault(tuple(row[field] for field in rank_key), []).append(row)
        dose_levels = {}
        for group_key, group in dose_groups.items():
            baseline = [row for row in group if float(row["spike_fraction_target"]) == 0]
            positive = sorted((row for row in group if float(row["spike_fraction_target"]) > 0),
                              key=lambda row: float(row["spike_fraction_target"]))
            if len(baseline) != 1 or len(positive) < 1:
                raise ValueError("each canonical dose series needs one baseline and positive doses: {}".format(group_key))
            dose_levels[id(baseline[0])] = "baseline"
            for index, row in enumerate(positive, start=1):
                dose_levels[id(row)] = "dose_{:02d}".format(index)
        for row in canonical:
            identity = tuple(row[field] for field in CONTEXT if field not in {"source_profile", "dose_level"})
            if identity in identities:
                raise ValueError("duplicate canonical context: {}".format(identity))
            identities.add(identity)
            target = alias_map.get((row["target_taxon"], row["profiler"]))
            if not target:
                raise ValueError("missing target alias for {} {}".format(
                    row["target_taxon"], row["profiler"]))
            record = {field: row[field] for field in CONTEXT if field != "dose_level"}
            record["dose_level"] = dose_levels[id(row)]
            record.update(target_taxon=row["target_taxon"], target_feature=target,
                          include=row["include"], exclusion_reason=row["exclusion_reason"])
            manifest.append(record)
            profile_key = (row["profiler"], str(Path(row["source_profile"]).resolve()))
            profiles[profile_key] = Path(row["source_profile"])
        abundance = []
        for (profiler, resolved), path in sorted(profiles.items()):
            if not path.is_file() or not path.stat().st_size:
                raise ValueError("missing native profile: {}".format(path))
            for feature, value in sorted(parse_profile(path, profiler).items()):
                abundance.append({"profiler": profiler, "source_profile": resolved,
                                  "feature": feature,
                                  "abundance_fraction": format(value, ".17g")})
        args.outdir.mkdir(parents=True, exist_ok=True)
        manifest_path = args.outdir / "biomarker_profile_manifest.tsv"
        abundance_path = args.outdir / "biomarker_abundance_long.tsv"
        for path, header, data in ((manifest_path, MANIFEST, manifest),
                                   (abundance_path, ABUNDANCE, abundance)):
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
                writer.writeheader(); writer.writerows(data)
        summary = args.outdir / "biomarker_input_summary.tsv"
        summary.write_text("metric\tvalue\nprofile_contexts\t{}\nunique_native_profiles\t{}\nobserved_species_rows\t{}\nstatus\tPASS\n".format(
            len(manifest), len(profiles), len(abundance)), encoding="utf-8")
        checksum = args.outdir / "biomarker_input.sha256"
        checksum.write_text("".join("{}  {}\n".format(sha256(path), path.resolve()) for path in
                                    (args.canonical, args.aliases, manifest_path, abundance_path, summary)),
                            encoding="utf-8")
        (args.outdir / "SUCCESS").write_text("status\tPASS\n", encoding="utf-8")
        print("[PASS] Expanded {} contexts from {} native profiles".format(len(manifest), len(profiles)))
    except (ValueError, OSError, KeyError) as error:
        raise SystemExit("[ERROR] {}".format(error)) from error


if __name__ == "__main__":
    main()
