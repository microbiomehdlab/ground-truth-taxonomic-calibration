#!/usr/bin/env python3
"""Build a supplementary audit of spike-genome representation in profiler databases.

The script deliberately separates three questions:
  1. Is the spike assembly identifier itself a UHGG species representative?
  2. Which UHGG species representative is closest by Mash distance?
  3. Which named MetaPhlAn 4 SGB(s) carry the target species label?

It never converts Mash distance into ANI. The reported "Mash similarity" is
100 * (1 - Mash distance) and is labelled as such in every output.
"""

from __future__ import annotations

import argparse
import bz2
import csv
import hashlib
import os
import pickle
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Iterable


TARGET_ORDER = ["Bfrag", "Csym", "Dpne", "Fnuc", "Hhat",
                "Pmic", "Pana", "Psto", "Porp", "Pint"]


def fail(message: str) -> "NoReturn":
    raise SystemExit(f"[ERROR] {message}")


def norm_accession(value: str) -> str:
    value = str(value).strip()
    accession = re.search(
        r"(MGYG\d+|GUT[_-]GENOME\d+|GC[AF]_\d+(?:\.\d+)?)",
        value,
        flags=re.I,
    )
    if accession:
        value = accession.group(1)
    else:
        value = Path(value).name
        while re.search(r"\.(fa|fna|fasta|gz)$", value, flags=re.I):
            value = re.sub(r"\.(fa|fna|fasta|gz)$", "", value, flags=re.I)
    return re.sub(r"[^A-Z0-9]", "", value.upper())


def norm_taxon(value: str) -> str:
    value = str(value).strip().lower().replace("_", " ")
    value = re.sub(r"^[a-z]__", "", value)
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def binomial(value: str) -> str:
    parts = norm_taxon(value).split()
    return " ".join(parts[:2])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_delimited(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        sample = handle.read(8192)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        except csv.Error:
            dialect = csv.excel_tab if "\t" in sample.splitlines()[0] else csv.excel
        return list(csv.DictReader(handle, dialect=dialect))


def require_columns(rows: list[dict[str, str]], columns: Iterable[str], label: str) -> None:
    if not rows:
        fail(f"{label} is empty")
    missing = [column for column in columns if column not in rows[0]]
    if missing:
        fail(f"{label} is missing columns: {', '.join(missing)}")


def load_panel(path: Path) -> list[dict[str, str]]:
    rows = read_delimited(path)
    require_columns(rows, ["label", "taxon_name", "assembly", "fasta"], "spike panel")
    by_label = {row["label"]: row for row in rows}
    missing = [label for label in TARGET_ORDER if label not in by_label]
    if missing:
        fail("spike panel is missing targets: " + ", ".join(missing))
    return [by_label[label] for label in TARGET_ORDER]


def load_metaphlan_names(path: Path, panel: list[dict[str, str]]) -> dict[str, set[str]]:
    names = {row["label"]: {row["taxon_name"], binomial(row["taxon_name"])}
             for row in panel}
    if not path.exists():
        return names
    rows = read_delimited(path)
    require_columns(rows, ["canonical", "alias", "tool"], "taxon-alias file")
    canonical_to_label = {norm_taxon(row["taxon_name"]): row["label"] for row in panel}
    for row in rows:
        if norm_taxon(row["tool"]) not in {"metaphlan4", "metaphlan 4"}:
            continue
        label = canonical_to_label.get(norm_taxon(row["canonical"]))
        if label:
            names[label].update({row["canonical"], row["alias"],
                                 binomial(row["canonical"]), binomial(row["alias"])})
    return names


def load_uhgg_metadata(
    path: Path,
) -> tuple[dict[str, dict[str, str]], set[str], set[str]]:
    by_genome: dict[str, dict[str, str]] = {}
    representative_ids: set[str] = set()
    representative_accessions: set[str] = set()
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = ["Genome", "Genome_accession", "Species_rep", "Lineage"]
        missing = [column for column in required if column not in (reader.fieldnames or [])]
        if missing:
            fail("UHGG metadata is missing columns: " + ", ".join(missing))
        for row in reader:
            genome = row["Genome"].strip()
            by_genome[norm_accession(genome)] = row
            if norm_accession(row["Species_rep"]) == norm_accession(genome):
                representative_ids.add(norm_accession(genome))
                representative_accessions.update({
                    norm_accession(genome),
                    norm_accession(row["Genome_accession"]),
                })
    if not representative_ids:
        fail("UHGG metadata contained no rows where Genome equals Species_rep")
    return by_genome, representative_ids, representative_accessions


def resolve_fasta(panel_path: Path, value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = panel_path.parent / path
    return path.resolve()


def closest_representatives(
    mash_bin: str,
    sketch: Path,
    queries: list[Path],
    representative_ids: set[str],
) -> dict[Path, dict[str, object]]:
    missing = [str(path) for path in queries if not path.is_file() or path.stat().st_size == 0]
    if missing:
        fail("missing or empty spike FASTA files: " + ", ".join(missing))
    command = [mash_bin, "dist", str(sketch), *map(str, queries)]
    try:
        process = subprocess.Popen(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        fail(f"Mash executable not found: {mash_bin}")
    if process.stdout is None:
        fail("Mash stdout was unavailable")
    best: dict[Path, dict[str, object]] = {}
    for line in process.stdout:
        line = line.rstrip("\n")
        fields = line.split("\t")
        if len(fields) < 5:
            process.kill()
            fail("unexpected Mash output line: " + line)
        reference, query, distance, pvalue, shared = fields[:5]
        ref_id = norm_accession(reference)
        if ref_id not in representative_ids:
            continue
        query_path = map_query(query, queries)
        if query_path is None:
            continue
        dist = float(distance)
        if query_path not in best or dist < float(best[query_path]["distance"]):
            best[query_path] = {
                "reference": ref_id,
                "reference_raw": reference,
                "distance": dist,
                "pvalue": pvalue,
                "shared_hashes": shared,
            }
    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code != 0:
        fail("Mash failed:\n" + stderr.strip())
    missing_queries = [path.name for path in queries if path not in best]
    if missing_queries:
        fail("Mash returned no UHGG species-representative match for: "
             + ", ".join(missing_queries))
    return best


def map_query(query_field: str, query_paths: list[Path]) -> Path | None:
    key = norm_accession(query_field)
    candidates = [path for path in query_paths
                  if key in {norm_accession(path.name), norm_accession(str(path))}]
    if len(candidates) == 1:
        return candidates[0]
    query_name = Path(query_field).name
    candidates = [path for path in query_paths if path.name == query_name]
    return candidates[0] if len(candidates) == 1 else None


def load_pickle(path: Path):
    with path.open("rb") as raw:
        magic = raw.read(3)
    opener = bz2.open if magic == b"BZh" else open
    with opener(path, "rb") as handle:
        return pickle.load(handle)


def metaphlan_sgb_map(path: Path, names: dict[str, set[str]]) -> dict[str, list[str]]:
    database = load_pickle(path)
    taxonomy = database.get("taxonomy")
    if not isinstance(taxonomy, dict):
        fail("MetaPhlAn pickle does not contain a dictionary named 'taxonomy'")
    result = {label: set() for label in names}
    for lineage in taxonomy:
        text = str(lineage)
        sgb_match = re.search(r"(?:^|[|;])t__(SGB\d+)", text)
        species_match = re.search(r"(?:^|[|;])s__([^|;]+)", text)
        if not sgb_match or not species_match:
            continue
        species = norm_taxon(species_match.group(1))
        species_binomial = binomial(species)
        for label, candidates in names.items():
            normalized = {norm_taxon(item) for item in candidates if item}
            binomials = {binomial(item) for item in candidates if item}
            if species in normalized or species_binomial in binomials:
                result[label].add(sgb_match.group(1))
    return {label: sorted(values, key=lambda x: int(re.sub(r"\D", "", x)))
            for label, values in result.items()}


def lineage_species(lineage: str) -> str:
    match = re.search(r"(?:^|;)s__([^;]+)", lineage)
    # Preserve GTDB/UHGG suffixes such as Clostridium_Q and hathewayi_A.
    # These underscores are part of the database taxon label; latex_escape()
    # handles them in the generated table.
    return match.group(1).strip() if match else "Unresolved"


def latex_escape(value: object) -> str:
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#",
        "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in text)


def write_outputs(rows: list[dict[str, object]], outdir: Path, provenance: list[str]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    for suffix, delimiter in [("tsv", "\t"), ("csv", ",")]:
        with (outdir / f"TableA7_target_database_representation.{suffix}").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter)
            writer.writeheader()
            writer.writerows(rows)

    tex = outdir / "TableA7_target_database_representation.tex"
    with tex.open("w", encoding="utf-8") as handle:
        handle.write(
            "\\begin{sidewaystable}[p]\n\\centering\n\\scriptsize\n"
            "\\setlength{\\tabcolsep}{4pt}\n"
            "\\caption{\\textbf{Representation of implanted target genomes in the profiler reference databases.} "
            "For each target genome, the closest UHGG~v2.0.2 species representative was identified using Mash. "
            "Mash-based similarity was calculated as $100\\times(1-D)$, where $D$ is the Mash distance, and is "
            "provided as an intuitive summary rather than as a direct measurement of average nucleotide identity. "
            "MetaPhlAn~4 SGBs were identified by matching the target species or its harmonized taxonomic synonym "
            "to the vJan25 database taxonomy. Multiple SGBs indicate that multiple SGBs carried the corresponding "
            "species-level label and do not represent multiple sequence-based assignments of the spike-in assembly.}\n"
            "\\label{tab:target_database_representation}\n"
            "\\begin{tabularx}{\\textheight}{"
            "l p{0.13\\textheight} p{0.14\\textheight} p{0.18\\textheight} "
            "p{0.08\\textheight} p{0.10\\textheight} Y}\n"
            "\\toprule\n"
            "\\textbf{Target} & \\textbf{Spike-in assembly} & "
            "\\textbf{Closest UHGG representative} & "
            "\\textbf{UHGG species label} & "
            "\\textbf{Mash distance} & "
            "\\textbf{Mash-based similarity} & "
            "\\textbf{MetaPhlAn vJan25 SGB(s)} \\\\\n"
            "\\midrule\n"
        )
        for row in rows:
            sgb = row["metaphlan4_sgb"] or "No named species-level SGB match"
            formatted_sgb = "; ".join(
                f"\\texttt{{{latex_escape(value.strip())}}}"
                for value in str(sgb).split(";")
            )
            handle.write(
                f"{latex_escape(row['label'])} & "
                f"\\texttt{{{latex_escape(row['spike_assembly'])}}} & "
                f"\\texttt{{{latex_escape(row['closest_uhgg_genome'])}}} & "
                f"\\textit{{{latex_escape(row['closest_uhgg_species'])}}} & "
                f"{float(row['closest_uhgg_mash_distance']):.4f} & "
                f"{float(row['closest_uhgg_mash_similarity_pct']):.2f}\\% & "
                f"{formatted_sgb} \\\\\n"
            )
        handle.write("\\bottomrule\n\\end{tabularx}\n\\end{sidewaystable}\n")
    (outdir / "TableA7_target_database_representation_PROVENANCE.txt").write_text(
        "\n".join(provenance) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spike-panel", required=True, type=Path)
    parser.add_argument("--aliases", required=True, type=Path)
    parser.add_argument("--uhgg-metadata", required=True, type=Path)
    parser.add_argument("--uhgg-mash-sketch", required=True, type=Path)
    parser.add_argument("--metaphlan-pkl", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--mash-bin", default="mash")
    args = parser.parse_args()

    for label, path in [
        ("spike panel", args.spike_panel), ("aliases", args.aliases),
        ("UHGG metadata", args.uhgg_metadata),
        ("UHGG Mash sketch", args.uhgg_mash_sketch),
        ("MetaPhlAn pickle", args.metaphlan_pkl),
    ]:
        if not path.is_file() or path.stat().st_size == 0:
            fail(f"{label} is missing or empty: {path}")

    panel = load_panel(args.spike_panel)
    names = load_metaphlan_names(args.aliases, panel)
    metadata, representative_ids, representative_accessions = load_uhgg_metadata(
        args.uhgg_metadata
    )
    query_paths = [resolve_fasta(args.spike_panel, row["fasta"]) for row in panel]
    closest = closest_representatives(
        args.mash_bin,
        args.uhgg_mash_sketch,
        query_paths,
        representative_ids,
    )
    sgb = metaphlan_sgb_map(args.metaphlan_pkl, names)

    rows: list[dict[str, object]] = []
    for target, fasta in zip(panel, query_paths):
        hit = closest[fasta]
        meta = metadata[str(hit["reference"])]
        assembly_key = norm_accession(target["assembly"])
        exact = assembly_key in representative_accessions
        rows.append({
            "label": target["label"],
            "canonical_target": target["taxon_name"],
            "spike_assembly": target["assembly"],
            "exact_uhgg_representative_id": "yes" if exact else "no",
            "closest_uhgg_genome": meta["Genome"],
            "closest_uhgg_original_accession": meta["Genome_accession"],
            "closest_uhgg_species": lineage_species(meta["Lineage"]),
            "closest_uhgg_mash_distance": f"{float(hit['distance']):.8f}",
            "closest_uhgg_mash_similarity_pct": f"{100 * (1 - float(hit['distance'])):.4f}",
            "mash_p_value": hit["pvalue"],
            "mash_shared_hashes": hit["shared_hashes"],
            "metaphlan4_sgb": "; ".join(sgb[target["label"]]),
        })

    provenance = [
        "Supplementary target-database representation audit",
        f"command={shlex.join(sys.argv)}",
        f"uhgg_release=UHGG v2.0.2",
        "uhgg_similarity_metric=100 * (1 - Mash distance); not ANI",
        "uhgg_search_space=species representatives identified by Genome == Species_rep",
        "metaphlan_release_expected=vJan25 / ChocoPhlAnSGB 202503",
        f"spike_panel_sha256={sha256(args.spike_panel)}",
        f"aliases_sha256={sha256(args.aliases)}",
        f"uhgg_metadata_sha256={sha256(args.uhgg_metadata)}",
        f"uhgg_mash_sketch_sha256={sha256(args.uhgg_mash_sketch)}",
        f"metaphlan_pkl_sha256={sha256(args.metaphlan_pkl)}",
    ]
    write_outputs(rows, args.outdir, provenance)
    print(f"[PASS] Wrote target-database representation table under: {args.outdir}")


if __name__ == "__main__":
    main()
