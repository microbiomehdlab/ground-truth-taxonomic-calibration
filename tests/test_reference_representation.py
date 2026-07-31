#!/usr/bin/env python3
"""Dependency-free smoke test for the reference-representation table builder."""

import bz2
import csv
import pickle
import subprocess
import tempfile
from pathlib import Path


TARGETS = [
    ("Bfrag", "Bacteroides fragilis", "GCF_000000001.1"),
    ("Csym", "Clostridium symbiosum", "GCF_000000002.1"),
    ("Dpne", "Dialister pneumosintes", "GCF_000000003.1"),
    ("Fnuc", "Fusobacterium nucleatum subsp. nucleatum", "GCF_000000004.1"),
    ("Hhat", "Hungatella hathewayi", "GCF_000000005.1"),
    ("Pmic", "Parvimonas micra", "GCF_000000006.1"),
    ("Pana", "Peptostreptococcus anaerobius", "GCF_000000007.1"),
    ("Psto", "Peptostreptococcus stomatis", "GCF_000000008.1"),
    ("Porp", "Porphyromonas asaccharolytica", "GCF_000000009.1"),
    ("Pint", "Prevotella intermedia", "GCF_000000010.1"),
]


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    builder = repo / "scripts" / "build_reference_representation_table.py"
    with tempfile.TemporaryDirectory(prefix="reference_audit_test.") as tmp_name:
        tmp = Path(tmp_name)
        panel = tmp / "spike_panel.tsv"
        aliases = tmp / "aliases.csv"
        metadata = tmp / "genomes-all_metadata.tsv"
        sketch = tmp / "all_genomes.msh"
        pkl = tmp / "mpa_vJan25.pkl"
        outdir = tmp / "out"

        with panel.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["label", "taxon_name", "assembly", "fasta", "weight", "url"])
            for index, (label, taxon, assembly) in enumerate(TARGETS, 1):
                fasta = tmp / f"{label}.fa"
                fasta.write_text(f">{label}\nACGTACGTACGTACGT\n")
                writer.writerow([label, taxon, assembly, fasta, 1, ""])

        with aliases.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["canonical", "alias", "tool", "spike_label"])
            for label, taxon, _ in TARGETS:
                alias = "Fusobacterium nucleatum" if label == "Fnuc" else taxon
                writer.writerow([taxon, alias, "metaphlan4", label])

        with metadata.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Genome", "Genome_type", "Genome_accession",
                             "Species_rep", "Lineage"])
            for index, (_, taxon, assembly) in enumerate(TARGETS, 1):
                genome = f"MGYG{index:09d}"
                uhgg_taxon = (
                    "Clostridium_Q symbiosum" if index == 2 else taxon
                )
                writer.writerow([
                    genome, "Isolate", assembly, genome,
                    f"d__Bacteria;p__Test;g__{uhgg_taxon.split()[0]};"
                    f"s__{uhgg_taxon}",
                ])

        taxonomy = {}
        for index, (_, taxon, _) in enumerate(TARGETS, 1):
            species = "Fusobacterium nucleatum" if taxon.startswith(
                "Fusobacterium nucleatum"
            ) else taxon
            taxonomy[
                f"k__Bacteria|p__Test|g__{species.split()[0]}|"
                f"s__{species.replace(' ', '_')}|t__SGB{index}"
            ] = 1000
        with bz2.open(pkl, "wb") as handle:
            pickle.dump({"taxonomy": taxonomy}, handle)
        sketch.write_bytes(b"synthetic sketch placeholder")

        fake_mash = tmp / "mash"
        fake_mash.write_text(
            "#!/usr/bin/env python3\n"
            "import pathlib, sys\n"
            "for i, query in enumerate(sys.argv[3:], 1):\n"
            " print(f'/db/MGYG{i:09d}.fna.gz\\t{query}\\t0.001\\t0\\t900/1000')\n"
        )
        fake_mash.chmod(0o755)

        subprocess.run([
            "python3", str(builder),
            "--spike-panel", str(panel),
            "--aliases", str(aliases),
            "--uhgg-metadata", str(metadata),
            "--uhgg-mash-sketch", str(sketch),
            "--metaphlan-pkl", str(pkl),
            "--mash-bin", str(fake_mash),
            "--outdir", str(outdir),
        ], check=True)

        table = list(csv.DictReader(
            (outdir / "TableA7_target_database_representation.tsv").open(),
            delimiter="\t",
        ))
        assert len(table) == 10
        assert [row["label"] for row in table] == [row[0] for row in TARGETS]
        assert all(row["exact_uhgg_representative_id"] == "yes" for row in table)
        assert all(row["closest_uhgg_mash_similarity_pct"] == "99.9000"
                   for row in table)
        assert table[1]["closest_uhgg_species"] == "Clostridium_Q symbiosum"
        assert table[3]["metaphlan4_sgb"] == "SGB4"
        tex = (outdir / "TableA7_target_database_representation.tex").read_text()
        assert "\\begin{sidewaystable}" in tex
        assert "Exact UHGG representative ID" not in tex
        assert "UHGG species label" in tex
        assert "Mash distance" in tex
        assert "Clostridium\\_Q symbiosum" in tex
        assert "\\textbackslash{}\\%" not in tex
        assert (outdir / "TableA7_target_database_representation_PROVENANCE.txt").stat().st_size > 0
        print("[PASS] Reference-representation table smoke test passed.")


if __name__ == "__main__":
    main()
