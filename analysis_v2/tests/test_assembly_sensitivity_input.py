#!/usr/bin/env python3
"""End-to-end fixture for the assembly-sensitivity canonical-input builder."""

from __future__ import annotations

import csv
import subprocess
import tempfile
from pathlib import Path


def write_tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def write_profiles(root: Path, profile_id: str, taxon: str, abundance: float) -> None:
    root.mkdir(parents=True, exist_ok=True)
    write_tsv(
        root / (profile_id + ".bracken.S.tsv"),
        ["name", "taxonomy_id", "taxonomy_lvl", "kraken_assigned_reads",
         "added_reads", "new_est_reads", "fraction_total_reads"],
        [[taxon, "1", "S", "1", "0", "1", abundance]],
    )
    metaphlan_taxon = taxon.replace(" ", "_")
    (root / (profile_id + ".metaphlan.tsv")).write_text(
        "#mpa_vJan25_CHOCOPhlAnSGB_202503\n"
        "#clade_name\tNCBI_tax_id\trelative_abundance\n"
        "k__Bacteria|p__Test|g__Test|s__{}\t1\t{}\n".format(
            metaphlan_taxon, abundance * 100
        ),
        encoding="utf-8",
    )


def main() -> None:
    repo = Path(__file__).resolve().parents[2]
    script = repo / "analysis_v2/scripts/build_assembly_sensitivity_input.py"
    with tempfile.TemporaryDirectory(prefix="assembly_sensitivity_input.") as name:
        root = Path(name)
        manifest = root / "manifest.tsv"
        arms = root / "arms.tsv"
        panel = root / "panel.tsv"
        aliases = root / "aliases.csv"
        sensitivity = root / "sensitivity"
        original = root / "original"
        baseline = root / "baseline"
        outdir = root / "output"

        write_tsv(manifest, ["sample_id", "Study", "Target_Condition"],
                  [["S1", "Study1", "CRC"]])
        write_tsv(arms, ["arm_label", "target_label", "assembly_accession"], [
            ["Pana_clean_A1", "Pana", "A1"],
            ["Pint_clean_A2", "Pint", "A2"],
        ])
        write_tsv(panel, ["label", "taxon_name"], [
            ["Pana", "Peptostreptococcus anaerobius"],
            ["Pint", "Prevotella intermedia"],
        ])
        with aliases.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["canonical", "alias", "tool"])
            for taxon in ("Peptostreptococcus anaerobius", "Prevotella intermedia"):
                writer.writerow([taxon, taxon, "kraken2_bracken"])
                writer.writerow([taxon, taxon, "metaphlan4"])

        write_profiles(baseline / "results/Study1/S1/profiles/baseline/S1", "S1",
                       "Peptostreptococcus anaerobius", 0.01)
        # Both target rows must coexist in each native baseline profile.
        bracken = baseline / "results/Study1/S1/profiles/baseline/S1/S1.bracken.S.tsv"
        with bracken.open("a", encoding="utf-8") as handle:
            handle.write("Prevotella intermedia\t2\tS\t1\t0\t1\t0.02\n")
        metaphlan = baseline / "results/Study1/S1/profiles/baseline/S1/S1.metaphlan.tsv"
        with metaphlan.open("a", encoding="utf-8") as handle:
            handle.write("k__Bacteria|p__Test|g__Test|s__Prevotella_intermedia\t2\t2.0\n")

        for arm, taxon in (("Pana_clean_A1", "Peptostreptococcus anaerobius"),
                           ("Pint_clean_A2", "Prevotella intermedia")):
            design_rows = []
            for index, fraction in enumerate((0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05), start=1):
                profile_id = "S1_{}_f{}".format(arm, index)
                design_rows.append(["S1", arm, fraction, index * 100, fraction,
                                    str(root / (profile_id + "_1.fastq.gz"))])
                write_profiles(
                    sensitivity / "results/Study1/S1/profiles/independent" / arm / profile_id,
                    profile_id, taxon, min(0.99, 0.01 + fraction),
                )
            write_tsv(
                sensitivity / "results/Study1/S1/spike_design/independent" / (arm + ".tsv"),
                ["sample_id", "label", "fraction", "N_inserted", "f_hat", "output_r1"],
                design_rows,
            )

            original_label = arm.split("_", 1)[0]
            original_rows = []
            for index, fraction in enumerate((0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05), start=1):
                profile_id = "S1_{}_f{}".format(original_label, index)
                original_rows.append(["S1", original_label, fraction, index * 100, fraction,
                                      str(root / (profile_id + "_1.fastq.gz"))])
                write_profiles(
                    original / "results/Study1/S1/profiles/independent" / original_label / profile_id,
                    profile_id, taxon, min(0.99, 0.008 + fraction),
                )
            write_tsv(
                original / "results/Study1/S1/spike_design/independent" / (original_label + ".tsv"),
                ["sample_id", "label", "fraction", "N_inserted", "f_hat", "output_r1"],
                original_rows,
            )

        result = subprocess.run([
            "python3", str(script), "--manifest", str(manifest), "--arms", str(arms),
            "--spike-panel", str(panel), "--aliases", str(aliases),
            "--sensitivity-root", str(sensitivity), "--original-root", str(original),
            "--baseline-root", str(baseline),
            "--outdir", str(outdir),
        ], text=True, capture_output=True)
        assert result.returncode == 0, result.stderr + result.stdout
        with (outdir / "canonical_input.tsv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == 56
        assert sum(float(row["spike_fraction_target"]) == 0 for row in rows) == 8
        assert {row["profiler"] for row in rows} == {"kraken2_bracken", "metaphlan4"}
        assert {row["target_label"] for row in rows} == {"Pana", "Pint"}
        assert {row["assembly_arm"] for row in rows} == {"original", "clean"}
        assert (outdir / "validation/SUCCESS").is_file()

    print("[PASS] assembly-sensitivity input-builder fixture")


if __name__ == "__main__":
    main()
