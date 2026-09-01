#!/usr/bin/env python3
"""Summarize native species assignments for the isolated pure-pool audit."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def norm(value: str) -> str:
    return " ".join(value.replace("_", " ").strip().lower().split())


def load_panel(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8") as handle:
        return {row["label"]: row["taxon_name"] for row in csv.DictReader(handle, delimiter="\t")}


def load_aliases(path: Path) -> dict[str, dict[str, set[str]]]:
    result: dict[str, dict[str, set[str]]] = {}
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            tool = norm(row["tool"])
            result.setdefault(tool, {}).setdefault(norm(row["canonical"]), set()).add(norm(row["alias"]))
    return result


def expected_names(canonical: str, tool: str, aliases: dict[str, dict[str, set[str]]]) -> set[str]:
    names = {norm(canonical)}
    names.update(aliases.get(norm(tool), {}).get(norm(canonical), set()))
    return names


def bracken(path: Path) -> list[tuple[str, float]]:
    with path.open(encoding="utf-8") as handle:
        rows = []
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append((row["name"].strip(), 100.0 * float(row["fraction_total_reads"])))
        return rows


def metaphlan(path: Path) -> list[tuple[str, float]]:
    rows = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            species = next((part[3:] for part in fields[0].split("|") if part.startswith("s__")), None)
            if species is not None:
                rows.append((species, float(fields[2])))
    return rows


def kraken_unclassified(path: Path) -> float:
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 6 and fields[3].strip() == "U":
                return float(fields[0].strip())
    raise ValueError(f"No unclassified row in {path}")


def summarize(rows: list[tuple[str, float]], expected: set[str]) -> tuple[float, float, int, str, float]:
    expected_pct = sum(value for name, value in rows if norm(name) in expected)
    off = [(name, value) for name, value in rows if norm(name) not in expected]
    top_name, top_value = max(off, key=lambda item: item[1], default=("", 0.0))
    return expected_pct, sum(value for _, value in off), len(rows), top_name, top_value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-root", required=True, type=Path)
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--aliases", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    panel = load_panel(args.panel)
    aliases = load_aliases(args.aliases)
    output_rows = []
    for label, canonical in panel.items():
        root = args.audit_root / label
        if not (root / "SUCCESS").is_file():
            raise SystemExit(f"[ERROR] Incomplete pure-pool audit: {label}")
        profile_root = root / "profiles" / f"pure_pool_{label}"
        kreport = profile_root / f"pure_pool_{label}.kraken2.report"
        inputs = {
            "kraken2_bracken": bracken(profile_root / f"pure_pool_{label}.bracken.S.tsv"),
            "metaphlan4": metaphlan(profile_root / f"pure_pool_{label}.metaphlan.tsv"),
        }
        for tool, assignments in inputs.items():
            expected_pct, off_pct, count, top_name, top_pct = summarize(
                assignments, expected_names(canonical, tool, aliases)
            )
            output_rows.append({
                "pool_label": label,
                "expected_taxon": canonical,
                "tool": tool,
                "expected_target_pct": f"{expected_pct:.12g}",
                "off_target_species_pct": f"{off_pct:.12g}",
                "reported_species_total_pct": f"{expected_pct + off_pct:.12g}",
                "reported_species_count": count,
                "top_off_target": top_name,
                "top_off_target_pct": f"{top_pct:.12g}",
                "kraken2_unclassified_pct": f"{kraken_unclassified(kreport):.12g}" if tool == "kraken2_bracken" else "NA",
            })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(output_rows[0])
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(output_rows)
    print(f"[PASS] Summarized {len(output_rows)} pure-pool profiler results: {args.output}")


if __name__ == "__main__":
    main()
