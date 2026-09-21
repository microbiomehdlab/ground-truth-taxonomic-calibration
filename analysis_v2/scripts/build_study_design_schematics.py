#!/usr/bin/env python3
"""Draw count-free, development-only vector schematics for Figure 1."""
import argparse
import hashlib
import html
from pathlib import Path


def box(x, y, width, title, detail, color):
    esc = html.escape
    return (f'<rect x="{x}" y="{y}" width="{width}" height="95" rx="12" '
            f'fill="{color}" stroke="#728494"/>'
            f'<text x="{x + 17}" y="{y + 32}" class="head">{esc(title)}</text>'
            f'<text x="{x + 17}" y="{y + 62}" class="body">{esc(detail)}</text>')


def svg(body, height=500):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="{height}" '
            f'viewBox="0 0 1400 {height}"><style>'
            '.title{font:bold 29px sans-serif;fill:#203540}'
            '.head{font:bold 19px sans-serif;fill:#203540}'
            '.body{font:16px sans-serif;fill:#36505b}'
            '.small{font:14px sans-serif;fill:#526873}'
            '.formula{font:23px sans-serif;fill:#203540}'
            '.arrow{stroke:#527687;stroke-width:3;fill:none;marker-end:url(#tip)}'
            '</style><defs><marker id="tip" viewBox="0 0 10 10" refX="9" refY="5" '
            'markerWidth="8" markerHeight="8" orient="auto-start-reverse">'
            '<path d="M0 0 L10 5 L0 10z" fill="#527687"/></marker></defs>'
            f'<rect width="1400" height="{height}" fill="white"/>{body}</svg>\n')


def build(outdir):
    if outdir.exists():
        raise ValueError(f"output directory exists: {outdir}")
    outdir.mkdir(parents=True)
    a = ['<text x="40" y="52" class="title">A. Paired perturbations in real CRC metagenomes</text>',
         box(40, 100, 290, "Unspiked baseline", "Control · Adenoma · CRC", "#e7f2f3"),
         box(415, 100, 250, "In silico spike-in", "Individual or 10-target mix", "#fff0d8"),
         box(750, 100, 260, "Two profilers", "Kraken2 + Bracken / MetaPhlAn 4", "#e8ebfb"),
         box(1095, 100, 260, "Paired output", "Same sample, before / after", "#e6f3e9"),
         '<path class="arrow" d="M335 148 H405"/><path class="arrow" d="M670 148 H740"/>'
         '<path class="arrow" d="M1015 148 H1085"/>',
         '<text x="40" y="275" class="head">Analytical layers</text>',
         box(40, 295, 290, "Detection", "Is the implanted taxon reported?", "#eff4fa"),
         box(385, 295, 290, "Quantitative recovery", "Observed versus expected", "#eff4fa"),
         box(730, 295, 290, "Non-implanted response", "Profiler cross-talk", "#eff4fa"),
         box(1075, 295, 280, "Biomarker stability", "CRC signal under perturbation", "#eff4fa"),
         '<text x="40" y="460" class="small">Conceptual design only · no sample counts or results shown · DEVELOPMENT ONLY</text>']
    b = ['<text x="40" y="52" class="title">B. Expected-abundance references follow the profiler scale</text>',
         box(40, 90, 620, "Kraken2 + Bracken", "Read-proportional abundance", "#e8f3ef"),
         box(740, 90, 620, "MetaPhlAn 4", "Genome-equivalent-like composition", "#e8ebfb"),
         '<text x="70" y="245" class="formula">expected = (1 - F) o + fₜ</text>',
         '<text x="770" y="245" class="formula">qₜ = fₜ Gₑff / Gₜ</text>',
         '<text x="770" y="300" class="formula">Q = Σ qₖ ; D = (1 - F) + Q</text>',
         '<text x="770" y="355" class="formula">expected = [(1 - F) o + qₜ] / D</text>',
         '<text x="70" y="420" class="small">o: unspiked target fraction · F: total implanted read fraction · fₜ: target read fraction</text>',
         '<text x="70" y="448" class="small">Gₜ: measured target genome size · Gₑff: baseline community effective genome size</text>',
         '<text x="70" y="476" class="small">No fitted genome-size constant · DEVELOPMENT ONLY</text>']
    outputs = {"figure_1A_study_design.svg": svg("".join(a)),
               "figure_1B_profiler_estimands.svg": svg("".join(b), 520)}
    for name, content in outputs.items():
        (outdir / name).write_text(content, encoding="utf-8")
    (outdir / "DEVELOPMENT_ONLY.txt").write_text("status\tDEVELOPMENT_ONLY\n")
    (outdir / "source.sha256").write_text("".join(
        f"{hashlib.sha256((outdir / name).read_bytes()).hexdigest()}  {name}\n"
        for name in outputs), encoding="utf-8")
    (outdir / "SUCCESS").write_text("status\tPASS\n")
    print(f"[PASS] two study-design SVGs: {outdir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    try:
        build(args.outdir)
    except (OSError, ValueError) as error:
        raise SystemExit(f"[ERROR] {error}") from error
