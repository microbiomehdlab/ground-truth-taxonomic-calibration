# Provisional paper figure book

This is a planning and review artifact, **not** a manuscript-ready figure set.
`PROVISIONAL_FIGURE_INVENTORY.tsv` specifies six main figures (two panels
each) and seven supplementary candidates. Its panel state is explicit:

- `checked`: the asset and seal must exist; this means only that the specified
  input gate was checked, not that the result is publication-ready.
- `provisional`: an available preview is copied with a visible label and hash;
  absent previews become `MISSING PREVIEW` placeholders.
- `placeholder`: no figure is supplied or silently substituted.

The book currently marks no scientific panel `checked`. Figure 2B and
supplementary S2–S5 point to the corrected three-cohort development figures;
all remain `provisional` pending review. Figure 1 needs a designed schematic
and frozen sample counts. Figures 3–6 need the corrected three-cohort analysis
and/or model gates described in the inventory. The four-taxon baseline figure
is held back until its Feng/Zeller denominator discrepancy is resolved.

André-only cluster command after pulling this code:

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
python3 analysis_v2/scripts/assemble_provisional_figure_book.py \
  --inventory analysis_v2/PROVISIONAL_FIGURE_INVENTORY.tsv \
  --project-root "$PWD" \
  --outdir "work/provisional_figure_book_${STAMP}"
```

Open `index.html` in the new directory. The builder copies previews into its
own `assets/`, writes per-asset SHA-256s to `FIGURE_MANIFEST.tsv`, records the
inventory checksum, and refuses to overwrite an existing output. It does not
run any analysis, profiling, model fitting, or figure-generation script.
Missing inputs remain plainly visible. Update the TSV deliberately as panels
are reviewed; do not promote a panel to `checked` merely because it renders.

Before freezing any publication figure, verify cohort completeness, the
reference scale, the statistical unit, sample/feature denominators, zero and
missing-data policy, uncertainty, and the exact plotted source rows. Rebuild
the book into a new directory after any source or state change.
