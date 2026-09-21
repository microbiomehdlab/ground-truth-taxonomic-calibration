# Provisional paper figure book

This is a planning and review artifact, **not** a manuscript-ready figure set.
`PROVISIONAL_FIGURE_INVENTORY.tsv` specifies six main figures (two panels
each) and seven supplementary candidates. Its panel state is explicit:

- `checked`: the asset and seal must exist; this means only that the specified
  input gate was checked, not that the result is publication-ready.
- `provisional`: an available preview is copied with a visible label and hash;
  absent previews become `MISSING PREVIEW` placeholders.
- `placeholder`: no figure is supplied or silently substituted.

The book currently marks no scientific panel `checked`. Figure 1A/B are
count-free design schematics, Figure 2A is a descriptive detection heatmap,
Figure 2B is quantitative recovery, and supplementary S2–S5 use the corrected
three-cohort development sources. All remain `provisional` pending review.
Figure 1 still needs final sample counts, and Figure 2A still needs model-based
contrasts and uncertainty for the publication version. Figures 3–6 need the corrected three-cohort analysis
and/or model gates described in the inventory. The four-taxon baseline figure
is held back until its Feng/Zeller denominator discrepancy is resolved.

André-only cluster command after pulling this code:

First produce the newly available Figure 1 and Figure 2A previews from the
already-derived corrected endpoint table (no profiling or fitting):

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
FIGURE_ROOT="$PWD/work/three_cohort_recoverability_20260921T154031Z"
ANALYSIS_SIF="/mnt/beegfs/apptainer/images/ground_truth_analysis_v2_1.sif"
python3 analysis_v2/scripts/build_study_design_schematics.py \
  --outdir "$FIGURE_ROOT/study_design"
apptainer exec --cleanenv --bind "$PWD:$PWD" --pwd "$PWD" "$ANALYSIS_SIF" \
  Rscript analysis_v2/scripts/plot_three_cohort_detection.R \
    --endpoints "$FIGURE_ROOT/endpoints/paired_endpoints.tsv" \
    --outdir "$FIGURE_ROOT/detection"
```

Figure 2A shows the fraction of biological samples with native nonzero target
abundance in each original-assembly community context. Its source TSV gives
`detected`, `samples`, and Wilson 95% intervals. These are descriptive
prevalences, not model-adjusted detection effects. The source gate requires
all cohort-condition-profiler-target-dose contexts, matching native abundance
flags, and corrected MetaPhlAn / unchanged Bracken reference types.

Then assemble a **new** book:

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
