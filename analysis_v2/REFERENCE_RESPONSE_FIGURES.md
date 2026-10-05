# Reference-response review figures

Verified caption addendum: REFERENCE_NE_REVIEW_20261005.md. All 291 MetaPhlAn
constant cases are all-zero native inputs at the lowest pilot dose; all 148
Bracken perfect-fit cases have zero controls and identical positive exposed
values. These are distinct states, not one shared missed-spike category. Existing
sealed figure exports are preserved; apply this caption clarification at paper
integration rather than modifying the archived bundle in place.

Descriptive review, not manuscript acceptance or production authorization.
Generate summaries from a checksum-verified completed reference report:

```bash
python3 -B analysis_v2/scripts/summarize_reference_response.py --report REPORT --out REVIEW
TZ=UTC Rscript analysis_v2/scripts/plot_reference_response.R REVIEW PATH_TO_SHARED_FIGURE_STYLE_R
```

The style file is supplied explicitly from the private manuscript project; no
manuscript files are committed. Profiler colors, font and dimensions follow that
shared theme. PDF/SVG/450-dpi PNG, a combined PDF, and grayscale/protanopia/
deuteranopia review PDFs are generated. Review figures visually before submission.

Each cohort/profiler has five pages, with all ten targets and n=10/20 retained:

1. Observed and reference positive-significance frequencies. Solid filled points
   denote observed, dashed open points reference. Full-family MaAsLin2 q≤0.05
   and positive coefficient are required. The denominator is every allocation,
   including non-estimable cases. Each arm's frequency is calculated from its
   own saved discovery flag, including when the other arm is non-estimable.
2. Separate non-estimability frequencies, distinguishing observed-only,
   reference-only and both NE. NE is not ordinary nonsignificance.
3. Observed-minus-reference coefficient median/IQR across allocations where
   both primary models are estimable. Zero marks agreement. Coefficients use
   log2(1+native abundance/1e-8), not untransformed abundance or a biological
   disease fold change. Negative differences do not alone prove attenuation.
4–5. Observed and reference standard-error median/IQR for the same matched
   estimable set. These are uncertainty diagnostics on the transformed scale.

Dose labels are nominal per-target final-library read-pair percentages, obtained
from the equal-weight ten-member nominal total community grid. Lines join tested
doses only. Actual reference computations use achieved integer-count fractions.
Each cell has 20 conditional allocations of existing people, not 20 independent
studies. IQRs describe allocation distributions, not confidence intervals.
The matched-estimable count is printed as m; excluded NE cases remain on page 2
and in exhaustive category tables. Median curves do not represent one fixed
allocation traced across dose. No new statistical tests or fits are run.

Source tables: frequencies.tsv (all seven primary/HC3 categories, including
zero-count states), significance.tsv (arm-specific primary flags), effects.tsv
(matched effects and standard errors). Input report checksums and summary script
hash are in provenance.json. R session information is stored alongside figures.
Do not treat pooled target/context rows as independent biological observations,
nor interpret the baseline-anchored reference as an unbiased true community.
