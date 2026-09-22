# Institute presentation: provisional plot set

These are talk-sized previews from the September 2026 finished-sample snapshot,
not final manuscript results. Slides 3–5 use the already-generated native
abundance and corrected paired endpoints; they do **not** wait for the
27.6-million-row response analyzer. The primary MetaPhlAn 4 recovery plot
uses the genome-equivalent expected abundance. A separate labelled sensitivity
plot compares that reference with the original read-proportional equation;
the latter is not a second measurement or the preferred estimand.

André runs this on lobo, preferably in a compute allocation:

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
git pull --ff-only origin revised-analysis-v2
FIGURE_ROOT="$PWD/work/three_cohort_recoverability_20260921T154031Z"
ANALYSIS_SIF="/mnt/beegfs/apptainer/images/ground_truth_analysis_v2_1.sif"
export FIGURE_ROOT ANALYSIS_SIF
OUTDIR="$FIGURE_ROOT/institute_presentation_$(date -u +%Y%m%dT%H%M%SZ)" \
  bash analysis_v2/run_institute_presentation_panels.sh
```

The command prints the new slides directory. It contains PNG, PDF, source TSV,
input/output checksum manifest, `DEVELOPMENT_ONLY.txt`, and `SUCCESS`. The
10-taxon original-reference atlas is in the sibling
`slide4_original_reference_atlas/` directory. Review plots visually before
putting them in a presentation.

| Slide | Asset | What it establishes |
| --- | --- | --- |
| 1 | No new asset; replace the unreadable literature heatmap only after verifying the underlying AUC values. | Adenoma detection is difficult to transfer across cohorts, not that no biomarker exists. |
| 2 | Existing study-design schematic. | Known in-silico input in real metagenomes; verify full-design versus completed-snapshot counts. |
| 3 | `slide3_baseline_visibility.png` | Unspiked baseline visibility for B. fragilis and F. nucleatum; sample-positive fractions with Wilson intervals. A reported zero is not proof of biological absence. |
| 4 | `slide4_focused_recovery.png` | B. fragilis, D. pneumoniae and F. nucleatum at three weak doses; points are positive-only medians and bars are sample IQR. Source table separately counts zero reports. |
| 4 original-reference atlas | `slide4_original_reference_atlas/three_cohort_quantitative_recovery_boxes.png` | All 10 implanted taxa, three doses, both profilers, three cohorts and clinical conditions. Both profilers use the original read-proportional expectation `E = (1-F)b + f`; boxes show positive-sample Q1–Q3. This is a labelled sensitivity representation, not the corrected MetaPhlAn primary estimand. |
| 4 sensitivity | `slide4_metaphlan_reference_comparison.png` | Same MetaPhlAn 4 observations with original and corrected expected-abundance equations, shown in separate rows. This isolates the effect of the reference scale; do not describe the original row as corrected recovery. |
| 5 | `slide5_fnuc_detection.png` | Yachida-only fraction reporting F. nucleatum with Wilson intervals by profiler, dose and clinical group. A positive report can include native baseline signal; a uniform spike is not evidence of a CRC association. |
| 6 | `biomarker_stress_test_feng_0p1/biomarker_stress_test.svg` | Feng-only technical stress test of baseline CRC calls; not biological truth or a three-cohort result. |

Slides 3–5 use finished-sample data and remain subject to the Feng/Zeller
denominator audit. The baseline bars and corrected recovery come from different
source tables by design. None of these figures changes the running response
analysis or the frozen disease-model calls.
