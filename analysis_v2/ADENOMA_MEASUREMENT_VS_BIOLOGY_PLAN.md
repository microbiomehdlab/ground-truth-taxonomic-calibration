# Adenoma non-calls: measurement limitation, smaller effect, or insufficient precision?

Status: analysis plan, not a completed inference. The focal directly spiked species are *Fusobacterium nucleatum*, *Parvimonas micra*, *Dialister pneumosintes* (reported as `Allisonella pneumosintes` by Kraken2 + Bracken), and *Peptostreptococcus stomatis*. Analyze each species in Feng, Yachida, and Zeller, separately for Kraken2 + Bracken and MetaPhlAn 4. Retain source feature names and the explicit panel aliases. Do not pool profiler abundance scales.

## Scientific question and interpretation rule

For each Adenoma-vs-Control call that is not BH-significant, distinguish (1) a demonstrably smaller stage effect, (2) an analytical limitation capable of hiding a meaningful effect, and (3) insufficient information. Nonsignificance alone supports none of these conclusions. A controlled spike establishes the response to a known addition in the tested material; it does not independently measure the true native abundance or prove biological absence. A biological-stage interpretation requires measurement adequacy *and* a precise, directly fitted CRC-vs-Adenoma contrast, ideally with orthogonal validation.

## Metric inventory for the paper

| Domain | Metric and operational definition | Why it matters |
| --- | --- | --- |
| Native detection | Positive / total and zero fraction, by species, cohort, profiler, and Control/Adenoma/CRC | Separates absence of reports from small positive values. |
| Native abundance | Full empirical distribution including zeros; positive-only median and IQR reported separately | Reveals group overlap and low-abundance tails without hiding zeros. |
| Native variability | IQR, robust spread, and prevalence heterogeneity by disease group; display individual points where feasible | A group may have a plausible effect but high sample-to-sample variation. |
| Association | Adjusted effect, SE, CI, model n, p, BH q, and feature-filter status for CRC-vs-Control and Adenoma-vs-Control | Distinguishes small estimates from imprecise estimates and filtered features. |
| Direct stage contrast | Fit CRC-vs-Adenoma from the original regression with its coefficient covariance, CI, and multiplicity policy | Do not subtract marginal CIs or infer stage differences from one significant and one nonsignificant comparison. |
| Paired detection rescue | Among samples with unspiked abundance zero, fraction with post-spike abundance positive at each dose | Tests whether a known low addition becomes observable in the same sample. |
| Dose-specific detection | Post-spike positive / tested at 0.01%, 0.05%, and 0.1%, by background condition | Tests whether adenoma background impairs detection. Avoid a precise LOD claim from only three doses and small n. |
| Quantitative recovery | Per-sample recovered added signal / expected added signal, with median and Q1–Q3 | Detection is not accurate quantification. Keep each profiler's reference scale explicit. |
| Recovery failure/tails | Fraction below 0.5, fraction above 1.5, and sample-level distribution | Finds severe under-/over-recovery concealed by a good median. |
| Dose response | Detection and recovery trend across tested doses; departures from monotonicity | Tests whether low-dose results extrapolate sensibly toward native abundance. |
| Native-to-spike gap | Compare native adenoma distribution with the lowest dose meeting a prespecified detection *and* quantitative-recovery criterion | Asks whether the likely native signal lies in a reliably measured range. Call this a tested-dose threshold, not a universal LOD/LOQ. |
| Off-target response | Change in every reported bystander taxon per unit known target addition after expected compositional change is accounted for | Tests signal redistribution, especially among Kraken *Fusobacterium* labels. |
| Combined-clade sensitivity | Prespecified genus/related-species aggregate, compared with the species-level result | Tests whether a coherent signal is split over labels; do not select a sum because it happens to be significant. |
| Sample-background interaction | Recovery/rescue versus native community composition, neighboring taxa, and Control/Adenoma/CRC | Tests sample-dependent interference rather than an intrinsic species limit. |
| Sequencing/pipeline QC | Usable microbial reads, host fraction, read length/quality, batch, failed profiles, database/marker coverage, and profiler settings by disease group | Detects unequal information or reference support across groups. |
| Statistical sensitivity | Simulation/injection of prespecified native-scale group effects through the observed response/error pattern, followed by the frozen biomarker model | Estimates power to recover a meaningful adenoma effect under measured technical variability. Treat as a sensitivity model, not observed truth. |
| Robustness/validation | Same-sample alternative defensible filters/models, held-out cohorts, and targeted independent assay if available | Checks whether the explanation transports beyond one analysis choice or platform. |

For detection criteria, specify in advance the required probability and acceptable recovery-error band; do not choose thresholds after seeing the disease calls. Report denominators and uncertainty intervals, especially where paired spike groups have roughly ten samples. Preserve the original expected-abundance equation and provenance for the chosen spike reference; Kraken read-proportional and MetaPhlAn genome-equivalent quantities are not numerically interchangeable. Distinguish independent single-species spikes from community mixtures. Spiking all clinical groups equally does not create a biological Adenoma-vs-Control contrast; any artificial group-enrichment experiment must state which group received the addition.

## Existing outputs and work still needed

- Already available: `plot_matched_crc_adenoma_species.py` generates matched association, full-baseline, and paired direct-spike sample/summary TSVs plus a four-species SVG. `plot_crc_candidate_mechanisms.py` generates focused Yachida *D. pneumosintes* and Zeller *P. micra* atlases. The response report contains `response_operator.tsv`; its cross-talk SVG displays only strongest bystander edges, so use the full TSV for focal taxa.
- Next targeted analyses: explicit paired zero-to-positive rescue, native-to-reliable-dose alignment, focal target-to-bystander routing, and background-by-recovery comparisons; direct covariance-aware CRC-vs-Adenoma contrasts from the original sample-level model; frozen-model spike-calibrated detectable-effect simulations.
- Requires extra data or experiments: orthogonal native-abundance validation, deeper/finer low-dose series if a precise detection/quantification limit is desired, and any unavailable read-level QC or marker-coverage audit.

## Decision framework

- **Technical limitation supported:** native adenoma values overlap a demonstrably unreliable range; paired low-dose spikes often fail or misallocate signal; and sensitivity analysis shows the observed error can erase a prespecified meaningful effect.
- **Smaller biological stage effect supported:** measurement at native-relevant abundance is adequate, cross-talk/QC do not account for the pattern, a precise direct CRC-vs-Adenoma contrast supports a smaller adenoma effect, and independent measurement agrees.
- **Inconclusive:** confidence intervals, spike sample sizes, dose resolution, or validation cannot distinguish the above. Never relabel this as biological absence.

## Presentation tomorrow: one defensible slide

Do not show all paper metrics. Use the updated four-species matched plot as a backup/detail figure, not a full-slide graphic: the text is too small for a room. For the main slide, use **two contrasting focused examples** from existing data: *F. nucleatum* (MetaPhlAn low-dose recovery can be weak in Zeller while Kraken recovers the spike) and *P. micra* or *D. pneumosintes* (some non-called adenoma comparisons have substantial spike recovery). Show three visual quantities per example: unspiked Adenoma detection/abundance, paired 0.01% recovery or rescue, and adjusted Adenoma effect with CI. If a new composite cannot be made and checked before the talk, show a readable crop of one existing mechanism atlas plus two large numeric callouts; keep the full matched plot in backup.

Suggested title: **Known spikes test—but do not fully explain—missing adenoma biomarker calls**.

Suggested spoken conclusion: **Some profiler-specific low-dose limitations are real, but several adenoma non-calls occur despite measurable spike recovery. The current data do not establish biological absence; direct stage contrasts and spike-calibrated sensitivity are the next tests.**

Avoid claiming that a nonsignificant adenoma q-value proves no adenoma association, that a spike was recovered at an absolute concentration equivalent across profilers, or that the spikes establish a precise universal limit of detection.
