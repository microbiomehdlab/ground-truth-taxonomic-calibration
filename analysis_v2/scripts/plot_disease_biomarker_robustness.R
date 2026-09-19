#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) stop("Missing ", flag)
  args[[hit + 1L]]
}
metrics_path <- value("--metrics")
ledger_path <- value("--ledger")
outdir <- value("--outdir")
analysis_status <- value("--analysis-status")
if (!analysis_status %in% c("DEVELOPMENT_ONLY", "DEFINITIVE")) {
  stop("Required: --analysis-status DEVELOPMENT_ONLY|DEFINITIVE")
}
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("OUTDIR must be new or empty: ", outdir)
}
dir.create(file.path(outdir, "figures"), recursive = TRUE)
dir.create(file.path(outdir, "figure_source"), recursive = TRUE)
dir.create(file.path(outdir, "provenance"), recursive = TRUE)

metrics <- read.delim(metrics_path, check.names = FALSE,
                      stringsAsFactors = FALSE, na.strings = "NA")
ledger <- read.delim(ledger_path, check.names = FALSE,
                     stringsAsFactors = FALSE, na.strings = "NA")
metrics_required <- c(
  "cohort", "analysis_population", "target_label", "profiler", "contrast",
  "spike_fraction_target", "q_threshold", "baseline_bystander_biomarkers",
  "gained_bystanders", "bystander_retention_rate",
  "bystander_induced_call_rate"
)
ledger_required <- c(
  "cohort", "analysis_population", "target_label", "profiler", "contrast",
  "spike_fraction_target", "q_threshold", "feature_role", "transition",
  "effect_change"
)
if (length(setdiff(metrics_required, names(metrics))))
  stop("Corrected bystander metrics are missing.")
if (length(setdiff(ledger_required, names(ledger))))
  stop("Disease transition ledger is incomplete.")

near <- function(x, y) abs(as.numeric(x) - y) <= 1e-12
metrics <- metrics[
  metrics$contrast == "CRC_vs_Control" & near(metrics$q_threshold, .05),
  , drop = FALSE
]
ledger <- ledger[
  ledger$contrast == "CRC_vs_Control" & near(ledger$q_threshold, .05) &
    ledger$feature_role == "bystander",
  , drop = FALSE
]
if (!nrow(metrics) || !nrow(ledger)) stop("No primary CRC q<=0.05 rows.")

profiler_labels <- c(
  kraken2_bracken = "Kraken2 + Bracken",
  metaphlan4 = "MetaPhlAn 4"
)
decorate <- function(data) {
  data$Cohort <- tools::toTitleCase(data$cohort)
  data$Experiment <- tools::toTitleCase(data$analysis_population)
  data$Profiler <- factor(profiler_labels[data$profiler],
                          levels = unname(profiler_labels))
  data$dose_percent <- 100 * as.numeric(data$spike_fraction_target)
  data
}
metrics <- decorate(metrics)
ledger <- decorate(ledger)

# Yachida records the achieved target fraction after integer read allocation,
# whereas legacy Feng/Zeller rows are usually serialized at the nominal dose.
# Collapse small allocation differences back onto the population-specific
# frozen grid so identical experimental doses form one plotting column.
community_doses <- c(.00001, .00005, .0001, .0005, .001, .005, .01)
independent_doses <- c(.0001, .0005, .001, .005, .01, .05)
map_nominal_dose <- function(value, population) {
  grid <- if (population == "community") community_doses else independent_doses
  distance <- abs(grid - as.numeric(value)) / grid
  nearest <- which.min(distance)
  if (!length(nearest) || distance[[nearest]] > .05) return(NA_real_)
  grid[[nearest]]
}
metrics$nominal_dose_fraction <- mapply(
  map_nominal_dose, metrics$spike_fraction_target,
  metrics$analysis_population, USE.NAMES = FALSE
)
if (any(!is.finite(metrics$nominal_dose_fraction))) {
  stop("A robustness row is outside its population-specific nominal dose grid.")
}
metrics$dose_percent <- 100 * metrics$nominal_dose_fraction
metrics$retention <- as.numeric(metrics$bystander_retention_rate)
metrics$induced_rate <- as.numeric(metrics$bystander_induced_call_rate)
metrics$gained_bystanders_numeric <- as.numeric(metrics$gained_bystanders)
ledger$absolute_effect_change <- abs(as.numeric(ledger$effect_change))
ledger$Transition <- factor(
  c(retained = "Retained", retained_sign_flipped = "Retained, sign-flipped",
    lost = "Lost", gained = "Gained")[ledger$transition],
  levels = c("Retained", "Lost", "Gained", "Retained, sign-flipped")
)

theme_paper <- theme_bw(10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = .25, color = "grey88"),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey94"),
    strip.text = element_text(face = "bold"),
    plot.title.position = "plot"
  )
colors <- c("Kraken2 + Bracken" = "#D55E00", "MetaPhlAn 4" = "#0072B2")
color_scale <- scale_color_manual(values = colors, drop = FALSE)
fill_scale <- scale_fill_manual(values = colors, drop = FALSE)
dose_scale <- scale_x_log10(
  breaks = c(.001, .005, .01, .05, .1, .5, 1, 5),
  labels = label_number(accuracy = .001)
)
save_plot <- function(plot, stem, width = 10, height = 6.3) {
  ggsave(file.path(outdir, "figures", paste0(stem, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf)
  ggsave(file.path(outdir, "figures", paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 320, bg = "white")
}
write_source <- function(data, name) {
  write.table(data, file.path(outdir, "figure_source", name), sep = "\t",
              quote = FALSE, row.names = FALSE, na = "NA")
}
summarize_targets <- function(data, field) {
  keys <- c("cohort", "Cohort", "analysis_population", "Experiment",
            "profiler", "Profiler", "dose_percent")
  valid <- data[is.finite(data[[field]]), , drop = FALSE]
  groups <- split(valid, interaction(valid[keys], drop = TRUE,
                                     lex.order = TRUE))
  do.call(rbind, lapply(groups, function(group) {
    values <- group[[field]]
    data.frame(group[1, keys, drop = FALSE], n_targets = nrow(group),
               median = median(values), q25 = unname(quantile(values, .25)),
               q75 = unname(quantile(values, .75)), check.names = FALSE)
  }))
}

# Individual target trajectories remain visible; the thick line and ribbon are
# the target-level median and IQR, so correlated doses are not mistaken for
# independent biological samples.
retention_summary <- summarize_targets(metrics, "retention")
denominator_groups <- split(
  metrics,
  interaction(metrics[c("cohort", "Cohort", "analysis_population",
                        "Experiment")], drop = TRUE, lex.order = TRUE)
)
denominator_annotations <- do.call(rbind, lapply(denominator_groups,
  function(group) {
    profiler_groups <- split(group, group$Profiler)
    labels <- vapply(names(profiler_groups), function(profiler) {
      values <- profiler_groups[[profiler]]$baseline_bystander_biomarkers
      span <- range(as.numeric(values))
      short <- if (profiler == "Kraken2 + Bracken") "K+B" else "MP4"
      value <- if (span[1] == span[2]) span[1] else paste(span, collapse = "-")
      paste0(short, " n=", value)
    }, character(1))
    data.frame(group[1, c("cohort", "Cohort", "analysis_population",
                          "Experiment"), drop = FALSE],
               annotation_x = min(group$dose_percent),
               label = paste(labels, collapse = "; "),
               stringsAsFactors = FALSE)
  }))
write_source(metrics, "bystander_metrics_q005.tsv")
write_source(retention_summary, "bystander_retention_summary.tsv")
write_source(denominator_annotations, "baseline_bystander_denominators.tsv")
p_retention <- ggplot(
  metrics[is.finite(metrics$retention), ],
  aes(dose_percent, retention, color = Profiler,
      group = interaction(Profiler, target_label))
) +
  geom_line(alpha = .22, linewidth = .45) +
  geom_point(alpha = .28, size = 1) +
  geom_ribbon(
    data = retention_summary,
    aes(x = dose_percent, ymin = q25, ymax = q75, fill = Profiler,
        group = Profiler),
    inherit.aes = FALSE, alpha = .13, color = NA
  ) +
  geom_line(
    data = retention_summary,
    aes(dose_percent, median, color = Profiler, group = Profiler),
    inherit.aes = FALSE, linewidth = 1.05
  ) +
  geom_point(
    data = retention_summary,
    aes(dose_percent, median, color = Profiler),
    inherit.aes = FALSE, size = 2
  ) +
  geom_text(
    data = denominator_annotations,
    aes(x = annotation_x, y = .025, label = label), inherit.aes = FALSE,
    hjust = -.02, vjust = 0, size = 2.55, color = "grey30"
  ) +
  facet_grid(Experiment ~ Cohort, scales = "free_x") + dose_scale +
  scale_y_continuous(limits = c(0, 1), labels = label_percent()) +
  color_scale + fill_scale +
  labs(
    x = "Implanted target fraction (%)",
    y = "Baseline bystander CRC-associated taxa retained",
    title = "Native CRC-associated call sets show context-specific fragility",
    subtitle = "Thin lines are implanted-species stress tests; thick lines and IQR ribbons summarize implanted species",
    caption = "The directly implanted taxon is excluded. IQR is not a confidence interval; loss denotes perturbation sensitivity, not proof of a false positive."
  ) + theme_paper
save_plot(p_retention, "bystander_biomarker_retention")

induced_summary <- summarize_targets(metrics, "induced_rate")
write_source(induced_summary, "bystander_induced_call_rate_summary.tsv")
p_induced <- ggplot(
  induced_summary,
  aes(dose_percent, 100 * median, color = Profiler, fill = Profiler,
      group = Profiler)
) +
  geom_ribbon(aes(ymin = 100 * q25, ymax = 100 * q75), alpha = .13,
              color = NA) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  facet_grid(Experiment ~ Cohort, scales = "free") + dose_scale +
  scale_y_continuous(
    trans = pseudo_log_trans(base = 10, sigma = .002),
    labels = label_number(suffix = "%", accuracy = .001)
  ) + color_scale + fill_scale +
  labs(
    x = "Implanted target fraction (%)",
    y = "Perturbation-induced significant bystander call rate",
    title = "Perturbation-induced bystander call rates across contexts",
    subtitle = "Gained non-target calls / eligible non-baseline features; target-level median and IQR",
    caption = "Panels use separate y scales. A low feature-normalized rate can still represent multiple calls; absolute counts are reported separately."
  ) + theme_paper
save_plot(p_induced, "bystander_induced_call_rate")

gained_summary <- summarize_targets(metrics, "gained_bystanders_numeric")
write_source(gained_summary, "bystander_gained_call_count_summary.tsv")
p_gained <- ggplot(
  gained_summary,
  aes(dose_percent, median, color = Profiler, fill = Profiler,
      group = Profiler)
) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = .13, color = NA) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  facet_grid(Experiment ~ Cohort, scales = "free") + dose_scale +
  scale_y_continuous(breaks = pretty_breaks()) + color_scale + fill_scale +
  labs(
    x = "Implanted target fraction (%)",
    y = "Gained significant bystander calls",
    title = "Absolute perturbation-induced call burden complements normalized rates",
    subtitle = "Target-level median and IQR; directly implanted taxa excluded",
    caption = "Panels use separate y scales. Counts are threshold transitions, not adjudicated false positives."
  ) + theme_paper
save_plot(p_gained, "bystander_gained_call_burden")

# A compact stress-test atlas preserves the identity of the implanted species.
heat <- metrics[is.finite(metrics$retention), ]
heat$Dose <- factor(
  format(heat$dose_percent, trim = TRUE, scientific = FALSE),
  levels = format(sort(unique(heat$dose_percent)), trim = TRUE,
                  scientific = FALSE)
)
cohort_order <- c("Feng", "Zeller", "Yachida")
cohort_order <- cohort_order[cohort_order %in% unique(heat$Cohort)]
cohort_order <- c(cohort_order, sort(setdiff(unique(heat$Cohort), cohort_order)))
panel_order <- unlist(lapply(c("Community", "Independent"), function(experiment) {
  unlist(lapply(unname(profiler_labels), function(profiler) {
    paste(cohort_order, experiment, profiler, sep = " | ")
  }))
}))
heat$Panel <- factor(paste(heat$Cohort, heat$Experiment, heat$Profiler,
                           sep = " | "), levels = panel_order)
panel_labels <- setNames(gsub(" \\| ", "\n", panel_order), panel_order)
missing_panel_names <- setdiff(panel_order, as.character(unique(heat$Panel)))
missing_panels <- data.frame(
  Panel = factor(missing_panel_names, levels = panel_order),
  label = rep("Undefined\n(no baseline bystander biomarkers)",
              length(missing_panel_names))
)
write_source(heat, "bystander_retention_atlas.tsv")
p_heat <- ggplot(heat, aes(Dose, target_label, fill = retention)) +
  geom_tile(color = "white", linewidth = .25) +
  geom_text(
    data = missing_panels,
    aes(x = 1, y = 1, label = label), inherit.aes = FALSE,
    color = "grey35", size = 3.2
  ) +
  facet_wrap(vars(Panel), ncol = length(cohort_order), scales = "free_x",
             drop = FALSE,
             labeller = as_labeller(panel_labels)) +
  scale_fill_viridis_c(option = "C", limits = c(0, 1),
                       labels = label_percent()) +
  labs(
    x = "Implanted target fraction (%)", y = "Implanted species",
    fill = "Bystanders\nretained",
    title = "A perturbation atlas identifies context-specific call-set fragility",
    subtitle = "CRC versus Control, BH q <= 0.05; directly implanted taxa excluded",
    caption = "Each cell is one controlled stress test, not an independent cohort estimate."
  ) + theme_paper +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")
save_plot(p_heat, "bystander_retention_atlas",
          max(10.5, 4.2 * length(cohort_order)), 8.2)

# Effect changes are descriptive until the paired ideal-counterfactual model
# supplies an uncertainty-aware test of coefficient distortion.
effect_source <- ledger[
  is.finite(ledger$absolute_effect_change) & !is.na(ledger$Transition),
  , drop = FALSE
]
write_source(effect_source, "bystander_transition_effect_changes.tsv")
p_effect <- ggplot(
  effect_source,
  aes(Transition, absolute_effect_change, fill = Transition)
) +
  geom_boxplot(outlier.shape = NA, width = .68) +
  facet_grid(Experiment ~ Cohort + Profiler, scales = "free_y") +
  scale_y_continuous(trans = pseudo_log_trans(base = 10, sigma = .01)) +
  scale_fill_manual(values = c(
    "Retained" = "#009E73", "Lost" = "#D55E00", "Gained" = "#0072B2",
    "Retained, sign-flipped" = "#CC79A7"
  ), drop = FALSE) +
  labs(
    x = NULL, y = "Absolute change in CRC-v-Control coefficient (log2 scale)",
    title = "Significance transitions mix threshold crossing with effect distortion",
    subtitle = "Distributions are descriptive; the paired counterfactual model supplies the formal distortion test",
    caption = "Panels use separate y scales. Rows are repeated feature-by-target-by-dose transitions, not independent replicates."
  ) + theme_paper + guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p_effect, "bystander_transition_effect_change", 12, 6.5)

if (analysis_status == "DEVELOPMENT_ONLY") {
  writeLines(c(
    "status=DEVELOPMENT_ONLY",
    "use_for_manuscript=NO",
    "direct_implanted_target_excluded=YES",
    "false_positive_interpretation=NOT_ESTABLISHED"
  ), file.path(outdir, "DEVELOPMENT_ONLY.txt"))
}
writeLines(c(
  "# Calibration-aware native biomarker robustness figures", "",
  paste0("Analysis status: `", analysis_status, "`."),
  if (analysis_status == "DEVELOPMENT_ONLY")
    "These figures are engineering outputs and must be regenerated from sealed definitive inputs."
  else
    "These figures inherit definitive status from a sealed non-development disease run.",
  "", "- Retention is target-excluded bystander retention.",
  "- Induced calls use the eligible non-baseline feature denominator.",
  "- Absolute gained-call burden is reported beside the normalized induced-call rate.",
  "- Community and independent baselines use different frozen sample panels; their denominators are not interchangeable.",
  "- IQR ribbons summarize implanted species and are not confidence intervals.",
  "- Lost/gained states are descriptive until paired counterfactual distortion is fitted.",
  "- No panel labels a disappearing native disease biomarker as a false positive."
), file.path(outdir, "FIGURE_GUIDE.md"))
files <- c(normalizePath(metrics_path), normalizePath(ledger_path),
           list.files(outdir, recursive = TRUE, full.names = TRUE))
status <- system2("sha256sum", files,
                  stdout = file.path(outdir, "provenance", "robustness_figures.sha256"))
if (!identical(status, 0L)) stop("Could not checksum robustness figures.")
writeLines(c(
  "analysis\tdisease_biomarker_robustness_figures",
  paste0("analysis_status\t", analysis_status),
  "status\tPASS"
), file.path(outdir, "SUCCESS"))
message("[PASS] Disease-biomarker robustness figures: ", normalizePath(outdir))
