#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  hit <- match(flag, args); if (is.na(hit)) return(default)
  if (hit == length(args)) stop("Missing value for ", flag)
  args[[hit + 1L]]
}
run_root <- value("--disease-run")
outdir <- value("--outdir")
report_status <- value("--report-status")
if (is.null(run_root) || is.null(outdir) || !report_status %in% c("DEVELOPMENT_ONLY", "DEFINITIVE"))
  stop("Required: --disease-run DIR --outdir DIR --report-status DEVELOPMENT_ONLY|DEFINITIVE")
required <- c(file.path(run_root, "SUCCESS"), file.path(run_root, "models", "SUCCESS"),
              file.path(run_root, "evaluation", "SUCCESS"),
              file.path(run_root, "models", "primary_disease_da_results.tsv"),
              file.path(run_root, "evaluation", "disease_biomarker_propagation_metrics.tsv"))
if (any(!file.exists(required)) || any(file.info(required)$size <= 0)) stop("Disease run is incomplete.")
source_development <- file.exists(file.path(run_root, "DEVELOPMENT_ONLY.txt"))
if (report_status == "DEFINITIVE" && source_development)
  stop("A development analysis cannot generate a definitive report.")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
for (subdir in c("tables", "figure_source", "figures", "diagnostics", "provenance"))
  dir.create(file.path(outdir, subdir), showWarnings = FALSE)

metrics <- read.delim(required[5], check.names = FALSE, stringsAsFactors = FALSE, na.strings = "NA")
models <- read.delim(required[4], check.names = FALSE, stringsAsFactors = FALSE)
metric_required <- c("cohort", "analysis_population", "target_label", "assembly_arm", "profiler",
  "contrast", "spike_fraction_target", "q_threshold", "baseline_biomarkers", "dose_biomarkers",
  "retained_biomarkers", "lost_biomarkers", "gained_biomarkers", "baseline_retention_rate",
  "target_effect_change_from_baseline", "biomarker_set_jaccard_vs_baseline")
if (length(setdiff(metric_required, names(metrics)))) stop("Propagation metrics lack report columns.")
for (field in c("spike_fraction_target", "q_threshold", "baseline_biomarkers", "dose_biomarkers",
                "retained_biomarkers", "lost_biomarkers", "gained_biomarkers", "baseline_retention_rate",
                "target_effect_change_from_baseline", "biomarker_set_jaccard_vs_baseline"))
  metrics[[field]] <- as.numeric(metrics[[field]])
if (!nrow(metrics) || anyNA(metrics[c("spike_fraction_target", "q_threshold")]))
  stop("Invalid propagation metrics.")
models$effect <- as.numeric(models$effect); models$standard_error <- as.numeric(models$standard_error)
models$lower_95 <- as.numeric(models$lower_95); models$upper_95 <- as.numeric(models$upper_95)

order_fields <- c("cohort", "analysis_population", "target_label", "assembly_arm", "profiler",
                  "contrast", "q_threshold", "spike_fraction_target")
metrics <- metrics[do.call(order, metrics[order_fields]), ]
write.table(metrics, file.path(outdir, "figure_source", "disease_biomarker_propagation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

summary_key <- c("cohort", "analysis_population", "profiler", "contrast", "q_threshold", "spike_fraction_target")
summary_groups <- split(metrics, interaction(metrics[summary_key], drop = TRUE, lex.order = TRUE))
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
summaries <- do.call(rbind, lapply(summary_groups, function(x) data.frame(
  x[1, summary_key, drop = FALSE], contexts = nrow(x),
  total_baseline_biomarkers = sum(x$baseline_biomarkers),
  total_dose_biomarkers = sum(x$dose_biomarkers),
  total_retained = sum(x$retained_biomarkers),
  mean_baseline_retention = safe_mean(x$baseline_retention_rate),
  overall_baseline_retention = if (sum(x$baseline_biomarkers) > 0)
    sum(x$retained_biomarkers) / sum(x$baseline_biomarkers) else NA_real_,
  total_lost = sum(x$lost_biomarkers), total_gained = sum(x$gained_biomarkers),
  median_jaccard = if (all(is.na(x$biomarker_set_jaccard_vs_baseline))) NA_real_ else
    median(x$biomarker_set_jaccard_vs_baseline, na.rm = TRUE),
  median_target_effect_change = median(x$target_effect_change_from_baseline))))
summaries <- summaries[do.call(order, summaries[summary_key]), ]
write.table(summaries, file.path(outdir, "tables", "disease_biomarker_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

covariate_fields <- intersect(c("cohort", "analysis_population", "target_label", "assembly_arm",
  "profiler", "model_spec", "n_samples", "covariates_used", "covariates_omitted"), names(models))
covariates <- unique(models[covariate_fields])
write.table(covariates, file.path(outdir, "tables", "model_covariate_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

labels <- c(kraken2_bracken = "Kraken2 + Bracken", metaphlan4 = "MetaPhlAn 4")
metrics$profiler_display <- unname(labels[metrics$profiler])
metrics$dose_percent <- 100 * metrics$spike_fraction_target
primary <- metrics[abs(metrics$q_threshold - .05) < 1e-12, ]
theme_report <- theme_bw(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
p1 <- ggplot(primary, aes(dose_percent, baseline_retention_rate, color = profiler_display, group = profiler_display)) +
  stat_summary(fun = mean, geom = "line", na.rm = TRUE) + stat_summary(fun = mean, geom = "point", na.rm = TRUE) +
  facet_grid(contrast ~ assembly_arm) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25)) +
  labs(x = "Implanted target fraction (%)", y = "Baseline disease biomarkers retained", color = "Profiler") + theme_report
p2 <- ggplot(primary, aes(dose_percent, biomarker_set_jaccard_vs_baseline,
                          color = profiler_display, group = profiler_display)) +
  stat_summary(fun = median, geom = "line", na.rm = TRUE) +
  stat_summary(fun = median, geom = "point", na.rm = TRUE) + facet_grid(contrast ~ assembly_arm) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25)) +
  labs(x = "Implanted target fraction (%)", y = "Biomarker-set Jaccard vs baseline", color = "Profiler") + theme_report
p3 <- ggplot(primary, aes(dose_percent, target_effect_change_from_baseline,
                          color = profiler_display, group = interaction(profiler_display, target_label))) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") + geom_line(alpha = .55) + geom_point(size = 1) +
  facet_grid(contrast ~ assembly_arm) +
  labs(x = "Implanted target fraction (%)", y = "Change in disease log2 effect vs baseline", color = "Profiler") + theme_report
for (item in list(list("baseline_biomarker_retention", p1), list("biomarker_set_stability", p2),
                  list("target_effect_change", p3))) {
  name <- item[[1]]; plot <- item[[2]]
  ggsave(file.path(outdir, "figures", paste0(name, ".pdf")), plot, width = 8, height = 5.5, units = "in")
  ggsave(file.path(outdir, "figures", paste0(name, ".png")), plot, width = 8, height = 5.5,
         units = "in", dpi = 300)
}

diagnostics <- data.frame(metric = c("source_rows", "summary_rows", "model_rows", "cohorts", "targets",
  "profilers", "contrasts", "q_thresholds", "missing_jaccard", "report_status"),
  value = c(nrow(metrics), nrow(summaries), nrow(models), length(unique(metrics$cohort)),
    length(unique(metrics$target_label)), length(unique(metrics$profiler)), length(unique(metrics$contrast)),
    length(unique(metrics$q_threshold)), sum(is.na(metrics$biomarker_set_jaccard_vs_baseline)), report_status))
write.table(diagnostics, file.path(outdir, "diagnostics", "report_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
captions <- c(
  "# Draft figure captions", "",
  "These captions inherit the analysis status recorded in `provenance/report_manifest.tsv`.", "",
  "## Baseline disease-biomarker retention", "Proportion of significant baseline disease biomarkers that remained significant after controlled read implantation at BH q <= 0.05. Contexts without baseline biomarkers are undefined and omitted from the summary.", "",
  "## Biomarker-set stability", "Median Jaccard similarity between each perturbed disease-biomarker set and its matched observed baseline call set at BH q <= 0.05. Empty-versus-empty call sets are undefined rather than treated as perfect stability.", "",
  "## Target effect change", "Change from baseline in the target species disease-contrast coefficient after controlled read implantation. Because targets are implanted across phenotype groups, target significance is a spurious-association diagnostic rather than recall. Lines connect dose levels within profiler and target; they are descriptive, not independent replicates.")
writeLines(captions, file.path(outdir, "captions.md"))
manifest <- data.frame(field = c("status", "source_analysis", "source_analysis_status", "created_at"),
  value = c(report_status, normalizePath(run_root), if (source_development) "DEVELOPMENT_ONLY" else "DEFINITIVE",
            format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
write.table(manifest, file.path(outdir, "provenance", "report_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
files <- c(required, list.files(outdir, recursive = TRUE, full.names = TRUE))
status <- system2("sha256sum", files, stdout = file.path(outdir, "provenance", "report.sha256"))
if (!identical(status, 0L)) stop("Could not seal report.")
writeLines(c("report\tdisease_biomarker", paste0("report_status\t", report_status), "status\tPASS"),
           file.path(outdir, "SUCCESS"))
message("[PASS] Disease-biomarker report completed: ", outdir)
