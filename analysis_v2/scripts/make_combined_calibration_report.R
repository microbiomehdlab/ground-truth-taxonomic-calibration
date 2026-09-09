#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Missing ", flag)
  args[[i + 1L]]
}
models_root <- value("--models-root")
outdir <- value("--outdir")
status <- value("--report-status")
if (!status %in% c("DEVELOPMENT_ONLY", "DEFINITIVE")) stop("Invalid report status.")
if (dir.exists(outdir) && length(list.files(outdir))) stop("OUTDIR must be new or empty.")

cohorts <- c("feng", "zeller")
populations <- c("independent", "community")
read_outputs <- function(cohort, population) {
  detection_dir <- file.path(models_root, paste0("detection_", cohort, "_", population))
  continuous_dir <- file.path(models_root, paste0("continuous_", cohort, "_", population))
  required <- c(
    file.path(detection_dir, "SUCCESS"),
    file.path(detection_dir, "detection_probability_curve.tsv"),
    file.path(detection_dir, "profiler_contrasts_by_dose.tsv"),
    file.path(detection_dir, "detection_model_diagnostics.tsv"),
    file.path(continuous_dir, "SUCCESS"),
    file.path(continuous_dir, "profiler_response_slopes.tsv"),
    file.path(continuous_dir, "primary_profiler_slope_contrast.tsv"),
    file.path(continuous_dir, "continuous_residuals_by_dose.tsv"),
    file.path(continuous_dir, "continuous_model_diagnostics.tsv")
  )
  if (any(!file.exists(required)) || any(file.info(required)$size <= 0)) {
    stop("Incomplete calibration model inputs for ", cohort, "/", population)
  }
  add_context <- function(path) {
    x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    x$cohort <- cohort
    x$analysis_population <- population
    x
  }
  list(
    curve = add_context(required[2]), detection_contrast = add_context(required[3]),
    detection_diagnostic = add_context(required[4]), slopes = add_context(required[6]),
    slope_contrast = add_context(required[7]), residual = add_context(required[8]),
    continuous_diagnostic = add_context(required[9]), required = required
  )
}

parts <- list()
for (cohort in cohorts) for (population in populations) {
  parts[[paste(cohort, population, sep = "_")]] <- read_outputs(cohort, population)
}
bind <- function(name) do.call(rbind, lapply(parts, `[[`, name))
curve <- bind("curve")
detection_contrast <- bind("detection_contrast")
slopes <- bind("slopes")
slope_contrast <- bind("slope_contrast")
residual <- bind("residual")
detection_diagnostic <- bind("detection_diagnostic")
continuous_diagnostic <- bind("continuous_diagnostic")

for (directory in c("tables", "figure_source", "figures", "diagnostics", "provenance")) {
  dir.create(file.path(outdir, directory), recursive = TRUE, showWarnings = FALSE)
}
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
write_tsv(detection_contrast, file.path(outdir, "tables", "detection_profiler_contrasts.tsv"))
write_tsv(slopes, file.path(outdir, "tables", "quantitative_response_slopes.tsv"))
write_tsv(slope_contrast, file.path(outdir, "tables", "quantitative_profiler_contrasts.tsv"))
write_tsv(detection_diagnostic, file.path(outdir, "diagnostics", "detection_model_diagnostics.tsv"))
write_tsv(continuous_diagnostic, file.path(outdir, "diagnostics", "continuous_model_diagnostics.tsv"))

curve$dose_percent <- 100 * as.numeric(curve$spike_fraction_target)
residual$dose_percent <- 100 * as.numeric(residual$spike_fraction_target)
labels <- c(kraken2_bracken = "Kraken2 + Bracken", metaphlan4 = "MetaPhlAn 4")
curve$profiler_display <- unname(labels[curve$profiler])
slopes$profiler_display <- unname(labels[slopes$profiler])
residual$profiler_display <- unname(labels[residual$profiler])
write_tsv(curve, file.path(outdir, "figure_source", "detection_probability_curves.tsv"))
write_tsv(slopes, file.path(outdir, "figure_source", "quantitative_response_slopes.tsv"))
write_tsv(residual, file.path(outdir, "figure_source", "quantitative_residuals_by_dose.tsv"))

theme_report <- theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
p_detection <- ggplot(curve, aes(dose_percent, predicted_detection_probability,
                                  color = profiler_display, group = profiler_display)) +
  geom_ribbon(aes(ymin = lower_95_descriptive, ymax = upper_95_descriptive,
                  fill = profiler_display), alpha = .12, color = NA) +
  geom_line() + geom_point() + facet_grid(analysis_population ~ cohort) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Implanted target fraction (%)", y = "Predicted detection probability",
       color = "Profiler", fill = "Profiler") + theme_report
p_slopes <- ggplot(slopes, aes(profiler_display, response_slope, color = profiler_display)) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  geom_pointrange(aes(ymin = lower_95, ymax = upper_95), position = position_dodge(width = .4)) +
  facet_grid(analysis_population ~ cohort) +
  labs(x = NULL, y = "Observed/expected response slope", color = "Profiler") + theme_report
p_residual <- ggplot(residual, aes(dose_percent, residual_median, color = profiler_display,
                                    group = profiler_display)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") + geom_line() + geom_point() +
  facet_grid(analysis_population ~ cohort) +
  labs(x = "Implanted target fraction (%)", y = "Median quantitative-model residual",
       color = "Profiler") + theme_report
plots <- list(detection_probability = p_detection, quantitative_response_slopes = p_slopes,
              quantitative_residuals = p_residual)
for (name in names(plots)) {
  ggsave(file.path(outdir, "figures", paste0(name, ".pdf")), plots[[name]], width = 7.5, height = 5)
  ggsave(file.path(outdir, "figures", paste0(name, ".png")), plots[[name]], width = 7.5, height = 5, dpi = 300)
}

writeLines(c(
  "# Draft combined calibration figure captions", "",
  "Detection probability is estimated separately by cohort and spike population from categorical-dose models that include the unspiked baseline.", "",
  "Quantitative response slopes summarize observed relative recovery per unit expected relative perturbation; the dashed line at one represents proportional recovery.", "",
  "Residual summaries diagnose dose-dependent departures from the prespecified linear quantitative model. Development outputs are not manuscript results."
), file.path(outdir, "captions.md"))
manifest <- data.frame(
  field = c("status", "analysis", "models_root", "cohorts", "populations", "created_at"),
  value = c(status, "combined_calibration", normalizePath(models_root),
            paste(cohorts, collapse = ","), paste(populations, collapse = ","),
            format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
write_tsv(manifest, file.path(outdir, "provenance", "report_manifest.tsv"))
inputs <- unique(unlist(lapply(parts, `[[`, "required")))
outputs <- list.files(outdir, recursive = TRUE, full.names = TRUE)
rc <- system2("sha256sum", c(inputs, outputs), stdout = file.path(outdir, "provenance", "report.sha256"))
if (rc != 0) stop("Could not seal combined calibration report.")
writeLines(c("report\tcombined_calibration", paste0("report_status\t", status), "status\tPASS"),
           file.path(outdir, "SUCCESS"))
message("[PASS] Combined calibration report completed: ", outdir)
