#!/usr/bin/env Rscript

if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("scales", quietly = TRUE)) {
  message("[SKIP] plotting dependencies unavailable")
  quit(status = 0)
}
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg)),
                                "../.."))
root <- tempfile("disease_robustness_figures.")
dir.create(root, recursive = TRUE)
metrics_path <- file.path(root, "metrics.tsv")
ledger_path <- file.path(root, "ledger.tsv")
outdir <- file.path(root, "out")

metrics <- expand.grid(
  cohort = c("feng", "zeller"),
  analysis_population = c("community", "independent"),
  target_label = c("Bfrag", "Fnuc"),
  profiler = c("kraken2_bracken", "metaphlan4"),
  spike_fraction_target = c(0.0001, 0.001),
  stringsAsFactors = FALSE
)
metrics$contrast <- "CRC_vs_Control"
metrics$q_threshold <- 0.05
metrics$baseline_bystander_biomarkers <- 20
metrics$gained_bystanders <- rep(c(0, 2), length.out = nrow(metrics))
metrics$bystander_retention_rate <- rep(c(0.9, 0.7),
                                        length.out = nrow(metrics))
metrics$bystander_induced_call_rate <- rep(c(0, 0.002),
                                           length.out = nrow(metrics))
write.table(metrics, metrics_path, sep = "\t", quote = FALSE,
            row.names = FALSE)

ledger <- expand.grid(
  cohort = c("feng", "zeller"),
  analysis_population = c("community", "independent"),
  target_label = c("Bfrag", "Fnuc"),
  profiler = c("kraken2_bracken", "metaphlan4"),
  spike_fraction_target = c(0.0001, 0.001),
  transition = c("retained", "lost", "gained"),
  stringsAsFactors = FALSE
)
ledger$contrast <- "CRC_vs_Control"
ledger$q_threshold <- 0.05
ledger$feature_role <- "bystander"
ledger$effect_change <- rep(c(0.01, -0.2, 0.3),
                            length.out = nrow(ledger))
write.table(ledger, ledger_path, sep = "\t", quote = FALSE,
            row.names = FALSE)

script <- file.path(repo,
                    "analysis_v2/scripts/plot_disease_biomarker_robustness.R")
result <- system2("Rscript", c(
  script, "--metrics", metrics_path, "--ledger", ledger_path,
  "--outdir", outdir, "--analysis-status", "DEVELOPMENT_ONLY"
), stdout = TRUE, stderr = TRUE)
status <- attr(result, "status")
if (is.null(status)) status <- 0L
if (status != 0L) stop(paste(result, collapse = "\n"))

required <- c(
  "SUCCESS", "DEVELOPMENT_ONLY.txt", "FIGURE_GUIDE.md",
  "figures/bystander_biomarker_retention.pdf",
  "figures/bystander_induced_call_rate.png",
  "figures/bystander_gained_call_burden.pdf",
  "figures/bystander_retention_atlas.pdf",
  "figures/bystander_transition_effect_change.png",
  "provenance/robustness_figures.sha256"
)
stopifnot(all(file.exists(file.path(outdir, required))))
atlas <- read.delim(file.path(outdir,
                              "figure_source/bystander_retention_atlas.tsv"),
                    check.names = FALSE)
stopifnot(nrow(atlas) == nrow(metrics), !any(grepl("\n", atlas$Panel)))
cat("[PASS] disease-biomarker robustness-figure fixture\n")
