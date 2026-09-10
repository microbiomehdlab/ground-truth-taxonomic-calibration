#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit)) return(default)
  if (hit == length(args)) stop("Missing value for ", flag)
  args[[hit + 1L]]
}
paths <- function(flag) strsplit(value(flag, ""), ",", fixed = TRUE)[[1]]
artificial_paths <- paths("--artificial-reports")
disease_paths <- paths("--disease-reports")
linkage_paths <- paths("--linkage-reports")
expected <- strsplit(value("--expected-cohorts", "yachida,feng,zeller"), ",", fixed = TRUE)[[1]]
status <- value("--report-status")
outdir <- value("--outdir")
synthesis <- value("--synthesis")
assembly_report <- value("--assembly-report")
if (is.null(outdir) || !status %in% c("DEVELOPMENT_ONLY", "DEFINITIVE"))
  stop("Required report inputs, --outdir, and valid --report-status.")
if (any(!nzchar(c(artificial_paths, disease_paths, linkage_paths)))) stop("Empty report input path.")
if (dir.exists(outdir) && length(list.files(outdir))) stop("OUTDIR must be new or empty.")

rbind_fill <- function(parts) {
  fields <- unique(unlist(lapply(parts, names), use.names = FALSE))
  parts <- lapply(parts, function(x) {
    for (field in setdiff(fields, names(x))) x[[field]] <- NA
    x[fields]
  })
  do.call(rbind, parts)
}
read_many <- function(roots, relative) rbind_fill(lapply(roots, function(root) {
  path <- file.path(root, relative)
  if (!file.exists(path) || file.info(path)$size <= 0) stop("Missing report table: ", path)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = "NA")
}))

artificial <- read_many(artificial_paths, "tables/artificial_biomarker_summary.tsv")
disease <- read_many(disease_paths, "tables/disease_biomarker_summary.tsv")
linkage <- read_many(linkage_paths, "figure_source/calibration_biomarker_linkage.tsv")
needed_artificial <- c("cohort", "analysis_population", "profiler", "q_threshold",
  "dose_percent_nominal", "target_recall", "mean_precision", "context_sum_off_target_calls")
needed_disease <- c("cohort", "analysis_population", "profiler", "contrast", "q_threshold",
  "spike_fraction_target", "overall_baseline_retention", "median_jaccard")
needed_linkage <- c("cohort", "analysis_population", "profiler", "q_threshold",
  "dose_percent_nominal", "median_response_ratio", "off_target_enriched_calls")
if (length(setdiff(needed_artificial, names(artificial)))) stop("Artificial summary schema mismatch.")
if (length(setdiff(needed_disease, names(disease)))) stop("Disease summary schema mismatch.")
if (length(setdiff(needed_linkage, names(linkage)))) stop("Calibration-linkage schema mismatch.")
observed <- sort(unique(c(artificial$cohort, disease$cohort, linkage$cohort)))
if (!setequal(observed, expected))
  stop("Observed cohorts do not match expected cohorts: ", paste(observed, collapse = ","))

for (directory in c("tables", "figure_source", "figures", "diagnostics", "provenance"))
  dir.create(file.path(outdir, directory), recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write_tsv(artificial, file.path(outdir, "tables", "three_cohort_artificial_biomarkers.tsv"))
write_tsv(disease, file.path(outdir, "tables", "three_cohort_disease_biomarkers.tsv"))
write_tsv(linkage, file.path(outdir, "tables", "three_cohort_calibration_linkage.tsv"))

primary_artificial <- artificial[abs(as.numeric(artificial$q_threshold) - .05) < 1e-12, , drop = FALSE]
if ("analysis_scope" %in% names(primary_artificial) && any(primary_artificial$analysis_scope == "pooled_primary"))
  primary_artificial <- primary_artificial[primary_artificial$analysis_scope == "pooled_primary", , drop = FALSE]
primary_disease <- disease[abs(as.numeric(disease$q_threshold) - .05) < 1e-12, , drop = FALSE]
primary_linkage <- linkage[abs(as.numeric(linkage$q_threshold) - .05) < 1e-12, , drop = FALSE]
if ("analysis_scope" %in% names(primary_linkage) && any(primary_linkage$analysis_scope == "pooled_primary"))
  primary_linkage <- primary_linkage[primary_linkage$analysis_scope == "pooled_primary", , drop = FALSE]

labels <- c(kraken2_bracken = "Kraken2 + Bracken", metaphlan4 = "MetaPhlAn 4")
for (name in c("primary_artificial", "primary_disease", "primary_linkage")) {
  frame <- get(name)
  frame$profiler_display <- unname(labels[frame$profiler])
  assign(name, frame)
}
write_tsv(primary_artificial, file.path(outdir, "figure_source", "artificial_biomarkers_q005.tsv"))
write_tsv(primary_disease, file.path(outdir, "figure_source", "disease_biomarkers_q005.tsv"))
write_tsv(primary_linkage, file.path(outdir, "figure_source", "calibration_linkage_q005.tsv"))

theme_publication <- theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
p1 <- ggplot(primary_artificial, aes(as.numeric(dose_percent_nominal), as.numeric(target_recall),
  color = profiler_display, group = profiler_display)) + geom_line() + geom_point() +
  facet_grid(analysis_population ~ cohort) + scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Implanted target fraction (%)", y = "Artificial-target recall", color = "Profiler") + theme_publication
p2 <- ggplot(primary_disease, aes(100 * as.numeric(spike_fraction_target),
  as.numeric(overall_baseline_retention), color = profiler_display,
  group = interaction(profiler_display, contrast))) + geom_line(alpha = .6) + geom_point() +
  facet_grid(analysis_population ~ cohort) + scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Implanted target fraction (%)", y = "Baseline disease biomarkers retained", color = "Profiler") + theme_publication
p3 <- ggplot(primary_linkage, aes(as.numeric(dose_percent_nominal), as.numeric(median_response_ratio),
  color = profiler_display, group = interaction(profiler_display, target_label))) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") + geom_line(alpha = .55) + geom_point() +
  facet_grid(analysis_population ~ cohort) +
  labs(x = "Implanted target fraction (%)", y = "Median read-perturbation response ratio", color = "Profiler") + theme_publication
p4 <- ggplot(primary_artificial, aes(as.numeric(dose_percent_nominal),
  as.numeric(context_sum_off_target_calls), color = profiler_display,
  group = profiler_display)) + geom_line() + geom_point() + facet_grid(analysis_population ~ cohort) +
  labs(x = "Implanted target fraction (%)", y = "Off-target enriched calls", color = "Profiler") + theme_publication
plots <- list(artificial_target_recall = p1, disease_biomarker_retention = p2,
              calibration_response_ratio = p3, off_target_discovery_burden = p4)
for (name in names(plots)) {
  ggsave(file.path(outdir, "figures", paste0(name, ".pdf")), plots[[name]], width = 9, height = 6)
  ggsave(file.path(outdir, "figures", paste0(name, ".png")), plots[[name]], width = 9, height = 6, dpi = 300)
}

extra_inputs <- character()
if (!is.null(synthesis)) {
  source <- file.path(synthesis, "results/random_effects_meta_analysis.tsv")
  if (!file.exists(source)) stop("Missing cross-cohort synthesis table.")
  file.copy(source, file.path(outdir, "tables", "cross_cohort_disease_meta_analysis.tsv"), overwrite = TRUE)
  extra_inputs <- c(extra_inputs, source)
}
if (!is.null(assembly_report)) {
  marker <- file.path(assembly_report, "SUCCESS")
  if (!file.exists(marker)) stop("Missing assembly-sensitivity report seal.")
  writeLines(normalizePath(assembly_report), file.path(outdir, "provenance", "assembly_sensitivity_report.txt"))
  extra_inputs <- c(extra_inputs, marker)
}
diagnostics <- data.frame(metric = c("cohorts", "artificial_rows", "disease_rows", "linkage_rows", "status"),
  value = c(length(observed), nrow(artificial), nrow(disease), nrow(linkage), status))
write_tsv(diagnostics, file.path(outdir, "diagnostics", "publication_report_diagnostics.tsv"))
writeLines(c("# Draft overview captions", "",
  "All panels preserve cohort-specific estimates and separate independent from community spike experiments.",
  "Artificial-target panels compare spiked profiles with matched unspiked profiles; disease-marker panels measure stability of native phenotype associations after the same controlled perturbation.",
  "Response ratios quantify recovery relative to exact read implantation on each profiler's native abundance scale and are not cellular-abundance estimates."),
  file.path(outdir, "captions.md"))
manifest <- data.frame(field = c("status", "expected_cohorts", "created_at"),
  value = c(status, paste(expected, collapse = ","), format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
write_tsv(manifest, file.path(outdir, "provenance", "report_manifest.tsv"))
if (status == "DEVELOPMENT_ONLY") writeLines(c("status\tDEVELOPMENT_ONLY", "use_for_manuscript\tNO"), file.path(outdir, "DEVELOPMENT_ONLY.txt"))
inputs <- c(unlist(lapply(artificial_paths, function(x) file.path(x, "tables/artificial_biomarker_summary.tsv"))),
  unlist(lapply(disease_paths, function(x) file.path(x, "tables/disease_biomarker_summary.tsv"))),
  unlist(lapply(linkage_paths, function(x) file.path(x, "figure_source/calibration_biomarker_linkage.tsv"))), extra_inputs)
outputs <- list.files(outdir, recursive = TRUE, full.names = TRUE)
rc <- system2("sha256sum", c(inputs, outputs), stdout = file.path(outdir, "provenance", "publication_report.sha256"))
if (!identical(rc, 0L)) stop("Could not seal publication report.")
writeLines(c("report\tthree_cohort_publication_overview", paste0("report_status\t", status), "status\tPASS"), file.path(outdir, "SUCCESS"))
message("[PASS] Three-cohort publication report completed: ", outdir)
