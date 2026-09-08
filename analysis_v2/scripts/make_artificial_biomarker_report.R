#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  hit <- match(flag, args); if (is.na(hit)) return(default)
  if (hit == length(args)) stop("Missing value for ", flag)
  args[[hit + 1L]]
}
run_root <- value("--paired-run")
secondary_root <- value("--secondary-paired-run")
outdir <- value("--outdir")
report_status <- value("--report-status")
if (is.null(run_root) || is.null(outdir) || !report_status %in% c("DEVELOPMENT_ONLY", "DEFINITIVE"))
  stop("Required: --paired-run DIR --outdir DIR --report-status DEVELOPMENT_ONLY|DEFINITIVE")
run_roots <- c(run_root, if (!is.null(secondary_root)) secondary_root)
required_for <- function(root) c(file.path(root, "SUCCESS"), file.path(root, "models", "SUCCESS"),
              file.path(run_root, "evaluation", "SUCCESS"),
              file.path(run_root, "models", "paired_da_results.tsv"),
              file.path(run_root, "evaluation", "biomarker_propagation_metrics.tsv"))
required <- unlist(lapply(run_roots, function(root) {
  x <- required_for(root); sub(run_root, root, x, fixed=TRUE)
}), use.names=FALSE)
if (any(!file.exists(required)) || any(file.info(required)$size <= 0)) stop("Paired run is incomplete.")
source_development <- any(file.exists(file.path(run_roots, "DEVELOPMENT_ONLY.txt")))
if (report_status == "DEFINITIVE" && source_development)
  stop("A development analysis cannot generate a definitive report.")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
for (subdir in c("tables", "figure_source", "figures", "diagnostics", "provenance"))
  dir.create(file.path(outdir, subdir), showWarnings = FALSE)

scope_for <- function(x) ifelse(x$contrast == "spiked_vs_matched_baseline__pooled",
  "pooled_primary", "phenotype_stratified_secondary")
metrics <- do.call(rbind, lapply(run_roots, function(root) {
  x <- read.delim(file.path(root,"evaluation","biomarker_propagation_metrics.tsv"), check.names=FALSE, stringsAsFactors=FALSE, na.strings="NA")
  x$analysis_scope <- scope_for(x); x
}))
calls <- do.call(rbind, lapply(run_roots, function(root) {
  x <- read.delim(file.path(root,"models","paired_da_results.tsv"), check.names=FALSE, stringsAsFactors=FALSE, na.strings="NA")
  x$analysis_scope <- if ("contrast" %in% names(x)) scope_for(x) else rep("unknown", nrow(x)); x
}))
needed <- c("cohort", "study", "analysis_population", "target_label", "assembly_arm", "profiler",
  "contrast", "spike_fraction_target", "q_threshold", "target_called", "enriched_calls",
  "off_target_enriched_calls", "precision", "recall", "f1", "target_effect", "target_q_value")
if (length(setdiff(needed, names(metrics)))) stop("Propagation metrics lack artificial-biomarker report columns.")
numeric_fields <- c("spike_fraction_target", "q_threshold", "target_called", "enriched_calls",
  "off_target_enriched_calls", "precision", "recall", "f1", "target_effect", "target_q_value")
for (field in numeric_fields) metrics[[field]] <- as.numeric(metrics[[field]])
if (!nrow(metrics) || anyNA(metrics[c("spike_fraction_target", "q_threshold", "target_called")]) ||
    any(metrics$spike_fraction_target <= 0)) stop("Invalid artificial-biomarker metrics.")

# Achieved fractions differ slightly because paired-read counts are integers and
# library sizes differ. Aggregate only on the frozen experimental dose grid,
# while preserving the exact achieved fraction in the detailed evidence table.
nominal_doses <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
nearest_nominal <- vapply(metrics$spike_fraction_target, function(dose) {
  nominal_doses[which.min(abs(nominal_doses - dose))]
}, numeric(1))
tolerance <- pmax(1e-10, nearest_nominal * 0.001)
if (any(abs(metrics$spike_fraction_target - nearest_nominal) > tolerance))
  stop("An achieved fraction does not match the frozen six-dose grid within tolerance.")
metrics$dose_fraction_nominal <- nearest_nominal
metrics$dose_percent_nominal <- 100 * nearest_nominal
metrics$dose_rank <- match(nearest_nominal, nominal_doses)

context_key <- c("analysis_scope", "cohort", "study", "analysis_population", "target_label", "assembly_arm",
                 "profiler", "contrast", "q_threshold")
metrics <- metrics[do.call(order, metrics[c(context_key, "dose_rank")]), ]
write.table(metrics, file.path(outdir, "figure_source", "artificial_biomarker_metrics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

contexts <- split(metrics, interaction(metrics[context_key], drop = TRUE, lex.order = TRUE))
minimum_dose <- do.call(rbind, lapply(contexts, function(x) {
  if (length(unique(x$dose_rank)) != nrow(x)) stop("A report context contains duplicate doses.")
  detected <- x[x$target_called == 1, , drop = FALSE]
  first_detected <- if (nrow(detected)) detected[which.min(detected$dose_rank), , drop = FALSE] else NULL
  sustained <- x[vapply(x$dose_rank, function(rank) all(x$target_called[x$dose_rank >= rank] == 1), logical(1)),,drop=FALSE]
  first_sustained <- if (nrow(sustained)) sustained[which.min(sustained$dose_rank),,drop=FALSE] else NULL
  data.frame(x[1, context_key, drop = FALSE], doses_tested = nrow(x),
    target_detected_any_dose = as.integer(nrow(detected) > 0),
    minimum_detected_nominal_fraction = if (nrow(detected)) first_detected$dose_fraction_nominal else NA_real_,
    minimum_detected_nominal_percent = if (nrow(detected)) first_detected$dose_percent_nominal else NA_real_,
    achieved_fraction_at_first_detection = if (nrow(detected)) first_detected$spike_fraction_target else NA_real_,
    minimum_sustained_detected_nominal_fraction = if (nrow(sustained)) first_sustained$dose_fraction_nominal else NA_real_,
    minimum_sustained_detected_nominal_percent = if (nrow(sustained)) first_sustained$dose_percent_nominal else NA_real_,
    target_detected_all_doses = as.integer(all(x$target_called == 1)),
    maximum_off_target_calls = max(x$off_target_enriched_calls),
    precision_at_first_detection = if (nrow(detected)) first_detected$precision else NA_real_,
    off_target_calls_at_first_detection = if (nrow(detected)) first_detected$off_target_enriched_calls else NA_real_,
    stringsAsFactors = FALSE)
}))
minimum_dose <- minimum_dose[do.call(order, minimum_dose[context_key]), ]
write.table(minimum_dose, file.path(outdir, "tables", "minimum_detectable_dose.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

summary_key <- c("analysis_scope", "cohort", "analysis_population", "profiler", "assembly_arm",
                 "q_threshold", "dose_rank", "dose_fraction_nominal", "dose_percent_nominal")
groups <- split(metrics, interaction(metrics[summary_key], drop = TRUE, lex.order = TRUE))
summary <- do.call(rbind, lapply(groups, function(x) data.frame(
  x[1, summary_key, drop = FALSE], contexts = nrow(x),
  median_achieved_fraction = median(x$spike_fraction_target),
  target_recall = mean(x$target_called), mean_precision = mean(x$precision),
  median_f1 = median(x$f1), context_sum_enriched_calls = sum(x$enriched_calls),
  context_sum_off_target_calls = sum(x$off_target_enriched_calls),
  median_target_effect = median(x$target_effect), stringsAsFactors = FALSE)))
summary <- summary[do.call(order, summary[summary_key]), ]
write.table(summary, file.path(outdir, "tables", "artificial_biomarker_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

labels <- c(kraken2_bracken = "Kraken2 + Bracken", metaphlan4 = "MetaPhlAn 4")
metrics$profiler_display <- unname(labels[metrics$profiler])
primary <- metrics[abs(metrics$q_threshold - .05) < 1e-12, ]
if (any(primary$analysis_scope == "pooled_primary")) primary <- primary[primary$analysis_scope == "pooled_primary",]
if (!nrow(primary)) stop("Primary q <= 0.05 metrics are absent.")
theme_report <- theme_bw(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
p1 <- ggplot(primary, aes(dose_percent_nominal, target_called, color = profiler_display, group = profiler_display)) +
  stat_summary(fun = mean, geom = "line") + stat_summary(fun = mean, geom = "point") +
  facet_grid(target_label ~ assembly_arm) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25)) +
  labs(x = "Implanted target fraction (%)", y = "Artificial-target recall", color = "Profiler") + theme_report
p2 <- ggplot(primary, aes(dose_percent_nominal, precision, color = profiler_display, group = profiler_display)) +
  stat_summary(fun = mean, geom = "line") + stat_summary(fun = mean, geom = "point") +
  facet_grid(target_label ~ assembly_arm) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25)) +
  labs(x = "Implanted target fraction (%)", y = "Mean precision of enriched calls", color = "Profiler") + theme_report
p3 <- ggplot(primary, aes(dose_percent_nominal, off_target_enriched_calls, color = profiler_display,
                          group = interaction(profiler_display, contrast))) +
  geom_line(alpha = .45) + geom_point(size = 1) + facet_grid(target_label ~ assembly_arm) +
  labs(x = "Implanted target fraction (%)", y = "Off-target enriched calls", color = "Profiler") + theme_report
p4 <- ggplot(primary, aes(dose_percent_nominal, target_effect, color = profiler_display,
                          group = interaction(profiler_display, contrast))) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") + geom_line(alpha = .45) +
  geom_point(size = 1) + facet_grid(target_label ~ assembly_arm) +
  labs(x = "Implanted target fraction (%)", y = "Paired target log2 abundance change", color = "Profiler") + theme_report
for (item in list(list("artificial_target_recall", p1), list("artificial_target_precision", p2),
                  list("off_target_discovery_burden", p3), list("artificial_target_effect", p4))) {
  name <- item[[1]]; plot <- item[[2]]
  ggsave(file.path(outdir, "figures", paste0(name, ".pdf")), plot, width = 8, height = 6, units = "in")
  ggsave(file.path(outdir, "figures", paste0(name, ".png")), plot, width = 8, height = 6,
         units = "in", dpi = 300)
}

diagnostics <- data.frame(metric = c("metric_rows", "model_rows", "contexts", "cohorts", "targets",
  "assembly_arms", "profilers", "q_thresholds", "report_status"), value = c(nrow(metrics), nrow(calls),
  nrow(minimum_dose), length(unique(metrics$cohort)), length(unique(metrics$target_label)),
  length(unique(metrics$assembly_arm)), length(unique(metrics$profiler)),
  length(unique(metrics$q_threshold)), report_status))
write.table(diagnostics, file.path(outdir, "diagnostics", "report_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Preserve the identities of significant non-target taxa, not only their count.
call_required <- c("cohort","study","analysis_population","target_label","assembly_arm","profiler",
                   "contrast","spike_fraction_target","feature","effect","q_value")
if (!length(setdiff(call_required,names(calls)))) {
  calls$spike_fraction_target <- as.numeric(calls$spike_fraction_target)
  calls$q_value <- as.numeric(calls$q_value); calls$effect <- as.numeric(calls$effect)
  calls$dose_fraction_nominal <- map <- vapply(calls$spike_fraction_target, function(d) nominal_doses[which.min(abs(nominal_doses-d))], numeric(1))
  alias_key <- unique(metrics[c("analysis_scope","cohort","study","analysis_population","target_label","assembly_arm","profiler","contrast","target_alias")])
  calls <- merge(calls,alias_key,by=c("analysis_scope","cohort","study","analysis_population","target_label","assembly_arm","profiler","contrast"),all.x=TRUE)
  ledger <- do.call(rbind,lapply(sort(unique(metrics$q_threshold)),function(q) {
    x <- calls[!is.na(calls$q_value)&calls$q_value<=q&calls$effect>0&calls$feature!=calls$target_alias,,drop=FALSE]
    x$q_threshold <- q; x
  }))
  keep <- c("analysis_scope","cohort","study","analysis_population","target_label","assembly_arm","profiler","contrast","dose_fraction_nominal","q_threshold","feature","effect","q_value")
  ledger <- ledger[keep]
  write.table(ledger,file.path(outdir,"tables","off_target_call_ledger.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
  if (nrow(ledger)) {
    rk <- c("analysis_scope","cohort","analysis_population","target_label","assembly_arm","profiler","q_threshold","feature")
    recurrence <- do.call(rbind,lapply(split(ledger,interaction(ledger[rk],drop=TRUE,lex.order=TRUE)),function(x)
      data.frame(x[1,rk,drop=FALSE],contexts_called=nrow(x),contrasts_called=length(unique(x$contrast)),
                 minimum_nominal_fraction=min(x$dose_fraction_nominal),maximum_nominal_fraction=max(x$dose_fraction_nominal))))
    recurrence <- recurrence[order(-recurrence$contexts_called,recurrence$feature),]
  } else recurrence <- data.frame()
  write.table(recurrence,file.path(outdir,"tables","recurrent_off_target_taxa.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
}
writeLines(c("# Draft figure captions", "",
  "These captions inherit the status recorded in `provenance/report_manifest.tsv`.", "",
  "## Artificial-target recall", "Proportion of matched phenotype-background contexts in which the implanted target was significantly enriched in spiked versus unmodified profiles at BH q <= 0.05.", "",
  "## Artificial-target precision", "Mean fraction of enriched differential-abundance calls corresponding to the prespecified implanted target. Zero-call contexts have precision zero.", "",
  "## Off-target discovery burden", "Number of non-target species called significantly enriched after controlled read implantation. Lines show phenotype-background contexts and are descriptive.", "",
  "## Artificial-target effect", "Mean within-sample log2 change in profiler-native target abundance after controlled read implantation relative to the matched unmodified library. This is sequencing-perturbation recovery, not cellular-abundance recovery."),
  file.path(outdir, "captions.md"))
manifest <- data.frame(field = c("status", "analysis", "source_analysis", "secondary_source_analysis", "source_analysis_status", "created_at"),
  value = c(report_status, "paired_artificial_biomarker_recovery", normalizePath(run_root),
            if(is.null(secondary_root)) "NONE" else normalizePath(secondary_root),
            if (source_development) "DEVELOPMENT_ONLY" else "DEFINITIVE",
            format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
write.table(manifest, file.path(outdir, "provenance", "report_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
files <- c(required, list.files(outdir, recursive = TRUE, full.names = TRUE))
status <- system2("sha256sum", files, stdout = file.path(outdir, "provenance", "report.sha256"))
if (!identical(status, 0L)) stop("Could not seal report.")
writeLines(c("report\tpaired_artificial_biomarker_recovery", paste0("report_status\t", report_status),
             "status\tPASS"), file.path(outdir, "SUCCESS"))
message("[PASS] Artificial-biomarker report completed: ", outdir)
