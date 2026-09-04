#!/usr/bin/env Rscript
set.seed(73)
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "../.."))
root <- tempfile("paired_biomarker."); dir.create(root, recursive = TRUE)
manifest_path <- file.path(root, "manifest.tsv")
abundance_path <- file.path(root, "abundance.tsv")
outdir <- file.path(root, "model")
samples <- paste0("S", 1:6)
levels <- c("baseline", "dose_01", "dose_02", "dose_03")
manifest <- expand.grid(sample_id = samples, profiler = c("kraken2_bracken", "metaphlan4"),
                        dose_level = levels, stringsAsFactors = FALSE)
manifest$cohort <- "yachida"; manifest$study <- "Study1"
manifest$analysis_population <- "independent"; manifest$condition <- "CRC"
manifest$target_label <- "Pana"; manifest$assembly_arm <- "original"
manifest$profile_id <- paste(manifest$sample_id, manifest$profiler, manifest$dose_level, sep = "_")
manifest$baseline_profile_id <- manifest$sample_id
manifest$spike_fraction_target <- c(baseline = 0, dose_01 = .001,
                                    dose_02 = .005, dose_03 = .01)[manifest$dose_level] *
  (1 + as.integer(sub("S", "", manifest$sample_id)) * 1e-5)
manifest$source_profile <- file.path(root, paste0(manifest$profile_id, ".tsv"))
manifest$target_taxon <- "Peptostreptococcus anaerobius"
manifest$target_feature <- "Peptostreptococcus anaerobius"
manifest$include <- 1; manifest$exclusion_reason <- ""
write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)

features <- c("Peptostreptococcus anaerobius", "Stable species", "Rare species")
abundance <- do.call(rbind, lapply(seq_len(nrow(manifest)), function(i) {
  row <- manifest[i, ]
  dose <- as.numeric(row$spike_fraction_target)
  sample_noise <- as.integer(sub("S", "", row$sample_id)) * 1e-5
  values <- c(0.0001 + dose * 2 + sample_noise, 0.02 + sample_noise,
              ifelse(row$sample_id == "S1" && row$dose_level == "dose_03", .001, 0))
  data.frame(profiler = row$profiler, source_profile = row$source_profile,
             feature = features, abundance_fraction = values)
}))
write.table(abundance, abundance_path, sep = "\t", quote = FALSE, row.names = FALSE)
script <- file.path(repo, "analysis_v2/scripts/fit_paired_biomarker_models.R")
result <- system2("Rscript", c(script, "--profile-manifest", manifest_path,
                               "--abundance-long", abundance_path, "--outdir", outdir),
                  stdout = TRUE, stderr = TRUE)
status <- attr(result, "status"); if (is.null(status)) status <- 0L
if (status != 0L) stop(paste(result, collapse = "\n"))
stopifnot(file.exists(file.path(outdir, "SUCCESS")))
calls <- read.delim(file.path(outdir, "paired_da_results.tsv"))
included <- calls[calls$include == 1, ]
stopifnot(all(c(0, .001000035, .005000175, .01000035) %in%
                round(unique(included$spike_fraction_target), 9)))
stopifnot(all(included$q_value >= 0 & included$q_value <= 1))
stopifnot(all(included$feature != "Rare species"))
stopifnot(any(calls$feature == "Rare species" & calls$include == 0))
stopifnot(nrow(included[included$feature == "Peptostreptococcus anaerobius" &
                          included$spike_fraction_target > 0, ]) == 6)
evaluation <- file.path(root, "evaluation")
eval_script <- file.path(repo, "analysis_v2/scripts/evaluate_biomarker_propagation.py")
eval_result <- system2("python3", c(eval_script, "--calls",
                                    file.path(outdir, "paired_da_results.tsv"),
                                    "--aliases", file.path(repo, "examples/spike_taxon_aliases.csv"),
                                    "--spike-panel", file.path(repo, "spikes/spike_panel.tsv"),
                                    "--outdir", evaluation), stdout = TRUE, stderr = TRUE)
eval_status <- attr(eval_result, "status"); if (is.null(eval_status)) eval_status <- 0L
if (eval_status != 0L) stop(paste(eval_result, collapse = "\n"))
stopifnot(file.exists(file.path(evaluation, "SUCCESS")))
metrics <- read.delim(file.path(evaluation, "biomarker_propagation_metrics.tsv"))
stopifnot(nrow(metrics) == 12L)
cat("[PASS] paired biomarker-model fixture\n")
