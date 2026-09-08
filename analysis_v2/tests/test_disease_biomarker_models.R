#!/usr/bin/env Rscript
set.seed(97)
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "../.."))
root <- tempfile("disease_biomarker."); dir.create(root, recursive = TRUE)
manifest_path <- file.path(root, "manifest.tsv"); abundance_path <- file.path(root, "abundance.tsv")
outdir <- file.path(root, "model")
samples <- paste0("S", 1:30)
condition <- rep(c("Control", "Adenoma", "CRC"), each = 10)
manifest <- expand.grid(sample_id = samples, profiler = c("kraken2_bracken", "metaphlan4"),
                        dose_level = c("baseline", "dose_01", "dose_02"), stringsAsFactors = FALSE)
index <- match(manifest$sample_id, samples)
manifest$cohort <- "synthetic"; manifest$study <- "Study1"; manifest$analysis_population <- "independent"
manifest$condition <- condition[index]; manifest$target_label <- "Pana"; manifest$assembly_arm <- "original"
manifest$profile_id <- paste(manifest$sample_id, manifest$profiler, manifest$dose_level, sep = "_")
manifest$baseline_profile_id <- manifest$sample_id
manifest$spike_fraction_target <- c(baseline = 0, dose_01 = .001, dose_02 = .01)[manifest$dose_level] *
  (1 + index * 1e-6)
manifest$source_profile <- file.path(root, paste0(manifest$profile_id, ".tsv"))
manifest$target_taxon <- "Peptostreptococcus anaerobius"
manifest$target_feature <- "Peptostreptococcus anaerobius"
manifest$age <- 45 + index; manifest$sex <- "Male"
manifest$bmi <- 20 + (index %% 10); manifest$bmi[index %in% c(2, 12)] <- NA
manifest$include <- 1; manifest$exclusion_reason <- ""
second <- manifest
second$target_label <- "Pint"; second$assembly_arm <- "clean"
second$target_taxon <- "Prevotella intermedia"; second$target_feature <- "Second target"
manifest <- rbind(manifest, second)
write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
features <- c("Peptostreptococcus anaerobius", "Second target", "CRC marker", "Constant", "Rare")
profile_rows <- manifest[!duplicated(paste(manifest$profiler, manifest$source_profile)), ]
abundance <- do.call(rbind, lapply(seq_len(nrow(profile_rows)), function(i) {
  row <- profile_rows[i, ]; j <- match(row$sample_id, samples); dose <- as.numeric(row$spike_fraction_target)
  values <- c(.0001 + dose, .0002 + dose,
              .001 + ifelse(condition[j] == "CRC", .01, 0) + runif(1, 0, 1e-4),
              .02, ifelse(j == 1 && row$dose_level == "dose_02", .001, 0))
  data.frame(profiler = row$profiler, source_profile = row$source_profile,
             feature = features, abundance_fraction = values)
}))
write.table(abundance, abundance_path, sep = "\t", quote = FALSE, row.names = FALSE)
script <- file.path(repo, "analysis_v2/scripts/fit_disease_biomarker_models.R")
result <- system2("Rscript", c(script, "--profile-manifest", manifest_path,
                               "--abundance-long", abundance_path, "--outdir", outdir),
                  stdout = TRUE, stderr = TRUE)
status <- attr(result, "status"); if (is.null(status)) status <- 0L
if (status != 0L) stop(paste(result, collapse = "\n"))
stopifnot(file.exists(file.path(outdir, "SUCCESS")))
primary <- read.delim(file.path(outdir, "primary_disease_da_results.tsv"))
sensitivity <- read.delim(file.path(outdir, "sensitivity_bmi_disease_da_results.tsv"))
stopifnot(all(c("CRC_vs_Control", "Adenoma_vs_Control") %in% primary$contrast))
stopifnot(all(primary$model_spec == "primary_age_sex"), all(primary$n_samples == 30))
stopifnot(all(sensitivity$n_samples == 28), all(primary$q_value >= 0 & primary$q_value <= 1))
stopifnot(any(primary$feature == "Constant" & primary$p_value == 1))
stopifnot(!any(primary$feature == "Rare"))
stopifnot(all(grepl("sex_invariant", primary$covariates_omitted)))
baseline <- primary[primary$spike_fraction_target == 0 & primary$contrast == "CRC_vs_Control" &
                      primary$feature == "CRC marker", ]
for (profiler in unique(baseline$profiler)) {
  rows <- baseline[baseline$profiler == profiler, ]
  stopifnot(length(unique(signif(rows$effect, 14))) == 1,
            length(unique(signif(rows$q_value, 14))) == 1)
}
cat("[PASS] disease biomarker-model fixture\n")
