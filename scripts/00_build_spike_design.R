#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))

get_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", ca, value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[[1]]), winslash = "/", mustWork = FALSE)))
  getwd()
}
script_dir <- get_script_dir()
source(file.path(script_dir, "..", "R", "common_utils.R"))
source(file.path(script_dir, "..", "R", "spike_design_utils.R"))

option_list <- list(
  make_option("--independent_manifest", type = "character", default = NULL),
  make_option("--community_manifest", type = "character", default = NULL),
  make_option("--spike_panel", type = "character"),
  make_option("--community_memberships", type = "character", default = NULL),
  make_option("--outdir", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$spike_panel) || is.null(opt$outdir)) {
  stop("Required: --spike_panel and --outdir", call. = FALSE)
}
if (is.null(opt$independent_manifest) && is.null(opt$community_manifest)) {
  stop("Provide at least one of --independent_manifest or --community_manifest", call. = FALSE)
}

ensure_dir(opt$outdir)
res <- build_spike_design(
  independent_manifest = opt$independent_manifest,
  community_manifest = opt$community_manifest,
  spike_panel = opt$spike_panel,
  community_memberships = opt$community_memberships
)

write_csv_safe(res$design, file.path(opt$outdir, "spike_design.tsv"))
write_csv_safe(res$meta, file.path(opt$outdir, "metadata_spiked.tsv"))
write_csv_safe(res$community_memberships, file.path(opt$outdir, "community_memberships.used.tsv"))

cat(sprintf("[OK] Wrote %s\n", file.path(opt$outdir, "spike_design.tsv")))
cat(sprintf("[OK] Wrote %s\n", file.path(opt$outdir, "metadata_spiked.tsv")))
cat(sprintf("[OK] Wrote %s\n", file.path(opt$outdir, "community_memberships.used.tsv")))
