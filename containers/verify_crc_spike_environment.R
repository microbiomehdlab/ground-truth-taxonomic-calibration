required <- c(
  "Maaslin2", "optparse", "data.table", "dplyr", "tibble", "readr",
  "tidyr", "stringr", "purrr", "rlang", "ggplot2", "forcats", "scales",
  "patchwork", "cowplot", "gridExtra", "pheatmap", "ggrepel",
  "RColorBrewer", "mgcv", "jsonlite", "sessioninfo", "png"
)
status <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
if (!all(status)) stop("Missing packages: ", paste(required[!status], collapse = ", "))

r_version <- paste(R.version$major, R.version$minor, sep = ".")
maaslin_version <- as.character(utils::packageVersion("Maaslin2"))
if (r_version != "4.3.3") stop("Expected R 4.3.3; found ", r_version)
if (maaslin_version != "1.18.0") {
  stop("Expected MaAsLin2 1.18.0; found ", maaslin_version)
}

required_args <- c(
  "input_data", "input_metadata", "output", "fixed_effects",
  "normalization", "transform", "analysis_method", "min_abundance",
  "min_prevalence", "min_variance", "max_significance", "correction",
  "standardize", "cores"
)
missing_args <- setdiff(required_args, names(formals(Maaslin2::Maaslin2)))
if (length(missing_args)) {
  stop("MaAsLin2 function is missing arguments: ", paste(missing_args, collapse = ", "))
}
cat("Environment verification: PASSED\n")
cat("R version: ", r_version, "\nMaAsLin2 version: ", maaslin_version, "\n", sep = "")
