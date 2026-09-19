required <- c(
  "Maaslin2", "optparse", "data.table", "dplyr", "tibble", "readr",
  "tidyr", "stringr", "purrr", "rlang", "ggplot2", "forcats", "scales",
  "patchwork", "cowplot", "gridExtra", "pheatmap", "ggrepel",
  "RColorBrewer", "mgcv", "sandwich", "arrow", "jsonlite",
  "sessioninfo", "png"
)
status <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
if (!all(status)) stop("Missing packages: ", paste(required[!status], collapse = ", "))

expected <- c(
  Maaslin2 = "1.18.0",
  sandwich = "3.1.1",
  arrow = "17.0.0"
)
observed <- vapply(names(expected), function(package) {
  as.character(utils::packageVersion(package))
}, character(1))
wrong <- names(expected)[observed != expected]
if (length(wrong)) {
  details <- paste0(wrong, " expected=", expected[wrong], " observed=", observed[wrong])
  stop("Unexpected package versions: ", paste(details, collapse = "; "))
}

r_version <- paste(R.version$major, R.version$minor, sep = ".")
if (r_version != "4.3.3") stop("Expected R 4.3.3; found ", r_version)

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

cat("R environment verification: PASSED\n")
cat("R version: ", r_version, "\n", sep = "")
for (package in names(expected)) {
  cat(package, " version: ", observed[[package]], "\n", sep = "")
}
