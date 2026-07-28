#!/usr/bin/env Rscript

# Validate the complete profiler-specific target mapping before analysis.
# Uses base R only so it can run early during preflight.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  pos <- match(flag, args)
  if (is.na(pos)) return(default)
  if (pos == length(args)) stop("Missing value after ", flag, call. = FALSE)
  args[[pos + 1L]]
}

panel_file <- get_arg("--spike-panel", "spike_panel.tsv")
alias_file <- get_arg("--aliases", "spike_taxon_aliases.csv")
kraken_file <- get_arg("--kraken", "kraken2_bracken_merged_unspecified.csv")
metaphlan_file <- get_arg("--metaphlan", "metaphlan4_merged_unspecified.csv")

for (path in c(panel_file, alias_file, kraken_file, metaphlan_file)) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing or empty required file: ", path, call. = FALSE)
  }
}

panel <- utils::read.delim(
  panel_file, stringsAsFactors = FALSE, check.names = FALSE,
  quote = "", comment.char = ""
)
aliases <- utils::read.csv(
  alias_file, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("NA")
)

panel_required <- c("label", "taxon_name")
alias_required <- c("canonical", "alias", "tool", "spike_label")
missing_panel <- setdiff(panel_required, names(panel))
missing_alias <- setdiff(alias_required, names(aliases))
if (length(missing_panel)) {
  stop("spike_panel.tsv is missing columns: ", paste(missing_panel, collapse = ", "), call. = FALSE)
}
if (length(missing_alias)) {
  stop("spike_taxon_aliases.csv is missing columns: ", paste(missing_alias, collapse = ", "), call. = FALSE)
}

trim_or_na <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}
panel$label <- trim_or_na(panel$label)
panel$taxon_name <- trim_or_na(panel$taxon_name)
aliases$canonical <- trim_or_na(aliases$canonical)
aliases$alias <- trim_or_na(aliases$alias)
aliases$tool <- tolower(trim_or_na(aliases$tool))
aliases$spike_label <- trim_or_na(aliases$spike_label)

if (anyNA(panel$label) || anyNA(panel$taxon_name)) {
  stop("spike_panel.tsv contains blank label or taxon_name values.", call. = FALSE)
}
if (anyDuplicated(panel$label)) stop("spike_panel.tsv contains duplicate labels.", call. = FALSE)
if (anyDuplicated(panel$taxon_name)) stop("spike_panel.tsv contains duplicate taxon_name values.", call. = FALSE)
if (nrow(panel) != 10L) {
  stop("Expected exactly 10 spike-panel targets; found ", nrow(panel), ".", call. = FALSE)
}

if (anyNA(aliases$canonical) || anyNA(aliases$alias) || anyNA(aliases$tool)) {
  stop("Alias canonical, alias, and tool values must be non-empty.", call. = FALSE)
}
if (any(!is.na(aliases$spike_label))) {
  stop(
    "All spike_label values must be blank so mappings apply to independent and community spike-ins.",
    call. = FALSE
  )
}

expected_tools <- c("kraken2_bracken", "metaphlan4")
unexpected_tools <- setdiff(unique(aliases$tool), expected_tools)
if (length(unexpected_tools)) {
  stop("Unexpected alias tool value(s): ", paste(unexpected_tools, collapse = ", "), call. = FALSE)
}

panel_taxa <- sort(panel$taxon_name)
alias_taxa <- sort(unique(aliases$canonical))
missing_taxa <- setdiff(panel_taxa, alias_taxa)
extra_taxa <- setdiff(alias_taxa, panel_taxa)
if (length(missing_taxa)) {
  stop("Canonical targets missing from alias file: ", paste(missing_taxa, collapse = "; "), call. = FALSE)
}
if (length(extra_taxa)) {
  stop("Alias canonicals absent from spike panel: ", paste(extra_taxa, collapse = "; "), call. = FALSE)
}

pair_key <- paste(aliases$canonical, aliases$tool, sep = "||")
duplicate_pairs <- unique(pair_key[duplicated(pair_key) | duplicated(pair_key, fromLast = TRUE)])
if (length(duplicate_pairs)) {
  stop(
    "Each canonical/tool pair must occur exactly once; duplicates: ",
    paste(duplicate_pairs, collapse = "; "),
    call. = FALSE
  )
}

expected_pairs <- as.vector(outer(panel$taxon_name, expected_tools, paste, sep = "||"))
missing_pairs <- setdiff(expected_pairs, pair_key)
extra_pairs <- setdiff(pair_key, expected_pairs)
if (length(missing_pairs) || length(extra_pairs) || nrow(aliases) != 20L) {
  msg <- c(
    if (length(missing_pairs)) paste0("missing: ", paste(missing_pairs, collapse = "; ")),
    if (length(extra_pairs)) paste0("unexpected: ", paste(extra_pairs, collapse = "; ")),
    if (nrow(aliases) != 20L) paste0("row count: expected 20, found ", nrow(aliases))
  )
  stop("Incomplete canonical/tool mapping (", paste(msg, collapse = " | "), ").", call. = FALSE)
}

ambiguous_key <- paste(aliases$tool, aliases$alias, sep = "||")
ambiguous <- split(aliases$canonical, ambiguous_key)
ambiguous <- ambiguous[vapply(ambiguous, function(x) length(unique(x)) > 1L, logical(1))]
if (length(ambiguous)) {
  stop(
    "A profiler alias maps to multiple canonical targets: ",
    paste(names(ambiguous), collapse = "; "),
    call. = FALSE
  )
}

headers <- list(
  kraken2_bracken = names(utils::read.csv(kraken_file, nrows = 0, check.names = FALSE)),
  metaphlan4 = names(utils::read.csv(metaphlan_file, nrows = 0, check.names = FALSE))
)
missing_reported <- character()
for (i in seq_len(nrow(aliases))) {
  tool <- aliases$tool[[i]]
  alias <- aliases$alias[[i]]
  if (!(alias %in% headers[[tool]])) {
    missing_reported <- c(
      missing_reported,
      paste0(tool, ": ", aliases$canonical[[i]], " -> ", alias)
    )
  }
}
if (length(missing_reported)) {
  stop(
    "Declared aliases absent from the corresponding unspiked profiler table: ",
    paste(missing_reported, collapse = "; "),
    call. = FALSE
  )
}

aliases <- aliases[order(aliases$tool, match(aliases$canonical, panel$taxon_name)), ]
cat("[PASS] Complete taxon-alias mapping validated.\n")
cat("       Canonical targets: ", length(panel_taxa), "\n", sep = "")
cat("       Profiler mappings: ", nrow(aliases), " (10 per profiler)\n", sep = "")
cat("       All declared aliases occur in their profiler abundance-table headers.\n")
