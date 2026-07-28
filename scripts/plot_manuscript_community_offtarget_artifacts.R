#!/usr/bin/env Rscript
# Polished v3: endpoint axes expanded for 0%/100% points; intended for vertical manuscript assembly.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(grid)
})

`%||%` <- function(x, y) {
  if (!is.null(x) && length(x) > 0 && !all(is.na(x)) && !(is.character(x) && !nzchar(x[1]))) x else y
}

option_list <- list(
  make_option("--indir", type = "character", default = "RUNS/maaslin_spike",
              help = "Input directory containing MaAsLin CSVs and run_design_manifest_ALLFILTERS.csv [default %default]"),
  make_option("--manifest-file", type = "character", default = NULL,
              help = "Optional direct path to run_design_manifest_ALLFILTERS.csv"),
  make_option("--member-file", type = "character", default = NULL,
              help = "Optional direct path to maaslin_member_detection_ALLFILTERS.csv"),
  make_option("--significant-file", type = "character", default = NULL,
              help = "Optional direct path to maaslin_significant_features_ALLFILTERS.csv"),
  make_option("--alias-file", type = "character", default = NULL,
              help = "Optional direct path to taxon_aliases.resolved.csv; used to exclude target aliases such as Fusobacterium nucleatum/F. nucleatum subsp. nucleatum"),
  make_option("--outdir", type = "character", default = "RUNS/plots_manuscript_community_offtarget_artifacts",
              help = "Output directory [default %default]"),
  make_option("--community-size", type = "integer", default = 10,
              help = "Number of taxa in the community spike. Effective per-species fraction = total community fraction / this value [default %default]"),
  make_option("--main_effective_fractions", type = "character", default = "0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01",
              help = "Effective per-species fractions used in Panels A/B/C [default %default = 0.001%%,0.005%%,0.01%%,0.05%%,0.10%%,0.50%%,1.00%%]"),
  make_option("--focus_effective_fractions", type = "character", default = "0.0001,0.01",
              help = "Comma-separated effective per-species fractions used in Panels E/F [default %default = 0.01%%,1.00%%]"),
  make_option("--filter-mode", type = "character", default = "original",
              help = "Filter mode to use if filter_mode column exists [default %default]"),
  make_option("--artifact-min-shift", type = "double", default = 0,
              help = "Minimum positive spiked-minus-original abundance shift counted as positive [default %default]"),
  make_option("--artifact-min-positive-rate", type = "double", default = 0.50,
              help = "Minimum fraction of samples with positive shift for a non-target taxon to count as abundance artefact [default %default]"),
  make_option("--error_threshold", type = "double", default = 0.05,
              help = "Threshold for high-error taxa in Panels C/D, as absolute relative error versus expected dilution [default %default = 5%%]"),
  make_option("--top-labels", type = "integer", default = 4,
              help = "Number of top taxa to label per profiler in Panels B-D [default %default]"),
  make_option("--max-rank", type = "integer", default = 200,
              help = "Maximum top-N taxa to display in ranked-overlap Panel D [default %default]"),
  make_option("--main-width", type = "double", default = 16.5,
              help = "Main figure width in inches [default %default]"),
  make_option("--main-height", type = "double", default = 20.8,
              help = "Main figure height in inches [default %default]"),
  make_option("--dpi", type = "integer", default = 320,
              help = "PNG resolution [default %default]")
)

raw_args <- commandArgs(trailingOnly = TRUE)
raw_args <- sub("^--main-effective-fractions$", "--main_effective_fractions", raw_args)
raw_args <- sub("^--focus-effective-fractions$", "--focus_effective_fractions", raw_args)
raw_args <- sub("^--focus-effective-fraction$", "--focus_effective_fractions", raw_args)
raw_args <- sub("^--artifact-min-shift$", "--artifact-min-shift", raw_args)
raw_args <- sub("^--artifact-min-positive-rate$", "--artifact-min-positive-rate", raw_args)
raw_args <- sub("^--error-threshold$", "--error_threshold", raw_args)
opt <- parse_args(OptionParser(option_list = option_list), args = raw_args)

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

manifest_file <- opt$`manifest-file` %||% file.path(opt$indir, "run_design_manifest_ALLFILTERS.csv")
member_file <- opt$`member-file` %||% file.path(opt$indir, "maaslin_member_detection_ALLFILTERS.csv")
significant_file <- opt$`significant-file` %||% file.path(opt$indir, "maaslin_significant_features_ALLFILTERS.csv")
alias_file <- opt$`alias-file` %||% file.path(opt$indir, "taxon_aliases.resolved.csv")

if (!file.exists(manifest_file)) stop("Could not find manifest file: ", manifest_file, call. = FALSE)
if (!file.exists(member_file)) stop("Could not find member detection file: ", member_file, call. = FALSE)
if (!file.exists(significant_file)) stop("Could not find significant-features file: ", significant_file, call. = FALSE)

parse_num_vector <- function(x) {
  x %>% strsplit(",") %>% unlist() %>% trimws() %>% .[nzchar(.)] %>% as.numeric()
}
fmt_fraction <- function(x) {
  p <- as.numeric(x) * 100
  ifelse(p < 0.01, sprintf("%.3f%%", p), sprintf("%.2f%%", p))
}
norm_key <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}
binomial_key <- function(x) {
  x <- norm_key(x)
  vapply(strsplit(x, " "), function(parts) {
    parts <- parts[nzchar(parts)]
    if (length(parts) >= 2) paste(parts[1:2], collapse = " ") else paste(parts, collapse = " ")
  }, character(1))
}
is_target_or_alias_key <- function(keys, target_keys, target_binomial_keys) {
  keys <- norm_key(keys)
  exact <- keys %in% target_keys
  binom <- binomial_key(keys)
  binom_match <- binom %in% target_binomial_keys
  exact | binom_match
}
short_taxon <- function(x, max_len = 34) {
  x <- as.character(x)
  ifelse(nchar(x) > max_len, paste0(substr(x, 1, max_len - 1), "\u2026"), x)
}
to_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  lx <- tolower(trimws(as.character(x)))
  lx %in% c("true", "t", "1", "yes", "y")
}

main_effective_fractions <- parse_num_vector(opt$main_effective_fractions)
focus_effective_fractions <- sort(unique(parse_num_vector(opt$focus_effective_fractions)))
if (length(focus_effective_fractions) == 0) stop("No values supplied to --focus-effective-fractions", call. = FALSE)
focus_labels <- fmt_fraction(focus_effective_fractions)
focus_label <- paste(focus_labels, collapse = "_")
error_threshold <- as.numeric(opt$error_threshold)
error_threshold_label <- fmt_fraction(error_threshold)

tool_labeller <- c(
  kraken2_bracken = "Kraken2 + Bracken",
  metaphlan4 = "MetaPhlAn 4"
)
tool_levels <- unname(tool_labeller)
profiler_cols <- c(
  "Kraken2 + Bracken" = "#2A9D8F",
  "MetaPhlAn 4" = "#7B61D1"
)

fraction_cols <- c(
  "0.001%" = "#244D8F",
  "0.005%" = "#3B6EA8",
  "0.01%"  = "#2A9D8F",
  "0.05%"  = "#63B179",
  "0.10%"  = "#D99B2B",
  "0.50%"  = "#B65B84",
  "1.00%"  = "#666666"
)

pub_theme <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid.major.y = element_line(colour = "#ECECEC", linewidth = 0.28),
      panel.grid.major.x = element_line(colour = "#F3F3F3", linewidth = 0.22),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "#4D4D4D", linewidth = 0.34),
      axis.ticks = element_line(colour = "#4D4D4D", linewidth = 0.30),
      axis.ticks.length = unit(2.0, "pt"),
      strip.background = element_rect(fill = "#F8F8F8", colour = "#D9D9D9", linewidth = 0.42),
      strip.text = element_text(face = "bold", colour = "#222222", margin = margin(3, 3, 3, 3)),
      axis.title = element_text(face = "bold", colour = "#222222"),
      axis.text = element_text(colour = "#444444"),
      legend.title = element_text(face = "bold", colour = "#222222", size = rel(0.92)),
      legend.text = element_text(colour = "#333333", size = rel(0.86)),
      legend.key = element_rect(fill = "white", colour = NA),
      plot.title = element_text(face = "bold", size = base_size + 2.4, hjust = 0, colour = "#111111"),
      plot.subtitle = element_blank(),
      plot.margin = margin(6, 6, 6, 6)
    )
}

save_plot_set <- function(plot_obj, stem, width, height, dpi = 320) {
  pdf_path <- file.path(opt$outdir, paste0(stem, ".pdf"))
  png_path <- file.path(opt$outdir, paste0(stem, ".png"))
  if (capabilities("cairo")) {
    ggsave(pdf_path, plot_obj, width = width, height = height, device = cairo_pdf, bg = "white")
  } else {
    ggsave(pdf_path, plot_obj, width = width, height = height, bg = "white")
  }
  ggsave(png_path, plot_obj, width = width, height = height, dpi = dpi, bg = "white")
}

# ----------------------------
# Read inputs
# ----------------------------
message("[INFO] Reading manifest: ", manifest_file)
manifest <- read_csv(manifest_file, show_col_types = FALSE)

required_manifest <- c("spike_label", "tool", "spike_fraction", "spiked_table", "original_table", "spike_mode")
missing_manifest <- setdiff(required_manifest, names(manifest))
if (length(missing_manifest) > 0) {
  stop("Manifest file is missing required columns: ", paste(missing_manifest, collapse = ", "), call. = FALSE)
}

message("[INFO] Reading member detection: ", member_file)
member <- read_csv(member_file, show_col_types = FALSE)
required_member <- c("tool", "spike_mode", "spike_label", "member_taxon")
missing_member <- setdiff(required_member, names(member))
if (length(missing_member) > 0) {
  stop("Member detection file is missing required columns: ", paste(missing_member, collapse = ", "), call. = FALSE)
}
if ("filter_mode" %in% names(member)) {
  member <- member %>% filter(filter_mode == opt$`filter-mode`)
}

target_taxa <- member %>%
  filter(spike_mode == "community", spike_label == "CRCpanel") %>%
  distinct(member_taxon) %>%
  mutate(target_key = norm_key(member_taxon))

if (!nrow(target_taxa)) {
  target_taxa <- member %>%
    filter(spike_mode == "independent", spike_label != "CRCpanel") %>%
    distinct(member_taxon) %>%
    mutate(target_key = norm_key(member_taxon))
}

target_keys_initial <- unique(target_taxa$target_key)

# Expand target keys with the alias table. This is critical because the spiked
# target may be named "Fusobacterium nucleatum subsp. nucleatum" in the spike
# manifest but "Fusobacterium nucleatum" in the profiler table or MaAsLin output.
target_spike_labels <- member %>%
  filter(spike_mode == "independent", member_taxon %in% target_taxa$member_taxon) %>%
  distinct(spike_label) %>%
  pull(spike_label)

target_keys <- target_keys_initial

if (file.exists(alias_file)) {
  message("[INFO] Reading aliases: ", alias_file)
  alias_tbl <- read_csv(alias_file, show_col_types = FALSE) %>%
    mutate(
      canonical_key2 = if ("canonical_key" %in% names(.)) canonical_key else norm_key(canonical),
      alias_key2 = if ("alias_key" %in% names(.)) alias_key else norm_key(alias),
      spike_label = if ("spike_label" %in% names(.)) as.character(spike_label) else NA_character_
    )
  
  alias_edges <- alias_tbl %>%
    select(canonical_key2, alias_key2, spike_label) %>%
    distinct()
  
  # Start with aliases explicitly associated with the spiked labels, if present,
  # and any canonical/alias entries overlapping the target member names.
  alias_target_keys <- alias_edges %>%
    filter(
      canonical_key2 %in% target_keys |
        alias_key2 %in% target_keys |
        (!is.na(spike_label) & spike_label %in% target_spike_labels)
    ) %>%
    select(canonical_key2, alias_key2) %>%
    unlist(use.names = FALSE) %>%
    unique()
  
  target_keys <- unique(c(target_keys, alias_target_keys))
  
  # A few alias tables contain chains. Expand until stable.
  repeat {
    new_keys <- alias_edges %>%
      filter(canonical_key2 %in% target_keys | alias_key2 %in% target_keys) %>%
      select(canonical_key2, alias_key2) %>%
      unlist(use.names = FALSE) %>%
      unique()
    expanded <- unique(c(target_keys, new_keys))
    if (length(expanded) == length(target_keys)) break
    target_keys <- expanded
  }
} else {
  warning("Alias file not found: ", alias_file,
          ". Target exclusion will use exact normalized target names only.")
}

target_keys <- unique(target_keys[!is.na(target_keys) & nzchar(target_keys)])

# Also exclude genus-species binomial keys for targets. This catches cases such
# as "Fusobacterium nucleatum" when the spike target is recorded as
# "Fusobacterium nucleatum subsp. nucleatum".  We use the binomial, not only the
# genus, so unrelated species from the same genus remain eligible off-targets.
target_binomial_keys <- unique(binomial_key(target_keys))
target_binomial_keys <- target_binomial_keys[!is.na(target_binomial_keys) & nzchar(target_binomial_keys)]
target_keys <- unique(c(target_keys, target_binomial_keys))

message("[INFO] Target taxa in community spike: ", length(target_keys_initial),
        "; target/alias/binomial keys excluded from off-target analyses: ", length(target_keys),
        " (binomials: ", paste(target_binomial_keys, collapse = "; "), ")")
write_csv(
  tibble(
    target_key = sort(unique(target_keys)),
    binomial_key = binomial_key(sort(unique(target_keys))),
    is_binomial_key = sort(unique(target_keys)) %in% target_binomial_keys
  ),
  file.path(opt$outdir, "target_and_alias_keys_excluded_from_offtarget_analysis.csv")
)

# ----------------------------
# Abundance artefact summaries from profile tables
# ----------------------------
detect_sample_col <- function(df) {
  candidates <- c("sample_id", "sample", "SampleID", "Sample", "id", "ID", "X", "...1")
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) return(hit[1])
  nonnum <- names(df)[vapply(df, function(z) !is.numeric(z), logical(1))]
  if (length(nonnum)) return(nonnum[1])
  names(df)[1]
}

clean_sample_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("^X(?=[0-9])", "", x, perl = TRUE)
  
  # Spiked profile rows often append the spike design/fraction to the original
  # sample id, e.g. ERR478961_CRCpanel_f0p0001.  Strip that suffix so the
  # spiked subset can be matched back to the larger original table.
  x <- gsub("_[A-Za-z0-9]+_f[0-9]+p[0-9]+$", "", x, perl = TRUE)
  x <- gsub("_[A-Za-z0-9]+_f[0-9]+$", "", x, perl = TRUE)
  x <- gsub("_f[0-9]+p[0-9]+$", "", x, perl = TRUE)
  x <- gsub("_f[0-9]+$", "", x, perl = TRUE)
  
  x <- gsub("\\.spiked$|_spiked$|\\.unspiked$|_unspiked$", "", x, ignore.case = TRUE)
  x <- gsub("\\s+", "", x)
  x
}

resolve_path <- function(path_value, manifest_dir = dirname(manifest_file)) {
  p <- as.character(path_value)
  if (file.exists(p)) return(p)
  p2 <- file.path(manifest_dir, basename(p))
  if (file.exists(p2)) return(p2)
  p3 <- file.path(getwd(), basename(p))
  if (file.exists(p3)) return(p3)
  return(p)
}

profile_cache <- new.env(parent = emptyenv())

make_row_oriented_candidate <- function(df) {
  sample_col <- detect_sample_col(df)
  df[[sample_col]] <- as.character(df[[sample_col]])
  taxa_cols <- setdiff(names(df), sample_col)
  df[taxa_cols] <- lapply(df[taxa_cols], function(z) suppressWarnings(as.numeric(z)))
  names(df)[names(df) == sample_col] <- ".sample_id"
  list(
    df = df,
    sample_col = ".sample_id",
    taxa_cols = taxa_cols,
    orientation = "samples_as_rows",
    target_overlap = sum(norm_key(taxa_cols) %in% target_keys, na.rm = TRUE)
  )
}

make_transposed_candidate <- function(df) {
  feature_col <- detect_sample_col(df)
  feature_names <- as.character(df[[feature_col]])
  sample_cols <- setdiff(names(df), feature_col)
  
  mat <- df %>%
    select(all_of(sample_cols)) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
    as.matrix()
  
  # Rows are features/taxa and columns are samples; transpose to samples x taxa.
  mat_t <- t(mat)
  taxa_names <- make.unique(ifelse(is.na(feature_names) | !nzchar(feature_names),
                                   paste0("feature_", seq_along(feature_names)),
                                   feature_names))
  
  out <- as.data.frame(mat_t, check.names = FALSE)
  colnames(out) <- taxa_names
  out$.sample_id <- rownames(out)
  out <- out %>% select(.sample_id, everything())
  
  list(
    df = out,
    sample_col = ".sample_id",
    taxa_cols = taxa_names,
    orientation = "samples_as_columns",
    target_overlap = sum(norm_key(taxa_names) %in% target_keys, na.rm = TRUE)
  )
}

read_profile <- function(path_value) {
  path <- resolve_path(path_value)
  if (!file.exists(path)) {
    stop("Could not find profile table referenced in manifest: ", path_value, call. = FALSE)
  }
  if (exists(path, envir = profile_cache, inherits = FALSE)) {
    return(get(path, envir = profile_cache, inherits = FALSE))
  }
  
  df_raw <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  
  row_candidate <- make_row_oriented_candidate(df_raw)
  transposed_candidate <- make_transposed_candidate(df_raw)
  
  # Pick the orientation with more target taxa in the taxa dimension.
  # This makes the script robust to files stored as samples x taxa or taxa x samples.
  if (transposed_candidate$target_overlap > row_candidate$target_overlap) {
    obj <- transposed_candidate
  } else {
    obj <- row_candidate
  }
  
  obj$df$.sample_id <- as.character(obj$df$.sample_id)
  obj$df$.sample_key <- clean_sample_id(obj$df$.sample_id)
  obj$taxa_cols <- setdiff(obj$taxa_cols, ".sample_key")
  
  message("[INFO]    loaded ", basename(path), " as ", obj$orientation,
          " (", nrow(obj$df), " samples, ", length(obj$taxa_cols),
          " taxa, target-overlap=", obj$target_overlap, ")")
  
  assign(path, obj, envir = profile_cache)
  obj
}

summarise_abundance_shift <- function(spiked_path, original_path, tool_pretty, spike_fraction_total, effective_fraction) {
  sp <- read_profile(spiked_path)
  og <- read_profile(original_path)
  
  sp_df <- sp$df
  og_df <- og$df
  
  if (!(".sample_key" %in% names(sp_df))) sp_df$.sample_key <- clean_sample_id(sp_df$.sample_id)
  if (!(".sample_key" %in% names(og_df))) og_df$.sample_key <- clean_sample_id(og_df$.sample_id)
  
  common_keys <- intersect(sp_df$.sample_key, og_df$.sample_key)
  
  if (!length(common_keys)) {
    if (nrow(sp_df) == nrow(og_df)) {
      warning("No shared sample IDs after normalization; pairing spiked/original tables by row order for: ",
              basename(resolve_path(spiked_path)), " / ", basename(resolve_path(original_path)))
      sp_df$.pair_key <- seq_len(nrow(sp_df))
      og_df$.pair_key <- seq_len(nrow(og_df))
      common_keys <- sp_df$.pair_key
      key_col <- ".pair_key"
    } else {
      stop(
        "No shared samples between spiked and original tables after normalization: ",
        spiked_path, " / ", original_path,
        "\nFirst spiked IDs: ", paste(head(sp_df$.sample_id, 5), collapse = ", "),
        "\nFirst original IDs: ", paste(head(og_df$.sample_id, 5), collapse = ", "),
        call. = FALSE
      )
    }
  } else {
    key_col <- ".sample_key"
  }
  
  message("[INFO]    matched ", length(common_keys), " spiked samples to original samples; ",
          "unmatched original samples are ignored for this comparison.")
  
  common_taxa <- intersect(sp$taxa_cols, og$taxa_cols)
  common_taxa <- setdiff(common_taxa, c(".sample_id", ".sample_key", ".pair_key"))
  if (!length(common_taxa)) {
    stop("No shared taxon columns between spiked and original tables: ", spiked_path, " / ", original_path, call. = FALSE)
  }
  
  taxon_keys <- norm_key(common_taxa)
  nontarget_idx <- !is_target_or_alias_key(taxon_keys, target_keys, target_binomial_keys)
  taxa_use <- common_taxa[nontarget_idx]
  keys_use <- taxon_keys[nontarget_idx]
  
  if (!length(taxa_use)) {
    stop("No non-target taxa available after excluding community target taxa.", call. = FALSE)
  }
  
  sp_mat <- sp_df %>%
    filter(.data[[key_col]] %in% common_keys) %>%
    arrange(match(.data[[key_col]], common_keys)) %>%
    select(all_of(taxa_use)) %>%
    as.matrix()
  
  og_mat <- og_df %>%
    filter(.data[[key_col]] %in% common_keys) %>%
    arrange(match(.data[[key_col]], common_keys)) %>%
    select(all_of(taxa_use)) %>%
    as.matrix()
  
  # Absolute abundance shift relative to the original table.
  delta <- sp_mat - og_mat
  delta[!is.finite(delta)] <- NA_real_
  
  # Relative error versus the expected dilution of non-target taxa.
  # For a community spike with total spike fraction f, a non-target taxon is
  # expected to be diluted to original_abundance * (1 - f).  Relative error is:
  #   (observed_spiked - expected_diluted) / expected_diluted
  # Rows with expected abundance <= 0 cannot define a relative error and are NA.
  expected_diluted <- og_mat * (1 - spike_fraction_total)
  rel_error <- (sp_mat - expected_diluted) / expected_diluted
  rel_error[!is.finite(rel_error) | expected_diluted <= 0] <- NA_real_
  
  positive_matrix <- delta > opt$`artifact-min-shift`
  positive_shift <- pmax(delta, 0)
  
  tibble(
    tool = tool_pretty,
    spike_fraction_total = spike_fraction_total,
    effective_fraction = effective_fraction,
    effective_fraction_label = fmt_fraction(effective_fraction),
    taxon = taxa_use,
    taxon_key = keys_use,
    n_samples = nrow(delta),
    median_shift = apply(delta, 2, median, na.rm = TRUE),
    mean_shift = colMeans(delta, na.rm = TRUE),
    mean_abs_shift = colMeans(abs(delta), na.rm = TRUE),
    median_abs_shift = apply(abs(delta), 2, median, na.rm = TRUE),
    sd_shift = apply(delta, 2, sd, na.rm = TRUE),
    iqr_shift = apply(delta, 2, IQR, na.rm = TRUE),
    mean_abs_relative_error = colMeans(abs(rel_error), na.rm = TRUE),
    median_abs_relative_error = apply(abs(rel_error), 2, median, na.rm = TRUE),
    sd_relative_error = apply(rel_error, 2, sd, na.rm = TRUE),
    iqr_relative_error = apply(rel_error, 2, IQR, na.rm = TRUE),
    positive_shift_rate = colMeans(positive_matrix, na.rm = TRUE),
    median_positive_shift = apply(positive_shift, 2, median, na.rm = TRUE),
    max_positive_shift = apply(positive_shift, 2, max, na.rm = TRUE)
  ) %>%
    mutate(
      artifact_score = positive_shift_rate * median_positive_shift,
      positive_abundance_artifact = positive_shift_rate >= opt$`artifact-min-positive-rate` &
        median_positive_shift > opt$`artifact-min-shift`
    )
}

community_manifest <- manifest %>%
  filter(spike_mode == "community", spike_label == "CRCpanel") %>%
  mutate(
    tool = recode(as.character(tool), !!!tool_labeller, .default = as.character(tool)),
    spike_fraction = as.numeric(spike_fraction),
    effective_fraction = spike_fraction / opt$`community-size`
  )

if (!nrow(community_manifest)) {
  stop("No community CRCpanel rows found in manifest.", call. = FALSE)
}

message("[INFO] Computing non-target abundance shifts from ", nrow(community_manifest), " community profile tables.")
abundance_list <- vector("list", nrow(community_manifest))
for (i in seq_len(nrow(community_manifest))) {
  row <- community_manifest[i, ]
  message("[INFO]  - ", row$tool, " total=", fmt_fraction(row$spike_fraction),
          " effective=", fmt_fraction(row$effective_fraction))
  abundance_list[[i]] <- summarise_abundance_shift(
    spiked_path = row$spiked_table,
    original_path = row$original_table,
    tool_pretty = row$tool,
    spike_fraction_total = row$spike_fraction,
    effective_fraction = row$effective_fraction
  )
}
abundance_artifacts <- bind_rows(abundance_list) %>%
  mutate(
    tool = factor(tool, levels = tool_levels),
    effective_fraction_label = factor(effective_fraction_label, levels = fmt_fraction(sort(unique(effective_fraction))))
  )

write_csv(abundance_artifacts, file.path(opt$outdir, "nontarget_abundance_artifacts_all_effective_fractions.csv"))

focus_abundance <- abundance_artifacts %>%
  filter(effective_fraction %in% focus_effective_fractions) %>%
  mutate(tool = factor(tool, levels = tool_levels))

if (!nrow(focus_abundance)) {
  stop("No abundance-artifact rows matched --focus-effective-fractions: ", paste(focus_effective_fractions, collapse = ","), call. = FALSE)
}

# ----------------------------
# DA off-target calls from MaAsLin significant features
# ----------------------------
message("[INFO] Reading significant features: ", significant_file)
sig_cols <- c(
  "feature_norm", "member_taxon", "is_target", "is_significant", "is_positive",
  "qval", "coef", "background_condition", "background_study", "tool",
  "spike_mode", "spike_label", "spike_fraction", "filter_mode"
)
sig_raw <- read_csv(
  significant_file,
  show_col_types = FALSE,
  col_select = any_of(sig_cols),
  progress = FALSE
)

if ("filter_mode" %in% names(sig_raw)) {
  sig_raw <- sig_raw %>% filter(filter_mode == opt$`filter-mode`)
}

required_sig <- c("tool", "spike_mode", "spike_label", "spike_fraction", "is_target", "is_significant", "is_positive", "qval", "coef")
missing_sig <- setdiff(required_sig, names(sig_raw))
if (length(missing_sig) > 0) {
  stop("Significant-features file is missing required columns: ", paste(missing_sig, collapse = ", "), call. = FALSE)
}

sig_tbl <- sig_raw %>%
  mutate(
    tool = recode(as.character(tool), !!!tool_labeller, .default = as.character(tool)),
    tool = factor(tool, levels = tool_levels),
    spike_fraction = as.numeric(spike_fraction),
    effective_fraction = if_else(spike_mode == "community", spike_fraction / opt$`community-size`, spike_fraction),
    effective_fraction_label = fmt_fraction(effective_fraction),
    is_target = to_logical_safe(is_target),
    is_significant = to_logical_safe(is_significant),
    is_positive = to_logical_safe(is_positive),
    qval = as.numeric(qval),
    coef = as.numeric(coef),
    feature_key = norm_key(coalesce(member_taxon, feature_norm)),
    # Force target/alias/binomial-matched taxa to be targets even if the upstream
    # significant-feature table did not flag the alias as target.
    is_target = is_target | is_target_or_alias_key(feature_key, target_keys, target_binomial_keys)
  ) %>%
  filter(spike_mode == "community", spike_label == "CRCpanel")

da_offtarget_summary <- sig_tbl %>%
  filter(!is_target) %>%
  group_by(tool, effective_fraction, effective_fraction_label, feature_key) %>%
  summarise(
    feature_label = first(coalesce(member_taxon, feature_norm)),
    any_offtarget_enriched_DA = any(is_significant & is_positive, na.rm = TRUE),
    min_q_offtarget_enriched = suppressWarnings(min(qval[is_significant & is_positive], na.rm = TRUE)),
    max_coef_offtarget_enriched = suppressWarnings(max(coef[is_significant & is_positive], na.rm = TRUE)),
    n_offtarget_enriched_contexts = sum(is_significant & is_positive, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    min_q_offtarget_enriched = if_else(is.infinite(min_q_offtarget_enriched), NA_real_, min_q_offtarget_enriched),
    max_coef_offtarget_enriched = if_else(is.infinite(max_coef_offtarget_enriched), NA_real_, max_coef_offtarget_enriched),
    neglog10_q_offtarget = if_else(any_offtarget_enriched_DA, -log10(pmax(min_q_offtarget_enriched, 1e-300)), 0)
  )

write_csv(da_offtarget_summary, file.path(opt$outdir, "DA_offtarget_enriched_summary_all_effective_fractions.csv"))

focus_da <- da_offtarget_summary %>%
  filter(effective_fraction %in% focus_effective_fractions) %>%
  select(
    tool, effective_fraction, feature_key, feature_label, any_offtarget_enriched_DA,
    min_q_offtarget_enriched, max_coef_offtarget_enriched,
    n_offtarget_enriched_contexts, neglog10_q_offtarget
  )

# ----------------------------
# Link abundance artefacts to off-target DA calls at the focus fraction
# ----------------------------
link_tbl <- focus_abundance %>%
  left_join(focus_da, by = c("tool", "effective_fraction", "taxon_key" = "feature_key")) %>%
  mutate(
    any_offtarget_enriched_DA = if_else(is.na(any_offtarget_enriched_DA), FALSE, any_offtarget_enriched_DA),
    neglog10_q_offtarget = if_else(is.na(neglog10_q_offtarget), 0, neglog10_q_offtarget),
    n_offtarget_enriched_contexts = if_else(is.na(n_offtarget_enriched_contexts), 0L, as.integer(n_offtarget_enriched_contexts)),
    display_label = coalesce(feature_label, taxon)
  ) %>%
  group_by(tool) %>%
  arrange(artifact_score, .by_group = TRUE) %>%
  mutate(
    abundance_artifact_percentile = 100 * percent_rank(artifact_score),
    abundance_artifact_rank_desc = min_rank(desc(artifact_score))
  ) %>%
  ungroup()

write_csv(link_tbl, file.path(opt$outdir, paste0("abundance_artifact_DA_overlap_focus_", gsub("%", "pct", focus_label), ".csv")))

# Sanity check: no known target/alias key should remain in the off-target table.
target_leak_check <- link_tbl %>%
  filter(
    is_target_or_alias_key(taxon_key, target_keys, target_binomial_keys) |
      is_target_or_alias_key(norm_key(display_label), target_keys, target_binomial_keys)
  ) %>%
  distinct(tool, taxon, taxon_key, display_label)
write_csv(target_leak_check, file.path(opt$outdir, "SANITY_CHECK_target_aliases_remaining_in_offtarget_table.csv"))
if (nrow(target_leak_check) > 0) {
  warning("Some target aliases still remained in the off-target table. See SANITY_CHECK_target_aliases_remaining_in_offtarget_table.csv")
}

# ----------------------------
# Link abundance artefacts to off-target DA calls across all requested fractions
# ----------------------------
link_tbl_all <- abundance_artifacts %>%
  filter(effective_fraction %in% main_effective_fractions) %>%
  left_join(
    da_offtarget_summary %>%
      select(
        tool, effective_fraction, feature_key, feature_label,
        any_offtarget_enriched_DA, min_q_offtarget_enriched,
        max_coef_offtarget_enriched, n_offtarget_enriched_contexts,
        neglog10_q_offtarget
      ),
    by = c("tool", "effective_fraction", "taxon_key" = "feature_key")
  ) %>%
  mutate(
    any_offtarget_enriched_DA = if_else(is.na(any_offtarget_enriched_DA), FALSE, any_offtarget_enriched_DA),
    neglog10_q_offtarget = if_else(is.na(neglog10_q_offtarget), 0, neglog10_q_offtarget),
    n_offtarget_enriched_contexts = if_else(is.na(n_offtarget_enriched_contexts), 0L, as.integer(n_offtarget_enriched_contexts)),
    display_label = coalesce(feature_label, taxon),
    effective_fraction_label = factor(effective_fraction_label, levels = fmt_fraction(main_effective_fractions)),
    artifact_score_for_outlier = artifact_score
  ) %>%
  group_by(tool, effective_fraction, effective_fraction_label) %>%
  mutate(
    q1_artifact = quantile(artifact_score_for_outlier, 0.25, na.rm = TRUE, names = FALSE),
    q3_artifact = quantile(artifact_score_for_outlier, 0.75, na.rm = TRUE, names = FALSE),
    iqr_artifact = IQR(artifact_score_for_outlier, na.rm = TRUE),
    upper_iqr_bound = q3_artifact + 1.5 * iqr_artifact,
    is_artifact_outlier = artifact_score_for_outlier > upper_iqr_bound & artifact_score_for_outlier > 0,
    abundance_artifact_percentile = 100 * percent_rank(artifact_score_for_outlier),
    abundance_artifact_rank_desc = min_rank(desc(artifact_score_for_outlier))
  ) %>%
  ungroup()

write_csv(link_tbl_all, file.path(opt$outdir, "abundance_artifact_DA_overlap_all_requested_fractions.csv"))

link_tbl <- link_tbl_all %>%
  filter(effective_fraction %in% focus_effective_fractions)

write_csv(link_tbl, file.path(opt$outdir, paste0("abundance_artifact_DA_overlap_focus_", gsub("%", "pct", focus_label), ".csv")))

target_leak_check <- link_tbl_all %>%
  filter(
    is_target_or_alias_key(taxon_key, target_keys, target_binomial_keys) |
      is_target_or_alias_key(norm_key(display_label), target_keys, target_binomial_keys)
  ) %>%
  distinct(tool, effective_fraction_label, taxon, taxon_key, display_label)
write_csv(target_leak_check, file.path(opt$outdir, "SANITY_CHECK_target_aliases_remaining_in_offtarget_table.csv"))
if (nrow(target_leak_check) > 0) {
  warning("Some target aliases still remained in the off-target table. See SANITY_CHECK_target_aliases_remaining_in_offtarget_table.csv")
}


# ----------------------------
# ----------------------------
# Main figure panels: combined artefact burden, enrichment, and error-variance landscapes
# ----------------------------
error_rate_summary <- link_tbl_all %>%
  group_by(tool, effective_fraction, effective_fraction_label) %>%
  summarise(
    n_non_target = n(),
    n_high_rel_error = sum(mean_abs_relative_error > error_threshold, na.rm = TRUE),
    high_rel_error_fraction_all = mean(mean_abs_relative_error > error_threshold, na.rm = TRUE),
    n_offtarget_DA = sum(any_offtarget_enriched_DA, na.rm = TRUE),
    n_DA_high_rel_error = sum(any_offtarget_enriched_DA & mean_abs_relative_error > error_threshold, na.rm = TRUE),
    high_rel_error_fraction_DA = if_else(n_offtarget_DA > 0, n_DA_high_rel_error / n_offtarget_DA, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    effective_fraction_label = factor(effective_fraction_label, levels = fmt_fraction(main_effective_fractions)),
    high_rel_error_fraction_all_pct = 100 * high_rel_error_fraction_all,
    high_rel_error_fraction_DA_pct = 100 * high_rel_error_fraction_DA
  )
write_csv(error_rate_summary, file.path(opt$outdir, "panel_A_high_relative_error_summary.csv"))

sd_rate_summary <- link_tbl_all %>%
  group_by(tool, effective_fraction, effective_fraction_label) %>%
  summarise(
    n_non_target = n(),
    n_high_rel_sd = sum(sd_relative_error > error_threshold, na.rm = TRUE),
    high_rel_sd_fraction_all = mean(sd_relative_error > error_threshold, na.rm = TRUE),
    n_offtarget_DA = sum(any_offtarget_enriched_DA, na.rm = TRUE),
    n_DA_high_rel_sd = sum(any_offtarget_enriched_DA & sd_relative_error > error_threshold, na.rm = TRUE),
    high_rel_sd_fraction_DA = if_else(n_offtarget_DA > 0, n_DA_high_rel_sd / n_offtarget_DA, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    effective_fraction_label = factor(effective_fraction_label, levels = fmt_fraction(main_effective_fractions)),
    high_rel_sd_fraction_all_pct = 100 * high_rel_sd_fraction_all,
    high_rel_sd_fraction_DA_pct = 100 * high_rel_sd_fraction_DA
  )
write_csv(sd_rate_summary, file.path(opt$outdir, "panel_B_high_relative_sd_summary.csv"))

burden_tbl <- bind_rows(
  error_rate_summary %>%
    transmute(
      tool, effective_fraction, effective_fraction_label,
      metric = paste0("Mean absolute relative error > ", error_threshold_label),
      burden_fraction = high_rel_error_fraction_all
    ),
  sd_rate_summary %>%
    transmute(
      tool, effective_fraction, effective_fraction_label,
      metric = paste0("SD of relative error > ", error_threshold_label),
      burden_fraction = high_rel_sd_fraction_all
    )
) %>%
  mutate(metric = factor(metric, levels = c(
    paste0("Mean absolute relative error > ", error_threshold_label),
    paste0("SD of relative error > ", error_threshold_label)
  )))
write_csv(burden_tbl, file.path(opt$outdir, "panel_A_combined_artefact_burden.csv"))

pA <- ggplot(
  burden_tbl,
  aes(x = effective_fraction_label, y = burden_fraction, colour = tool, group = tool)
) +
  geom_line(linewidth = 0.95, alpha = 0.96) +
  geom_point(size = 2.8, alpha = 0.98) +
  facet_wrap(~ metric, ncol = 1) +
  scale_color_manual(values = profiler_cols, name = "Profiler") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, NA),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "Profiler-induced artefact burden across spike fractions",
    x = "Effective per-species fraction",
    y = "Non-target taxa exceeding threshold (%)"
  ) +
  pub_theme(10.7) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12.3, face = "bold", margin = margin(b = 5)),
    axis.title = element_text(size = 10.4, face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8.1),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#E8E8E8", linewidth = 0.28),
    strip.text = element_text(size = 9.1, face = "bold"),
    plot.margin = margin(8, 10, 8, 8)
  )

enrichment_tbl <- bind_rows(
  error_rate_summary %>%
    transmute(
      tool, effective_fraction, effective_fraction_label,
      metric = paste0("Mean absolute relative error > ", error_threshold_label),
      all_non_target = high_rel_error_fraction_all,
      DA_offtargets = high_rel_error_fraction_DA
    ),
  sd_rate_summary %>%
    transmute(
      tool, effective_fraction, effective_fraction_label,
      metric = paste0("SD of relative error > ", error_threshold_label),
      all_non_target = high_rel_sd_fraction_all,
      DA_offtargets = high_rel_sd_fraction_DA
    )
) %>%
  mutate(metric = factor(metric, levels = levels(burden_tbl$metric)))
write_csv(enrichment_tbl, file.path(opt$outdir, "panel_B_combined_enrichment_summary.csv"))

enrichment_long <- enrichment_tbl %>%
  pivot_longer(
    cols = c(all_non_target, DA_offtargets),
    names_to = "group",
    values_to = "fraction_value"
  ) %>%
  mutate(
    group = recode(group,
                   all_non_target = "All non-target taxa",
                   DA_offtargets = "Off-target DA taxa"
    ),
    group = factor(group, levels = c("All non-target taxa", "Off-target DA taxa"))
  )

pB <- ggplot() +
  geom_segment(
    data = enrichment_tbl,
    aes(
      y = effective_fraction_label, yend = effective_fraction_label,
      x = all_non_target, xend = DA_offtargets
    ),
    colour = "#CFCFCF",
    linewidth = 0.62,
    lineend = "round"
  ) +
  geom_point(
    data = enrichment_long %>% filter(group == "All non-target taxa"),
    aes(x = fraction_value, y = effective_fraction_label, colour = group),
    size = 2.6, alpha = 0.98, na.rm = TRUE
  ) +
  geom_point(
    data = enrichment_long %>% filter(group == "Off-target DA taxa"),
    aes(x = fraction_value, y = effective_fraction_label, colour = group),
    size = 2.95, alpha = 0.99, na.rm = TRUE
  ) +
  facet_grid(metric ~ tool) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    breaks = c(0, 0.5, 1),
    limits = c(-0.035, 1.035),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_colour_manual(
    values = c("All non-target taxa" = "#9E9E9E", "Off-target DA taxa" = "#B2182B"),
    name = NULL
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3.0, alpha = 1))) +
  labs(
    title = "DA off-targets are enriched among artefactual taxa",
    x = "Fraction of taxa exceeding threshold",
    y = "Effective per-species fraction"
  ) +
  pub_theme(10.7) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.title = element_text(size = 12.3, face = "bold", margin = margin(b = 5)),
    axis.title = element_text(size = 10.4, face = "bold"),
    axis.text = element_text(size = 8.0),
    panel.grid.major = element_line(colour = "#E9E9E9", linewidth = 0.25),
    strip.text = element_text(size = 8.9, face = "bold"),
    plot.margin = margin(8, 10, 8, 8)
  )

panel_scatter_tbl <- link_tbl %>%
  mutate(
    effective_fraction_label = factor(effective_fraction_label, levels = focus_labels),
    mean_abs_relative_error_pct = 100 * mean_abs_relative_error,
    sd_relative_error_pct = 100 * sd_relative_error,
    da_status = if_else(any_offtarget_enriched_DA, "Off-target enriched DA call", "No enriched DA call")
  ) %>%
  group_by(tool, effective_fraction_label) %>%
  arrange(desc(any_offtarget_enriched_DA), desc(neglog10_q_offtarget), desc(mean_abs_relative_error + sd_relative_error), .by_group = TRUE) %>%
  mutate(
    label_rank = row_number(),
    label_me = any_offtarget_enriched_DA & label_rank <= opt$`top-labels`,
    label_text = if_else(label_me, short_taxon(display_label, 26), NA_character_)
  ) %>%
  ungroup()
write_csv(panel_scatter_tbl, file.path(opt$outdir, paste0("panel_C_D_focus_", gsub("%", "pct", focus_label), "_relative_error_variance_scatter.csv")))

pC <- ggplot() +
  geom_point(
    data = panel_scatter_tbl %>% filter(!any_offtarget_enriched_DA),
    aes(x = mean_abs_relative_error_pct, y = sd_relative_error_pct),
    colour = "#C9C9C9", fill = "#C9C9C9", alpha = 0.22, size = 0.80, stroke = 0
  ) +
  geom_point(
    data = panel_scatter_tbl %>% filter(any_offtarget_enriched_DA),
    aes(x = mean_abs_relative_error_pct, y = sd_relative_error_pct),
    colour = "#B2182B", fill = "#B2182B", alpha = 0.94, size = 1.85, stroke = 0
  ) +
  geom_vline(xintercept = 100 * error_threshold, linetype = 2, linewidth = 0.45, colour = "#7A7A7A") +
  geom_hline(yintercept = 100 * error_threshold, linetype = 2, linewidth = 0.45, colour = "#7A7A7A") +
  geom_text(
    data = panel_scatter_tbl %>% filter(label_me),
    aes(x = mean_abs_relative_error_pct, y = sd_relative_error_pct, label = label_text),
    colour = "#B2182B", size = 2.00, hjust = -0.02, vjust = -0.18,
    check_overlap = TRUE, show.legend = FALSE
  ) +
  facet_grid(rows = vars(effective_fraction_label), cols = vars(tool)) +
  scale_x_continuous(labels = label_number(accuracy = 1)) +
  scale_y_continuous(labels = label_number(accuracy = 1)) +
  labs(
    title = "Error–variability landscape at representative spike fractions",
    x = "Mean absolute relative error (%)",
    y = "SD of relative error (%)"
  ) +
  pub_theme(10.3) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12.1, face = "bold", margin = margin(b = 5)),
    axis.title = element_text(size = 10.2, face = "bold"),
    axis.text = element_text(size = 7.8),
    panel.grid.major = element_line(colour = "#EAEAEA", linewidth = 0.24),
    strip.text = element_text(size = 8.6, face = "bold"),
    plot.margin = margin(8, 10, 8, 8)
  )

pD <- ggplot() +
  geom_point(
    data = panel_scatter_tbl %>% filter(!any_offtarget_enriched_DA),
    aes(x = mean_abs_relative_error_pct, y = sd_relative_error_pct),
    colour = "#C9C9C9", fill = "#C9C9C9", alpha = 0.20, size = 0.80, stroke = 0
  ) +
  geom_point(
    data = panel_scatter_tbl %>% filter(any_offtarget_enriched_DA),
    aes(x = mean_abs_relative_error_pct, y = sd_relative_error_pct),
    colour = "#B2182B", fill = "#B2182B", alpha = 0.95, size = 1.85, stroke = 0
  ) +
  geom_vline(xintercept = 100 * error_threshold, linetype = 2, linewidth = 0.45, colour = "#7A7A7A") +
  geom_hline(yintercept = 100 * error_threshold, linetype = 2, linewidth = 0.45, colour = "#7A7A7A") +
  geom_text(
    data = panel_scatter_tbl %>% filter(label_me, mean_abs_relative_error_pct <= 20, sd_relative_error_pct <= 20),
    aes(x = mean_abs_relative_error_pct, y = sd_relative_error_pct, label = label_text),
    colour = "#B2182B", size = 2.00, hjust = -0.02, vjust = -0.18,
    check_overlap = TRUE, show.legend = FALSE
  ) +
  facet_grid(rows = vars(effective_fraction_label), cols = vars(tool)) +
  coord_cartesian(xlim = c(0, 20), ylim = c(0, 20), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20), labels = label_number(accuracy = 1)) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20), labels = label_number(accuracy = 1)) +
  labs(
    title = "Zoomed error–variability cloud (0–20%)",
    x = "Mean absolute relative error (%)",
    y = "SD of relative error (%)"
  ) +
  pub_theme(10.3) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12.1, face = "bold", margin = margin(b = 5)),
    axis.title = element_text(size = 10.2, face = "bold"),
    axis.text = element_text(size = 7.8),
    panel.grid.major = element_line(colour = "#EAEAEA", linewidth = 0.24),
    strip.text = element_text(size = 8.6, face = "bold"),
    plot.margin = margin(8, 10, 8, 8)
  )

# Supplementary: ranked-capture curves retained for quantitative detail.
rank_capture_all <- link_tbl_all %>%
  group_by(tool, effective_fraction, effective_fraction_label) %>%
  arrange(desc(artifact_score_for_outlier), .by_group = TRUE) %>%
  mutate(
    rank = row_number(),
    hit = any_offtarget_enriched_DA,
    total_hits = sum(hit, na.rm = TRUE),
    n_taxa = n(),
    cumulative_hits = cumsum(hit),
    cumulative_fraction_DA = if_else(total_hits > 0, cumulative_hits / total_hits, NA_real_),
    top_fraction_taxa = rank / n_taxa
  ) %>%
  ungroup()
write_csv(rank_capture_all, file.path(opt$outdir, "supp_ranked_capture_all_fractions.csv"))

selected_capture_labels <- fmt_fraction(c(0.0001, 0.001, 0.01))
selected_capture_labels <- selected_capture_labels[selected_capture_labels %in% levels(link_tbl_all$effective_fraction_label)]
capture_plot_tbl <- rank_capture_all %>%
  filter(effective_fraction_label %in% selected_capture_labels) %>%
  group_by(tool, effective_fraction_label) %>%
  filter(rank <= min(opt$`max-rank`, max(rank, na.rm = TRUE))) %>%
  ungroup()

pS_rank <- ggplot(capture_plot_tbl, aes(x = top_fraction_taxa, y = cumulative_fraction_DA, colour = effective_fraction_label)) +
  geom_abline(slope = 1, intercept = 0, colour = "#8A8A8A", linetype = 2, linewidth = 0.45) +
  geom_line(linewidth = 0.80, alpha = 0.94, na.rm = TRUE) +
  facet_wrap(~ tool, nrow = 1) +
  scale_colour_manual(values = fraction_cols[names(fraction_cols) %in% selected_capture_labels], name = "Effective\nfraction") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0.02, 0.03))) +
  labs(
    title = "Supplementary. High artefact-score taxa capture DA off-targets",
    x = "Top fraction of non-target taxa ranked by abundance artefact score",
    y = "Cumulative fraction of\noff-target DA calls"
  ) +
  pub_theme(10.2) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(colour = "#EFEFEF", linewidth = 0.25),
    strip.text = element_text(size = 9.2, face = "bold"),
    plot.margin = margin(4, 8, 4, 4)
  )

# Supplementary top tables/plots.
top_artifacts <- abundance_artifacts %>%
  group_by(tool, effective_fraction, effective_fraction_label) %>%
  arrange(desc(artifact_score), .by_group = TRUE) %>%
  slice_head(n = 50) %>%
  ungroup()
write_csv(top_artifacts, file.path(opt$outdir, "supp_top50_nontarget_abundance_artifacts_by_fraction.csv"))

top_da_overlap <- link_tbl %>%
  filter(any_offtarget_enriched_DA) %>%
  arrange(tool, desc(neglog10_q_offtarget), desc(artifact_score))
write_csv(top_da_overlap, file.path(opt$outdir, paste0("supp_DA_offtargets_with_abundance_artifact_metrics_focus_", gsub("%", "pct", focus_label), ".csv")))


# ----------------------------
# ----------------------------
# Main figure
# ----------------------------
main_title <- "Community spike-ins reveal that off-target DA calls are linked to measurable profiler-induced abundance artefacts"

main_fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(0.95, 1.28), widths = c(1.0, 1.0)) +
  plot_annotation(
    title = main_title,
    tag_levels = "A",
    tag_suffix = ".",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16.0, hjust = 0, colour = "#111111", margin = margin(b = 8)),
      plot.tag = element_text(face = "bold", size = 13.5, colour = "#111111")
    )
  ) &
  theme(
    plot.tag.position = c(0.01, 0.995),
    plot.margin = margin(10, 12, 10, 12)
  )

save_plot_set(pA, "panel_A_combined_artefact_burden", width = 8.0, height = 6.2, dpi = opt$dpi)
save_plot_set(pB, "panel_B_combined_DA_enrichment_among_artefactual_taxa", width = 9.3, height = 6.4, dpi = opt$dpi)
save_plot_set(pC, "panel_C_error_variability_landscape_focus_fractions", width = 9.0, height = 8.0, dpi = opt$dpi)
save_plot_set(pD, "panel_D_error_variability_landscape_focus_fractions_zoom20", width = 9.0, height = 8.0, dpi = opt$dpi)
save_plot_set(pS_rank, "supp_ranked_capture_DA_offtargets_by_artifact_score", width = 8.4, height = 5.0, dpi = opt$dpi)
save_plot_set(main_fig, "manuscript_community_offtarget_artifacts_overview", width = opt$`main-width`, height = opt$`main-height`, dpi = opt$dpi)

readme <- c(
  "Run example:",
  "Rscript scripts/05_plot_manuscript_community_offtarget_artifacts.R \\",
  "  --indir RUNS/maaslin_spike \\",
  "  --outdir RUNS/plots_manuscript_community_offtarget_artifacts \\",
  "  --community-size 10 \\",
  "  --main-effective-fractions 0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01 \\",
  "  --focus-effective-fraction 0.01 \\",
  "  --error-threshold 0.05",
  "",
  "Required inputs:",
  "- run_design_manifest_ALLFILTERS.csv, with spiked_table and original_table paths",
  "- maaslin_member_detection_ALLFILTERS.csv",
  "- maaslin_significant_features_ALLFILTERS.csv",
  "- taxon_aliases.resolved.csv is strongly recommended so target aliases are excluded from off-target analyses",
  "",
  "Main outputs:",
  "- manuscript_community_offtarget_artifacts_overview.(pdf|png)",
  "- panel_A_high_relative_error_rate.(pdf|png)",
  "- panel_B_DA_offtargets_enriched_among_high_relative_error_taxa.(pdf|png)",
  "- panel_C_high_relative_sd_rate.(pdf|png)",
  "- panel_D_DA_offtargets_enriched_among_high_relative_sd_taxa.(pdf|png)",
  "- panel_E_focus_relative_error_variance_scatter.(pdf|png)",
  "- panel_F_focus_relative_error_variance_scatter_zoom20.(pdf|png)",
  "- supp_ranked_capture_DA_offtargets_by_artifact_score.(pdf|png)",
  "",
  "Key tables:",
  "- nontarget_abundance_artifacts_all_effective_fractions.csv",
  "- DA_offtarget_enriched_summary_all_effective_fractions.csv",
  "- panel_A_B_high_relative_error_summary.csv",
  "- panel_C_D_high_relative_sd_summary.csv",
  "- panel_E_F_focus_<fraction1>_<fraction2>_relative_error_variance_scatter.csv",
  "- abundance_artifact_DA_overlap_focus_<fraction>.csv",
  "- supp_ranked_capture_all_fractions.csv",
  "",
  "Definitions:",
  "- Off-target DA call: non-target taxon significant and enriched in MaAsLin.",
  "- High relative error taxon: mean absolute relative error versus expected dilution > --error-threshold.",
  "- High relative SD taxon: SD of relative error versus expected dilution > --error-threshold.",
  "- Panels A/B summarise the burden of high-relative-error taxa and their enrichment among DA off-targets.",
  "- Panels C/D repeat the same logic for high relative SD.",
  "- Panels E/F show the joint relationship between relative error and relative SD across the selected focus effective fractions, including a 0–20% zoom.",
  "- Community total spike fractions are divided by --community-size to obtain effective per-species fractions."
)
writeLines(readme, con = file.path(opt$outdir, "README_community_offtarget_artifacts.txt"))

message("[OK] Wrote community off-target artefact figure set under: ", opt$outdir)
