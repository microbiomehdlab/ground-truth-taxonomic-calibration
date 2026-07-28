#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tibble)
})

get_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", ca, value = TRUE)
  if (length(hit) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", hit[[1]]), winslash = "/", mustWork = FALSE)))
  }
  getwd()
}

script_dir <- get_script_dir()

source(file.path(script_dir, "..", "R", "common_utils.R"))
source(file.path(script_dir, "..", "R", "io_utils.R"))
source(file.path(script_dir, "..", "R", "spike_metrics_utils.R"))
source(file.path(script_dir, "..", "R", "spike_plot_utils.R"))

# ------------------------------------------------------------------------------
# Taxon alias handling
# ------------------------------------------------------------------------------

normalize_taxon_key <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[[:space:]]+", " ", x)
  x
}

default_spike_taxon_aliases <- function() {
  tibble::tibble(
    canonical = character(),
    alias = character(),
    tool = character(),
    spike_label = character()
  )
}

read_taxon_aliases <- function(path = NULL) {
  aliases <- default_spike_taxon_aliases()
  
  if (!is.null(path) && nzchar(path)) {
    ext <- tolower(tools::file_ext(path))
    
    if (ext %in% c("tsv", "txt")) {
      user_tbl <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
      user_tbl <- tibble::as_tibble(user_tbl)
    } else {
      user_tbl <- safe_read_csv(path)
    }
    
    required_cols <- c("canonical", "alias")
    missing_cols <- setdiff(required_cols, names(user_tbl))
    if (length(missing_cols) > 0) {
      stop(
        sprintf("Alias file is missing required columns: %s", paste(missing_cols, collapse = ", ")),
        call. = FALSE
      )
    }
    
    if (!"tool" %in% names(user_tbl)) user_tbl$tool <- NA_character_
    if (!"spike_label" %in% names(user_tbl)) user_tbl$spike_label <- NA_character_
    
    aliases <- bind_rows(aliases, user_tbl)
  }
  
  aliases %>%
    mutate(
      canonical = trimws(.data$canonical),
      alias = trimws(.data$alias),
      tool = ifelse(is.na(.data$tool) | .data$tool == "", NA_character_, tolower(trimws(.data$tool))),
      spike_label = ifelse(is.na(.data$spike_label) | .data$spike_label == "", NA_character_, trimws(.data$spike_label))
    ) %>%
    distinct()
}

prepare_alias_table <- function(aliases, design = NULL) {
  out <- aliases %>%
    mutate(
      canonical_key = normalize_taxon_key(.data$canonical),
      alias_key = normalize_taxon_key(.data$alias)
    )
  
  if (nrow(out) > 0) {
    self_maps <- out %>%
      distinct(.data$canonical, .data$tool, .data$spike_label) %>%
      transmute(
        canonical = .data$canonical,
        alias = .data$canonical,
        tool = .data$tool,
        spike_label = .data$spike_label,
        canonical_key = normalize_taxon_key(.data$canonical),
        alias_key = normalize_taxon_key(.data$canonical)
      )
    
    out <- bind_rows(out, self_maps) %>% distinct()
  }
  
  if (!is.null(design) && "member_taxon" %in% names(design)) {
    design_self <- design %>%
      distinct(.data$member_taxon, .data$spike_label) %>%
      transmute(
        canonical = .data$member_taxon,
        alias = .data$member_taxon,
        tool = NA_character_,
        spike_label = .data$spike_label,
        canonical_key = normalize_taxon_key(.data$member_taxon),
        alias_key = normalize_taxon_key(.data$member_taxon)
      )
    out <- bind_rows(out, design_self) %>% distinct()
  }
  
  out
}

build_alias_lookup <- function(aliases) {
  if (is.null(aliases) || nrow(aliases) == 0) {
    return(setNames(character(), character()))
  }
  
  lut_df <- aliases %>%
    mutate(
      tool_key = ifelse(is.na(.data$tool), "*", .data$tool),
      spike_key = ifelse(is.na(.data$spike_label), "*", .data$spike_label),
      specificity =
        ifelse(.data$tool_key != "*", 1L, 0L) +
        ifelse(.data$spike_key != "*", 1L, 0L),
      alias_pref = ifelse(.data$alias_key != .data$canonical_key, 1L, 0L),
      lookup_key = paste(.data$tool_key, .data$spike_key, .data$alias_key, sep = "||")
    ) %>%
    arrange(desc(.data$specificity), desc(.data$alias_pref)) %>%
    distinct(.data$lookup_key, .keep_all = TRUE)
  
  setNames(lut_df$canonical, lut_df$lookup_key)
}

canonicalize_taxa_vector <- function(x, alias_lookup, tool = NULL, spike_label = NULL) {
  if (is.null(x)) return(x)
  
  x_chr <- as.character(x)
  n <- length(x_chr)
  if (n == 0L) return(x_chr)
  if (length(alias_lookup) == 0L) return(x_chr)
  
  tool_vec <- if (is.null(tool)) rep(NA_character_, n) else as.character(tool)
  spike_vec <- if (is.null(spike_label)) rep(NA_character_, n) else as.character(spike_label)
  
  if (length(tool_vec) == 1L && n > 1L) tool_vec <- rep(tool_vec, n)
  if (length(spike_vec) == 1L && n > 1L) spike_vec <- rep(spike_vec, n)
  
  if (length(tool_vec) != n) stop("tool must have length 1 or same length as x", call. = FALSE)
  if (length(spike_vec) != n) stop("spike_label must have length 1 or same length as x", call. = FALSE)
  
  missing_x <- is.na(x_chr) | !nzchar(trimws(x_chr))
  out <- x_chr
  
  x_key <- normalize_taxon_key(ifelse(missing_x, "", x_chr))
  tool_key <- ifelse(is.na(tool_vec) | !nzchar(trimws(tool_vec)), "*", tolower(trimws(tool_vec)))
  spike_key <- ifelse(is.na(spike_vec) | !nzchar(trimws(spike_vec)), "*", trimws(spike_vec))
  
  k1 <- paste(tool_key, spike_key, x_key, sep = "||")
  k2 <- paste(tool_key, "*", x_key, sep = "||")
  k3 <- paste("*", spike_key, x_key, sep = "||")
  k4 <- paste("*", "*", x_key, sep = "||")
  
  hit1 <- unname(alias_lookup[k1])
  hit2 <- unname(alias_lookup[k2])
  hit3 <- unname(alias_lookup[k3])
  hit4 <- unname(alias_lookup[k4])
  
  replace_idx <- !missing_x & !is.na(hit1)
  out[replace_idx] <- hit1[replace_idx]
  
  replace_idx <- !missing_x & is.na(hit1) & !is.na(hit2)
  out[replace_idx] <- hit2[replace_idx]
  
  replace_idx <- !missing_x & is.na(hit1) & is.na(hit2) & !is.na(hit3)
  out[replace_idx] <- hit3[replace_idx]
  
  replace_idx <- !missing_x & is.na(hit1) & is.na(hit2) & is.na(hit3) & !is.na(hit4)
  out[replace_idx] <- hit4[replace_idx]
  
  out
}

canonicalize_taxonomy_columns <- function(df, alias_lookup) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (length(alias_lookup) == 0L) return(df)
  
  taxon_cols <- intersect(
    c(
      "member_taxon",
      "taxon",
      "species",
      "species_name",
      "organism",
      "expected_taxon",
      "observed_taxon",
      "baseline_taxon",
      "target_taxon",
      "detected_taxon"
    ),
    names(df)
  )
  
  if (length(taxon_cols) == 0) return(df)
  
  tool_vec <- if ("tool" %in% names(df)) df$tool else rep(NA_character_, nrow(df))
  spike_vec <- if ("spike_label" %in% names(df)) df$spike_label else rep(NA_character_, nrow(df))
  
  for (col in taxon_cols) {
    raw_col <- paste0(col, "_raw")
    if (!raw_col %in% names(df)) df[[raw_col]] <- df[[col]]
    
    df[[col]] <- canonicalize_taxa_vector(
      x = df[[col]],
      alias_lookup = alias_lookup,
      tool = tool_vec,
      spike_label = spike_vec
    )
  }
  
  df
}

choose_context_alias_rows <- function(aliases, tool_i, spike_i) {
  if (is.null(aliases) || nrow(aliases) == 0) return(aliases[0, , drop = FALSE])
  
  aliases %>%
    mutate(
      tool_match = dplyr::case_when(
        is.na(.data$tool) ~ 0L,
        .data$tool == tool_i ~ 1L,
        TRUE ~ -1L
      ),
      spike_match = dplyr::case_when(
        is.na(.data$spike_label) ~ 0L,
        .data$spike_label == spike_i ~ 1L,
        TRUE ~ -1L
      )
    ) %>%
    filter(.data$tool_match >= 0L, .data$spike_match >= 0L) %>%
    mutate(
      specificity = .data$tool_match + .data$spike_match,
      alias_pref = ifelse(.data$alias_key != .data$canonical_key, 1L, 0L)
    ) %>%
    arrange(desc(.data$specificity), desc(.data$alias_pref)) %>%
    distinct(.data$canonical_key, .keep_all = TRUE)
}

expand_design_for_row <- function(design_df, row_df, aliases) {
  if (is.null(design_df) || nrow(design_df) == 0) return(design_df)
  if (!"member_taxon" %in% names(design_df)) return(design_df)
  if (is.null(aliases) || nrow(aliases) == 0) return(design_df)
  
  tool_i <- tolower(trimws(as.character(row_df$tool[[1]])))
  spike_i <- trimws(as.character(row_df$spike_label[[1]]))
  
  a_sub <- choose_context_alias_rows(aliases, tool_i, spike_i)
  if (nrow(a_sub) == 0) return(design_df)
  
  alias_map <- setNames(a_sub$alias, a_sub$canonical_key)
  
  design_aug <- design_df %>%
    mutate(
      member_taxon_canonical = .data$member_taxon,
      member_taxon_key = normalize_taxon_key(.data$member_taxon),
      member_taxon = dplyr::coalesce(unname(alias_map[.data$member_taxon_key]), .data$member_taxon)
    )
  
  design_aug
}

expand_meta_for_row <- function(meta_df, row_df, aliases) {
  if (is.null(meta_df) || nrow(meta_df) == 0) return(meta_df)
  if (!"member_taxon" %in% names(meta_df)) return(meta_df)
  if (is.null(aliases) || nrow(aliases) == 0) return(meta_df)
  
  tool_i <- tolower(trimws(as.character(row_df$tool[[1]])))
  spike_i <- trimws(as.character(row_df$spike_label[[1]]))
  
  a_sub <- choose_context_alias_rows(aliases, tool_i, spike_i)
  if (nrow(a_sub) == 0) return(meta_df)
  
  alias_map <- setNames(a_sub$alias, a_sub$canonical_key)
  
  meta_aug <- meta_df %>%
    mutate(
      member_taxon_canonical = .data$member_taxon,
      member_taxon_key = normalize_taxon_key(.data$member_taxon),
      member_taxon = dplyr::coalesce(unname(alias_map[.data$member_taxon_key]), .data$member_taxon)
    )
  
  meta_aug
}

compute_spike_metrics_for_manifest_row_aliased <- function(row, design, meta, aliases) {
  design_expanded <- expand_design_for_row(design, row, aliases)
  meta_expanded <- expand_meta_for_row(meta, row, aliases)
  
  if ("member_taxon" %in% names(design_expanded)) {
    target_taxa <- unique(design_expanded$member_taxon[design_expanded$spike_label == row$spike_label[[1]]])
    message(sprintf(
      "[ALIAS] %s | %s -> %s",
      row$spike_label[[1]],
      row$tool[[1]],
      paste(target_taxa, collapse = "; ")
    ))
  }
  
  res <- compute_spike_metrics_for_manifest_row(
    row,
    design = design_expanded,
    meta = meta_expanded
  )
  
  if (is.null(res)) return(NULL)
  res
}

# ------------------------------------------------------------------------------
# CLI options
# ------------------------------------------------------------------------------

option_list <- list(
  make_option("--design", type = "character"),
  make_option("--meta", type = "character", default = NULL),
  make_option("--run_manifest", type = "character", default = NULL),
  make_option("--profile_root", type = "character", default = NULL),
  make_option("--original_tables", type = "character", default = NULL),
  make_option("--auto_manifest_out", type = "character", default = NULL),
  make_option("--analysis_mode", type = "character", default = "all"),
  make_option("--taxon_aliases", type = "character", default = NULL,
              help = "Optional CSV/TSV with columns: canonical, alias, tool, spike_label"),
  make_option("--outdir", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$design) || is.null(opt$outdir)) {
  stop("Required: --design and --outdir", call. = FALSE)
}

ensure_dir(opt$outdir)

# ------------------------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------------------------

design <- safe_read_csv(opt$design)
meta <- if (!is.null(opt$meta) && nzchar(opt$meta)) read_metadata_table(opt$meta) else NULL

run_manifest <- resolve_run_manifest(
  run_manifest = opt$run_manifest,
  profile_root = opt$profile_root,
  original_tables_arg = opt$original_tables,
  design = design
)

mode_req <- tolower(trimws(opt$analysis_mode))
if (!mode_req %in% c("all", "independent", "community")) {
  stop("--analysis_mode must be one of: all, independent, community", call. = FALSE)
}

# ------------------------------------------------------------------------------
# Taxon alias table
# ------------------------------------------------------------------------------

aliases <- read_taxon_aliases(opt$taxon_aliases)
aliases <- prepare_alias_table(aliases, design = design)
alias_lookup <- build_alias_lookup(aliases)

write_csv_safe(aliases, file.path(opt$outdir, "taxon_aliases.resolved.csv"))

# ------------------------------------------------------------------------------
# Filter analysis mode
# ------------------------------------------------------------------------------

if (mode_req != "all") {
  mode_labels <- design %>%
    dplyr::filter(.data$spike_mode == !!mode_req) %>%
    dplyr::pull(.data$spike_label) %>%
    unique()
  
  design <- design %>% dplyr::filter(.data$spike_mode == !!mode_req)
  
  if (!is.null(meta) && "spike_mode" %in% names(meta)) {
    meta <- meta %>% dplyr::filter(.data$spike_mode == !!mode_req)
  }
  
  run_manifest <- run_manifest %>% dplyr::filter(.data$spike_label %in% mode_labels)
}

if (nrow(run_manifest) == 0) {
  stop("No run manifest rows were resolved", call. = FALSE)
}

write_csv_safe(run_manifest, file.path(opt$outdir, "run_manifest.resolved.csv"))
if (!is.null(opt$auto_manifest_out) && nzchar(opt$auto_manifest_out)) {
  write_csv_safe(run_manifest, opt$auto_manifest_out)
}

# ------------------------------------------------------------------------------
# Main loop
# ------------------------------------------------------------------------------

all_baseline <- list()
all_trace <- list()
all_combined <- list()
all_targets <- list()
idx <- 1L

for (i in seq_len(nrow(run_manifest))) {
  message(sprintf(
    "[INFO] Processing %d/%d: %s | %s | %s",
    i,
    nrow(run_manifest),
    run_manifest$spike_label[[i]],
    run_manifest$tool[[i]],
    run_manifest$spike_fraction[[i]]
  ))
  
  row <- run_manifest[i, , drop = FALSE]
  
  res <- compute_spike_metrics_for_manifest_row_aliased(
    row = row,
    design = design,
    meta = meta,
    aliases = aliases
  )
  
  if (is.null(res)) next
  
  all_baseline[[idx]] <- res$baseline
  all_trace[[idx]] <- res$trace
  all_combined[[idx]] <- res$combined
  all_targets[[idx]] <- res$target_member_errors
  idx <- idx + 1L
}

baseline <- dplyr::bind_rows(all_baseline)
trace <- dplyr::bind_rows(all_trace)
combined <- dplyr::bind_rows(all_combined)
targets <- dplyr::bind_rows(all_targets)

baseline <- canonicalize_taxonomy_columns(baseline, alias_lookup)
trace <- canonicalize_taxonomy_columns(trace, alias_lookup)
combined <- canonicalize_taxonomy_columns(combined, alias_lookup)
targets <- canonicalize_taxonomy_columns(targets, alias_lookup)

species_errors <- if (nrow(trace) > 0) summarise_species_errors(trace) else tibble::tibble()
species_errors <- canonicalize_taxonomy_columns(species_errors, alias_lookup)

design_long <- design

design_summary <- design %>%
  dplyr::group_by(.data$spike_mode, .data$spike_label, .data$member_taxon) %>%
  dplyr::summarise(
    n_samples = dplyr::n_distinct(.data$sample_id),
    min_fraction = min(.data$spike_fraction_total, na.rm = TRUE),
    max_fraction = max(.data$spike_fraction_total, na.rm = TRUE),
    .groups = "drop"
  )

lod_summary <- if (nrow(targets) > 0) {
  compute_lod_summary(targets, detection_threshold = 0.95)
} else {
  tibble::tibble()
}
lod_summary <- canonicalize_taxonomy_columns(lod_summary, alias_lookup)

analysis_overview <- tibble::tibble(
  table_name = c(
    "design_rows",
    "manifest_rows",
    "combined_rows",
    "trace_rows",
    "target_rows",
    "baseline_rows",
    "alias_rows"
  ),
  n = c(
    nrow(design_long),
    nrow(run_manifest),
    nrow(combined),
    nrow(trace),
    nrow(targets),
    nrow(baseline),
    nrow(aliases)
  )
)

write_csv_safe(analysis_overview, file.path(opt$outdir, "analysis_overview.csv"))
write_csv_safe(combined, file.path(opt$outdir, "combined_with_condition.csv"))
write_csv_safe(baseline, file.path(opt$outdir, "baseline_with_condition.csv"))
write_csv_safe(species_errors, file.path(opt$outdir, "species_errors_with_condition.csv"))
write_csv_safe(trace, file.path(opt$outdir, "species_trace_with_condition.csv"))
write_csv_safe(targets, file.path(opt$outdir, "target_member_errors_with_condition.csv"))
write_csv_safe(design_long, file.path(opt$outdir, "spike_design_long.csv"))
write_csv_safe(design_summary, file.path(opt$outdir, "spike_design_summary.csv"))
write_csv_safe(lod_summary, file.path(opt$outdir, "lod_summary.csv"))

write_bundle <- function(bundle_dir, design_df, run_manifest_df, combined_df, baseline_df, trace_df,
                         targets_df, species_errors_df, lod_df, aliases_df) {
  ensure_dir(bundle_dir)
  ensure_dir(file.path(bundle_dir, "plots"))
  
  write_csv_safe(run_manifest_df, file.path(bundle_dir, "run_manifest.resolved.csv"))
  write_csv_safe(combined_df, file.path(bundle_dir, "combined_with_condition.csv"))
  write_csv_safe(baseline_df, file.path(bundle_dir, "baseline_with_condition.csv"))
  write_csv_safe(species_errors_df, file.path(bundle_dir, "species_errors_with_condition.csv"))
  write_csv_safe(trace_df, file.path(bundle_dir, "species_trace_with_condition.csv"))
  write_csv_safe(targets_df, file.path(bundle_dir, "target_member_errors_with_condition.csv"))
  write_csv_safe(design_df, file.path(bundle_dir, "spike_design_long.csv"))
  write_csv_safe(lod_df, file.path(bundle_dir, "lod_summary.csv"))
  write_csv_safe(aliases_df, file.path(bundle_dir, "taxon_aliases.resolved.csv"))
  
  generate_all_spike_plots(combined_df, baseline_df, trace_df, targets_df, file.path(bundle_dir, "plots"))
  
  community_labels <- design_df %>%
    dplyr::group_by(.data$spike_label) %>%
    dplyr::summarise(n_members = dplyr::n_distinct(.data$member_taxon), .groups = "drop") %>%
    dplyr::filter(.data$n_members > 1) %>%
    dplyr::pull(.data$spike_label)
  
  if (length(community_labels) > 0) {
    c_trace <- trace_df %>% dplyr::filter(.data$spike_label %in% community_labels)
    c_targ <- targets_df %>% dplyr::filter(.data$spike_label %in% community_labels)
    generate_community_member_plots(c_trace, c_targ, file.path(bundle_dir, "plots", "community_member_resolved"))
  }
  
  if ("Target_Condition" %in% names(combined_df)) {
    conditions <- sort(unique(stats::na.omit(combined_df$Target_Condition)))
    for (cond in conditions) {
      safe_cond <- sanitize_slug(cond)
      out_c <- file.path(bundle_dir, "plots", "by_condition", safe_cond)
      
      c_comb <- combined_df %>% dplyr::filter(.data$Target_Condition == !!cond)
      c_base <- if ("Target_Condition" %in% names(baseline_df)) baseline_df %>% dplyr::filter(.data$Target_Condition == !!cond) else baseline_df[0, , drop = FALSE]
      c_trace <- if ("Target_Condition" %in% names(trace_df)) trace_df %>% dplyr::filter(.data$Target_Condition == !!cond) else trace_df[0, , drop = FALSE]
      c_targ <- if ("Target_Condition" %in% names(targets_df)) targets_df %>% dplyr::filter(.data$Target_Condition == !!cond) else targets_df[0, , drop = FALSE]
      generate_all_spike_plots(c_comb, c_base, c_trace, c_targ, out_c)
    }
  }
  
  labels <- sort(unique(stats::na.omit(design_df$spike_label)))
  for (lab in labels) {
    out_l <- file.path(bundle_dir, "plots", "by_spike_label", sanitize_slug(lab))
    
    d_lab <- design_df %>% dplyr::filter(.data$spike_label == !!lab)
    m_lab <- run_manifest_df %>% dplyr::filter(.data$spike_label == !!lab)
    c_lab <- combined_df %>% dplyr::filter(.data$spike_label == !!lab)
    b_lab <- baseline_df %>% dplyr::filter(.data$spike_label == !!lab)
    t_lab <- trace_df %>% dplyr::filter(.data$spike_label == !!lab)
    g_lab <- targets_df %>% dplyr::filter(.data$spike_label == !!lab)
    s_lab <- species_errors_df %>% dplyr::filter(.data$spike_label == !!lab)
    l_lab <- lod_df %>% dplyr::filter(.data$spike_label == !!lab)
    a_lab <- aliases_df %>% dplyr::filter(is.na(.data$spike_label) | .data$spike_label == !!lab)
    
    ensure_dir(out_l)
    write_csv_safe(m_lab, file.path(out_l, "run_manifest.resolved.csv"))
    write_csv_safe(c_lab, file.path(out_l, "combined_with_condition.csv"))
    write_csv_safe(b_lab, file.path(out_l, "baseline_with_condition.csv"))
    write_csv_safe(s_lab, file.path(out_l, "species_errors_with_condition.csv"))
    write_csv_safe(t_lab, file.path(out_l, "species_trace_with_condition.csv"))
    write_csv_safe(g_lab, file.path(out_l, "target_member_errors_with_condition.csv"))
    write_csv_safe(d_lab, file.path(out_l, "spike_design_long.csv"))
    write_csv_safe(l_lab, file.path(out_l, "lod_summary.csv"))
    write_csv_safe(a_lab, file.path(out_l, "taxon_aliases.resolved.csv"))
    
    generate_all_spike_plots(c_lab, b_lab, t_lab, g_lab, out_l)
    
    if ("member_taxon" %in% names(d_lab) && dplyr::n_distinct(d_lab$member_taxon) > 1) {
      generate_community_member_plots(t_lab, g_lab, file.path(out_l, "member_resolved"))
    }
  }
}

if (mode_req == "all") {
  modes_present <- sort(unique(stats::na.omit(design$spike_mode)))
  
  for (mode in modes_present) {
    labels_mode <- design %>%
      dplyr::filter(.data$spike_mode == !!mode) %>%
      dplyr::pull(.data$spike_label) %>%
      unique()
    
    design_m <- design %>% dplyr::filter(.data$spike_mode == !!mode)
    run_manifest_m <- run_manifest %>% dplyr::filter(.data$spike_label %in% labels_mode)
    combined_m <- combined %>% dplyr::filter(.data$spike_mode == !!mode)
    baseline_m <- baseline %>% dplyr::filter(.data$spike_mode == !!mode)
    trace_m <- trace %>% dplyr::filter(.data$spike_mode == !!mode)
    targets_m <- targets %>% dplyr::filter(.data$spike_mode == !!mode)
    species_errors_m <- species_errors %>% dplyr::filter(.data$spike_mode == !!mode)
    lod_m <- lod_summary %>% dplyr::filter(.data$spike_label %in% labels_mode)
    aliases_m <- aliases %>% dplyr::filter(is.na(.data$spike_label) | .data$spike_label %in% labels_mode)
    
    write_bundle(
      file.path(opt$outdir, mode),
      design_m,
      run_manifest_m,
      combined_m,
      baseline_m,
      trace_m,
      targets_m,
      species_errors_m,
      lod_m,
      aliases_m
    )
  }
} else {
  labels_mode <- design %>% dplyr::pull(.data$spike_label) %>% unique()
  lod_m <- lod_summary %>% dplyr::filter(.data$spike_label %in% labels_mode)
  aliases_m <- aliases %>% dplyr::filter(is.na(.data$spike_label) | .data$spike_label %in% labels_mode)
  
  write_bundle(
    opt$outdir,
    design,
    run_manifest,
    combined,
    baseline,
    trace,
    targets,
    species_errors,
    lod_m,
    aliases_m
  )
}

cat(sprintf("[OK] Wrote outputs under %s\n", opt$outdir))
