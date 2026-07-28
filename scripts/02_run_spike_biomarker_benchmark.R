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
source(file.path(script_dir, "..", "R", "io_utils.R"))
source(file.path(script_dir, "..", "R", "maaslin_utils.R"))

option_list <- list(
  make_option("--design", type = "character"),
  make_option("--meta", type = "character"),
  make_option("--metrics_dir", type = "character"),
  make_option("--run_manifest", type = "character", default = NULL),
  make_option("--profile_root", type = "character", default = NULL),
  make_option("--original_tables", type = "character", default = NULL),
  make_option("--taxon_aliases", type = "character", default = NULL,
              help = "Optional CSV/TSV with columns: canonical, alias, tool, spike_label. If omitted, metrics_dir/taxon_aliases.resolved.csv is used when present."),
  make_option("--filter_modes", type = "character", default = "original,std,errorvar"),
  make_option("--background_conditions", type = "character", default = "ALL"),
  make_option("--background_condition", type = "character", default = NULL),
  make_option("--background_studies", type = "character", default = "ALL"),
  make_option("--background_study", type = "character", default = NULL),
  make_option("--fixed_effect", type = "character", default = "spike_status"),
  make_option("--reference_level", type = "character", default = "unspiked"),
  make_option("--q_threshold", type = "double", default = 0.25),
  make_option("--min_samples_per_group", type = "integer", default = 5L),
  make_option("--outdir", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

req <- c("design", "meta", "metrics_dir", "outdir")
for (nm in req) if (is.null(opt[[nm]])) stop(sprintf("Required: --%s", nm), call. = FALSE)

ensure_dir(opt$outdir)

spike_alias_normalize_key <- function(x) {
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
  
  if (!is.null(path) && nzchar(path) && file.exists(path)) {
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
    
    aliases <- dplyr::bind_rows(aliases, user_tbl)
  }
  
  aliases %>%
    dplyr::mutate(
      canonical = trimws(as.character(.data$canonical)),
      alias = trimws(as.character(.data$alias)),
      tool = dplyr::if_else(is.na(.data$tool) | .data$tool == "", NA_character_, tolower(trimws(as.character(.data$tool)))),
      spike_label = dplyr::if_else(is.na(.data$spike_label) | .data$spike_label == "", NA_character_, trimws(as.character(.data$spike_label)))
    ) %>%
    dplyr::distinct()
}

prepare_alias_table <- function(aliases, design = NULL) {
  out <- aliases %>%
    dplyr::mutate(
      canonical_key = spike_alias_normalize_key(.data$canonical),
      alias_key = spike_alias_normalize_key(.data$alias)
    )
  
  if (nrow(out) > 0) {
    self_maps <- out %>%
      dplyr::distinct(.data$canonical, .data$tool, .data$spike_label) %>%
      dplyr::transmute(
        canonical = .data$canonical,
        alias = .data$canonical,
        tool = .data$tool,
        spike_label = .data$spike_label,
        canonical_key = spike_alias_normalize_key(.data$canonical),
        alias_key = spike_alias_normalize_key(.data$canonical)
      )
    out <- dplyr::bind_rows(out, self_maps) %>% dplyr::distinct()
  }
  
  if (!is.null(design) && "member_taxon" %in% names(design)) {
    design_self <- design %>%
      dplyr::distinct(.data$member_taxon, .data$spike_label) %>%
      dplyr::transmute(
        canonical = .data$member_taxon,
        alias = .data$member_taxon,
        tool = NA_character_,
        spike_label = .data$spike_label,
        canonical_key = spike_alias_normalize_key(.data$member_taxon),
        alias_key = spike_alias_normalize_key(.data$member_taxon)
      )
    out <- dplyr::bind_rows(out, design_self) %>% dplyr::distinct()
  }
  
  out
}

build_alias_lookup <- function(aliases) {
  if (is.null(aliases) || nrow(aliases) == 0) {
    return(stats::setNames(character(), character()))
  }
  
  lut_df <- aliases %>%
    dplyr::mutate(
      tool_key = dplyr::if_else(is.na(.data$tool), "*", .data$tool),
      spike_key = dplyr::if_else(is.na(.data$spike_label), "*", .data$spike_label),
      specificity =
        dplyr::if_else(.data$tool_key != "*", 1L, 0L) +
        dplyr::if_else(.data$spike_key != "*", 1L, 0L),
      alias_pref = dplyr::if_else(.data$alias_key != .data$canonical_key, 1L, 0L),
      lookup_key = paste(.data$tool_key, .data$spike_key, .data$alias_key, sep = "||")
    ) %>%
    dplyr::arrange(dplyr::desc(.data$specificity), dplyr::desc(.data$alias_pref)) %>%
    dplyr::distinct(.data$lookup_key, .keep_all = TRUE)
  
  stats::setNames(lut_df$canonical, lut_df$lookup_key)
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
  
  x_key <- spike_alias_normalize_key(ifelse(missing_x, "", x_chr))
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

choose_context_alias_rows <- function(aliases, tool_i, spike_i) {
  if (is.null(aliases) || nrow(aliases) == 0) return(aliases[0, , drop = FALSE])
  
  aliases %>%
    dplyr::mutate(
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
    dplyr::filter(.data$tool_match >= 0L, .data$spike_match >= 0L) %>%
    dplyr::mutate(
      specificity = .data$tool_match + .data$spike_match,
      alias_pref = dplyr::if_else(.data$alias_key != .data$canonical_key, 1L, 0L)
    ) %>%
    dplyr::arrange(dplyr::desc(.data$specificity), dplyr::desc(.data$alias_pref)) %>%
    dplyr::distinct(.data$canonical_key, .keep_all = TRUE)
}

resolve_alias_path <- function(opt_path, metrics_dir) {
  if (!is.null(opt_path) && nzchar(opt_path)) return(opt_path)
  candidate <- file.path(metrics_dir, "taxon_aliases.resolved.csv")
  if (file.exists(candidate)) return(candidate)
  NULL
}

build_target_reference <- function(target_taxa, aliases, tool_i, spike_label_i) {
  target_taxa <- unique(as.character(target_taxa))
  target_taxa <- target_taxa[!is.na(target_taxa) & nzchar(target_taxa)]
  if (!length(target_taxa)) return(tibble::tibble(member_taxon = character(), target_norm = character()))
  
  base_ref <- tibble::tibble(
    member_taxon = target_taxa,
    target_name = target_taxa
  )
  
  a_sub <- choose_context_alias_rows(aliases, tolower(trimws(tool_i)), trimws(spike_label_i))
  if (nrow(a_sub) > 0) {
    alias_ref <- a_sub %>%
      dplyr::filter(.data$canonical %in% target_taxa) %>%
      dplyr::transmute(
        member_taxon = .data$canonical,
        target_name = .data$alias
      )
    base_ref <- dplyr::bind_rows(base_ref, alias_ref)
  }
  
  base_ref %>%
    dplyr::mutate(target_norm = normalize_taxon_core(.data$target_name)) %>%
    dplyr::filter(!is.na(.data$target_norm), nzchar(.data$target_norm)) %>%
    dplyr::distinct(.data$member_taxon, .data$target_norm)
}

design <- safe_read_csv(opt$design)
meta <- read_metadata_table(opt$meta)
run_manifest <- resolve_run_manifest(
  run_manifest = opt$run_manifest,
  profile_root = opt$profile_root,
  original_tables_arg = opt$original_tables,
  design = design
)

alias_path <- resolve_alias_path(opt$taxon_aliases, opt$metrics_dir)
aliases <- read_taxon_aliases(alias_path)
aliases <- prepare_alias_table(aliases, design = design)
alias_lookup <- build_alias_lookup(aliases)

if (nrow(aliases) > 0) {
  write_csv_safe(aliases, file.path(opt$outdir, "taxon_aliases.resolved.csv"))
}

if (!"base_id" %in% names(meta)) meta$base_id <- meta$sample_id
if (!"Target_Condition" %in% names(meta)) stop("Metadata must contain Target_Condition", call. = FALSE)
if (!"Study" %in% names(meta)) meta$Study <- NA_character_
meta$Study <- as.character(meta$Study)
meta$Study[is.na(meta$Study) | !nzchar(meta$Study)] <- "ALL"

filter_modes <- strsplit(opt$filter_modes, ",", fixed = TRUE)[[1]]
filter_modes <- trimws(filter_modes)
filter_modes <- filter_modes[nzchar(filter_modes)]
if (length(filter_modes) == 0) stop("No valid --filter_modes supplied", call. = FALSE)

selected_conditions_arg <- opt$background_conditions
if (!is.null(opt$background_condition) && nzchar(opt$background_condition)) {
  selected_conditions_arg <- opt$background_condition
}
selected_studies_arg <- opt$background_studies
if (!is.null(opt$background_study) && nzchar(opt$background_study)) {
  selected_studies_arg <- opt$background_study
}

available_conditions <- meta %>%
  dplyr::filter(!is.na(.data$Target_Condition), nzchar(.data$Target_Condition)) %>%
  dplyr::distinct(.data$Target_Condition) %>%
  dplyr::pull(.data$Target_Condition) %>%
  sort()

if (!length(available_conditions)) stop("No non-missing Target_Condition values found in metadata", call. = FALSE)

if (toupper(selected_conditions_arg) == "ALL") {
  selected_conditions <- available_conditions
} else {
  selected_conditions <- strsplit(selected_conditions_arg, ",", fixed = TRUE)[[1]]
  selected_conditions <- trimws(selected_conditions)
  selected_conditions <- selected_conditions[nzchar(selected_conditions)]
}

missing_conditions <- setdiff(selected_conditions, available_conditions)
if (length(missing_conditions) > 0) {
  stop(sprintf(
    "Unknown background condition(s): %s. Available: %s",
    paste(missing_conditions, collapse = ", "),
    paste(available_conditions, collapse = ", ")
  ), call. = FALSE)
}

available_studies <- meta %>%
  dplyr::filter(!is.na(.data$Study), nzchar(.data$Study)) %>%
  dplyr::distinct(.data$Study) %>%
  dplyr::pull(.data$Study) %>%
  sort()

if (!length(available_studies)) available_studies <- "ALL"

if (toupper(selected_studies_arg) == "ALL") {
  selected_studies <- available_studies
} else {
  selected_studies <- strsplit(selected_studies_arg, ",", fixed = TRUE)[[1]]
  selected_studies <- trimws(selected_studies)
  selected_studies <- selected_studies[nzchar(selected_studies)]
}

missing_studies <- setdiff(selected_studies, available_studies)
if (length(missing_studies) > 0) {
  stop(sprintf(
    "Unknown background study/studies: %s. Available: %s",
    paste(missing_studies, collapse = ", "),
    paste(available_studies, collapse = ", ")
  ), call. = FALSE)
}

bind_abundance_tables <- function(orig_tab, spiked_tab) {
  taxa <- union(setdiff(names(orig_tab), "sample_id"), setdiff(names(spiked_tab), "sample_id"))
  if (length(taxa) == 0) stop("No taxa columns available after combining abundance tables", call. = FALSE)
  
  add_missing <- function(df, taxa) {
    miss <- setdiff(taxa, names(df))
    for (nm in miss) df[[nm]] <- 0
    df %>% dplyr::select(sample_id, dplyr::all_of(taxa))
  }
  
  dplyr::bind_rows(add_missing(orig_tab, taxa), add_missing(spiked_tab, taxa))
}

guess_feature_columns <- function(df) {
  list(
    feature = first_existing(names(df), c("feature", "Feature", "taxon", "Taxon")),
    qval = first_existing(names(df), c("qval", "q_value", "qvalue", "adj_pval", "padj")),
    pval = first_existing(names(df), c("pval", "p_value", "pvalue")),
    coef = first_existing(names(df), c("coef", "coefficient", "estimate", "beta")),
    metadata = first_existing(names(df), c("metadata", "Metadata")),
    value = first_existing(names(df), c("value", "Value", "level", "Level"))
  )
}

run_one_maaslin_with_reference <- function(abundance_df, metadata_df, fixed_effect, reference_level, outdir) {
  require_maaslin2()
  ensure_dir(outdir)
  
  input_data <- make_table_matrix(abundance_df)
  input_meta <- metadata_df %>% as.data.frame()
  rownames(input_meta) <- input_meta$sample_id
  input_meta$sample_id <- NULL
  
  if (!(fixed_effect %in% names(input_meta))) {
    stop(sprintf("Fixed effect '%s' not present in metadata", fixed_effect), call. = FALSE)
  }
  
  levs <- unique(as.character(input_meta[[fixed_effect]]))
  levs <- c(reference_level, setdiff(levs, reference_level))
  input_meta[[fixed_effect]] <- factor(as.character(input_meta[[fixed_effect]]), levels = unique(levs))
  
  maaslin_fun <- getExportedValue("Maaslin2", "Maaslin2")
  args <- list(
    input_data = input_data,
    input_metadata = input_meta,
    output = outdir,
    fixed_effects = fixed_effect,
    normalization = "NONE",
    transform = "LOG",
    analysis_method = "LM",
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    reference = paste(fixed_effect, reference_level, sep = ",")
  )
  do.call(maaslin_fun, args)
}

annotate_maaslin_results <- function(res_df, target_taxa, fixed_effect, q_threshold,
                                     aliases = NULL, alias_lookup = NULL,
                                     tool = NULL, spike_label = NULL) {
  if (!nrow(res_df)) return(tibble::tibble())
  
  cols <- guess_feature_columns(res_df)
  if (is.na(cols$feature)) return(tibble::tibble())
  
  out <- res_df
  if (is.na(cols$qval)) out$qval_proxy <- NA_real_
  if (is.na(cols$pval)) out$pval_proxy <- NA_real_
  if (is.na(cols$coef)) out$coef_proxy <- NA_real_
  if (is.na(cols$metadata)) out$metadata_proxy <- fixed_effect
  if (is.na(cols$value)) out$value_proxy <- NA_character_
  
  cols <- guess_feature_columns(out)
  if (is.na(cols$qval)) cols$qval <- "qval_proxy"
  if (is.na(cols$pval)) cols$pval <- "pval_proxy"
  if (is.na(cols$coef)) cols$coef <- "coef_proxy"
  if (is.na(cols$metadata)) cols$metadata <- "metadata_proxy"
  if (is.na(cols$value)) cols$value <- "value_proxy"
  
  target_ref <- build_target_reference(target_taxa, aliases, tool_i = tool %||% NA_character_, spike_label_i = spike_label %||% NA_character_)
  
  out <- out %>%
    dplyr::mutate(
      feature_raw = as.character(.data[[cols$feature]]),
      qval = suppressWarnings(as.numeric(.data[[cols$qval]])),
      pval = suppressWarnings(as.numeric(.data[[cols$pval]])),
      coef = suppressWarnings(as.numeric(.data[[cols$coef]])),
      metadata_name = as.character(.data[[cols$metadata]]),
      metadata_value = as.character(.data[[cols$value]])
    )
  
  if (length(alias_lookup) > 0L) {
    out$feature <- canonicalize_taxa_vector(out$feature_raw, alias_lookup, tool = tool, spike_label = spike_label)
  } else {
    out$feature <- out$feature_raw
  }
  
  out %>%
    dplyr::mutate(feature_norm = normalize_taxon_core(.data$feature)) %>%
    dplyr::left_join(target_ref, by = c("feature_norm" = "target_norm")) %>%
    dplyr::mutate(
      is_target = !is.na(.data$member_taxon),
      is_significant = is.finite(.data$qval) & .data$qval < q_threshold,
      effect_direction = dplyr::case_when(
        is.na(.data$coef) ~ NA_character_,
        .data$coef > 0 ~ "positive",
        .data$coef < 0 ~ "negative",
        TRUE ~ "zero"
      ),
      is_positive = is.finite(.data$coef) & .data$coef > 0,
      is_negative = is.finite(.data$coef) & .data$coef < 0
    )
}

summarise_member_detection <- function(annotated_df, row_meta, target_taxa) {
  if (length(target_taxa) == 0) return(tibble::tibble())
  target_tbl <- tibble::tibble(
    member_taxon = unique(as.character(target_taxa)),
    feature_norm = normalize_taxon_core(unique(as.character(target_taxa)))
  )
  
  if (!nrow(annotated_df)) {
    return(target_tbl %>%
             dplyr::transmute(
               background_condition = row_meta$background_condition[[1]],
               background_study = row_meta$background_study[[1]] %||% "ALL",
               tool = row_meta$tool[[1]],
               spike_mode = row_meta$spike_mode[[1]] %||% NA_character_,
               spike_label = row_meta$spike_label[[1]],
               spike_fraction = row_meta$spike_fraction[[1]],
               filter_mode = row_meta$filter_mode[[1]],
               member_taxon = .data$member_taxon,
               member_detected_any = FALSE,
               member_detected_positive = FALSE,
               min_q_target = NA_real_,
               max_coef_target = NA_real_
             ))
  }
  
  target_rows <- annotated_df %>% dplyr::filter(.data$is_target)
  target_tbl %>%
    dplyr::left_join(
      target_rows %>%
        dplyr::group_by(.data$member_taxon) %>%
        dplyr::summarise(
          member_detected_any = any(.data$is_significant, na.rm = TRUE),
          member_detected_positive = any(.data$is_significant & .data$is_positive, na.rm = TRUE),
          min_q_target = suppressWarnings(min(.data$qval[is.finite(.data$qval)], na.rm = TRUE)),
          max_coef_target = suppressWarnings(max(.data$coef[is.finite(.data$coef)], na.rm = TRUE)),
          .groups = "drop"
        ),
      by = "member_taxon"
    ) %>%
    dplyr::mutate(
      member_detected_any = dplyr::coalesce(.data$member_detected_any, FALSE),
      member_detected_positive = dplyr::coalesce(.data$member_detected_positive, FALSE),
      min_q_target = dplyr::if_else(is.infinite(.data$min_q_target), NA_real_, .data$min_q_target),
      background_condition = row_meta$background_condition[[1]],
      background_study = row_meta$background_study[[1]] %||% "ALL",
      tool = row_meta$tool[[1]],
      spike_mode = row_meta$spike_mode[[1]] %||% NA_character_,
      spike_label = row_meta$spike_label[[1]],
      spike_fraction = row_meta$spike_fraction[[1]],
      filter_mode = row_meta$filter_mode[[1]]
    ) %>%
    dplyr::select("background_condition", "background_study", "tool", "spike_mode", "spike_label", "spike_fraction", "filter_mode", dplyr::everything())
}

base_meta <- meta %>%
  dplyr::distinct(.data$base_id, .data$original_id, .data$Target_Condition, .data$Study)

summary_rows <- list()
member_rows <- list()
sig_feature_rows <- list()
idx <- 1L
midx <- 1L
sidx <- 1L

for (cond in selected_conditions) {
  eligible_base_meta_cond <- base_meta %>% dplyr::filter(.data$Target_Condition == cond)
  if (nrow(eligible_base_meta_cond) == 0) {
    message(sprintf("[WARN] No base samples found for background condition '%s'; skipping", cond))
    next
  }

  cond_studies <- intersect(selected_studies, sort(unique(eligible_base_meta_cond$Study)))
  if (!length(cond_studies)) {
    message(sprintf("[WARN] No matching studies found for background condition '%s'; skipping", cond))
    next
  }

  for (bg_study in cond_studies) {
    eligible_base_meta <- eligible_base_meta_cond %>% dplyr::filter(.data$Study == bg_study)
    if (nrow(eligible_base_meta) == 0) next

    message(sprintf("[INFO] Running background condition: %s | study: %s", cond, bg_study))

    for (i in seq_len(nrow(run_manifest))) {
      man <- run_manifest[i, , drop = FALSE]
      tool_i <- as.character(man$tool[[1]])
      spike_label_i <- as.character(man$spike_label[[1]])
      spike_mode_i <- if ("spike_mode" %in% names(man)) as.character(man$spike_mode[[1]]) else NA_character_
      spike_fraction_i <- as.numeric(man$spike_fraction[[1]])

      if (!"original_table" %in% names(man) || is.na(man$original_table[[1]]) || !nzchar(man$original_table[[1]])) {
        message(sprintf("[WARN] Skipping %s/%s/%s because original_table is missing", tool_i, spike_label_i, fraction_to_tag(spike_fraction_i)))
        next
      }

      tab_spiked <- read_abundance_table(man$spiked_table[[1]])
      tab_orig <- read_abundance_table(man$original_table[[1]])

      spiked_meta <- meta %>%
        dplyr::filter(.data$sample_id %in% tab_spiked$sample_id, .data$Target_Condition == cond, .data$Study == bg_study)

      if (!"base_id" %in% names(spiked_meta)) spiked_meta$base_id <- spiked_meta$sample_id
      if (!is.na(spike_mode_i) && "spike_mode" %in% names(spiked_meta)) {
        spiked_meta <- spiked_meta %>% dplyr::filter(.data$spike_mode == spike_mode_i)
      }
      if ("spike_label" %in% names(spiked_meta)) {
        spiked_meta <- spiked_meta %>% dplyr::filter(.data$spike_label == spike_label_i)
      }
      if ("spike_fraction_total" %in% names(spiked_meta)) {
        spiked_meta <- spiked_meta %>% dplyr::filter(abs(.data$spike_fraction_total - spike_fraction_i) < 1e-12)
      }

      candidate_base_ids <- intersect(unique(spiked_meta$base_id), unique(tab_orig$sample_id))
      if (length(candidate_base_ids) == 0) {
        message(sprintf("[WARN] No overlapping base IDs for %s | %s | %s | %s | %s", cond, bg_study, tool_i, spike_label_i, fraction_to_tag(spike_fraction_i)))
        next
      }

      spiked_meta <- spiked_meta %>% dplyr::filter(.data$base_id %in% candidate_base_ids)
      tab_spiked_sub <- tab_spiked %>% dplyr::filter(.data$sample_id %in% spiked_meta$sample_id)
      tab_orig_sub <- tab_orig %>% dplyr::filter(.data$sample_id %in% candidate_base_ids)

      orig_meta <- eligible_base_meta %>%
        dplyr::filter(.data$base_id %in% candidate_base_ids) %>%
        dplyr::transmute(
          sample_id = .data$base_id,
          base_id = .data$base_id,
          original_id = .data$original_id,
          Target_Condition = .data$Target_Condition,
          Study = .data$Study,
          spike_status = "unspiked"
        )

      spiked_meta2 <- spiked_meta %>%
        dplyr::transmute(
          sample_id = .data$sample_id,
          base_id = .data$base_id,
          original_id = .data$original_id,
          Target_Condition = .data$Target_Condition,
          Study = .data$Study,
          spike_status = "spiked"
        )

      meta_sub <- dplyr::bind_rows(orig_meta, spiked_meta2) %>%
        dplyr::mutate(spike_status = factor(.data$spike_status, levels = c(opt$reference_level, setdiff(c("unspiked", "spiked"), opt$reference_level))))

      if (nrow(tab_orig_sub) != nrow(orig_meta) || nrow(tab_spiked_sub) != nrow(spiked_meta2)) {
        message(sprintf("[WARN] Sample/metadata mismatch for %s | %s | %s | %s | %s; skipping", cond, bg_study, tool_i, spike_label_i, fraction_to_tag(spike_fraction_i)))
        next
      }

      if (sum(meta_sub$spike_status == "unspiked") < opt$min_samples_per_group || sum(meta_sub$spike_status == "spiked") < opt$min_samples_per_group) {
        message(sprintf("[WARN] Skipping %s | %s | %s | %s | %s due to insufficient samples per group", cond, bg_study, tool_i, spike_label_i, fraction_to_tag(spike_fraction_i)))
        next
      }

      target_taxa <- design %>%
        dplyr::filter(.data$spike_label == spike_label_i, abs(.data$spike_fraction_total - spike_fraction_i) < 1e-12) %>%
        dplyr::pull(.data$member_taxon) %>%
        unique()

      combined_tab <- bind_abundance_tables(tab_orig_sub, tab_spiked_sub)

      for (mode in filter_modes) {
        filt <- filter_abundance_table(
          combined_tab,
          mode = mode,
          tool = tool_i,
          metrics_dir = opt$metrics_dir,
          alias_table = aliases,
          spike_label = spike_label_i
        )

        run_dir <- file.path(
          opt$outdir,
          "by_condition_and_study",
          sanitize_slug(cond),
          sanitize_slug(bg_study),
          tool_i,
          spike_label_i,
          fraction_to_tag(spike_fraction_i),
          mode
        )
        ensure_dir(run_dir)
        write_csv_safe(filt$audit, file.path(run_dir, "filter_audit.csv"))
        writeLines(filt$audit$taxon[filt$audit$keep], file.path(run_dir, "kept_taxa.txt"))
        writeLines(filt$audit$taxon[!filt$audit$keep], file.path(run_dir, "removed_taxa.txt"))

        run_one_maaslin_with_reference(
          abundance_df = filt$table,
          metadata_df = meta_sub,
          fixed_effect = opt$fixed_effect,
          reference_level = opt$reference_level,
          outdir = run_dir
        )

        res <- read_maaslin_results(run_dir)
        anno <- annotate_maaslin_results(
          res,
          target_taxa = target_taxa,
          fixed_effect = opt$fixed_effect,
          q_threshold = opt$q_threshold,
          aliases = aliases,
          alias_lookup = alias_lookup,
          tool = tool_i,
          spike_label = spike_label_i
        )

        enriched_sig <- anno %>% dplyr::filter(.data$is_significant, .data$is_positive)
        depleted_sig <- anno %>% dplyr::filter(.data$is_significant, .data$is_negative)
        target_sig <- anno %>% dplyr::filter(.data$is_significant, .data$is_target)
        target_sig_enriched <- enriched_sig %>% dplyr::filter(.data$is_target)
        target_sig_depleted <- depleted_sig %>% dplyr::filter(.data$is_target)
        offtarget_sig <- anno %>% dplyr::filter(.data$is_significant, !.data$is_target)
        offtarget_sig_enriched <- enriched_sig %>% dplyr::filter(!.data$is_target)
        offtarget_sig_depleted <- depleted_sig %>% dplyr::filter(!.data$is_target)

        summary_rows[[idx]] <- tibble::tibble(
          background_condition = cond,
          background_study = bg_study,
          tool = tool_i,
          spike_mode = spike_mode_i,
          spike_label = spike_label_i,
          spike_fraction = spike_fraction_i,
          filter_mode = mode,
          output_dir = run_dir,
          n_unspiked_samples = sum(meta_sub$spike_status == "unspiked"),
          n_spiked_samples = sum(meta_sub$spike_status == "spiked"),
          n_taxa_input = ncol(combined_tab) - 1,
          n_taxa_kept = ncol(filt$table) - 1,
          n_target_members = length(unique(target_taxa)),
          n_sig_total = nrow(anno %>% dplyr::filter(.data$is_significant)),
          n_sig_enriched = nrow(enriched_sig),
          n_sig_depleted = nrow(depleted_sig),
          n_sig_target = nrow(target_sig),
          n_sig_target_enriched = nrow(target_sig_enriched),
          n_sig_target_depleted = nrow(target_sig_depleted),
          n_sig_offtarget = nrow(offtarget_sig),
          n_sig_offtarget_enriched = nrow(offtarget_sig_enriched),
          n_sig_offtarget_depleted = nrow(offtarget_sig_depleted),
          any_target_sig = nrow(target_sig) > 0,
          any_target_sig_enriched = nrow(target_sig_enriched) > 0
        )
        idx <- idx + 1L

        member_meta <- tibble::tibble(
          background_condition = cond,
          background_study = bg_study,
          tool = tool_i,
          spike_mode = spike_mode_i,
          spike_label = spike_label_i,
          spike_fraction = spike_fraction_i,
          filter_mode = mode
        )
        member_out <- summarise_member_detection(anno, member_meta, target_taxa)
        if (nrow(member_out) > 0) {
          member_rows[[midx]] <- member_out
          midx <- midx + 1L
        }

        if (nrow(anno) > 0) {
          sig_feature_rows[[sidx]] <- anno %>%
            dplyr::mutate(
              background_condition = cond,
              background_study = bg_study,
              tool = tool_i,
              spike_mode = spike_mode_i,
              spike_label = spike_label_i,
              spike_fraction = spike_fraction_i,
              filter_mode = mode
            )
          sidx <- sidx + 1L
        }
      }
    }
  }
}

summary_df <- dplyr::bind_rows(summary_rows)
member_df <- dplyr::bind_rows(member_rows)
sig_df <- dplyr::bind_rows(sig_feature_rows)

write_csv_safe(summary_df, file.path(opt$outdir, "maaslin_run_summary_ALLFILTERS.csv"))
write_csv_safe(run_manifest, file.path(opt$outdir, "run_design_manifest_ALLFILTERS.csv"))
write_csv_safe(tibble::tibble(background_condition = selected_conditions), file.path(opt$outdir, "background_conditions_used.csv"))
write_csv_safe(tibble::tibble(background_study = selected_studies), file.path(opt$outdir, "background_studies_used.csv"))
if (nrow(aliases) > 0) write_csv_safe(aliases, file.path(opt$outdir, "taxon_aliases.resolved.csv"))
if (nrow(member_df) > 0) write_csv_safe(member_df, file.path(opt$outdir, "maaslin_member_detection_ALLFILTERS.csv"))
if (nrow(sig_df) > 0) write_csv_safe(sig_df, file.path(opt$outdir, "maaslin_significant_features_ALLFILTERS.csv"))

cat(sprintf("[OK] Wrote outputs under %s\n", opt$outdir))
