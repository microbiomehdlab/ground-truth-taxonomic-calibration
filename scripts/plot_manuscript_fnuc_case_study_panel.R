#!/usr/bin/env Rscript
# Publication-polished Fig. 2 case-study panel: compact recovery-class panel with stronger fraction blocks.
# v3 publication polish for the revised Fig. 2 recovery-class panel.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(scales)
  library(patchwork)
  library(forcats)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

# Robust CLI parser so dashed flags like --metadata-input always work
args <- commandArgs(trailingOnly = TRUE)
arg_map <- list()
i <- 1
while (i <= length(args)) {
  a <- args[[i]]
  if (startsWith(a, '--')) {
    if (grepl('=', a, fixed = TRUE)) {
      key <- sub('^--', '', strsplit(a, '=', fixed = TRUE)[[1]][1])
      val <- sub('^[^=]*=', '', a)
      arg_map[[key]] <- val
      arg_map[[gsub('-', '_', key)]] <- val
    } else {
      key <- sub('^--', '', a)
      if (i < length(args) && !startsWith(args[[i+1]], '--')) {
        val <- args[[i+1]]
        arg_map[[key]] <- val
        arg_map[[gsub('-', '_', key)]] <- val
        i <- i + 1
      } else {
        arg_map[[key]] <- TRUE
        arg_map[[gsub('-', '_', key)]] <- TRUE
      }
    }
  }
  i <- i + 1
}
get_arg <- function(name, default = NULL) {
  arg_map[[name]] %||% arg_map[[gsub('-', '_', name)]] %||% default
}
opt <- list(
  indir = get_arg('indir'),
  outdir = get_arg('outdir', 'RUNS/plots_manuscript_fnuc_case_study_boxplotC'),
  spike_label = get_arg('spike-label', 'Fnuc'),
  manifest_input = get_arg('manifest-input'),
  design_input = get_arg('design-input'),
  aliases_input = get_arg('aliases-input'),
  target_input = get_arg('target-input'),
  kraken_input = get_arg('kraken-input'),
  metaphlan_input = get_arg('metaphlan-input'),
  metadata_input = get_arg('metadata-input'),
  metadata_sample_col = get_arg('metadata-sample-col'),
  metadata_condition_col = get_arg('metadata-condition-col'),
  metadata_study_col = get_arg('metadata-study-col'),
  abundance_scale = get_arg('abundance-scale', 'auto'),
  prevalence_threshold = as.numeric(get_arg('prevalence-threshold', 0)),
  max_weak_fraction = as.numeric(get_arg('max-weak-fraction', 0.001)),
  good_relative_error = as.numeric(get_arg('good-relative-error', 0.10)),
  avg_relative_error = as.numeric(get_arg('avg-relative-error', 0.50)),
  allow_metadata_subset = isTRUE(get_arg('allow-metadata-subset', FALSE)),
  panelc_ymax = as.numeric(get_arg('panelc-ymax', NA)),
  width = as.numeric(get_arg('width', 13.4)),
  height = as.numeric(get_arg('height', 10.0)),
  dpi = as.integer(get_arg('dpi', 400))
)
message('[DEBUG] Parsed --indir: ', opt$indir %||% '<NULL>')
message('[DEBUG] Parsed --metadata-input: ', opt$metadata_input %||% '<NULL>')
message('[DEBUG] Parsed --kraken-input: ', opt$kraken_input %||% '<NULL>')
message('[DEBUG] Parsed --metaphlan-input: ', opt$metaphlan_input %||% '<NULL>')
stop_if_missing_file <- function(path, label) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) stop(sprintf("Missing %s: %s", label, path %||% "NULL"), call. = FALSE)
}
read_any_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) read_tsv(path, show_col_types = FALSE) else read_csv(path, show_col_types = FALSE)
}
first_existing <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit)) hit[[1]] else NA_character_
}
clean_key <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\|", " ")
  x <- str_replace_all(x, "(^|[ ;])s__", " ")
  x <- str_replace_all(x, "(^|[ ;])t__", " ")
  x <- str_replace_all(x, "_", " ")
  x <- str_to_lower(x)
  x <- str_replace_all(x, "[^a-z0-9]+", " ")
  str_squish(x)
}
sample_key <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\\\", "/")
  x <- basename(x)
  x <- str_replace(x, "\\.gz$", "")
  x <- str_replace(x, "\\.bz2$", "")
  x <- str_replace(x, "\\.xz$", "")
  x <- str_replace(x, "\\.(fastq|fq|fasta|fa|sam|bam|csv|tsv|txt)$", "")
  x <- str_replace(x, "(_R[12]|_[12])$", "")
  x <- str_replace(x, "(_kneaddata|_metaphlan.*|_profile.*|_bracken.*|_kraken.*)$", "")
  str_trim(x)
}
species_fallback <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "_", " ")
  x <- str_replace_all(x, "^(s__|t__)", "")
  x <- str_squish(x)
  vapply(str_split(x, "\\s+"), function(tok) {
    tok <- tok[nzchar(tok)]
    if (length(tok) >= 2) paste(tok[1:2], collapse = " ") else paste(tok, collapse = " ")
  }, character(1))
}
fmt_fraction <- function(x) sprintf("%.2f%%", 100 * x)
pretty_tool <- function(x) dplyr::case_when(
  str_detect(x, regex("kraken|bracken", ignore_case = TRUE)) ~ "Kraken2 + Bracken",
  str_detect(x, regex("metaphlan", ignore_case = TRUE)) ~ "MetaPhlAn 4",
  TRUE ~ as.character(x)
)

if (is.null(opt$indir) || !nzchar(opt$indir)) stop("Required: --indir", call. = FALSE)
if (is.null(opt$metadata_input) || !nzchar(opt$metadata_input)) stop("Required: --metadata-input", call. = FALSE)
if (opt$good_relative_error >= opt$avg_relative_error) stop("--good-relative-error must be smaller than --avg-relative-error", call. = FALSE)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

indir <- opt$indir
manifest_path <- opt$manifest_input %||% file.path(indir, "run_manifest.resolved.csv")
design_path <- opt$design_input %||% file.path(indir, "spike_design_long.csv")
aliases_path <- opt$aliases_input %||% file.path(indir, "taxon_aliases.resolved.csv")
target_path <- opt$target_input %||% file.path(indir, "target_member_errors_with_condition.csv")
stop_if_missing_file(design_path, "spike_design_long.csv")
stop_if_missing_file(target_path, "target_member_errors_with_condition.csv")
stop_if_missing_file(opt$metadata_input, "metadata table")

design <- read_any_table(design_path)
targets <- read_any_table(target_path)
manifest <- if (file.exists(manifest_path)) read_any_table(manifest_path) else NULL

selected_label <- opt$spike_label
label_design <- design %>% filter(.data$spike_label == selected_label)
if (nrow(label_design) == 0) stop(sprintf("spike_label '%s' not found in spike_design_long.csv", selected_label), call. = FALSE)
label_taxa <- label_design %>% distinct(spike_label, member_taxon)
taxon_title <- unique(species_fallback(label_taxa$member_taxon))[[1]]

read_abundance_wide <- function(path, abundance_scale = "auto") {
  stop_if_missing_file(path, "abundance table")
  dat <- read_any_table(path)
  sid <- first_existing(names(dat), c("sample_id", "sample", "Sample", "SampleID", "ID", "base_id"))
  if (is.na(sid)) sid <- names(dat)[[1]]
  dat <- dat %>% rename(sample_id = all_of(sid))
  dat$sample_id <- as.character(dat$sample_id)
  taxa <- setdiff(names(dat), "sample_id")
  dat[taxa] <- lapply(dat[taxa], function(x) { x <- suppressWarnings(as.numeric(x)); x[is.na(x)] <- 0; x })
  mat <- as.matrix(dat[taxa]); storage.mode(mat) <- "double"
  finite_vals <- as.numeric(mat[is.finite(mat)])
  row_sums <- rowSums(mat, na.rm = TRUE)
  max_value <- if (length(finite_vals)) max(finite_vals, na.rm = TRUE) else NA_real_
  med_sum <- if (length(row_sums)) median(row_sums, na.rm = TRUE) else NA_real_
  inferred <- if (identical(abundance_scale, "auto")) {
    if ((is.finite(max_value) && max_value > 1.0000001) || (is.finite(med_sum) && med_sum > 1.5)) "percent" else "fraction"
  } else abundance_scale
  if (identical(inferred, "percent")) dat[taxa] <- lapply(dat[taxa], function(x) x / 100)
  dat
}
make_candidate_map <- function(label_taxa, tool_name, aliases) {
  base <- label_taxa %>% transmute(spike_label, canonical = member_taxon, candidate = member_taxon)
  fallback <- label_taxa %>% transmute(spike_label, canonical = member_taxon, candidate = species_fallback(member_taxon))
  al <- tibble(canonical = character(), alias = character(), tool = character(), spike_label = character(), candidate = character())
  if (nrow(aliases) > 0) {
    alias_tool_regex <- ifelse(str_detect(tool_name, "Kraken"), "kraken|bracken", "metaphlan")
    al_exact <- aliases %>%
      mutate(tool_chr = coalesce(as.character(tool), "")) %>%
      filter(tool_chr == "" | str_detect(tool_chr, regex(alias_tool_regex, ignore_case = TRUE))) %>%
      inner_join(label_taxa %>% distinct(spike_label, member_taxon), by = c("canonical" = "member_taxon"), relationship = "many-to-many") %>%
      transmute(spike_label, canonical, candidate = alias)
    al_species <- aliases %>%
      mutate(tool_chr = coalesce(as.character(tool), ""), canonical_species = species_fallback(canonical)) %>%
      filter(tool_chr == "" | str_detect(tool_chr, regex(alias_tool_regex, ignore_case = TRUE))) %>%
      inner_join(label_taxa %>% distinct(spike_label, member_taxon) %>% mutate(canonical_species = species_fallback(member_taxon)), by = "canonical_species", relationship = "many-to-many") %>%
      transmute(spike_label, canonical = member_taxon, candidate = alias)
    al <- bind_rows(al_exact, al_species)
  }
  bind_rows(base, fallback, al) %>%
    filter(!is.na(candidate), nzchar(candidate)) %>%
    distinct(spike_label, canonical, candidate) %>%
    mutate(candidate_key = clean_key(candidate))
}
extract_label_abundances <- function(wide, tool_display, label_taxa, aliases) {
  taxa_cols <- setdiff(names(wide), "sample_id")
  col_map <- tibble(column = taxa_cols, column_key = clean_key(taxa_cols))
  candidates <- make_candidate_map(label_taxa, tool_display, aliases)
  matches <- candidates %>%
    inner_join(col_map, by = c("candidate_key" = "column_key"), relationship = "many-to-many") %>%
    distinct(spike_label, canonical, candidate, column)
  cols <- matches %>% filter(spike_label == selected_label) %>% pull(column) %>% unique()
  abundance <- if (length(cols) == 0) rep(0, nrow(wide)) else rowSums(wide[, cols, drop = FALSE], na.rm = TRUE)
  list(data = tibble(sample_id = wide$sample_id, tool = tool_display, abundance = abundance), matches = matches)
}
choose_best_sample_col <- function(md, full_sample_ids, explicit_col = NULL) {
  if (!is.null(explicit_col) && nzchar(explicit_col)) {
    if (!explicit_col %in% names(md)) stop("--metadata-sample-col not found in metadata.", call. = FALSE)
    return(explicit_col)
  }
  full_keys <- unique(sample_key(full_sample_ids))
  candidate_cols <- names(md)[vapply(md, function(z) is.character(z) || is.factor(z) || is.numeric(z), logical(1))]
  if (!length(candidate_cols)) candidate_cols <- names(md)
  scores <- bind_rows(lapply(candidate_cols, function(cc) {
    vals <- as.character(md[[cc]])
    keys <- unique(sample_key(vals[!is.na(vals) & nzchar(vals)]))
    overlap <- length(intersect(keys, full_keys))
    tibble(column = cc, overlap = overlap, frac_of_full = overlap / length(full_keys))
  })) %>% arrange(desc(overlap), desc(frac_of_full))
  scores$column[[1]]
}
read_metadata <- function(path, full_sample_ids) {
  md <- read_any_table(path)
  sid <- choose_best_sample_col(md, full_sample_ids, opt$metadata_sample_col)
  cond <- if (!is.null(opt$metadata_condition_col) && nzchar(opt$metadata_condition_col)) opt$metadata_condition_col else first_existing(names(md), c("Target_Condition","condition","Condition","group","Group","diagnosis","Diagnosis","disease","Disease","status","Status"))
  study <- if (!is.null(opt$metadata_study_col) && nzchar(opt$metadata_study_col)) opt$metadata_study_col else first_existing(names(md), c("Study","study","Dataset","dataset","Cohort","cohort","dataset_name","source"))
  if (is.na(cond) || !cond %in% names(md)) stop("Metadata needs a condition column.", call. = FALSE)
  if (is.na(study) || !study %in% names(md)) stop("Metadata needs a study column.", call. = FALSE)
  md %>%
    transmute(sample_key = sample_key(.data[[sid]]),
              Target_Condition = as.character(.data[[cond]]),
              Study = as.character(.data[[study]])) %>%
    filter(!is.na(sample_key), nzchar(sample_key), !is.na(Target_Condition), nzchar(Target_Condition), !is.na(Study), nzchar(Study)) %>%
    distinct(sample_key, Target_Condition, Study, .keep_all = TRUE)
}

if (is.null(opt$kraken_input) || is.null(opt$metaphlan_input)) {
  if (is.null(manifest) || !all(c("tool", "original_table") %in% names(manifest))) {
    stop("Provide --kraken-input and --metaphlan-input, or give --indir with run_manifest.resolved.csv.", call. = FALSE)
  }
  originals <- manifest %>% distinct(tool, original_table) %>% filter(!is.na(original_table), nzchar(original_table))
  if (is.null(opt$kraken_input)) opt$kraken_input <- originals %>% filter(str_detect(tool, regex("kraken|bracken", ignore_case = TRUE))) %>% pull(original_table) %>% unique() %>% .[[1]]
  if (is.null(opt$metaphlan_input)) opt$metaphlan_input <- originals %>% filter(str_detect(tool, regex("metaphlan", ignore_case = TRUE))) %>% pull(original_table) %>% unique() %>% .[[1]]
}
stop_if_missing_file(opt$kraken_input, "Kraken2+Bracken abundance table")
stop_if_missing_file(opt$metaphlan_input, "MetaPhlAn 4 abundance table")

aliases <- tibble(canonical = character(), alias = character(), tool = character())
if (!is.null(aliases_path) && file.exists(aliases_path)) {
  aliases <- read_any_table(aliases_path) %>%
    transmute(canonical = as.character(canonical), alias = as.character(alias), tool = as.character(tool)) %>%
    filter(!is.na(canonical), !is.na(alias))
}

kraken_wide <- read_abundance_wide(opt$kraken_input, opt$abundance_scale)
metaphlan_wide <- read_abundance_wide(opt$metaphlan_input, opt$abundance_scale)
kr <- extract_label_abundances(kraken_wide, "Kraken2 + Bracken", label_taxa, aliases)
mp <- extract_label_abundances(metaphlan_wide, "MetaPhlAn 4", label_taxa, aliases)
full_baseline <- bind_rows(kr$data, mp$data) %>% mutate(sample_key = sample_key(sample_id))
metadata <- read_metadata(opt$metadata_input, full_baseline$sample_id)
coverage_join <- full_baseline %>% distinct(tool, sample_id, sample_key) %>% left_join(metadata, by = "sample_key")
coverage <- coverage_join %>% summarise(total_tool_samples = n(), with_metadata = sum(!is.na(Target_Condition) & !is.na(Study)), metadata_fraction = with_metadata / total_tool_samples)
if (!isTRUE(opt$allow_metadata_subset) && coverage$metadata_fraction < 0.95) {
  stop(sprintf("Metadata cover only %.1f%% of tool-samples. Use a fuller metadata table or rerun with --allow-metadata-subset if intentional.", 100 * coverage$metadata_fraction), call. = FALSE)
}

condition_levels <- c("Control", "Adenoma", "colorectal carcinoma")
study_levels <- c("FengQ_2015", "ZellerG_2014")
tool_levels <- c("Kraken2 + Bracken", "MetaPhlAn 4")

baseline_df <- full_baseline %>%
  left_join(metadata, by = "sample_key") %>%
  filter(!is.na(Target_Condition), !is.na(Study)) %>%
  mutate(Target_Condition = factor(Target_Condition, levels = condition_levels),
         Study = factor(Study, levels = study_levels),
         tool = factor(tool, levels = tool_levels),
         positive = abundance > opt$prevalence_threshold) %>%
  filter(!is.na(Target_Condition), !is.na(Study), !is.na(tool))

wilson <- function(x, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2)
  p <- x / n
  den <- 1 + z^2 / n
  centre <- p + z^2 / (2 * n)
  half <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n)
  tibble(ymin = pmax(0, (centre - half) / den), ymax = pmin(1, (centre + half) / den))
}

baseline_prev <- baseline_df %>%
  group_by(tool, Study, Target_Condition) %>%
  summarise(n = n_distinct(sample_id), positive_n = sum(positive, na.rm = TRUE), prevalence = positive_n / n, .groups = "drop") %>%
  rowwise() %>% mutate(ci = list(wilson(positive_n, n))) %>% unnest(ci) %>% ungroup()

recovery_df <- targets %>%
  filter(spike_label == selected_label) %>%
  { if ("target_flag" %in% names(.)) dplyr::filter(., !is.na(target_flag) & target_flag) else . } %>%
  mutate(
    tool = factor(pretty_tool(tool), levels = tool_levels),
    Target_Condition = factor(Target_Condition, levels = condition_levels),
    Study = factor(Study, levels = study_levels),
    rel_dev = abs(observed_over_expected - 1),
    detected = as.logical(detected)
  ) %>%
  filter(!is.na(tool), !is.na(Target_Condition), !is.na(Study))
if (nrow(recovery_df) == 0) stop(sprintf("No rows found in target_member_errors_with_condition.csv for spike_label '%s'", selected_label), call. = FALSE)

all_spikes_df <- recovery_df %>%
  mutate(
    spike_fraction_total = as.numeric(spike_fraction_total),
    spike_fraction_label = factor(fmt_fraction(spike_fraction_total), levels = fmt_fraction(sort(unique(spike_fraction_total)))),
    recovery_class = case_when(
      !detected ~ "Poor / missed",
      rel_dev <= opt$good_relative_error ~ "Good",
      rel_dev <= opt$avg_relative_error ~ "Intermediate",
      TRUE ~ "Poor / missed"
    ),
    recovery_class = factor(recovery_class, levels = c("Good", "Intermediate", "Poor / missed"))
  )

weak_df <- all_spikes_df %>% filter(spike_fraction_total <= opt$max_weak_fraction)
if (nrow(weak_df) == 0) stop("No spike rows remained after applying --max-weak-fraction.", call. = FALSE)

# Revised Fig. 2 panel C order: show the highest weak spike fraction at the top
# and the lowest at the bottom, i.e. 0.10% -> 0.05% -> 0.01%.
weak_fraction_levels_desc <- weak_df %>%
  distinct(spike_fraction_total, spike_fraction_label) %>%
  arrange(desc(spike_fraction_total)) %>%
  pull(spike_fraction_label) %>%
  as.character()

class_sum <- weak_df %>%
  group_by(tool, spike_fraction_total, spike_fraction_label, Target_Condition, recovery_class) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(total_n = sum(n), pct = n / total_n) %>%
  ungroup()
class_labels <- class_sum %>% mutate(label = ifelse(pct >= 0.08, paste0(round(100 * pct), "%"), ""))

class_sum_all <- all_spikes_df %>%
  group_by(tool, spike_fraction_total, spike_fraction_label, Target_Condition, recovery_class) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(total_n = sum(n), pct = n / total_n) %>%
  ungroup()
class_labels_all <- class_sum_all %>% mutate(label = ifelse(pct >= 0.08, paste0(round(100 * pct), "%"), ""))

# Panel C data: pooled across studies, show all samples via boxplots
box_df <- weak_df %>%
  mutate(
    spike_fraction_label = factor(spike_fraction_label, levels = levels(spike_fraction_label)),
    Target_Condition = factor(Target_Condition, levels = condition_levels),
    tool = factor(tool, levels = tool_levels)
  )

box_df_all <- all_spikes_df %>%
  mutate(
    spike_fraction_label = factor(spike_fraction_label, levels = levels(spike_fraction_label)),
    Target_Condition = factor(Target_Condition, levels = condition_levels),
    tool = factor(tool, levels = tool_levels)
  )

box_stats <- box_df %>%
  group_by(tool, spike_fraction_label, Target_Condition) %>%
  summarise(
    n = n(),
    median_oe = median(observed_over_expected, na.rm = TRUE),
    q1 = quantile(observed_over_expected, 0.25, na.rm = TRUE),
    q3 = quantile(observed_over_expected, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

box_stats_all <- box_df_all %>%
  group_by(tool, spike_fraction_label, Target_Condition) %>%
  summarise(
    n = n(),
    median_oe = median(observed_over_expected, na.rm = TRUE),
    q1 = quantile(observed_over_expected, 0.25, na.rm = TRUE),
    q3 = quantile(observed_over_expected, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# y-axis limit for panel C: make shared across tools and robust to extreme outliers
if (is.na(opt$panelc_ymax)) {
  upper <- quantile(box_df$observed_over_expected[is.finite(box_df$observed_over_expected)], 0.98, na.rm = TRUE)
  upper <- max(upper, 1.4)
  panelc_ymax <- upper * 1.08
} else {
  panelc_ymax <- opt$panelc_ymax
}
panelc_ymax <- max(panelc_ymax, 1.2)

if (is.na(opt$panelc_ymax)) {
  upper_all <- quantile(box_df_all$observed_over_expected[is.finite(box_df_all$observed_over_expected)], 0.98, na.rm = TRUE)
  upper_all <- max(upper_all, 1.4)
  panelc_ymax_all <- upper_all * 1.08
} else {
  panelc_ymax_all <- opt$panelc_ymax
}
panelc_ymax_all <- max(panelc_ymax_all, 1.2)

cond_cols <- c("Control" = "#4C78A8", "Adenoma" = "#D9A441", "colorectal carcinoma" = "#D55E00")
class_cols <- c("Good" = "#009E73", "Intermediate" = "#E69F00", "Poor / missed" = "#D55E00")

theme_pub <- function(base_size = 10.6) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "#6F6F6F", fill = NA, linewidth = 0.46),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.25, colour = "#ECECEC"),
      axis.line = element_blank(),
      axis.ticks = element_line(colour = "#4D4D4D", linewidth = 0.30),
      strip.background = element_rect(fill = "#EDEDED", colour = "#6F6F6F", linewidth = 0.50),
      strip.text = element_text(face = "bold", colour = "#222222", margin = margin(4.8, 6, 4.8, 6)),
      legend.position = "bottom",
      legend.justification = "center",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold", size = rel(0.92)),
      legend.text = element_text(size = rel(0.88)),
      legend.key.width = unit(12, "pt"),
      legend.key.height = unit(10, "pt"),
      plot.title = element_text(face = "bold", size = rel(1.14), hjust = 0, margin = margin(b = 4)),
      plot.subtitle = element_text(size = rel(0.90), colour = "#444444", margin = margin(b = 4)),
      axis.title = element_text(face = "bold", colour = "#222222"),
      axis.text = element_text(colour = "#333333"),
      panel.spacing = unit(0.82, "lines"),
      plot.margin = margin(6, 7, 6, 6)
    )
}

dodge_w <- 0.78

make_panel_B <- function(class_tbl, label_tbl, title_text) {
  class_tbl2 <- class_tbl %>%
    mutate(
      Target_Condition = factor(Target_Condition, levels = condition_levels),
      condition_short = factor(
        recode(as.character(Target_Condition),
               "Control" = "Control",
               "Adenoma" = "Adenoma",
               "colorectal carcinoma" = "CRC",
               .default = as.character(Target_Condition)),
        levels = c("Control", "Adenoma", "CRC")
      )
    )
  label_tbl2 <- label_tbl %>%
    mutate(
      Target_Condition = factor(Target_Condition, levels = condition_levels),
      condition_short = factor(
        recode(as.character(Target_Condition),
               "Control" = "Control",
               "Adenoma" = "Adenoma",
               "colorectal carcinoma" = "CRC",
               .default = as.character(Target_Condition)),
        levels = c("Control", "Adenoma", "CRC")
      ),
      show_label = !is.na(label) & pct >= 0.075
    ) %>%
    filter(show_label)
  
  ggplot(class_tbl2, aes(x = condition_short, y = pct, fill = recovery_class)) +
    geom_col(width = 0.82, colour = "white", linewidth = 0.38) +
    geom_text(data = label_tbl2, aes(label = label),
              position = position_stack(vjust = 0.5), size = 2.85,
              colour = "black", lineheight = 0.90, show.legend = FALSE) +
    facet_grid(tool ~ spike_fraction_label) +
    scale_fill_manual(values = class_cols, drop = FALSE, name = "Recovery class") +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.01),
                       breaks = c(0, 0.25, 0.50, 0.75, 1),
                       expand = expansion(mult = c(0, 0.01))) +
    coord_cartesian(clip = "off") +
    labs(
      title = title_text,
      subtitle = paste0("Good ≤ ", scales::percent(opt$good_relative_error, accuracy = 1),
                        "; Intermediate > ", scales::percent(opt$good_relative_error, accuracy = 1),
                        " and ≤ ", scales::percent(opt$avg_relative_error, accuracy = 1),
                        "; Poor / missed > ", scales::percent(opt$avg_relative_error, accuracy = 1),
                        " or undetected"),
      x = NULL,
      y = "Samples"
    ) +
    theme_pub(10.6) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 8.8),
      axis.text.y = element_text(size = 8.9),
      strip.text = element_text(size = 9.0, face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.margin = margin(6, 8, 6, 6)
    )
}

make_panel_B_compact <- function(class_tbl, label_tbl, title_text) {
  class_tbl2 <- class_tbl %>%
    mutate(
      spike_fraction_label = factor(as.character(spike_fraction_label), levels = weak_fraction_levels_desc),
      Target_Condition = factor(Target_Condition, levels = condition_levels),
      condition_short = factor(
        recode(as.character(Target_Condition),
               "Control" = "Control",
               "Adenoma" = "Adenoma",
               "colorectal carcinoma" = "CRC",
               .default = as.character(Target_Condition)),
        levels = c("Control", "Adenoma", "CRC")
      )
    )
  label_tbl2 <- label_tbl %>%
    mutate(
      spike_fraction_label = factor(as.character(spike_fraction_label), levels = weak_fraction_levels_desc),
      Target_Condition = factor(Target_Condition, levels = condition_levels),
      condition_short = factor(
        recode(as.character(Target_Condition),
               "Control" = "Control",
               "Adenoma" = "Adenoma",
               "colorectal carcinoma" = "CRC",
               .default = as.character(Target_Condition)),
        levels = c("Control", "Adenoma", "CRC")
      ),
      show_label = !is.na(label) & pct >= 0.14
    ) %>%
    filter(show_label)
  
  ggplot(class_tbl2, aes(x = condition_short, y = pct, fill = recovery_class)) +
    geom_col(width = 0.78, colour = "white", linewidth = 0.42) +
    geom_text(data = label_tbl2, aes(label = label),
              position = position_stack(vjust = 0.5), size = 2.55,
              colour = "black", lineheight = 0.90, show.legend = FALSE) +
    facet_grid(spike_fraction_label ~ tool) +
    scale_fill_manual(values = class_cols, drop = FALSE, name = "Recovery class") +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.01),
                       breaks = c(0, 0.5, 1),
                       expand = expansion(mult = c(0, 0.01))) +
    coord_cartesian(clip = "off") +
    labs(
      title = title_text,
      subtitle = paste0(
        "Good ≤", scales::percent(opt$good_relative_error, accuracy = 1),
        "\nIntermediate >", scales::percent(opt$good_relative_error, accuracy = 1),
        " and ≤", scales::percent(opt$avg_relative_error, accuracy = 1),
        "\nPoor / missed >", scales::percent(opt$avg_relative_error, accuracy = 1),
        " or undetected"
      ),
      x = NULL,
      y = "Samples"
    ) +
    theme_pub(9.4) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 7.9),
      axis.text.y = element_text(size = 8.1),
      strip.text = element_text(size = 8.4, face = "bold"),
      strip.background = element_rect(fill = "#EDEDED", colour = "#5F5F5F", linewidth = 0.55),
      panel.border = element_rect(colour = "#6F6F6F", fill = NA, linewidth = 0.48),
      panel.spacing.x = unit(0.80, "lines"),
      panel.spacing.y = unit(0.95, "lines"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 8.1, face = "bold"),
      legend.text = element_text(size = 7.8),
      plot.title = element_text(face = "bold", size = 12.5, hjust = 0, margin = margin(b = 3)),
      # Extra subtitle-to-strip spacing aligns the tool facet boxes with the
      # Bias/Variability facet boxes in Fig. 2 panel B after assembly.
      # Use a larger internal gap rather than adding a cowplot spacer,
      # because this keeps the C-panel title aligned with the B-panel title.
      plot.subtitle = element_text(size = 8.2, lineheight = 1.05, colour = "#333333", margin = margin(b = 10)),
      plot.margin = margin(5, 7, 5, 6)
    )
}


make_panel_C <- function(box_tbl, ymax_value, title_text) {
  ggplot(box_tbl, aes(x = spike_fraction_label, y = observed_over_expected, fill = Target_Condition)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.90, ymax = 1.10, fill = "#1B9E77", alpha = 0.09) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.50, ymax = 0.90, fill = "#D9A441", alpha = 0.08) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.10, ymax = 1.50, fill = "#D9A441", alpha = 0.08) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.42, colour = "black") +
    geom_boxplot(
      aes(group = interaction(spike_fraction_label, Target_Condition)),
      position = position_dodge2(width = dodge_w, preserve = "single"),
      width = 0.64,
      outlier.shape = 16,
      outlier.size = 0.9,
      outlier.alpha = 0.22,
      linewidth = 0.30,
      colour = "grey25"
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      aes(group = Target_Condition),
      position = position_dodge(width = dodge_w),
      size = 1.5,
      shape = 21,
      stroke = 0.20,
      colour = "black",
      fill = "white",
      show.legend = FALSE
    ) +
    facet_grid(tool ~ .) +
    scale_fill_manual(values = cond_cols, drop = FALSE, name = "Condition") +
    coord_cartesian(ylim = c(0, ymax_value)) +
    labs(title = title_text, x = "Spike fraction", y = "Observed / expected target mass") +
    theme_pub() +
    theme(legend.position = "bottom")
}

save_png_pdf <- function(plot_obj, stem, width_in, height_in) {
  ggsave(file.path(opt$outdir, paste0(stem, ".png")), plot_obj, width = width_in, height = height_in, dpi = opt$dpi, units = "in")
  ggsave(file.path(opt$outdir, paste0(stem, ".pdf")), plot_obj, width = width_in, height = height_in, units = "in", device = cairo_pdf)
  saveRDS(plot_obj, file.path(opt$outdir, paste0(stem, ".rds")))
}

pA <- ggplot(baseline_prev, aes(x = Target_Condition, y = prevalence, fill = Target_Condition)) +
  geom_col(width = 0.66, colour = "grey25", linewidth = 0.22, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.15, linewidth = 0.28) +
  facet_grid(tool ~ Study) +
  scale_fill_manual(values = cond_cols, drop = FALSE, name = "Condition") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = "A. Baseline detectability",
    x = NULL,
    y = "Positive baseline samples"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1), legend.position = "none")

pB <- make_panel_B(class_sum, class_labels, "B. Recovery quality classes at low spike fractions")
# Revised manuscript Fig. 2 uses recovery quality classes as panel C, replacing
# the former baseline-detectability inset.
pB_revised_C <- make_panel_B_compact(class_sum, class_labels, "C. F. nucleatum recovery quality")
pC <- make_panel_C(box_df, panelc_ymax, "C. Quantitative recovery across all samples")

# Separate supplementary versions across all spike fractions
pB_all <- make_panel_B(class_sum_all, class_labels_all, "Recovery quality classes across all spike fractions")
pC_all <- make_panel_C(box_df_all, panelc_ymax_all, "Quantitative recovery across all spike fractions")

final_plot <- ((pA | pB) / pC) +
  plot_layout(widths = c(1.02, 1.18), heights = c(1.0, 0.96), guides = "collect") +
  plot_annotation(
    title = paste0("Figure 2. ", taxon_title, " illustrates how profiler choice can obscure weak adenoma-associated signal"),
    theme = theme(plot.title = element_text(face = "bold", size = 14.2, hjust = 0), legend.position = "bottom")
  ) &
  theme(legend.position = "bottom", legend.justification = "center", legend.box = "horizontal")

caption_text <- paste0(
  "Figure 2. ", taxon_title, " illustrates how profiler choice can obscure weak adenoma-associated signal. ",
  "(A) Full-cohort prevalence of positive baseline ", taxon_title, " in unspiked samples, stratified by study and profiler. ",
  "(B) Recovery quality classes for the weak-signal regime (spike fractions ≤ ", fmt_fraction(opt$max_weak_fraction), "). Samples were classified as Good when |observed/expected - 1| ≤ ",
  opt$good_relative_error, ", Intermediate when the deviation was > ", opt$good_relative_error, " and ≤ ", opt$avg_relative_error,
  ", and Poor / missed otherwise; non-detections were counted as Poor / missed. ",
  "(C) Boxplots of observed/expected target mass across all samples for each low spike fraction and condition, pooled across studies within each profiler. Box centers and spreads summarize the median and interquartile range directly, enabling visual comparison of underestimation, overestimation, and variability."
)

writeLines(caption_text, con = file.path(opt$outdir, paste0("caption_", selected_label, ".txt")))
write_csv(baseline_prev, file.path(opt$outdir, paste0(selected_label, "_baseline_prevalence_summary.csv")))
write_csv(class_sum, file.path(opt$outdir, paste0(selected_label, "_weak_signal_recovery_classes.csv")))
write_csv(class_sum_all, file.path(opt$outdir, paste0(selected_label, "_all_spikes_recovery_classes.csv")))
write_csv(box_stats, file.path(opt$outdir, paste0(selected_label, "_weak_signal_boxplot_summary.csv")))
write_csv(box_stats_all, file.path(opt$outdir, paste0(selected_label, "_all_spikes_boxplot_summary.csv")))
write_csv(tibble(panelc_ymax = panelc_ymax, panelc_ymax_all = panelc_ymax_all), file.path(opt$outdir, paste0(selected_label, "_panelC_axis.csv")))

# Main combined figure
save_png_pdf(final_plot, paste0("manuscript_case_study_boxplotC_", selected_label), opt$width, opt$height)
# Individual main-text panels
save_png_pdf(pA, paste0("panelA_baseline_detectability_", selected_label), 7.2, 5.3)
save_png_pdf(pB, paste0("panelB_recovery_classes_lowspikes_", selected_label), 11.2, 4.85)
save_png_pdf(pB_revised_C, paste0("panelC_recovery_classes_lowspikes_", selected_label), 4.35, 6.85)
save_png_pdf(pC, paste0("panelC_quant_recovery_lowspikes_", selected_label), 9.2, 5.2)
# Supplementary all-spike versions
save_png_pdf(pB_all, paste0("supplement_panelB_recovery_classes_allspikes_", selected_label), 11.0, 5.8)
save_png_pdf(pC_all, paste0("supplement_panelC_quant_recovery_allspikes_", selected_label), 9.6, 5.8)
save_png_pdf((pB_all / pC_all) + plot_layout(heights = c(1.0, 0.95), guides = "collect") & theme(legend.position = "bottom"),
             paste0("supplement_allspikes_BC_", selected_label), 11.0, 10.0)

message("[OK] Wrote outputs under: ", normalizePath(opt$outdir, mustWork = FALSE))
