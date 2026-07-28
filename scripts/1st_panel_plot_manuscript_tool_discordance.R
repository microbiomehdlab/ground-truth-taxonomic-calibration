#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(scales)
  library(patchwork)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

option_list <- list(
  make_option("--indir", type = "character", default = NULL),
  make_option("--manifest-input", type = "character", default = NULL),
  make_option("--design-input", type = "character", default = NULL),
  make_option("--aliases-input", type = "character", default = NULL),
  make_option("--kraken-input", type = "character", default = NULL),
  make_option("--metaphlan-input", type = "character", default = NULL),
  make_option("--metadata-input", type = "character", default = NULL),
  make_option("--metadata-sample-col", type = "character", default = NULL),
  make_option("--metadata-condition-col", type = "character", default = NULL),
  make_option("--metadata-study-col", type = "character", default = NULL),
  make_option("--outdir", type = "character", default = "RUNS/plots_manuscript_full_baseline"),
  make_option("--spike-labels", type = "character", default = NULL),
  make_option("--prevalence-threshold", type = "double", default = 0),
  make_option("--abundance-scale", type = "character", default = "auto"),
  make_option("--allow-metadata-subset", action = "store_true", default = FALSE),
  make_option("--width", type = "double", default = 14),
  make_option("--height", type = "double", default = 9),
  make_option("--dpi", type = "integer", default = 400)
)
opt <- parse_args(OptionParser(option_list = option_list))

cli_value <- function(flag) {
  args <- commandArgs(trailingOnly = TRUE)
  eq_prefix <- paste0(flag, "=")
  hit <- grep(paste0("^", gsub("-", "\\-", eq_prefix)), args, value = TRUE)
  if (length(hit) > 0) return(sub(eq_prefix, "", hit[[1]], fixed = TRUE))
  idx <- which(args == flag)
  if (length(idx) > 0 && idx[[1]] < length(args)) return(args[[idx[[1]] + 1]])
  NULL
}
opt_value <- function(opt, stem) {
  keys <- unique(c(stem, gsub("-", "_", stem), gsub("-", ".", stem), gsub("_", "-", stem), gsub("_", ".", stem)))
  for (k in keys) {
    val <- opt[[k]]
    if (!is.null(val) && length(val) > 0 && !all(is.na(val))) return(val)
  }
  NULL
}
for (nm in c("indir","manifest-input","design-input","aliases-input","kraken-input","metaphlan-input","metadata-input","metadata-sample-col","metadata-condition-col","metadata-study-col","outdir","spike-labels","abundance-scale")) {
  opt[[gsub("-", "_", nm)]] <- opt_value(opt, nm) %||% cli_value(paste0("--", nm)) %||% opt[[gsub("-", "_", nm)]]
}
opt$prevalence_threshold <- as.numeric(opt_value(opt, "prevalence-threshold") %||% cli_value("--prevalence-threshold") %||% 0)
opt$allow_metadata_subset <- isTRUE(opt_value(opt, "allow-metadata-subset")) || "--allow-metadata-subset" %in% commandArgs(trailingOnly = TRUE)
opt$width <- as.numeric(opt_value(opt, "width") %||% cli_value("--width") %||% 14)
opt$height <- as.numeric(opt_value(opt, "height") %||% cli_value("--height") %||% 9)
opt$dpi <- as.integer(opt_value(opt, "dpi") %||% cli_value("--dpi") %||% 400)

message("[DEBUG] Parsed --spike-labels: ", opt$spike_labels %||% "<NULL>")
message("[DEBUG] Parsed --metadata-input: ", opt$metadata_input %||% "<NULL>")
message("[DEBUG] Parsed --metadata-sample-col: ", opt$metadata_sample_col %||% "<auto>")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing_file <- function(path, label) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) stop(sprintf("Missing %s: %s", label, path %||% "NULL"), call. = FALSE)
}
read_any_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) readr::read_tsv(path, show_col_types = FALSE) else readr::read_csv(path, show_col_types = FALSE)
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

choose_best_sample_col <- function(md, full_sample_ids, explicit_col = NULL) {
  if (!is.null(explicit_col) && nzchar(explicit_col)) {
    if (!explicit_col %in% names(md)) stop("--metadata-sample-col not found. Available columns: ", paste(names(md), collapse = ", "), call. = FALSE)
    return(explicit_col)
  }
  full_keys <- unique(sample_key(full_sample_ids))
  candidate_cols <- names(md)[vapply(md, function(z) is.character(z) || is.factor(z) || is.numeric(z), logical(1))]
  if (!length(candidate_cols)) candidate_cols <- names(md)
  scores <- bind_rows(lapply(candidate_cols, function(cc) {
    vals <- as.character(md[[cc]])
    keys <- unique(sample_key(vals[!is.na(vals) & nzchar(vals)]))
    overlap <- length(intersect(keys, full_keys))
    tibble(column = cc, overlap = overlap, n_values = length(keys), frac_of_full = overlap / length(full_keys))
  })) %>% arrange(desc(overlap), desc(frac_of_full))
  write_csv(scores, file.path(opt$outdir, "metadata_sample_column_match_scores.csv"))
  best <- scores$column[[1]]
  message("[INFO] Auto-selected metadata sample column: ", best, " (", scores$overlap[[1]], "/", length(full_keys), " abundance sample IDs matched after normalization)")
  best
}

read_abundance_wide <- function(path, abundance_scale = "auto") {
  stop_if_missing_file(path, "abundance table")
  dat <- read_any_table(path)
  sid <- first_existing(names(dat), c("sample_id", "sample", "Sample", "SampleID", "ID", "base_id"))
  if (is.na(sid)) sid <- names(dat)[[1]]
  dat <- dat %>% rename(sample_id = all_of(sid))
  dat$sample_id <- as.character(dat$sample_id)
  taxa <- setdiff(names(dat), "sample_id")
  dat[taxa] <- lapply(dat[taxa], function(x) {x <- suppressWarnings(as.numeric(x)); x[is.na(x)] <- 0; x})
  mat <- as.matrix(dat[taxa]); storage.mode(mat) <- "double"
  finite_vals <- as.numeric(mat[is.finite(mat)]); row_sums <- rowSums(mat, na.rm = TRUE)
  max_value <- if (length(finite_vals)) max(finite_vals, na.rm = TRUE) else NA_real_
  med_sum <- if (length(row_sums)) median(row_sums, na.rm = TRUE) else NA_real_
  inferred <- if (identical(abundance_scale, "auto")) {
    if ((is.finite(max_value) && max_value > 1.0000001) || (is.finite(med_sum) && med_sum > 1.5)) "percent" else "fraction"
  } else abundance_scale
  if (identical(inferred, "percent")) dat[taxa] <- lapply(dat[taxa], function(x) x / 100)
  message(sprintf("[INFO] %s: %d samples, %d taxa, scale=%s (max=%.4g, median row sum=%.4g)", basename(path), nrow(dat), length(taxa), inferred, max_value, med_sum))
  dat
}

indir <- opt$indir
manifest_path <- opt$manifest_input %||% if (!is.null(indir)) file.path(indir, "run_manifest.resolved.csv") else NULL
design_path <- opt$design_input %||% if (!is.null(indir)) file.path(indir, "spike_design_long.csv") else NULL
aliases_path <- opt$aliases_input %||% if (!is.null(indir)) file.path(indir, "taxon_aliases.resolved.csv") else NULL
stop_if_missing_file(design_path, "spike design")
design <- read_any_table(design_path)
manifest <- if (!is.null(manifest_path) && file.exists(manifest_path)) read_any_table(manifest_path) else NULL

if (is.null(opt$kraken_input) || is.null(opt$metaphlan_input)) {
  if (is.null(manifest) || !all(c("tool", "original_table") %in% names(manifest))) stop("Provide --kraken-input and --metaphlan-input, or --indir with run_manifest.resolved.csv.", call. = FALSE)
  originals <- manifest %>% distinct(tool, original_table) %>% filter(!is.na(original_table), nzchar(original_table))
  if (is.null(opt$kraken_input)) opt$kraken_input <- originals %>% filter(str_detect(tool, regex("kraken|bracken", ignore_case = TRUE))) %>% pull(original_table) %>% unique() %>% .[[1]]
  if (is.null(opt$metaphlan_input)) opt$metaphlan_input <- originals %>% filter(str_detect(tool, regex("metaphlan", ignore_case = TRUE))) %>% pull(original_table) %>% unique() %>% .[[1]]
}
stop_if_missing_file(opt$kraken_input, "Kraken2+Bracken original/full baseline table")
stop_if_missing_file(opt$metaphlan_input, "MetaPhlAn4 original/full baseline table")
message("[INFO] Kraken2+Bracken full baseline: ", opt$kraken_input)
message("[INFO] MetaPhlAn4 full baseline: ", opt$metaphlan_input)

if (is.null(opt$spike_labels) || !nzchar(opt$spike_labels)) {
  selected_labels <- design %>% distinct(spike_label) %>% pull(spike_label) %>% as.character()
  message("[INFO] No --spike-labels provided; using all labels in design (", length(selected_labels), ").")
} else {
  selected_labels <- str_split(opt$spike_labels, ",", simplify = FALSE)[[1]] %>% str_trim()
  selected_labels <- selected_labels[nzchar(selected_labels)]
  message("[INFO] Requested labels before sorting (", length(selected_labels), "): ", paste(selected_labels, collapse = ", "))
}
selected_labels <- sort(unique(selected_labels))
message("[INFO] Final spike-label order used in the figure (alphabetical): ", paste(selected_labels, collapse = ", "))

label_taxa <- design %>% filter(.data$spike_label %in% selected_labels) %>% distinct(spike_label, member_taxon) %>% filter(!is.na(member_taxon), nzchar(member_taxon))
missing_design <- setdiff(selected_labels, unique(label_taxa$spike_label))
if (length(missing_design)) stop("Requested labels not found in spike_design_long.csv: ", paste(missing_design, collapse = ", "), call. = FALSE)

aliases <- tibble(canonical = character(), alias = character(), tool = character())
if (!is.null(aliases_path) && file.exists(aliases_path)) {
  aliases <- read_any_table(aliases_path) %>% transmute(canonical = as.character(.data$canonical), alias = as.character(.data$alias), tool = as.character(.data$tool)) %>% filter(!is.na(canonical), !is.na(alias))
}
# Return a species-level fallback from a possibly subspecies/strain-level name.
# Example: "Fusobacterium nucleatum subsp. vincentii" -> "Fusobacterium nucleatum".
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

make_candidate_map <- function(label_taxa, tool_name, aliases) {
  base <- label_taxa %>%
    transmute(spike_label, canonical = member_taxon, candidate = member_taxon)
  
  # Add binomial/species-level fallbacks. This is crucial when the spike design uses
  # a subspecies/strain but the full baseline table reports the parent species.
  fallback <- label_taxa %>%
    transmute(spike_label, canonical = member_taxon, candidate = species_fallback(member_taxon))
  
  al <- tibble(spike_label = character(), canonical = character(), candidate = character())
  if (nrow(aliases) > 0) {
    alias_tool_regex <- ifelse(str_detect(tool_name, "Kraken"), "kraken|bracken", "metaphlan")
    
    # IMPORTANT: aliases with missing/blank tool are generic and should be allowed
    # for both tools. The previous version dropped these rows, which made Fnuc fail
    # because taxon_aliases.resolved.csv has tool = NA for Fnuc.
    al_exact <- aliases %>%
      mutate(tool_chr = coalesce(as.character(tool), "")) %>%
      filter(tool_chr == "" | str_detect(tool_chr, regex(alias_tool_regex, ignore_case = TRUE))) %>%
      inner_join(label_taxa %>% distinct(spike_label, member_taxon),
                 by = c("canonical" = "member_taxon"),
                 relationship = "many-to-many") %>%
      transmute(spike_label, canonical, candidate = alias)
    
    # Also allow aliases whose canonical collapses to the same species-level binomial
    # as the design member_taxon. This catches canonical subspecies/strain strings
    # that do not exactly equal the design value.
    al_species <- aliases %>%
      mutate(tool_chr = coalesce(as.character(tool), ""),
             canonical_species = species_fallback(canonical)) %>%
      filter(tool_chr == "" | str_detect(tool_chr, regex(alias_tool_regex, ignore_case = TRUE))) %>%
      inner_join(label_taxa %>%
                   distinct(spike_label, member_taxon) %>%
                   mutate(canonical_species = species_fallback(member_taxon)),
                 by = "canonical_species",
                 relationship = "many-to-many") %>%
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
  matches <- candidates %>% inner_join(col_map, by = c("candidate_key" = "column_key"), relationship = "many-to-many") %>% distinct(spike_label, canonical, candidate, column)
  rows <- vector("list", length(selected_labels))
  for (i in seq_along(selected_labels)) {
    lab <- selected_labels[[i]]
    cols <- matches %>% filter(spike_label == lab) %>% pull(column) %>% unique()
    abund <- if (length(cols) == 0) rep(0, nrow(wide)) else rowSums(wide[, cols, drop = FALSE], na.rm = TRUE)
    rows[[i]] <- tibble(sample_id = wide$sample_id, tool = tool_display, spike_label = lab, abundance = abund)
  }
  list(data = bind_rows(rows), matches = matches %>% mutate(tool = tool_display))
}

kraken_wide <- read_abundance_wide(opt$kraken_input, opt$abundance_scale)
metaphlan_wide <- read_abundance_wide(opt$metaphlan_input, opt$abundance_scale)
kr <- extract_label_abundances(kraken_wide, "Kraken2 + Bracken", label_taxa, aliases)
mp <- extract_label_abundances(metaphlan_wide, "MetaPhlAn 4", label_taxa, aliases)
full_baseline <- bind_rows(kr$data, mp$data) %>% mutate(sample_key = sample_key(sample_id))
all_matches <- bind_rows(kr$matches, mp$matches)
write_csv(all_matches, file.path(opt$outdir, "full_baseline_taxon_column_matches.csv"))
unmatched <- expand.grid(spike_label = selected_labels, tool = c("Kraken2 + Bracken", "MetaPhlAn 4"), stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  anti_join(all_matches %>% distinct(spike_label, tool), by = c("spike_label", "tool"))
write_csv(unmatched, file.path(opt$outdir, "full_baseline_unmatched_requested_labels.csv"))
if (nrow(unmatched) > 0) {
  message("[WARN] Some requested labels had no exact/alias/species-level column match. See full_baseline_unmatched_requested_labels.csv")
}

read_metadata <- function(path, full_sample_ids) {
  md <- read_any_table(path)
  sid <- choose_best_sample_col(md, full_sample_ids, opt$metadata_sample_col)
  cond <- if (!is.null(opt$metadata_condition_col) && nzchar(opt$metadata_condition_col)) opt$metadata_condition_col else first_existing(names(md), c("Target_Condition", "condition", "Condition", "group", "Group", "diagnosis", "Diagnosis", "disease", "Disease", "status", "Status"))
  study <- if (!is.null(opt$metadata_study_col) && nzchar(opt$metadata_study_col)) opt$metadata_study_col else first_existing(names(md), c("Study", "study", "Dataset", "dataset", "Cohort", "cohort", "original_id", "dataset_name", "source"))
  if (is.na(cond) || !cond %in% names(md)) stop("Metadata must contain a condition column. Use --metadata-condition-col. Available columns: ", paste(names(md), collapse = ", "), call. = FALSE)
  if (is.na(study) || !study %in% names(md)) stop("Metadata must contain a Study/Dataset/Cohort column. Use --metadata-study-col. Available columns: ", paste(names(md), collapse = ", "), call. = FALSE)
  out <- md %>% transmute(sample_id_metadata_original = as.character(.data[[sid]]), sample_key = sample_key(.data[[sid]]), Target_Condition = as.character(.data[[cond]]), Study = as.character(.data[[study]])) %>%
    filter(!is.na(sample_key), nzchar(sample_key), !is.na(Target_Condition), nzchar(Target_Condition), !is.na(Study), nzchar(Study)) %>% distinct(sample_key, Target_Condition, Study, .keep_all = TRUE)
  write_csv(head(out, 50), file.path(opt$outdir, "metadata_normalized_preview.csv"))
  out
}

if (!is.null(opt$metadata_input) && nzchar(opt$metadata_input)) {
  stop_if_missing_file(opt$metadata_input, "full sample metadata")
  metadata <- read_metadata(opt$metadata_input, full_baseline %>% distinct(sample_id) %>% pull(sample_id))
  metadata_source <- opt$metadata_input
} else {
  fallback_md <- if (!is.null(indir)) file.path(indir, "baseline_with_condition.csv") else NULL
  if (!is.null(fallback_md) && file.exists(fallback_md)) {
    tmp <- read_any_table(fallback_md)
    sid <- first_existing(names(tmp), c("base_id", "sample_id", "sample", "Sample", "SampleID"))
    metadata <- tmp %>% transmute(sample_id_metadata_original = as.character(.data[[sid]]), sample_key = sample_key(sample_id_metadata_original), Target_Condition = as.character(.data$Target_Condition), Study = as.character(.data$Study)) %>%
      filter(!is.na(sample_key), nzchar(sample_key), !is.na(Target_Condition), !is.na(Study)) %>% distinct(sample_key, Target_Condition, Study, .keep_all = TRUE)
    metadata_source <- fallback_md
  } else stop("No --metadata-input supplied and no fallback metadata found.", call. = FALSE)
  message("[WARN] No --metadata-input supplied. Falling back to metadata from ", metadata_source)
}

all_samples <- full_baseline %>% distinct(tool, sample_id, sample_key)
coverage_join <- all_samples %>% left_join(metadata, by = "sample_key")
coverage <- coverage_join %>% summarise(total_tool_samples = n(), with_metadata = sum(!is.na(Target_Condition) & !is.na(Study)), metadata_fraction = with_metadata / total_tool_samples)
write_csv(coverage, file.path(opt$outdir, "full_baseline_metadata_coverage.csv"))
write_csv(coverage_join %>% filter(is.na(Target_Condition) | is.na(Study)) %>% head(100), file.path(opt$outdir, "unmatched_abundance_sample_ids_first100.csv"))
message(sprintf("[INFO] Metadata coverage: %d/%d tool-samples (%.1f%%)", coverage$with_metadata, coverage$total_tool_samples, 100 * coverage$metadata_fraction))
if (!isTRUE(opt$allow_metadata_subset) && coverage$metadata_fraction < 0.95) {
  stop(paste0("Metadata cover only ", round(100 * coverage$metadata_fraction, 1), "% of full abundance table tool-samples. Check metadata_sample_column_match_scores.csv and unmatched_abundance_sample_ids_first100.csv, or provide --metadata-sample-col / --metadata-condition-col / --metadata-study-col."), call. = FALSE)
}

plot_df <- full_baseline %>% left_join(metadata, by = "sample_key") %>% filter(!is.na(Target_Condition), !is.na(Study)) %>% mutate(
  spike_label = factor(spike_label, levels = selected_labels),
  Target_Condition = dplyr::recode(
    as.character(Target_Condition),
    "colorectal carcinoma" = "CRC",
    "Colorectal carcinoma" = "CRC",
    "CRC" = "CRC",
    .default = as.character(Target_Condition)
  ),
  Target_Condition = factor(Target_Condition, levels = c("Control", "Adenoma", "CRC")),
  positive = abundance > opt$prevalence_threshold
) %>% filter(!is.na(spike_label), !is.na(Target_Condition))
if (nrow(plot_df) == 0) stop("No rows left after joining metadata and filtering labels/conditions.", call. = FALSE)

wilson <- function(x, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2); p <- x / n; den <- 1 + z^2 / n
  centre <- p + z^2 / (2 * n); half <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n)
  tibble(ymin = pmax(0, (centre - half) / den), ymax = pmin(1, (centre + half) / den))
}
prev <- plot_df %>% group_by(Study, spike_label, Target_Condition, tool) %>% summarise(n = n_distinct(sample_id), positive_n = sum(positive, na.rm = TRUE), prevalence = positive_n / n, .groups = "drop") %>% rowwise() %>% mutate(ci = list(wilson(positive_n, n))) %>% unnest(ci) %>% ungroup()
abund <- plot_df %>% filter(positive)
write_csv(prev, file.path(opt$outdir, "full_baseline_prevalence_summary.csv"))
write_csv(abund, file.path(opt$outdir, "full_baseline_positive_abundance_long.csv"))

tool_cols <- c("Kraken2 + Bracken" = "#2A9D8F", "MetaPhlAn 4" = "#7B6DCC")
abundance_percent_labels <- function(x) {
  # Data remain proportions and the scale remains logarithmic. Only the
  # displayed labels are converted to percent units.
  out <- trimws(formatC(
    as.numeric(x) * 100,
    format = "fg",
    digits = 6,
    drop0trailing = TRUE
  ))
  out[is.na(x)] <- NA_character_
  paste0(out, "%")
}
theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size) + theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey88"),
    panel.border = element_rect(linewidth = 0.35, colour = "grey35"), strip.background = element_rect(fill = "grey96", colour = "grey80", linewidth = 0.3),
    strip.text = element_text(face = "bold"), legend.position = "top", legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = rel(1.15)), plot.subtitle = element_text(size = rel(0.95)), panel.spacing = unit(0.8, "lines")
  )
}
pA <- ggplot(prev, aes(x = Target_Condition, y = prevalence, fill = tool)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = "grey25", linewidth = 0.2, alpha = 0.95) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax, colour = tool), position = position_dodge(width = 0.72), width = 0.18, linewidth = 0.35) +
  facet_grid(Study ~ spike_label) + scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.04))) +
  scale_fill_manual(values = tool_cols, drop = FALSE) + scale_colour_manual(values = tool_cols, drop = FALSE) +
  labs(title = "A. Full-cohort baseline prevalence by condition", subtitle = paste0("Bars show fraction of samples with baseline abundance > ", opt$prevalence_threshold, "; error bars are 95% Wilson CIs."), x = NULL, y = "Positive baseline samples") + theme_pub()
pB <- ggplot(abund, aes(x = Target_Condition, y = abundance, colour = tool, fill = tool)) +
  geom_boxplot(aes(group = interaction(Target_Condition, tool)), position = position_dodge(width = 0.72), width = 0.55, outlier.shape = NA, linewidth = 0.3, alpha = 0.22) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.72), size = 0.8, alpha = 0.65, stroke = 0) +
  facet_grid(Study ~ spike_label, scales = "free_y") +
  scale_y_log10(labels = abundance_percent_labels) +
  scale_fill_manual(values = tool_cols, drop = FALSE) + scale_colour_manual(values = tool_cols, drop = FALSE) +
  labs(title = "B. Full-cohort baseline abundance among positive samples", subtitle = paste0("Only samples with baseline abundance > ", opt$prevalence_threshold, " are shown here."), x = NULL, y = "Baseline abundance") + theme_pub() + theme(legend.position = "none")
final <- pA / pB + plot_layout(heights = c(1, 1.15))

ggsave(file.path(opt$outdir, "manuscript_full_baseline_discordance_panel.png"), final, width = opt$width, height = opt$height, dpi = opt$dpi, units = "in")
ggsave(file.path(opt$outdir, "manuscript_full_baseline_discordance_panel.pdf"), final, width = opt$width, height = opt$height, units = "in", device = cairo_pdf)
message("[OK] Final spike labels in figure (", length(selected_labels), "): ", paste(selected_labels, collapse = ", "))
message("[OK] Wrote full-baseline manuscript panel under ", normalizePath(opt$outdir, mustWork = FALSE))
