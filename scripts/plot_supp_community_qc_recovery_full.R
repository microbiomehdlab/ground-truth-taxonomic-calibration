#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1]]
}

run_root <- get_arg("--run-root", "RUNS_publication_original_unpaired_q010_cluster")

c8_outdir <- get_arg(
  "--c8-outdir",
  file.path(
    run_root,
    "revised_manuscript_figures",
    "source_panels",
    "FigC8_community_depth_qc_recovery"
  )
)

input_file <- get_arg(
  "--community-target-file",
  file.path(c8_outdir, "community_target_level_depth_recovery.tsv")
)

outdir <- get_arg(
  "--outdir",
  file.path(
    run_root,
    "revised_manuscript_figures",
    "source_panels",
    "community_qc_recovery_FULL_reproducible"
  )
)
dpi <- as.integer(get_arg("--dpi", "450"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)

message("[INFO] Input: ", input_file)
message("[INFO] Outdir: ", outdir)

if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file, "\nRun the community FigC8 script first.")
}

x <- read_tsv(input_file, show_col_types = FALSE)

required <- c(
  "sample_id",
  "Study",
  "Target_Condition",
  "tool_label",
  "target_taxon",
  "effective_fraction_per_species",
  "implanted_read_pairs_target",
  "log10_implanted",
  "relative_error",
  "observed_over_expected",
  "good_recovery",
  "raw_read_pairs_million",
  "postqc_read_pairs_million",
  "retention_pct"
)

missing <- setdiff(required, names(x))
if (length(missing) > 0) {
  stop("Input table is missing required columns: ", paste(missing, collapse = ", "))
}

x <- x %>%
  mutate(
    Target_Condition = dplyr::recode(
      as.character(Target_Condition),
      "colorectal carcinoma" = "CRC",
      "Colorectal carcinoma" = "CRC",
      .default = as.character(Target_Condition)
    ),
    Target_Condition = factor(Target_Condition, levels = c("Control", "Adenoma", "CRC"))
  )

canonical_taxon <- function(z) {
  z <- as.character(z)
  z_low <- tolower(z)

  case_when(
    str_detect(z_low, "bacteroides.*fragilis") ~ "Bacteroides fragilis",
    str_detect(z_low, "clostridium.*symbiosum") ~ "Clostridium symbiosum",
    str_detect(z_low, "dialister.*pneumosintes|allisonella.*pneumosintes") ~ "Dialister pneumosintes",
    str_detect(z_low, "fusobacterium.*nucleatum") ~ "Fusobacterium nucleatum",
    str_detect(z_low, "hungatella.*hathewayi") ~ "Hungatella hathewayi",
    str_detect(z_low, "parvimonas.*micra") ~ "Parvimonas micra",
    str_detect(z_low, "peptostreptococcus.*anaerobius") ~ "Peptostreptococcus anaerobius",
    str_detect(z_low, "peptostreptococcus.*stomatis") ~ "Peptostreptococcus stomatis",
    str_detect(z_low, "porphyromonas.*asaccharolytica") ~ "Porphyromonas asaccharolytica",
    str_detect(z_low, "prevotella.*intermedia") ~ "Prevotella intermedia",
    TRUE ~ NA_character_
  )
}

taxon_levels <- c(
  "Bacteroides fragilis",
  "Clostridium symbiosum",
  "Dialister pneumosintes",
  "Fusobacterium nucleatum",
  "Hungatella hathewayi",
  "Parvimonas micra",
  "Peptostreptococcus anaerobius",
  "Peptostreptococcus stomatis",
  "Porphyromonas asaccharolytica",
  "Prevotella intermedia"
)
taxon_codes <- c(
  "Bacteroides fragilis" = "Bfrag",
  "Clostridium symbiosum" = "Csym",
  "Dialister pneumosintes" = "Dpne",
  "Fusobacterium nucleatum" = "Fnuc",
  "Hungatella hathewayi" = "Hhat",
  "Parvimonas micra" = "Pmic",
  "Peptostreptococcus anaerobius" = "Pana",
  "Peptostreptococcus stomatis" = "Psto",
  "Porphyromonas asaccharolytica" = "Porp",
  "Prevotella intermedia" = "Pint"
)

recovery_class_levels <- c(
  "Good",
  "Intermediate",
  "Poor / missed"
)

recovery_class_cols <- c(
  "Good" = "#009E73",
  "Intermediate" = "#E69F00",
  "Poor / missed" = "#D55E00"
)

profiler_cols <- c(
  "Kraken2 + Bracken" = "#009E73",
  "MetaPhlAn 4" = "#6F5BD3"
)

save_pdf <- function(filename, plot, width, height) {
  if (capabilities("cairo")) {
    ggsave(filename, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  } else {
    ggsave(filename, plot, width = width, height = height, bg = "white")
  }
}

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "grey20"),
      plot.title = element_text(face = "bold", hjust = 0, colour = "grey10"),
      plot.subtitle = element_text(size = rel(0.92), colour = "grey25"),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

pct_label_from_log <- function(v) {
  vals <- 10^v
  paste0(format(signif(vals, 3), trim = TRUE, scientific = FALSE), "%")
}

# ============================================================
# 1. Clean canonical target-level table
# ============================================================

target_level_raw <- x %>%
  mutate(
    effective_fraction_per_species = as.numeric(effective_fraction_per_species),
    effective_fraction_pct = 100 * effective_fraction_per_species,
    log10_effective_fraction_pct = log10(effective_fraction_pct),
    implanted_read_pairs_target = as.numeric(implanted_read_pairs_target),
    log10_implanted = as.numeric(log10_implanted),
    relative_error = suppressWarnings(as.numeric(relative_error)),
    observed_over_expected = suppressWarnings(as.numeric(observed_over_expected)),
    canonical_taxon = canonical_taxon(target_taxon),
    taxon_label = factor(
      unname(taxon_codes[canonical_taxon]),
      levels = unname(taxon_codes[taxon_levels])
    ),
    tool_label = factor(tool_label, levels = c("Kraken2 + Bracken", "MetaPhlAn 4")),
    Study = factor(dplyr::recode(
      as.character(Study),
      "FengQ_2015" = "FengQ 2015",
      "ZellerG_2014" = "Zeller 2014",
      "Zeller_2014" = "Zeller 2014",
      .default = as.character(Study)
    )),
    Target_Condition = factor(Target_Condition)
  ) %>%
  filter(!is.na(canonical_taxon)) %>%
  filter(!is.na(effective_fraction_pct), effective_fraction_pct > 0) %>%
  filter(!is.na(implanted_read_pairs_target), implanted_read_pairs_target > 0)

alias_diagnostics <- target_level_raw %>%
  count(target_taxon, canonical_taxon, sort = TRUE)

write_tsv(alias_diagnostics, file.path(outdir, "tables", "raw_to_canonical_taxon_diagnostics.tsv"))

target_level <- target_level_raw %>%
  group_by(
    sample_id,
    Study,
    Target_Condition,
    tool_label,
    effective_fraction_per_species,
    effective_fraction_pct,
    log10_effective_fraction_pct,
    canonical_taxon
  ) %>%
  summarise(
    implanted_read_pairs_target = median(implanted_read_pairs_target, na.rm = TRUE),
    log10_implanted = median(log10_implanted, na.rm = TRUE),
    relative_error = median(relative_error, na.rm = TRUE),
    observed_over_expected = median(observed_over_expected, na.rm = TRUE),
    raw_read_pairs_million = first(raw_read_pairs_million),
    postqc_read_pairs_million = first(postqc_read_pairs_million),
    retention_pct = first(retention_pct),
    .groups = "drop"
  ) %>%
  mutate(
    relative_error = ifelse(is.infinite(relative_error), NA_real_, relative_error),
    observed_over_expected = ifelse(is.infinite(observed_over_expected), NA_real_, observed_over_expected),
    recovery_class = case_when(
      is.na(observed_over_expected) | is.na(relative_error) ~ "Poor / missed",
      observed_over_expected <= 0 ~ "Poor / missed",
      abs(relative_error) <= 0.10 ~ "Good",
      abs(relative_error) <= 0.50 ~ "Intermediate",
      TRUE ~ "Poor / missed"
    ),
    recovery_class = factor(recovery_class, levels = recovery_class_levels),
    is_good = recovery_class == "Good",
    is_average = recovery_class == "Intermediate",
    is_poor_missed = recovery_class == "Poor / missed",
    canonical_taxon = factor(canonical_taxon, levels = taxon_levels)
  )

write_tsv(target_level, file.path(outdir, "tables", "community_target_level_clean_canonical.tsv"))

target_counts <- target_level %>%
  count(canonical_taxon, sort = TRUE)

write_tsv(target_counts, file.path(outdir, "tables", "canonical_taxon_counts.tsv"))

message("[INFO] Target-level rows after canonical filtering: ", nrow(target_level))
message("[INFO] Unique canonical taxa: ", n_distinct(target_level$canonical_taxon))
print(target_counts, n = Inf)

fraction_breaks <- sort(unique(target_level$effective_fraction_pct))
fraction_breaks_log <- log10(fraction_breaks)

# ============================================================
# 2. Sample-level recovery-class proportions
# ============================================================

sample_keys <- target_level %>%
  distinct(
    sample_id,
    Study,
    Target_Condition,
    tool_label,
    effective_fraction_per_species,
    effective_fraction_pct,
    log10_effective_fraction_pct,
    implanted_read_pairs_target,
    log10_implanted
  )

sample_class_counts <- target_level %>%
  count(
    sample_id,
    Study,
    Target_Condition,
    tool_label,
    effective_fraction_per_species,
    effective_fraction_pct,
    log10_effective_fraction_pct,
    implanted_read_pairs_target,
    log10_implanted,
    recovery_class,
    name = "n_class"
  )

sample_class <- tidyr::expand_grid(
  sample_keys,
  recovery_class = factor(recovery_class_levels, levels = recovery_class_levels)
) %>%
  left_join(
    sample_class_counts,
    by = c(
      "sample_id",
      "Study",
      "Target_Condition",
      "tool_label",
      "effective_fraction_per_species",
      "effective_fraction_pct",
      "log10_effective_fraction_pct",
      "implanted_read_pairs_target",
      "log10_implanted",
      "recovery_class"
    )
  ) %>%
  mutate(n_class = ifelse(is.na(n_class), 0L, n_class)) %>%
  group_by(
    sample_id,
    Study,
    Target_Condition,
    tool_label,
    effective_fraction_per_species,
    effective_fraction_pct,
    log10_effective_fraction_pct,
    implanted_read_pairs_target,
    log10_implanted
  ) %>%
  mutate(
    n_targets = sum(n_class),
    class_prop = ifelse(n_targets > 0, n_class / n_targets, NA_real_),
    class_pct = 100 * class_prop
  ) %>%
  ungroup() %>%
  filter(!is.na(class_prop), n_targets > 0)

write_tsv(sample_class, file.path(outdir, "tables", "community_sample_level_recovery_class_proportions.tsv"))

# ============================================================
# B8. QC only: A/B/C
# ============================================================

qc_plot <- target_level %>%
  distinct(
    sample_id,
    Study,
    Target_Condition,
    raw_read_pairs_million,
    postqc_read_pairs_million,
    retention_pct
  ) %>%
  filter(!is.na(postqc_read_pairs_million), !is.na(raw_read_pairs_million))

pA <- ggplot(qc_plot, aes(x = Target_Condition, y = postqc_read_pairs_million)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.35) +
  geom_jitter(width = 0.15, alpha = 0.45, size = 0.9) +
  facet_wrap(~ Study, scales = "free_x") +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  labs(
    title = "A. Post-QC read-pair depth",
    x = NULL,
    y = "Read pairs after QC (millions)"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

pB <- ggplot(qc_plot, aes(x = Target_Condition, y = retention_pct)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.35) +
  geom_jitter(width = 0.15, alpha = 0.45, size = 0.9) +
  facet_wrap(~ Study, scales = "free_x") +
  scale_y_continuous(labels = function(v) paste0(round(v, 1), "%")) +
  labs(
    title = "B. Read retention after preprocessing",
    x = NULL,
    y = "Post-QC / raw read pairs (%)"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

pC <- ggplot(qc_plot, aes(x = raw_read_pairs_million, y = postqc_read_pairs_million)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.4) +
  geom_point(aes(shape = Study), alpha = 0.55, size = 1.5) +
  labs(
    title = "C. Raw versus post-QC depth",
    x = "Raw read pairs before QC (millions)",
    y = "Read pairs after QC (millions)",
    shape = "Cohort"
  ) +
  theme_pub()

fig_b8 <- (pA | pB) / pC +
  plot_layout(heights = c(1, 1.05))

save_pdf(file.path(outdir, "FigB8_QC_only_ABC.pdf"), fig_b8, width = 12.5, height = 9.2)
ggsave(file.path(outdir, "FigB8_QC_only_ABC.png"), fig_b8, width = 12.5, height = 9.2, dpi = dpi, bg = "white")

# ============================================================
# B9. Recovery classes vs effective spike fraction
# ============================================================

plot_class_fraction <- sample_class %>%
  mutate(
    panel_label = case_when(
      recovery_class == "Good" ~ "A. Good recovery",
      recovery_class == "Intermediate" ~ "B. Intermediate recovery",
      recovery_class == "Poor / missed" ~ "C. Poor / missed recovery",
      TRUE ~ as.character(recovery_class)
    ),
    panel_label = factor(
      panel_label,
      levels = c("A. Good recovery", "B. Intermediate recovery", "C. Poor / missed recovery")
    )
  )

fig_b9 <- ggplot(
  plot_class_fraction,
  aes(
    x = log10_effective_fraction_pct,
    y = class_prop,
    colour = tool_label,
    shape = Study,
    weight = n_targets
  )
) +
  geom_point(alpha = 0.42, size = 1.35, position = position_jitter(width = 0.015, height = 0)) +
  geom_smooth(
    aes(fill = tool_label),
    method = "glm",
    method.args = list(family = quasibinomial()),
    se = TRUE,
    alpha = 0.10,
    linewidth = 0.8
  ) +
  facet_wrap(~ panel_label, ncol = 3) +
  scale_x_continuous(
    breaks = fraction_breaks_log,
    labels = pct_label_from_log
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1)
  ) +
  scale_colour_manual(values = profiler_cols, drop = FALSE) +
  scale_fill_manual(values = profiler_cols, drop = FALSE, guide = "none") +
  labs(
    title = "Community recovery classes across effective spike fractions",
    subtitle = "Good ≤10%; Intermediate >10% and ≤50%; Poor / missed >50% or undetected",
    x = "Effective spike fraction per community member",
    y = "Fraction of community members",
    colour = "Profiler",
    shape = "Cohort"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_pdf(file.path(outdir, "FigB9_recovery_classes_vs_spike_fraction.pdf"), fig_b9, width = 12, height = 5.8)
ggsave(file.path(outdir, "FigB9_recovery_classes_vs_spike_fraction.png"), fig_b9, width = 12, height = 5.8, dpi = dpi, bg = "white")

# ============================================================
# B10. Recovery-class composition by low/medium/high read support
# ============================================================

target_support <- target_level %>%
  mutate(
    support_bin_num = ntile(log10_implanted, 3),
    support_bin = factor(
      support_bin_num,
      levels = c(1, 2, 3),
      labels = c("Low read support", "Medium read support", "High read support")
    )
  )

binned <- target_support %>%
  count(tool_label, support_bin, recovery_class, name = "n") %>%
  group_by(tool_label, support_bin) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    recovery_class = factor(recovery_class, levels = recovery_class_levels),
    label_code = dplyr::recode(
      as.character(recovery_class),
      "Good" = "G",
      "Intermediate" = "I",
      "Poor / missed" = "P/M"
    ),
    label_colour = if_else(recovery_class == "Intermediate", "#1A1A1A", "white")
  )

write_tsv(binned, file.path(outdir, "tables", "recovery_class_composition_by_read_support.tsv"))

support_ranges <- target_support %>%
  group_by(support_bin) %>%
  summarise(
    min_implanted_pairs = min(implanted_read_pairs_target, na.rm = TRUE),
    median_implanted_pairs = median(implanted_read_pairs_target, na.rm = TRUE),
    max_implanted_pairs = max(implanted_read_pairs_target, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(support_ranges, file.path(outdir, "tables", "read_support_bin_ranges.tsv"))

fig_b10 <- ggplot(
  binned,
  aes(x = support_bin, y = prop, fill = recovery_class)
) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.3) +
  geom_text(
    aes(
      label = ifelse(
        prop >= 0.08,
        paste(label_code, percent(prop, accuracy = 1)),
        ""
      ),
      colour = label_colour
    ),
    position = position_stack(vjust = 0.5),
    size = 3.2
  ) +
  facet_wrap(~ tool_label, ncol = 1) +
  scale_fill_manual(
    values = recovery_class_cols,
    breaks = recovery_class_levels,
    drop = FALSE
  ) +
  scale_colour_identity() +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Recovery-class composition by implanted read-support tertile",
    subtitle = paste0(
      "Good (G) ≤10%; Intermediate (I) >10% and ≤50%; ",
      "Poor / missed (P/M) >50% or undetected"
    ),
    x = "Implanted read-support tertile",
    y = "Target-level observations (%)",
    fill = "Recovery class"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

save_pdf(file.path(outdir, "FigB10_recovery_class_composition_by_read_support.pdf"), fig_b10, width = 8.8, height = 7.4)
ggsave(file.path(outdir, "FigB10_recovery_class_composition_by_read_support.png"), fig_b10, width = 8.8, height = 7.4, dpi = dpi, bg = "white")

# ============================================================
# Taxon-level rates for B11-B14
# ============================================================

taxon_fraction_summary <- target_level %>%
  group_by(canonical_taxon, tool_label, Study, effective_fraction_pct, log10_effective_fraction_pct) %>%
  summarise(
    n = n(),
    good_rate = mean(recovery_class == "Good", na.rm = TRUE),
    average_rate = mean(recovery_class == "Intermediate", na.rm = TRUE),
    poor_missed_rate = mean(recovery_class == "Poor / missed", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    taxon_label = factor(
      unname(taxon_codes[as.character(canonical_taxon)]),
      levels = unname(taxon_codes[taxon_levels])
    )
  )

write_tsv(taxon_fraction_summary, file.path(outdir, "tables", "taxon_recovery_class_rates_by_fraction.tsv"))

taxon_fraction_long <- taxon_fraction_summary %>%
  pivot_longer(
    cols = c(good_rate, average_rate, poor_missed_rate),
    names_to = "rate_type",
    values_to = "rate"
  ) %>%
  mutate(
    recovery_panel = case_when(
      rate_type == "good_rate" ~ "Good recovery",
      rate_type == "average_rate" ~ "Intermediate recovery",
      rate_type == "poor_missed_rate" ~ "Poor / missed recovery",
      TRUE ~ rate_type
    ),
    recovery_panel = factor(
      recovery_panel,
      levels = c("Good recovery", "Intermediate recovery", "Poor / missed recovery")
    )
  )

write_tsv(taxon_fraction_long, file.path(outdir, "tables", "taxon_recovery_class_rates_long.tsv"))

make_taxon_plot <- function(dat, y_col, y_label, title_text, subtitle_text, outfile_prefix) {
  p <- ggplot(
    dat,
    aes(
      x = log10_effective_fraction_pct,
      y = .data[[y_col]],
      colour = tool_label,
      linetype = tool_label,
      shape = Study,
      group = interaction(tool_label, Study)
    )
  ) +
    geom_line(alpha = 0.65, linewidth = 0.55) +
    geom_point(size = 1.8, alpha = 0.8) +
    facet_wrap(~ taxon_label, ncol = 5) +
    scale_x_continuous(
      breaks = fraction_breaks_log,
      labels = pct_label_from_log
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_colour_manual(values = profiler_cols, drop = FALSE) +
    scale_linetype_manual(
      values = c("Kraken2 + Bracken" = "solid", "MetaPhlAn 4" = "22"),
      drop = FALSE
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Effective spike fraction per community member",
      y = y_label,
      colour = "Profiler",
      linetype = "Profiler",
      shape = "Cohort"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  ggsave(file.path(outdir, paste0(outfile_prefix, ".png")), p, width = 10, height = 7.5, dpi = dpi, bg = "white")
  save_pdf(file.path(outdir, paste0(outfile_prefix, ".pdf")), p, width = 10, height = 7.5)

  invisible(p)
}

# ============================================================
# B11. Taxon-specific good recovery
# ============================================================

make_taxon_plot(
  taxon_fraction_summary,
  "good_rate",
  "Good recovery rate",
  "Good recovery across community spike fractions",
  "Absolute relative recovery error ≤10%",
  "FigB11_taxon_good_recovery_vs_spike_fraction"
)

# ============================================================
# B12. Taxon-specific average/intermediate recovery
# ============================================================

make_taxon_plot(
  taxon_fraction_summary,
  "average_rate",
  "Intermediate recovery rate",
  "Intermediate recovery across community spike fractions",
  "Absolute relative recovery error >10% and ≤50%",
  "FigB12_taxon_average_intermediate_recovery_vs_spike_fraction"
)

# ============================================================
# B13. Taxon-specific poor/missed recovery
# ============================================================

make_taxon_plot(
  taxon_fraction_summary,
  "poor_missed_rate",
  "Poor / missed recovery rate",
  "Poor or missed recovery across community spike fractions",
  "Absolute relative recovery error >50% or target undetected",
  "FigB13_taxon_poor_missed_recovery_vs_spike_fraction"
)

# ============================================================
# B14. Optional combined taxon recovery-class overview
# ============================================================

fig_b14 <- ggplot(
  taxon_fraction_long,
  aes(
    x = log10_effective_fraction_pct,
    y = rate,
    colour = tool_label,
    shape = Study,
    group = interaction(tool_label, Study)
  )
) +
  geom_line(alpha = 0.65, linewidth = 0.45) +
  geom_point(size = 1.35, alpha = 0.8) +
  facet_grid(recovery_panel ~ canonical_taxon) +
  scale_x_continuous(
    breaks = fraction_breaks_log,
    labels = pct_label_from_log
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1)
  ) +
  labs(
    x = "Effective spike fraction per community member",
    y = "Recovery-class rate",
    colour = "Profiler",
    shape = "Cohort"
  ) +
  theme_pub(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
    strip.text.x = element_text(size = 7),
    strip.text.y = element_text(angle = 0)
  )

save_pdf(file.path(outdir, "FigB14_combined_taxon_recovery_classes_vs_spike_fraction.pdf"), fig_b14, width = 19, height = 10.5)
ggsave(file.path(outdir, "FigB14_combined_taxon_recovery_classes_vs_spike_fraction.png"), fig_b14, width = 19, height = 10.5, dpi = dpi, bg = "white")

# ============================================================
# Caption and manifest tables
# ============================================================

class_by_fraction <- target_level %>%
  count(effective_fraction_pct, tool_label, recovery_class, name = "n") %>%
  group_by(effective_fraction_pct, tool_label) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

write_tsv(class_by_fraction, file.path(outdir, "tables", "caption_recovery_class_by_fraction.tsv"))

taxon_caption <- taxon_fraction_summary %>%
  group_by(canonical_taxon, tool_label) %>%
  summarise(
    max_good_rate = max(good_rate, na.rm = TRUE),
    min_poor_missed_rate = min(poor_missed_rate, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(taxon_caption, file.path(outdir, "tables", "caption_taxon_recovery_summary.tsv"))

writeLines(
  c(
    "Generated outputs",
    "=================",
    "",
    file.path(outdir, "FigB8_QC_only_ABC.pdf"),
    file.path(outdir, "FigB9_recovery_classes_vs_spike_fraction.pdf"),
    file.path(outdir, "FigB10_recovery_class_composition_by_read_support.pdf"),
    file.path(outdir, "FigB11_taxon_good_recovery_vs_spike_fraction.pdf"),
    file.path(outdir, "FigB12_taxon_average_intermediate_recovery_vs_spike_fraction.pdf"),
    file.path(outdir, "FigB13_taxon_poor_missed_recovery_vs_spike_fraction.pdf"),
    file.path(outdir, "FigB14_combined_taxon_recovery_classes_vs_spike_fraction.pdf"),
    "",
    "Recovery-class thresholds",
    "=========================",
    "",
    "Good: abs(relative_error) <= 0.10",
    "Intermediate: 0.10 < abs(relative_error) <= 0.50",
    "Poor / missed: abs(relative_error) > 0.50, zero recovery, missing recovery, or non-finite recovery",
    "",
    "Recovery-class colours",
    "======================",
    "",
    "Good: green",
    "Intermediate: yellow",
    "Poor / missed: orange",
    "",
    "Spike-fraction convention",
    "=========================",
    "",
    "Effective spike fraction is expressed per community member.",
    "For a 10-member community spike-in, effective per-member fraction = total community spike fraction / 10.",
    "",
    "Canonical taxa",
    "==============",
    "",
    paste(taxon_levels, collapse = "\n")
  ),
  con = file.path(outdir, "MANIFEST.txt")
)

message("[OK] Wrote all reproducible supplementary figures to: ", outdir)
message("[OK] Main outputs:")
message("  FigB8_QC_only_ABC")
message("  FigB9_recovery_classes_vs_spike_fraction")
message("  FigB10_recovery_class_composition_by_read_support")
message("  FigB11_taxon_good_recovery_vs_spike_fraction")
message("  FigB12_taxon_average_intermediate_recovery_vs_spike_fraction")
message("  FigB13_taxon_poor_missed_recovery_vs_spike_fraction")
message("  FigB14_combined_taxon_recovery_classes_vs_spike_fraction")
