#!/usr/bin/env Rscript
# Polished v4: fixed manuscript-wide taxon order and improved label legibility.

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

manuscript_taxon_order <- c(
  "Bfrag", "Csym", "Dpne", "Fnuc", "Hhat",
  "Pmic", "Pana", "Psto", "Porp", "Pint"
)

option_list <- list(
  make_option("--indir", type = "character", default = "RUNS/spike_metrics",
              help = "Input directory containing target_member_errors_with_condition.csv [default %default]"),
  make_option("--input-file", type = "character", default = NULL,
              help = "Optional direct path to target_member_errors_with_condition.csv"),
  make_option("--outdir", type = "character", default = "RUNS/plots_manuscript_community_spike_concordance",
              help = "Output directory [default %default]"),
  make_option("--community-size", type = "integer", default = 10,
              help = "Number of species in the community spike. Effective per-species fraction = total community fraction / this value [default %default]"),
  make_option("--matched-effective-fractions", type = "character", default = "0.0001,0.0005,0.001",
              help = "Effective per-species fractions used for independent-vs-community concordance [default %default]"),
  make_option("--main_effective_fractions", type = "character", default = "0.00001,0.00005,0.0001,0.0005,0.001",
              help = "Effective per-species fractions shown in main community-only Panels B and C [default %default = 0.001%,0.005%,0.01%,0.05%,0.10%]"),
  make_option("--spike-label-order", type = "character", default = paste(manuscript_taxon_order, collapse = ","),
              help = "Comma-separated manuscript species order, displayed top-to-bottom [default %default]"),
  make_option("--good-threshold", type = "double", default = 0.10,
              help = "Absolute relative-error threshold for Good recovery [default %default]"),
  make_option("--average-threshold", type = "double", default = 0.50,
              help = "Absolute relative-error threshold for Intermediate recovery; larger values are Poor/missed [default %default]"),
  make_option("--main-width", type = "double", default = 15.5,
              help = "Main figure width in inches [default %default]"),
  make_option("--main-height", type = "double", default = 16.0,
              help = "Main figure height in inches [default %default]"),
  make_option("--dpi", type = "integer", default = 320,
              help = "PNG resolution [default %default]")
)
raw_args <- commandArgs(trailingOnly = TRUE)
raw_args <- sub("^--main-effective-fractions$", "--main_effective_fractions", raw_args)
opt <- parse_args(OptionParser(option_list = option_list), args = raw_args)

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

input_file <- opt$`input-file` %||% file.path(opt$indir, "target_member_errors_with_condition.csv")
if (!file.exists(input_file)) stop("Could not find target-member table at: ", input_file, call. = FALSE)

parse_num_vector <- function(x) {
  x %>% strsplit(",") %>% unlist() %>% trimws() %>% .[nzchar(.)] %>% as.numeric()
}
parse_chr_vector <- function(x) {
  x %>% strsplit(",") %>% unlist() %>% trimws() %>% .[nzchar(.)]
}

fmt_fraction <- function(x) {
  p <- as.numeric(x) * 100
  ifelse(p < 0.01, sprintf("%.3f%%", p), sprintf("%.2f%%", p))
}

matched_effective_fractions <- parse_num_vector(opt$`matched-effective-fractions`)
main_effective_fractions <- parse_num_vector(opt$main_effective_fractions)
preferred_label_order <- parse_chr_vector(opt$`spike-label-order`)
if (!identical(preferred_label_order, manuscript_taxon_order)) {
  stop(
    "--spike-label-order must match the manuscript-wide order exactly: ",
    paste(manuscript_taxon_order, collapse = ","),
    call. = FALSE
  )
}

tool_labeller <- c(
  kraken2_bracken = "Kraken2 + Bracken",
  metaphlan4 = "MetaPhlAn 4"
)
tool_levels <- unname(tool_labeller)

profiler_cols <- c(
  "Kraken2 + Bracken" = "#009E73",
  "MetaPhlAn 4" = "#6F5BD3"
)
condition_cols <- c(
  "Adenoma" = "#D8A03A",
  "CRC" = "#D95F02",
  "Control" = "#4C78A8"
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
  saveRDS(plot_obj, file.path(opt$outdir, paste0(stem, ".rds")))
}

classify_recovery <- function(abs_rel_err, ratio, good_thr, avg_thr) {
  case_when(
    is.na(ratio) | is.na(abs_rel_err) ~ "Poor / missed",
    abs_rel_err <= good_thr ~ "Good",
    abs_rel_err <= avg_thr ~ "Intermediate",
    TRUE ~ "Poor / missed"
  )
}

message("[INFO] Reading: ", input_file)
raw_tbl <- read_csv(input_file, show_col_types = FALSE)

required_cols <- c(
  "spike_mode", "spike_label", "tool", "spike_fraction_total",
  "relative_error", "observed_over_expected", "member_taxon",
  "Target_Condition", "Study", "base_id"
)
missing_cols <- setdiff(required_cols, names(raw_tbl))
if (length(missing_cols) > 0) {
  stop("Input table is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

# Map community-panel member taxa back onto the independent-spike short labels.
label_map <- raw_tbl %>%
  filter(spike_mode == "independent", !is.na(spike_label), spike_label != "CRCpanel") %>%
  distinct(member_taxon, spike_label) %>%
  group_by(member_taxon) %>%
  summarise(spike_label_short = first(spike_label), .groups = "drop")

tbl <- raw_tbl %>%
  filter(spike_mode %in% c("independent", "community")) %>%
  mutate(
    tool = recode(tool, !!!tool_labeller),
    tool = factor(tool, levels = tool_levels),
    spike_fraction_total = as.numeric(spike_fraction_total),
    observed_over_expected = as.numeric(observed_over_expected),
    abs_relative_error = abs(as.numeric(relative_error)),
    effective_fraction = if_else(
      spike_mode == "community",
      spike_fraction_total / opt$`community-size`,
      spike_fraction_total
    ),
    effective_fraction_label = fmt_fraction(effective_fraction),
    recovery_class = classify_recovery(
      abs_relative_error,
      observed_over_expected,
      opt$`good-threshold`,
      opt$`average-threshold`
    ),
    recovery_class = factor(recovery_class, levels = c("Good", "Intermediate", "Poor / missed"))
  ) %>%
  left_join(label_map, by = "member_taxon") %>%
  mutate(
    spike_label_clean = case_when(
      spike_mode == "independent" ~ as.character(spike_label),
      spike_mode == "community" ~ spike_label_short,
      TRUE ~ as.character(spike_label)
    )
  ) %>%
  filter(!is.na(spike_label_clean), spike_label_clean != "CRCpanel")

# Retain target rows only if the table has target_flag.
if ("target_flag" %in% names(tbl)) {
  tbl <- tbl %>% filter(is.na(target_flag) | target_flag)
}

if (!nrow(tbl)) stop("No usable target-member rows remained after filtering.", call. = FALSE)

message("[INFO] Rows after filtering: ", nrow(tbl))
message("[INFO] Independent base samples: ", n_distinct(tbl$base_id[tbl$spike_mode == "independent"]))
message("[INFO] Community base samples: ", n_distinct(tbl$base_id[tbl$spike_mode == "community"]))
message("[INFO] Community total fractions: ",
        paste(sort(unique(tbl$spike_fraction_total[tbl$spike_mode == "community"])), collapse = ", "))
message("[INFO] Community effective per-species fractions: ",
        paste(fmt_fraction(sort(unique(tbl$effective_fraction[tbl$spike_mode == "community"]))), collapse = ", "))
message("[INFO] Main community-only panels use effective fractions: ",
        paste(fmt_fraction(main_effective_fractions), collapse = ", "))
message("[INFO] Supplementary community-only panels retain all effective fractions.")

# Shared manuscript order. ggplot draws the first y-factor level at the bottom,
# so reverse once here to display Bfrag -> Pint from top to bottom.
missing_manuscript_taxa <- setdiff(manuscript_taxon_order, unique(tbl$spike_label_clean))
if (length(missing_manuscript_taxa)) {
  stop(
    "Input is missing manuscript taxa required for Figure 5: ",
    paste(missing_manuscript_taxa, collapse = ", "),
    call. = FALSE
  )
}
display_levels <- rev(manuscript_taxon_order)

# ----------------------------
# Summary functions
# ----------------------------
rate_summary <- function(dat, by_context = FALSE) {
  group_cols <- c("spike_mode", "tool", "spike_label_clean", "effective_fraction", "effective_fraction_label")
  if (by_context) group_cols <- c(group_cols, "Target_Condition", "Study")
  
  dat %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n = n(),
      good_rate = mean(recovery_class == "Good", na.rm = TRUE),
      avg_rate = mean(recovery_class == "Intermediate", na.rm = TRUE),
      bad_rate = mean(recovery_class == "Poor / missed", na.rm = TRUE),
      median_recovery = median(observed_over_expected, na.rm = TRUE),
      iqr_recovery = IQR(observed_over_expected, na.rm = TRUE),
      .groups = "drop"
    )
}

summary_context <- rate_summary(tbl, by_context = TRUE)
summary_overall <- rate_summary(tbl, by_context = FALSE)

write_csv(summary_context, file.path(opt$outdir, "community_independent_summary_by_context.csv"))
write_csv(summary_overall, file.path(opt$outdir, "community_independent_summary_overall.csv"))

# ----------------------------
# Panel A: independent-versus-community recovery at matched effective fractions
# ----------------------------
matched_labels <- fmt_fraction(matched_effective_fractions)

ind_ctx <- summary_context %>%
  filter(spike_mode == "independent", effective_fraction %in% matched_effective_fractions) %>%
  select(tool, spike_label_clean, effective_fraction, effective_fraction_label,
         Target_Condition, Study, ind_good_rate = good_rate)

comm_ctx <- summary_context %>%
  filter(spike_mode == "community", effective_fraction %in% matched_effective_fractions) %>%
  select(tool, spike_label_clean, effective_fraction, effective_fraction_label,
         Target_Condition, Study, comm_good_rate = good_rate)

concordance_tbl <- inner_join(
  ind_ctx, comm_ctx,
  by = c("tool", "spike_label_clean", "effective_fraction", "effective_fraction_label", "Target_Condition", "Study")
) %>%
  mutate(
    effective_fraction_label = factor(effective_fraction_label, levels = matched_labels),
    Target_Condition = factor(
      dplyr::recode(
        as.character(Target_Condition),
        "colorectal carcinoma" = "CRC",
        "Colorectal carcinoma" = "CRC",
        .default = as.character(Target_Condition)
      ),
      levels = c("Adenoma", "CRC", "Control")
    ),
    Study_short = recode(Study, FengQ_2015 = "FengQ", ZellerG_2014 = "Zeller", .default = Study)
  )

if (!nrow(concordance_tbl)) {
  warning("No matched independent/community rows found for the requested effective fractions. Supplementary concordance panel will be empty.")
}

rho_tbl <- concordance_tbl %>%
  group_by(tool, effective_fraction_label) %>%
  summarise(
    n = n(),
    rho = suppressWarnings(cor(ind_good_rate, comm_good_rate, method = "spearman", use = "complete.obs")),
    median_delta = median(comm_good_rate - ind_good_rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = ifelse(is.finite(rho),
                   paste0("rho=", sprintf("%.2f", rho), "\nΔ=", percent(median_delta, accuracy = 1)),
                   paste0("n=", n))
  )

write_csv(concordance_tbl, file.path(opt$outdir, "supp_concordance_matched_effective_fractions_by_context.csv"))
write_csv(rho_tbl, file.path(opt$outdir, "supp_concordance_spearman_summary_by_context.csv"))

# Main Panel A: aggregate across study/background to make the comparison directly readable.
mode_rate_tbl <- summary_overall %>%
  filter(spike_mode %in% c("independent", "community"),
         effective_fraction %in% matched_effective_fractions) %>%
  mutate(
    spike_mode_pretty = recode(
      spike_mode,
      independent = "Independent",
      community = "Community"
    ),
    effective_fraction_label = factor(effective_fraction_label, levels = matched_labels),
    spike_label_plot = factor(spike_label_clean, levels = display_levels)
  ) %>%
  group_by(tool, spike_label_clean, spike_label_plot, effective_fraction, effective_fraction_label, spike_mode_pretty) %>%
  summarise(good_rate = mean(good_rate, na.rm = TRUE), .groups = "drop")

mode_pair_tbl <- mode_rate_tbl %>%
  tidyr::pivot_wider(
    names_from = spike_mode_pretty,
    values_from = good_rate
  ) %>%
  mutate(
    spike_label_plot = factor(as.character(spike_label_plot), levels = display_levels)
  )

write_csv(mode_rate_tbl, file.path(opt$outdir, "panel_A_mode_good_recovery_matched_effective_fractions.csv"))

mode_cols <- c(
  "Independent" = "#4C78A8",
  "Community" = "#D95F02"
)

pA <- ggplot(mode_pair_tbl, aes(y = spike_label_plot)) +
  geom_segment(
    aes(x = Independent, xend = Community, yend = spike_label_plot),
    colour = "#C4C4C4",
    linewidth = 0.52,
    lineend = "round",
    na.rm = TRUE
  ) +
  geom_point(
    aes(x = Independent, colour = "Independent", shape = "Independent"),
    size = 2.9,
    alpha = 0.98,
    position = position_nudge(y = 0.10),
    na.rm = TRUE
  ) +
  geom_point(
    aes(x = Community, colour = "Community", shape = "Community"),
    size = 2.9,
    alpha = 0.98,
    position = position_nudge(y = -0.10),
    na.rm = TRUE
  ) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "#9B9B9B", linewidth = 0.40) +
  facet_grid(tool ~ effective_fraction_label) +
  scale_x_continuous(
    # Expand slightly beyond the biological 0–100% range so endpoint points are fully visible.
    limits = c(-0.035, 1.035),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_colour_manual(values = mode_cols, name = "Spike design") +
  scale_shape_manual(values = c("Independent" = 16, "Community" = 17), name = "Spike design") +
  labs(
    title = "A. Independent- and community-spike recovery at matched fractions",
    x = "Good recovery (absolute relative error ≤10%)",
    y = NULL
  ) +
  pub_theme(10.7) +
  theme(
    axis.text.y = element_text(size = 9.6, face = "bold"),
    axis.text.x = element_text(size = 9.1),
    strip.text = element_text(size = 9.8, face = "bold", margin = margin(4, 4, 4, 4)),
    panel.grid.major.y = element_line(colour = "#F1F1F1", linewidth = 0.24),
    panel.grid.major.x = element_line(colour = "#F4F4F4", linewidth = 0.22),
    panel.spacing.x = unit(0.95, "lines"),
    panel.spacing.y = unit(0.68, "lines"),
    legend.position = "right",
    plot.margin = margin(4, 8, 4, 4)
  )

# Supplementary version of the original background/study-resolved concordance scatter.
pA_concordance_supp <- ggplot(concordance_tbl, aes(x = ind_good_rate, y = comm_good_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#8A8A8A", linewidth = 0.45) +
  geom_point(
    aes(fill = Target_Condition, shape = Study_short),
    size = 2.45,
    stroke = 0.35,
    colour = "#333333",
    alpha = 0.85,
    na.rm = TRUE
  ) +
  geom_text(
    data = rho_tbl,
    aes(x = 0.04, y = 0.96, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.0,
    colour = "#333333"
  ) +
  facet_grid(tool ~ effective_fraction_label) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(labels = percent_format(accuracy = 1), breaks = c(0, 0.5, 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = c(0, 0.5, 1)) +
  scale_fill_manual(values = condition_cols, name = "Background") +
  scale_shape_manual(values = c("FengQ" = 21, "Zeller" = 24), name = "Study") +
  labs(
    title = "Supplementary. Background-resolved concordance between independent and community spikes",
    x = "Independent spikes: % good recovery\nbalanced subset",
    y = "Community spikes: % good recovery\nfull cohort"
  ) +
  pub_theme(10.7) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(colour = "#EEEEEE", linewidth = 0.25),
    strip.text = element_text(size = 9.2, face = "bold"),
    panel.spacing = unit(0.62, "lines"),
    plot.margin = margin(4, 8, 4, 4)
  )

# ----------------------------
# Panel B: community full-cohort recovery across all effective per-species fractions
# ----------------------------
community_overall_all <- summary_overall %>%
  filter(spike_mode == "community") %>%
  mutate(
    spike_label_plot = factor(spike_label_clean, levels = display_levels),
    effective_fraction_label = factor(
      effective_fraction_label,
      levels = fmt_fraction(sort(unique(effective_fraction)))
    )
  )

community_overall <- community_overall_all %>%
  filter(effective_fraction %in% main_effective_fractions) %>%
  mutate(
    effective_fraction_label = factor(
      effective_fraction_label,
      levels = fmt_fraction(main_effective_fractions)
    )
  )

if (!nrow(community_overall)) {
  stop("No community rows matched --main-effective-fractions: ",
       paste(main_effective_fractions, collapse = ", "), call. = FALSE)
}

make_rate_dumbbell <- function(dat, rate_col, title_text, x_label) {
  pair_tbl <- dat %>%
    mutate(
      tool_short = case_when(
        tool == "Kraken2 + Bracken" ~ "kraken",
        tool == "MetaPhlAn 4" ~ "metaphlan",
        TRUE ~ as.character(tool)
      )
    ) %>%
    group_by(spike_label_plot, effective_fraction_label, tool_short) %>%
    summarise(rate_value = mean(.data[[rate_col]], na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = tool_short,
      values_from = rate_value,
      values_fill = 0
    ) %>%
    mutate(spike_label_plot = factor(as.character(spike_label_plot), levels = display_levels))
  
  ggplot(pair_tbl, aes(y = spike_label_plot)) +
    geom_segment(
      aes(x = kraken, xend = metaphlan, yend = spike_label_plot),
      colour = "#C4C4C4",
      linewidth = 0.50,
      lineend = "round"
    ) +
    geom_point(
      aes(x = kraken, colour = "Kraken2 + Bracken", shape = "Kraken2 + Bracken"),
      size = 2.75, alpha = 0.98, position = position_nudge(y = 0.10)
    ) +
    geom_point(
      aes(x = metaphlan, colour = "MetaPhlAn 4", shape = "MetaPhlAn 4"),
      size = 2.75, alpha = 0.98, position = position_nudge(y = -0.10)
    ) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "#9B9B9B", linewidth = 0.38) +
    facet_wrap(~ effective_fraction_label, nrow = 1) +
    scale_x_continuous(
      limits = c(-0.035, 1.035),
      breaks = c(0, 0.5, 1),
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_colour_manual(values = profiler_cols, name = "Profiler") +
    scale_shape_manual(values = c("Kraken2 + Bracken" = 16, "MetaPhlAn 4" = 17), name = "Profiler") +
    labs(title = title_text, x = x_label, y = NULL) +
    pub_theme(10.6) +
    theme(
      axis.text.y = element_text(size = 9.4, face = "bold"),
      axis.text.x = element_text(size = 8.9),
      strip.text = element_text(size = 9.3, face = "bold", margin = margin(4, 3, 4, 3)),
      panel.grid.major.y = element_line(colour = "#F1F1F1", linewidth = 0.24),
      panel.grid.major.x = element_line(colour = "#F4F4F4", linewidth = 0.22),
      panel.spacing.x = unit(1.05, "lines"),
      legend.position = "right",
      plot.margin = margin(4, 8, 4, 4)
    )
}

pB <- make_rate_dumbbell(
  community_overall,
  rate_col = "good_rate",
  title_text = "B. Full-cohort community spike recovery across ultra-low and weak effective fractions",
  x_label = "Good recovery (absolute relative error ≤10%)"
)

pB_avg <- make_rate_dumbbell(
  community_overall,
  rate_col = "avg_rate",
  title_text = "Supplementary. Full-cohort community spike intermediate recovery",
  x_label = "% samples with intermediate recovery"
)

pB_bad <- make_rate_dumbbell(
  community_overall,
  rate_col = "bad_rate",
  title_text = "Supplementary. Full-cohort community spike poor or missed recovery",
  x_label = "% samples with poor or missed recovery"
)

pB_all <- make_rate_dumbbell(
  community_overall_all,
  rate_col = "good_rate",
  title_text = "Supplementary. Full-cohort community spike good recovery across all effective fractions",
  x_label = "% samples with good recovery"
)

pB_avg_all <- make_rate_dumbbell(
  community_overall_all,
  rate_col = "avg_rate",
  title_text = "Supplementary. Full-cohort community spike intermediate recovery across all effective fractions",
  x_label = "% samples with intermediate recovery"
)

pB_bad_all <- make_rate_dumbbell(
  community_overall_all,
  rate_col = "bad_rate",
  title_text = "Supplementary. Full-cohort community spike poor or missed recovery across all effective fractions",
  x_label = "% samples with poor or missed recovery"
)

# ----------------------------
# Panel C: community full-cohort bias and variability
# ----------------------------
fraction_levels_main <- levels(community_overall$effective_fraction_label)
fraction_levels_all <- levels(community_overall_all$effective_fraction_label)
fraction_cols <- c(
  "0.001%" = "#244D8F",
  "0.005%" = "#3B6EA8",
  "0.01%"  = "#2A9D8F",
  "0.05%"  = "#63B179",
  "0.10%"  = "#D99B2B",
  "0.50%"  = "#B65B84",
  "1.00%"  = "#666666"
)
fraction_cols_use <- fraction_cols[fraction_levels_main]
missing_cols <- setdiff(fraction_levels_main, names(fraction_cols_use))
if (length(missing_cols) > 0) {
  fraction_cols_use <- c(fraction_cols_use, setNames(hue_pal()(length(missing_cols)), missing_cols))
}

offset_map <- setNames(seq(-0.26, 0.26, length.out = length(fraction_levels_main)), fraction_levels_main)

bias_long <- community_overall %>%
  mutate(
    spike_label_chr = as.character(spike_label_clean),
    y_base = match(spike_label_chr, display_levels),
    y_plot = y_base + unname(offset_map[as.character(effective_fraction_label)])
  ) %>%
  select(tool, spike_label_chr, effective_fraction_label, median_recovery, iqr_recovery, y_plot) %>%
  pivot_longer(
    cols = c(median_recovery, iqr_recovery),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    median_recovery = "Bias\nmedian observed/expected",
                    iqr_recovery = "Variability\nIQR observed/expected"),
    metric = factor(metric, levels = c("Bias\nmedian observed/expected", "Variability\nIQR observed/expected"))
  )

ref_lines <- tibble::tibble(
  metric = factor(c("Bias\nmedian observed/expected", "Variability\nIQR observed/expected"),
                  levels = c("Bias\nmedian observed/expected", "Variability\nIQR observed/expected")),
  xintercept = c(1, 0)
)
ideal_bands <- tibble::tibble(
  metric = factor(c("Bias\nmedian observed/expected", "Variability\nIQR observed/expected"),
                  levels = c("Bias\nmedian observed/expected", "Variability\nIQR observed/expected")),
  xmin = c(0.9, 0.0),
  xmax = c(1.1, 0.12)
)

pC <- ggplot(bias_long, aes(x = value, y = y_plot, color = effective_fraction_label)) +
  geom_rect(
    data = ideal_bands,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#EAF4EE",
    alpha = 0.75,
    colour = NA
  ) +
  geom_vline(
    data = ref_lines,
    aes(xintercept = xintercept),
    inherit.aes = FALSE,
    linetype = 2,
    colour = "#6A6A6A",
    linewidth = 0.42
  ) +
  geom_hline(yintercept = seq_along(display_levels), colour = "#F0F0F0", linewidth = 0.23) +
  geom_line(
    aes(group = interaction(tool, metric, spike_label_chr)),
    colour = "#C7C7C7",
    linewidth = 0.40,
    alpha = 0.82,
    na.rm = TRUE
  ) +
  geom_point(size = 2.05, alpha = 0.94, na.rm = TRUE) +
  facet_grid(tool ~ metric, scales = "free_x") +
  scale_y_continuous(
    breaks = seq_along(display_levels),
    labels = display_levels,
    expand = expansion(mult = c(0.035, 0.035))
  ) +
  scale_color_manual(values = fraction_cols_use, name = "Effective\nfraction") +
  labs(
    title = "C. Full-cohort community spike bias and variability across ultra-low and weak effective fractions",
    x = NULL,
    y = NULL
  ) +
  pub_theme(10.6) +
  theme(
    axis.text.y = element_text(size = 8.7),
    axis.text.x = element_text(size = 8.6),
    strip.text = element_text(size = 9.3, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "#EFEFEF", linewidth = 0.24),
    legend.position = "right",
    panel.spacing.x = unit(1.0, "lines"),
    plot.margin = margin(4, 8, 4, 4)
  )


# Supplementary all-fraction bias/variability plot.
fraction_cols_all_use <- fraction_cols[fraction_levels_all]
missing_cols_all <- setdiff(fraction_levels_all, names(fraction_cols_all_use))
if (length(missing_cols_all) > 0) {
  fraction_cols_all_use <- c(fraction_cols_all_use, setNames(hue_pal()(length(missing_cols_all)), missing_cols_all))
}
offset_map_all <- setNames(seq(-0.30, 0.30, length.out = length(fraction_levels_all)), fraction_levels_all)

bias_long_all <- community_overall_all %>%
  mutate(
    spike_label_chr = as.character(spike_label_clean),
    y_base = match(spike_label_chr, display_levels),
    y_plot = y_base + unname(offset_map_all[as.character(effective_fraction_label)])
  ) %>%
  select(tool, spike_label_chr, effective_fraction_label, median_recovery, iqr_recovery, y_plot) %>%
  pivot_longer(
    cols = c(median_recovery, iqr_recovery),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    median_recovery = "Bias\nmedian observed/expected",
                    iqr_recovery = "Variability\nIQR observed/expected"),
    metric = factor(metric, levels = c("Bias\nmedian observed/expected", "Variability\nIQR observed/expected"))
  )

pC_all <- ggplot(bias_long_all, aes(x = value, y = y_plot, color = effective_fraction_label)) +
  geom_rect(
    data = ideal_bands,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#EAF4EE",
    alpha = 0.75,
    colour = NA
  ) +
  geom_vline(
    data = ref_lines,
    aes(xintercept = xintercept),
    inherit.aes = FALSE,
    linetype = 2,
    colour = "#6A6A6A",
    linewidth = 0.42
  ) +
  geom_hline(yintercept = seq_along(display_levels), colour = "#F0F0F0", linewidth = 0.23) +
  geom_line(
    aes(group = interaction(tool, metric, spike_label_chr)),
    colour = "#C7C7C7",
    linewidth = 0.40,
    alpha = 0.82,
    na.rm = TRUE
  ) +
  geom_point(size = 1.95, alpha = 0.94, na.rm = TRUE) +
  facet_grid(tool ~ metric, scales = "free_x") +
  scale_y_continuous(
    breaks = seq_along(display_levels),
    labels = display_levels,
    expand = expansion(mult = c(0.035, 0.035))
  ) +
  scale_color_manual(values = fraction_cols_all_use, name = "Effective\nfraction") +
  labs(
    title = "Supplementary. Full-cohort community spike bias and variability across all effective fractions",
    x = NULL,
    y = NULL
  ) +
  pub_theme(10.6) +
  theme(
    axis.text.y = element_text(size = 8.7),
    axis.text.x = element_text(size = 8.6),
    strip.text = element_text(size = 9.3, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "#EFEFEF", linewidth = 0.24),
    legend.position = "right",
    panel.spacing.x = unit(1.0, "lines"),
    plot.margin = margin(4, 8, 4, 4)
  )

main_title <- "Community spike-ins broadly reproduce independent-spike recovery trends"
main_fig <- pA / pB / pC +
  plot_layout(heights = c(0.95, 0.80, 1.0)) +
  plot_annotation(
    title = main_title,
    theme = theme(plot.title = element_text(face = "bold", size = 18.0, hjust = 0, colour = "#111111"))
  )

save_plot_set(pA, "panel_A_independent_vs_community_good_recovery_dumbbell", width = 12.0, height = 6.5, dpi = opt$dpi)
save_plot_set(pA_concordance_supp, "supp_panel_A_background_resolved_independent_vs_community_concordance", width = 11.5, height = 6.8, dpi = opt$dpi)
save_plot_set(pB, "panel_B_community_good_recovery_dumbbell_main_effective_fractions", width = 12.5, height = 5.5, dpi = opt$dpi)
save_plot_set(pC, "panel_C_community_bias_variability_main_effective_fractions", width = 12.0, height = 6.5, dpi = opt$dpi)
save_plot_set(main_fig, "manuscript_community_spike_concordance_overview", width = opt$`main-width`, height = opt$`main-height`, dpi = opt$dpi)

save_plot_set(pB_avg, "supp_panel_B_community_average_recovery_dumbbell_main_effective_fractions", width = 12.5, height = 5.5, dpi = opt$dpi)
save_plot_set(pB_bad, "supp_panel_B_community_bad_recovery_dumbbell_main_effective_fractions", width = 12.5, height = 5.5, dpi = opt$dpi)

save_plot_set(pB_all, "supp_panel_B_community_good_recovery_dumbbell_all_effective_fractions", width = 15.0, height = 5.5, dpi = opt$dpi)
save_plot_set(pB_avg_all, "supp_panel_B_community_average_recovery_dumbbell_all_effective_fractions", width = 15.0, height = 5.5, dpi = opt$dpi)
save_plot_set(pB_bad_all, "supp_panel_B_community_bad_recovery_dumbbell_all_effective_fractions", width = 15.0, height = 5.5, dpi = opt$dpi)
save_plot_set(pC_all, "supp_panel_C_community_bias_variability_all_effective_fractions", width = 12.5, height = 6.5, dpi = opt$dpi)

readme <- c(
  "Run example:",
  "Rscript scripts/03_plot_manuscript_community_spike_concordance_MAINFRACTIONS.R \\",
  "  --input-file RUNS/spike_metrics/target_member_errors_with_condition.csv \\",
  "  --outdir RUNS/plots_manuscript_community_spike_concordance \\",
  "  --community-size 10 \\",
  "  --main-effective-fractions 0.00001,0.00005,0.0001,0.0005,0.001",
  "",
  "Main outputs:",
  "- manuscript_community_spike_concordance_overview.(pdf|png)",
  "- panel_A_independent_vs_community_good_recovery_dumbbell.(pdf|png)",
  "- panel_B_community_good_recovery_dumbbell_main_effective_fractions.(pdf|png)",
  "- panel_C_community_bias_variability_main_effective_fractions.(pdf|png)",
  "",
  "Supplementary outputs:",
  "- supp_panel_A_background_resolved_independent_vs_community_concordance.(pdf|png)",
  "- supp_panel_B_community_average_recovery_dumbbell_main_effective_fractions.(pdf|png)",
  "- supp_panel_B_community_bad_recovery_dumbbell_main_effective_fractions.(pdf|png)",
  "- supp_panel_B_community_good_recovery_dumbbell_all_effective_fractions.(pdf|png)",
  "- supp_panel_B_community_average_recovery_dumbbell_all_effective_fractions.(pdf|png)",
  "- supp_panel_B_community_bad_recovery_dumbbell_all_effective_fractions.(pdf|png)",
  "- supp_panel_C_community_bias_variability_all_effective_fractions.(pdf|png)",
  "",
  "Notes:",
  "- Community total spike fractions are divided by --community-size to obtain effective per-species fractions.",
  "- Main Panel A compares independent spikes from the balanced subset with community spikes from the full cohort at matched effective fractions after aggregating across background/study.",
  "- The supplementary Panel A concordance scatter retains the background/study-resolved comparison.",
  "- Main Panel B and C use all full-cohort community samples but only --main-effective-fractions.",
  "- The script normalizes --main-effective-fractions to --main_effective_fractions internally for optparse compatibility.",
  "- Supplementary Panel B and C variants retain all effective per-species fractions."
)
writeLines(readme, con = file.path(opt$outdir, "README_community_spike_concordance.txt"))

message("[OK] Wrote community-spike concordance figure set under: ", opt$outdir)
