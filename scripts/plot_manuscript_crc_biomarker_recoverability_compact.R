#!/usr/bin/env Rscript
# v4: improves Figure 4 readability; Panel A is a labelled threshold heatmap.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(forcats)
  library(grid)
})

`%||%` <- function(x, y) {
  if (!is.null(x) && length(x) > 0 && !all(is.na(x)) && !(is.character(x) && !nzchar(x[1]))) x else y
}

option_list <- list(
  make_option("--summary-dir", type = "character", default = ".",
              help = "Directory containing detection_threshold_by_study_v6.csv and species_driver_summary_by_study_condition_fraction_v6.csv [default %default]"),
  make_option("--threshold-csv", type = "character", default = NULL,
              help = "Optional explicit path to detection_threshold_by_study_v6.csv"),
  make_option("--driver-csv", type = "character", default = NULL,
              help = "Optional explicit path to species_driver_summary_by_study_condition_fraction_v6.csv"),
  make_option("--outdir", type = "character", default = "RUNS/plots_manuscript_crc_biomarker_recoverability_assocC_intuitive_dotC_naturestyle",
              help = "Output directory [default %default]"),
  make_option("--focus-fraction", type = "double", default = 0.0001,
              help = "Spike fraction used in Panels B and C [default %default = 0.01%%]"),
  make_option("--spike-label-order", type = "character", default = "Bfrag,Csym,Pmic,Hhat,Pana,Psto,Dpne,Fnuc,Pint,Porp",
              help = "Comma-separated order for biomarkers [default %default]"),
  make_option("--tool-order", type = "character", default = "kraken2_bracken,metaphlan4",
              help = "Comma-separated tool order [default %default]"),
  make_option("--condition-order", type = "character", default = "Control,Adenoma,CRC",
              help = "Comma-separated background condition order [default %default]"),
  make_option("--main-width", type = "double", default = 18.0,
              help = "Main figure width in inches [default %default]"),
  make_option("--main-height", type = "double", default = 17.2,
              help = "Main figure height in inches [default %default]"),
  make_option("--dpi", type = "integer", default = 320,
              help = "PNG resolution [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

parse_csv_arg <- function(x) trimws(unlist(strsplit(x, ",", fixed = TRUE)))

# Standardize condition names before converting them to factors. This keeps the
# requested order stable even if an input table uses "Adenomas", "CRC", or the
# full "colorectal carcinoma" label.
normalize_condition <- function(x) {
  x_chr <- trimws(as.character(x))
  x_key <- tolower(x_chr)
  
  dplyr::case_when(
    x_key %in% c("control", "healthy", "normal") ~ "Control",
    x_key %in% c("adenoma", "adenomas") ~ "Adenoma",
    x_key %in% c("crc", "colorectal carcinoma", "colorectal cancer") ~ "CRC",
    TRUE ~ x_chr
  )
}

spike_order <- parse_csv_arg(opt$`spike-label-order`)
tool_order_raw <- parse_csv_arg(opt$`tool-order`)
condition_order <- normalize_condition(parse_csv_arg(opt$`condition-order`))
study_levels <- c("FengQ_2015", "ZellerG_2014")
study_short_map <- c("FengQ_2015" = "FengQ", "ZellerG_2014" = "Zeller")

threshold_csv <- opt$`threshold-csv` %||% file.path(opt$`summary-dir`, "detection_threshold_by_study_v6.csv")
driver_csv <- opt$`driver-csv` %||% file.path(opt$`summary-dir`, "species_driver_summary_by_study_condition_fraction_v6.csv")

if (!file.exists(threshold_csv)) stop("Could not find threshold CSV: ", threshold_csv, call. = FALSE)
if (!file.exists(driver_csv)) stop("Could not find driver CSV: ", driver_csv, call. = FALSE)

threshold_tbl <- read_csv(threshold_csv, show_col_types = FALSE)
driver_tbl <- read_csv(driver_csv, show_col_types = FALSE)

naturestyle_tool <- function(x) dplyr::recode(
  as.character(x),
  kraken2_bracken = "Kraken2 + Bracken",
  metaphlan4 = "MetaPhlAn 4",
  .default = as.character(x)
)

tool_order_naturestyle <- naturestyle_tool(tool_order_raw)

tool_cols <- c(
  "Kraken2 + Bracken" = "#1B9E77",
  "MetaPhlAn 4" = "#7570B3"
)
condition_cols <- c(
  "Control" = "#4C78A8",
  "Adenoma" = "#D8A03A",
  "CRC" = "#D95F02"
)

nature_theme <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid.major.y = element_line(colour = "#ECECEC", linewidth = 0.28),
      panel.grid.major.x = element_line(colour = "#F2F2F2", linewidth = 0.22),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "#4D4D4D", linewidth = 0.34),
      axis.ticks = element_line(colour = "#4D4D4D", linewidth = 0.30),
      axis.ticks.length = unit(2.0, "pt"),
      strip.background = element_rect(fill = "#F8F8F8", colour = "#D9D9D9", linewidth = 0.42),
      strip.text = element_text(face = "bold", colour = "#222222", margin = margin(3, 3, 3, 3)),
      axis.title = element_text(face = "bold", colour = "#222222"),
      axis.text = element_text(colour = "#3F3F3F"),
      legend.title = element_text(face = "bold", colour = "#222222", size = rel(0.92)),
      legend.text = element_text(colour = "#333333", size = rel(0.86)),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.key.height = unit(10.5, "pt"),
      legend.key.width = unit(12, "pt"),
      legend.spacing.y = unit(2, "pt"),
      legend.box.spacing = unit(4, "pt"),
      plot.title = element_text(face = "bold", size = base_size + 2.3, hjust = 0, colour = "#111111"),
      plot.subtitle = element_blank(),
      plot.margin = margin(6, 6, 6, 6)
    )
}

save_plot_pair <- function(p, stem, width, height, dpi = 320) {
  pdf_path <- file.path(opt$outdir, paste0(stem, ".pdf"))
  png_path <- file.path(opt$outdir, paste0(stem, ".png"))
  if (capabilities("cairo")) {
    ggsave(pdf_path, p, width = width, height = height, device = cairo_pdf, bg = "white")
  } else {
    ggsave(pdf_path, p, width = width, height = height, bg = "white")
  }
  ggsave(png_path, p, width = width, height = height, dpi = dpi, bg = "white")
}

# ----------------------------
# Panel A: threshold heatmap (publication default) plus conservative dumbbell
# ----------------------------
threshold_to_step <- function(x) {
  dplyr::case_when(
    !is.finite(x) ~ 7,
    abs(x - 0.0001) < 1e-12 ~ 1,
    abs(x - 0.0005) < 1e-12 ~ 2,
    abs(x - 0.001)  < 1e-12 ~ 3,
    abs(x - 0.005)  < 1e-12 ~ 4,
    abs(x - 0.01)   < 1e-12 ~ 5,
    abs(x - 0.05)   < 1e-12 ~ 6,
    TRUE ~ 7
  )
}

threshold_axis <- tibble::tibble(
  threshold_step = 1:7,
  threshold_y = 8 - threshold_step,
  threshold_label = c("0.01%", "0.05%", "0.10%", "0.50%", "1.00%", "5.00%", "NR")
)

panelA <- threshold_tbl %>%
  transmute(
    tool = factor(naturestyle_tool(tool), levels = tool_order_naturestyle),
    background_study = factor(background_study, levels = study_levels),
    study_label = factor(study_short_map[as.character(background_study)], levels = c("FengQ", "Zeller")),
    background_condition = factor(normalize_condition(background_condition), levels = condition_order),
    spike_label = factor(row_label, levels = spike_order),
    threshold_fraction = as.numeric(threshold_fraction),
    threshold_step = threshold_to_step(as.numeric(threshold_fraction))
  ) %>%
  mutate(
    threshold_y = 8 - threshold_step,
    threshold_label = factor(
      threshold_axis$threshold_label[match(threshold_step, threshold_axis$threshold_step)],
      levels = threshold_axis$threshold_label
    ),
    spike_x = as.numeric(spike_label),
    x_plot = spike_x + ifelse(tool == "Kraken2 + Bracken", -0.14, 0.14)
  )

write_csv(
  panelA %>%
    transmute(
      profiler = as.character(tool),
      study = as.character(study_label),
      background = as.character(background_condition),
      biomarker = as.character(spike_label),
      threshold_fraction,
      displayed_threshold = as.character(threshold_label)
    ),
  file.path(opt$outdir, "panel_A_minimum_spike_fraction_heatmap_data.csv")
)

panelA_segments <- panelA %>%
  select(study_label, background_condition, spike_label, spike_x, tool, threshold_y) %>%
  tidyr::pivot_wider(names_from = tool, values_from = threshold_y) %>%
  mutate(x_left = spike_x - 0.14, x_right = spike_x + 0.14)

pA_conservative <- ggplot() +
  geom_hline(
    data = threshold_axis,
    aes(yintercept = threshold_y),
    color = "#ECECEC",
    linewidth = 0.32
  ) +
  geom_segment(
    data = panelA_segments,
    aes(x = x_left, xend = x_right, y = `Kraken2 + Bracken`, yend = `MetaPhlAn 4`),
    color = "#C7C7C7",
    linewidth = 0.70,
    lineend = "round",
    na.rm = TRUE
  ) +
  geom_point(
    data = panelA,
    aes(x = x_plot, y = threshold_y, color = tool),
    size = 2.65,
    alpha = 0.96,
    na.rm = TRUE
  ) +
  facet_grid(study_label ~ background_condition) +
  scale_color_manual(values = tool_cols, name = "Profiler") +
  guides(color = guide_legend(order = 1, override.aes = list(size = 3.2, alpha = 1))) +
  scale_x_continuous(
    breaks = seq_along(spike_order),
    labels = spike_order,
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  scale_y_continuous(
    breaks = threshold_axis$threshold_y,
    labels = threshold_axis$threshold_label,
    limits = c(0.65, 7.35),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "A. Minimum spike fraction required for biomarker recovery",
    x = "Biomarker",
    y = "Minimum spike fraction"
  ) +
  nature_theme(11.3) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    panel.spacing = unit(0.7, "lines"),
    legend.position = "right"
  )

threshold_cols <- c(
  "0.01%" = "#006D2C",
  "0.05%" = "#31A354",
  "0.10%" = "#78C679",
  "0.50%" = "#FDD49E",
  "1.00%" = "#FC8D59",
  "5.00%" = "#D7301F",
  "NR" = "#BDBDBD"
)

panelA_heatmap <- panelA %>%
  mutate(
    cohort_background = factor(
      paste(as.character(study_label), as.character(background_condition), sep = "__"),
      levels = c(
        "FengQ__Control", "FengQ__Adenoma", "FengQ__CRC",
        "Zeller__Control", "Zeller__Adenoma", "Zeller__CRC"
      )
    ),
    spike_label = factor(as.character(spike_label), levels = rev(spike_order)),
    label_colour = ifelse(as.character(threshold_label) %in% c("0.01%", "0.05%", "5.00%"), "white", "#202020")
  )

pA <- ggplot(panelA_heatmap, aes(x = cohort_background, y = spike_label, fill = threshold_label)) +
  geom_tile(colour = "white", linewidth = 0.75, width = 0.98, height = 0.98) +
  geom_text(
    aes(label = threshold_label, colour = label_colour),
    size = 3.15,
    fontface = "bold",
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 3.5, colour = "#555555", linewidth = 0.55) +
  facet_grid(. ~ tool) +
  scale_fill_manual(
    values = threshold_cols,
    limits = threshold_axis$threshold_label,
    drop = FALSE,
    name = "Minimum fraction\n(lower = better)"
  ) +
  scale_colour_identity() +
  scale_x_discrete(
    labels = c(
      "FengQ__Control" = "Control\nFengQ",
      "FengQ__Adenoma" = "Adenoma\nFengQ",
      "FengQ__CRC" = "CRC\nFengQ",
      "Zeller__Control" = "Control\nZeller",
      "Zeller__Adenoma" = "Adenoma\nZeller",
      "Zeller__CRC" = "CRC\nZeller"
    ),
    drop = FALSE,
    expand = expansion(add = 0)
  ) +
  scale_y_discrete(drop = FALSE, expand = expansion(add = 0)) +
  guides(fill = guide_legend(
    order = 1,
    ncol = 1,
    byrow = TRUE,
    override.aes = list(colour = NA)
  )) +
  labs(
    title = "A. Minimum spike fraction required for biomarker recovery",
    subtitle = "Lower values indicate better recovery; NR, not recovered at any tested fraction",
    x = "Background and cohort",
    y = "Biomarker"
  ) +
  nature_theme(11.5) +
  theme(
    plot.subtitle = element_text(size = 10.2, colour = "#444444", margin = margin(0, 0, 7, 0)),
    strip.text.x = element_text(size = 11.0, margin = margin(5, 5, 5, 5)),
    axis.text.x = element_text(size = 9.2, lineheight = 0.92),
    axis.text.y = element_text(size = 9.5, face = "bold"),
    panel.grid = element_blank(),
    panel.spacing.x = unit(1.0, "lines"),
    legend.position = "right",
    legend.key.height = unit(15, "pt"),
    plot.margin = margin(7, 8, 7, 8)
  )

# ----------------------------
# Prepare driver data for Panels B and C
# ----------------------------
focus_fraction <- as.numeric(opt$`focus-fraction`)
focus_fraction_label <- percent(focus_fraction, accuracy = 0.01)

main_driver <- driver_tbl %>%
  filter(abs(as.numeric(spike_fraction) - focus_fraction) < 1e-12) %>%
  transmute(
    tool = factor(naturestyle_tool(tool), levels = tool_order_naturestyle),
    background_condition = factor(normalize_condition(background_condition), levels = condition_order),
    background_study = factor(background_study, levels = study_levels),
    study_label = factor(study_short_map[as.character(background_study)], levels = c("FengQ", "Zeller")),
    spike_label = factor(spike_label, levels = spike_order),
    baseline_log_raw = as.numeric(baseline_log_raw),
    tax_detect_raw = as.numeric(tax_detect_raw),
    tax_recovery_raw = as.numeric(tax_recovery_raw),
    recovery_var_raw = as.numeric(recovery_var_raw),
    fp_taxa_raw = as.numeric(fp_taxa_raw),
    q_strength_raw = as.numeric(q_strength_raw)
  )

if (!nrow(main_driver)) {
  stop(sprintf("No driver rows found at spike_fraction = %s", focus_fraction), call. = FALSE)
}

write_csv(main_driver, file.path(opt$outdir, paste0("driver_data_at_", gsub("%", "pct", focus_fraction_label), ".csv")))

# ----------------------------
# Panel B: biomarker significance versus baseline / taxonomy-side drivers
# ----------------------------
driver_levels <- c("Baseline abundance", "Tax detect", "Recovery error", "Recovery variability")

relationship_tbl <- bind_rows(
  main_driver %>%
    transmute(tool, background_condition, background_study, study_label, spike_label,
              driver = factor("Baseline abundance", levels = driver_levels),
              driver_value = baseline_log_raw,
              q_strength_raw = q_strength_raw),
  main_driver %>%
    transmute(tool, background_condition, background_study, study_label, spike_label,
              driver = factor("Tax detect", levels = driver_levels),
              driver_value = tax_detect_raw,
              q_strength_raw = q_strength_raw),
  main_driver %>%
    transmute(tool, background_condition, background_study, study_label, spike_label,
              driver = factor("Recovery error", levels = driver_levels),
              driver_value = abs(log2(tax_recovery_raw)),
              q_strength_raw = q_strength_raw),
  main_driver %>%
    transmute(tool, background_condition, background_study, study_label, spike_label,
              driver = factor("Recovery variability", levels = driver_levels),
              driver_value = recovery_var_raw,
              q_strength_raw = q_strength_raw)
) %>%
  filter(is.finite(driver_value), is.finite(q_strength_raw))

q_line <- -log10(0.05)
q_plot_max <- suppressWarnings(as.numeric(stats::quantile(
  relationship_tbl$q_strength_raw[is.finite(relationship_tbl$q_strength_raw)],
  probs = 0.99,
  na.rm = TRUE,
  names = FALSE
)))
if (!is.finite(q_plot_max) || q_plot_max <= 0) q_plot_max <- max(relationship_tbl$q_strength_raw, na.rm = TRUE)
if (!is.finite(q_plot_max) || q_plot_max <= 0) q_plot_max <- 10
q_plot_max <- max(q_line * 1.5, q_plot_max * 1.08)

pB_scatter <- ggplot(
  relationship_tbl,
  aes(x = driver_value, y = q_strength_raw, fill = background_condition, shape = study_label)
) +
  geom_hline(yintercept = q_line, linetype = 3, linewidth = 0.48, color = "#666666") +
  geom_point(size = 2.25, alpha = 0.82, stroke = 0.36, colour = "#3A3A3A", na.rm = TRUE) +
  facet_grid(tool ~ driver, scales = "free_x") +
  scale_fill_manual(values = condition_cols, name = "Background") +
  scale_shape_manual(values = c("FengQ" = 21, "Zeller" = 24), name = "Study") +
  guides(
    fill = guide_legend(order = 2, override.aes = list(shape = 21, size = 3.0, alpha = 1)),
    shape = guide_legend(order = 3, override.aes = list(fill = "grey70", size = 3.0, alpha = 1))
  ) +
  coord_cartesian(ylim = c(0, q_plot_max), clip = "off") +
  labs(
    x = NULL,
    y = expression(Biomarker~significance~(-log[10](q)))
  ) +
  nature_theme(10.6) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "#EEEEEE", colour = "#C7C7C7", linewidth = 0.5),
    strip.text.x = element_text(face = "bold", size = 10.4, margin = margin(5, 4, 5, 4)),
    strip.text.y = element_text(face = "bold", size = 10.0, margin = margin(5, 5, 5, 5)),
    axis.text.x = element_text(size = 8.4),
    axis.text.y = element_text(size = 8.6),
    panel.spacing.x = unit(1.05, "lines"),
    panel.spacing.y = unit(0.68, "lines"),
    plot.margin = margin(2, 8, 2, 8)
  )

trend_tbl <- relationship_tbl %>%
  group_by(driver, tool) %>%
  mutate(
    driver_bin = dplyr::ntile(driver_value, 3),
    driver_bin = factor(driver_bin, levels = 1:3, labels = c("Low", "Mid", "High"))
  ) %>%
  ungroup() %>%
  group_by(driver, tool, driver_bin) %>%
  summarise(
    n = n(),
    median_q = median(q_strength_raw, na.rm = TRUE),
    q25 = quantile(q_strength_raw, 0.25, na.rm = TRUE, names = FALSE),
    q75 = quantile(q_strength_raw, 0.75, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  )

write_csv(trend_tbl, file.path(opt$outdir, paste0("panel_B_binned_median_trends_at_", gsub("%", "pct", focus_fraction_label), ".csv")))

pB_trend <- ggplot(trend_tbl, aes(x = driver_bin, y = median_q, group = tool, color = tool)) +
  geom_hline(yintercept = q_line, linetype = 3, linewidth = 0.40, color = "#666666") +
  geom_line(linewidth = 0.88, alpha = 0.95) +
  geom_point(size = 2.2, alpha = 0.98) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.09, linewidth = 0.42, alpha = 0.75) +
  facet_wrap(~ driver, nrow = 1) +
  scale_color_manual(values = tool_cols, name = "Profiler") +
  coord_cartesian(ylim = c(0, q_plot_max), clip = "off") +
  labs(
    x = "Driver tertile",
    y = expression(Median~ -log[10](q))
  ) +
  nature_theme(9.6) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8.7, face = "bold", margin = margin(3, 2, 3, 2)),
    strip.background = element_rect(fill = "#F5F5F5", colour = "#D5D5D5", linewidth = 0.35),
    axis.text.x = element_text(size = 8.0),
    axis.text.y = element_text(size = 8.0),
    axis.title = element_text(size = 9.0, face = "bold"),
    panel.spacing.x = unit(1.05, "lines"),
    plot.margin = margin(0, 8, 4, 8)
  )

pB_body <- (pB_scatter / pB_trend) +
  plot_layout(heights = c(1.0, 0.34), guides = "collect") &
  theme(legend.position = "right")

pB_header <- ggplot() +
  theme_void() +
  labs(title = paste0("B. Biomarker significance versus baseline and recovery-quality drivers at ", focus_fraction_label)) +
  theme(
    plot.title = element_text(face = "bold", size = 14.5, hjust = 0, colour = "#111111",
                              margin = margin(0, 0, 3, 0)),
    plot.margin = margin(0, 4, 0, 4)
  )

pB <- pB_header / pB_body +
  plot_layout(heights = c(0.06, 1.0))

# ----------------------------
# Panel C: compact association summary matrix
# ----------------------------
assoc_tbl <- relationship_tbl %>%
  mutate(
    combo = factor(
      paste0(as.character(background_condition), "\n", as.character(study_label)),
      levels = c(
        "Control\nFengQ", "Control\nZeller",
        "Adenoma\nFengQ", "Adenoma\nZeller",
        "CRC\nFengQ", "CRC\nZeller"
      )
    )
  ) %>%
  group_by(tool, driver, background_condition, study_label, combo) %>%
  summarise(
    n = sum(is.finite(driver_value) & is.finite(q_strength_raw)),
    spearman_rho = suppressWarnings(cor(driver_value, q_strength_raw, method = "spearman", use = "complete.obs")),
    .groups = "drop"
  ) %>%
  mutate(driver = factor(driver, levels = rev(driver_levels)))

write_csv(assoc_tbl, file.path(opt$outdir, paste0("panel_C_association_summary_at_", gsub("%", "pct", focus_fraction_label), ".csv")))

assoc_tbl <- assoc_tbl %>%
  mutate(
    abs_rho = abs(spearman_rho),
    driver = fct_recode(
      driver,
      "Detectability" = "Tax detect"
    )
  )

pC <- ggplot(assoc_tbl, aes(x = combo, y = driver)) +
  geom_tile(width = 0.94, height = 0.88, fill = "#FAFAFA", colour = "#E2E2E2", linewidth = 0.38) +
  geom_vline(xintercept = c(2.5, 4.5), colour = "#CFCFCF", linewidth = 0.35) +
  geom_point(
    aes(size = abs_rho, fill = spearman_rho),
    shape = 21,
    colour = "#333333",
    stroke = 0.22,
    alpha = 0.96,
    na.rm = TRUE
  ) +
  facet_grid(tool ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#B94A48",
    mid = "#F7F7F7",
    high = "#3366A8",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    na.value = "#DDDDDD",
    name = "Spearman\nrho"
  ) +
  scale_size_continuous(
    range = c(1.8, 7.5),
    limits = c(0, 1),
    breaks = c(0.25, 0.50, 0.75, 1.00),
    name = "Association\nstrength |rho|"
  ) +
  guides(
    fill = guide_colorbar(order = 4, barheight = unit(48, "pt"), barwidth = unit(10, "pt")),
    size = guide_legend(order = 5, override.aes = list(fill = "#777777", colour = "#333333", alpha = 0.90))
  ) +
  labs(
    title = paste0("C. Summary associations with biomarker significance at ", focus_fraction_label),
    x = NULL,
    y = NULL
  ) +
  nature_theme(10.4) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.text.y = element_text(size = 10.2, margin = margin(6, 6, 6, 6)),
    axis.text.x = element_text(size = 9.5, lineheight = 0.92),
    axis.text.y = element_text(size = 10.2, face = "bold"),
    panel.spacing.y = unit(0.72, "lines"),
    plot.margin = margin(4, 8, 8, 8)
  )

# ----------------------------
# Main figure and outputs
# ----------------------------
main_fig <- pA / pB / pC +
  plot_layout(heights = c(0.90, 1.05, 0.46), guides = "keep") +
  plot_annotation(
    title = "Canonical CRC biomarkers differ in biomarker recoverability across profilers and backgrounds"
  ) &
  theme(
    plot.title = element_text(face = "bold", size = 17.0, hjust = 0, colour = "#111111"),
    legend.position = "right"
  )

save_plot_pair(pA, "panel_A_crc_biomarker_thresholds_dumbbell_naturestyle", width = 16.2, height = 7.8, dpi = opt$dpi)
save_plot_pair(pA, "panel_A_crc_biomarker_thresholds_heatmap_naturestyle", width = 16.2, height = 7.8, dpi = opt$dpi)
save_plot_pair(pA_conservative, "panel_A_crc_biomarker_thresholds_dumbbell_conservative_naturestyle", width = 16.2, height = 7.8, dpi = opt$dpi)
save_plot_pair(pB, "panel_B_biomarker_outcome_vs_intuitive_drivers_with_trends_naturestyle", width = 15.8, height = 8.5, dpi = opt$dpi)
save_plot_pair(pC, "panel_C_intuitive_association_dot_summary_naturestyle", width = 14.0, height = 5.5, dpi = opt$dpi)
save_plot_pair(main_fig, "manuscript_crc_biomarker_recoverability_with_assocC_intuitive_dotC_naturestyle_panelBtitle", width = opt$`main-width`, height = opt$`main-height`, dpi = opt$dpi)

readme <- c(
  "Run example:",
  "Rscript scripts/03_plot_manuscript_crc_biomarker_recoverability_with_assocC_intuitive_DOTC_NATURESTYLE_panelBtitle.R \\",
  "  --summary-dir RUNS/species_driver_and_thresholds \\",
  "  --outdir RUNS/plots_manuscript_crc_biomarker_recoverability_assocC_intuitive_dotC_naturestyle \\",
  "  --focus-fraction 0.0001",
  "",
  "Expected inputs:",
  "- detection_threshold_by_study_v6.csv",
  "- species_driver_summary_by_study_condition_fraction_v6.csv",
  "",
  "Main outputs:",
  "- manuscript_crc_biomarker_recoverability_with_assocC_intuitive_dotC_naturestyle_panelBtitle.pdf/png",
  "- panel_A_crc_biomarker_thresholds_heatmap_naturestyle.pdf/png (publication default)",
  "- panel_A_crc_biomarker_thresholds_dumbbell_naturestyle.pdf/png (same heatmap under the legacy filename used by the figure assembler)",
  "- panel_A_crc_biomarker_thresholds_dumbbell_conservative_naturestyle.pdf/png (comparison candidate)",
  "- panel_B_biomarker_outcome_vs_intuitive_drivers_with_trends_naturestyle.pdf/png",
  "- panel_C_intuitive_association_dot_summary_naturestyle.pdf/png",
  "",
  "Tabular outputs:",
  paste0("- driver_data_at_", gsub("%", "pct", focus_fraction_label), ".csv"),
  paste0("- panel_B_binned_median_trends_at_", gsub("%", "pct", focus_fraction_label), ".csv"),
  paste0("- panel_C_association_summary_at_", gsub("%", "pct", focus_fraction_label), ".csv"),
  "",
  "Notes:",
  "- Panel A is a labelled heatmap. Rows are biomarkers, columns are cohort/background combinations, and profilers are shown in separate facets.",
  "- Panel A encodes the minimum spike fraction required for recovery; lower is better and NR means not recovered at any tested fraction.",
  "- The conservative connected-point version is retained as a separate comparison output; the heatmap is the publication default.",
  "- Recovery error is abs(log2(observed/expected)); 0 is ideal, higher is worse.",
  "- Panels B/C intentionally omit the false-positive taxa driver; off-target false positives are handled in the community artefact/filtering figures.",
  "- Panel C dot color shows Spearman rho with biomarker significance (-log10(q)): blue = positive, red = negative.",
  "- Panel C dot size shows the absolute association strength |rho|.",
  "- These revisions change only aesthetics and layout, not the underlying data, statistical summaries, or conclusions.",
  "",
  "Suggested Figure 4 caption text for Panel A:",
  "(A) Minimum spike fraction required for each independently spiked taxon to be recovered as a significantly enriched biomarker. Rows show taxa; columns show cohort and diagnostic background; profilers are displayed in separate facets. Tile colour and text give the minimum recovered fraction (lower values indicate better recovery); NR indicates that the target was not recovered at any tested fraction."
)
writeLines(readme, file.path(opt$outdir, "README_crc_biomarker_recoverability_assocC_naturestyle.txt"))

message("[OK] Wrote naturestyle CRC biomarker recoverability figure set under ", opt$outdir)
