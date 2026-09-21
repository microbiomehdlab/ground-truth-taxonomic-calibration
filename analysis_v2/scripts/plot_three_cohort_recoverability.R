#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
value <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Missing ", flag)
  args[[i + 1L]]
}
source_dir <- value("--source-dir")
outdir <- value("--outdir")
if (dir.exists(outdir)) stop("Output directory already exists")
if (!file.exists(file.path(source_dir, "SUCCESS"))) stop("Unsealed figure source")
dir.create(outdir, recursive = TRUE)
thresholds <- read.delim(file.path(source_dir, "minimum_biomarker_fraction.tsv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
drivers <- read.delim(file.path(source_dir, "biomarker_recovery_drivers.tsv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(thresholds) != 180L || nrow(drivers) != 720L)
  stop("Incomplete three-cohort ten-target figure source")
cohort_levels <- c("feng", "yachida", "zeller")
cohort_labels <- c("Feng", "Yachida", "Zeller")
condition_levels <- c("Control", "Adenoma", "CRC")
profiler_levels <- c("kraken2_bracken", "metaphlan4")
profiler_labels <- c("Kraken2 + Bracken", "MetaPhlAn 4")
taxon_levels <- c("Bfrag", "Csym", "Dpne", "Fnuc", "Hhat", "Pmic",
                  "Pana", "Psto", "Porp", "Pint")
driver_levels <- c("baseline_log10", "detectability", "recovery_error", "recovery_iqr")
driver_labels <- c("Baseline abundance\nlog10(median + 1e-6)",
                   "Target detectability\nfraction observed", "Recovery error\nmedian absolute relative error",
                   "Recovery variability\nIQR recovered/implanted")
for (name in c("thresholds", "drivers")) {
  x <- get(name)
  x$cohort <- factor(x$cohort, levels = cohort_levels, labels = cohort_labels)
  x$condition <- factor(x$condition, levels = condition_levels)
  x$profiler <- factor(x$profiler, levels = profiler_levels, labels = profiler_labels)
  x$target_label <- factor(x$target_label, levels = rev(taxon_levels))
  if (anyNA(x$cohort) || anyNA(x$condition) || anyNA(x$profiler) || anyNA(x$target_label))
    stop("Unexpected figure facet or taxon")
  assign(name, x)
}
thresholds$minimum_fraction <- as.character(thresholds$minimum_fraction)
dose_values <- c(.0001, .0005, .001, .005, .01, .05)
dose_labels <- c("0.01%", "0.05%", "0.10%", "0.50%", "1.00%", "5.00%")
thresholds$dose_label <- ifelse(thresholds$minimum_fraction == "NR", "NR",
                                dose_labels[match(as.numeric(thresholds$minimum_fraction), dose_values)])
if (anyNA(thresholds$dose_label)) stop("Threshold outside frozen dose grid")
thresholds$dose_label <- factor(thresholds$dose_label, levels = c(dose_labels, "NR"))
thresholds$context <- factor(paste(thresholds$condition, thresholds$cohort, sep = "\n"),
                             levels = as.vector(outer(condition_levels, cohort_labels, paste, sep = "\n")))
drivers$driver <- factor(drivers$driver, levels = driver_levels, labels = driver_labels)
drivers$driver_value <- as.numeric(drivers$driver_value)
drivers$biomarker_strength <- as.numeric(drivers$biomarker_strength)
if (anyNA(drivers$driver) || any(!is.finite(drivers$driver_value)) ||
    any(!is.finite(drivers$biomarker_strength))) stop("Invalid driver data")
drivers$context <- factor(paste(drivers$condition, drivers$cohort, sep = "\n"),
                          levels = levels(thresholds$context))
theme_fig <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(),
  strip.background = element_rect(fill = "grey95"),
  axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
pA <- ggplot(thresholds, aes(context, target_label, fill = dose_label)) +
  geom_tile(color = "white", linewidth = .35) +
  geom_text(aes(label = dose_label), size = 2.5) +
  facet_grid(. ~ profiler) +
  scale_fill_manual(values = c("0.01%" = "#12345B", "0.05%" = "#344E73",
    "0.10%" = "#60718C", "0.50%" = "#9C946A", "1.00%" = "#D9BE55",
    "5.00%" = "#F4D842", "NR" = "#D9D9D9"), drop = FALSE) +
  labs(title = "A. Minimum spike fraction yielding a significant enriched target",
       subtitle = "Independent spikes; BH q ≤ 0.05; lower is better; NR = never recovered",
       x = NULL, y = "Implanted taxon", fill = "Minimum fraction") + theme_fig
colors <- c(Control = "#4C78A8", Adenoma = "#D8A03A", CRC = "#D95F02")
shapes <- c(Feng = 16, Yachida = 17, Zeller = 15)
pB <- ggplot(drivers, aes(driver_value, biomarker_strength,
                          color = condition, shape = cohort)) +
  geom_hline(yintercept = -log10(.05), linetype = 3, color = "grey55") +
  geom_point(size = 1.8, alpha = .75) +
  facet_grid(profiler ~ driver, scales = "free_x") +
  scale_color_manual(values = colors) + scale_shape_manual(values = shapes) +
  labs(title = "B. Biomarker significance versus baseline and recovery drivers at 0.01%",
       subtitle = "Each point is one target × cohort × condition; dashed line marks q = 0.05",
       x = "Driver value", y = expression(-log[10](q)), color = "Condition", shape = "Cohort") + theme_fig
correlation <- do.call(rbind, lapply(split(drivers,
  interaction(drivers$cohort, drivers$condition, drivers$profiler, drivers$driver, drop = TRUE)),
  function(part) {
    if (nrow(part) != 10L) stop("Correlation group does not contain ten targets")
    rho <- if (length(unique(part$driver_value)) < 2L ||
               length(unique(part$biomarker_strength)) < 2L) NA_real_ else
      suppressWarnings(cor(part$driver_value, part$biomarker_strength, method = "spearman"))
    data.frame(cohort = part$cohort[1], condition = part$condition[1],
               profiler = part$profiler[1], driver = part$driver[1],
               rho = rho, n_targets = nrow(part))
  }))
correlation$context <- factor(paste(correlation$condition, correlation$cohort, sep = "\n"),
                              levels = levels(thresholds$context))
correlation$driver <- factor(as.character(correlation$driver), levels = driver_labels)
write.table(correlation, file.path(outdir, "panel_C_spearman_source.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
pC <- ggplot(correlation, aes(context, driver, color = rho, size = abs(rho))) +
  geom_point(na.rm = TRUE) + facet_grid(profiler ~ .) +
  scale_color_gradient2(low = "#3B6FB6", mid = "white", high = "#C65353",
                        midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
  scale_size(range = c(1, 7), limits = c(0, 1)) +
  labs(title = "C. Descriptive target-level associations with biomarker significance",
       subtitle = "Spearman ρ across ten implanted taxa per cohort/condition; no inferential claim",
       x = NULL, y = NULL, color = "Spearman ρ", size = "|ρ|") + theme_fig
save_plot <- function(plot, name, width, height) {
  ggsave(file.path(outdir, paste0(name, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(outdir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 250)
}
save_plot(pA, "panel_A_minimum_spike_fraction", 17, 7)
save_plot(pB, "panel_B_recovery_drivers", 15, 8)
save_plot(pC, "panel_C_spearman_associations", 15, 6)
draw_combined <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(3, 1,
    heights = grid::unit(c(7, 8, 6), "null"))))
  print(pA, vp = grid::viewport(layout.pos.row = 1), newpage = FALSE)
  print(pB, vp = grid::viewport(layout.pos.row = 2), newpage = FALSE)
  print(pC, vp = grid::viewport(layout.pos.row = 3), newpage = FALSE)
  grid::popViewport()
}
pdf(file.path(outdir, "three_cohort_recoverability_combined.pdf"), width = 17, height = 21)
draw_combined()
dev.off()
png(file.path(outdir, "three_cohort_recoverability_combined.png"),
    width = 3400, height = 4200, res = 200)
draw_combined()
dev.off()
writeLines(c("DEVELOPMENT_ONLY", "Awaits corrected profiler-scale endpoints from all three cohorts."),
           file.path(outdir, "DEVELOPMENT_ONLY.txt"))
writeLines("status\tPASS", file.path(outdir, "SUCCESS"))
message("[PASS] Three-cohort recoverability figure: ", outdir)
