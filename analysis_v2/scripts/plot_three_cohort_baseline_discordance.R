#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
value <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Missing ", flag)
  args[[i + 1L]]
}
input <- value("--input")
outdir <- value("--outdir")
if (dir.exists(outdir)) stop("Output directory already exists: ", outdir)
dir.create(outdir, recursive = TRUE)
x <- read.delim(input, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("cohort", "sample_id", "condition", "profiler", "target_label", "abundance_fraction")
if (length(setdiff(required, names(x)))) stop("Missing baseline figure columns")
if (!setequal(unique(x$cohort), c("feng", "yachida", "zeller"))) stop("Three-cohort gate failed")
if (!setequal(unique(x$profiler), c("kraken2_bracken", "metaphlan4"))) stop("Two-profiler gate failed")
x$abundance_fraction <- as.numeric(x$abundance_fraction)
if (any(!is.finite(x$abundance_fraction) | x$abundance_fraction < 0 |
        x$abundance_fraction > 1)) stop("Invalid abundance")
x$cohort <- factor(x$cohort, levels = c("feng", "yachida", "zeller"),
                   labels = c("Feng", "Yachida", "Zeller"))
x$profiler <- factor(x$profiler, levels = c("kraken2_bracken", "metaphlan4"),
                     labels = c("Kraken2 + Bracken", "MetaPhlAn 4"))
x$target_label <- factor(x$target_label, levels = c("Bfrag", "Fnuc", "Pint", "Pmic"))
condition <- tolower(trimws(x$condition))
condition[condition %in% c("colorectal carcinoma", "cancer", "case")] <- "crc"
condition[condition %in% c("healthy", "normal")] <- "control"
if (any(!condition %in% c("control", "adenoma", "crc")))
  stop("Unknown clinical condition: ", paste(unique(condition[!condition %in% c("control", "adenoma", "crc")]), collapse = ", "))
x$condition <- factor(condition, levels = c("control", "adenoma", "crc"),
                      labels = c("Control", "Adenoma", "CRC"))
x$positive <- x$abundance_fraction > 0
group_fields <- c("cohort", "target_label", "condition", "profiler")
counts <- aggregate(x$positive, x[group_fields], function(v) c(n = length(v), positive = sum(v)))
counts$n <- counts$x[, "n"]
counts$positive_n <- counts$x[, "positive"]
counts$x <- NULL
counts$prevalence <- counts$positive_n / counts$n
z <- qnorm(.975)
denominator <- 1 + z^2 / counts$n
center <- (counts$prevalence + z^2 / (2 * counts$n)) / denominator
radius <- z * sqrt(counts$prevalence * (1 - counts$prevalence) / counts$n +
                   z^2 / (4 * counts$n^2)) / denominator
counts$lower_95 <- pmax(0, center - radius)
counts$upper_95 <- pmin(1, center + radius)
write.table(counts, file.path(outdir, "baseline_prevalence_wilson.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
positive <- x[x$positive, ]
write.table(positive, file.path(outdir, "baseline_positive_abundance.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
colors <- c("Kraken2 + Bracken" = "#009E73", "MetaPhlAn 4" = "#6F5BD3")
theme_plot <- theme_bw(base_size = 10) + theme(
  panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
  strip.background = element_rect(fill = "grey95"),
  axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")
p1 <- ggplot(counts, aes(condition, prevalence, fill = profiler)) +
  geom_col(position = position_dodge(width = .72), width = .62) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                position = position_dodge(width = .72), width = .16) +
  facet_grid(cohort ~ target_label) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(limits = c(0, 1), labels = function(v) paste0(round(100*v), "%")) +
  labs(title = "A. Full-cohort baseline prevalence by condition",
       subtitle = "Bars: fraction of included unspiked samples with positive abundance; 95% Wilson intervals",
       x = NULL, y = "Positive baseline samples", fill = NULL) + theme_plot
p2 <- ggplot(positive, aes(condition, abundance_fraction, color = profiler, fill = profiler)) +
  geom_boxplot(position = position_dodge(width = .72), width = .58,
               outlier.shape = NA, alpha = .22) +
  geom_point(position = position_jitterdodge(jitter.width = .10, dodge.width = .72),
             alpha = .35, size = .55) +
  facet_grid(cohort ~ target_label, scales = "free_y") +
  scale_color_manual(values = colors) + scale_fill_manual(values = colors) +
  scale_y_log10(labels = function(v) paste0(signif(100*v, 2), "%")) +
  labs(title = "B. Baseline abundance among positive samples",
       subtitle = "Only positive samples; abundance shown as percent on a log scale",
       x = NULL, y = "Baseline abundance (%)", color = NULL, fill = NULL) +
  theme_plot + theme(legend.position = "none")
for (name in c("baseline_prevalence", "baseline_positive_abundance")) {
  plot <- if (name == "baseline_prevalence") p1 else p2
  width <- 14
  height <- 6
  ggsave(file.path(outdir, paste0(name, ".pdf")), plot,
         width = width, height = height, units = "in")
  ggsave(file.path(outdir, paste0(name, ".png")), plot,
         width = width, height = height, units = "in", dpi = 250)
}
draw_combined <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1,
                                                  heights = grid::unit(c(1, 1.1), "null"))))
  print(p1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1), newpage = FALSE)
  print(p2, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1), newpage = FALSE)
  grid::popViewport()
}
pdf(file.path(outdir, "baseline_three_cohort_combined.pdf"), width = 15, height = 12)
draw_combined()
dev.off()
png(file.path(outdir, "baseline_three_cohort_combined.png"), width = 3750,
    height = 3000, res = 250)
draw_combined()
dev.off()
writeLines(c("DEVELOPMENT_ONLY", "Three-cohort input includes historical Feng/Zeller profiles."),
           file.path(outdir, "DEVELOPMENT_ONLY.txt"))
writeLines("status\tPASS", file.path(outdir, "SUCCESS"))
message("[PASS] Three-cohort baseline draft figures: ", outdir)
