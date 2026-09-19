#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop(paste("Missing", flag))
  args[[i + 1]]
}
input <- value("--scores")
outdir <- value("--outdir")
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE))) {
  stop("OUTDIR must be new or empty")
}
d <- read.delim(input, check.names = FALSE, stringsAsFactors = FALSE)
needed <- c("cohort", "analysis_population", "profiler", "feature",
            "perturbation_persistence", "effect_fidelity",
            "perturbation_reliability_score", "external_replication_status", "evidence_tier")
if (!all(needed %in% names(d)) || !nrow(d)) stop("Score table is incomplete")
for (x in c("perturbation_persistence", "effect_fidelity", "perturbation_reliability_score")) {
  if (any(!is.finite(d[[x]]) | d[[x]] < 0 | d[[x]] > 1)) stop(paste("Invalid", x))
}

display_profiler <- c(kraken2_bracken = "Kraken2 + Bracken", metaphlan4 = "MetaPhlAn 4")
d$Profiler <- unname(display_profiler[d$profiler])
d$Profiler[is.na(d$Profiler)] <- d$profiler[is.na(d$Profiler)]
d$Population <- tools::toTitleCase(d$analysis_population)
d$Cohort <- tools::toTitleCase(d$cohort)
d$Replication <- factor(d$external_replication_status,
  levels = c("DIRECTIONALLY_REPLICATED", "NOT_SIGNIFICANT_ELSEWHERE",
             "DIRECTIONALLY_DISCORDANT", "NOT_IN_OTHER_COHORT_UNIVERSE",
             "NOT_AVAILABLE_IN_SPARSE_LEDGER"),
  labels = c("Directionally replicated", "Not significant elsewhere",
             "Directionally discordant", "Absent from other feature universe",
             "Other-cohort status unavailable"))
pal <- c("Directionally replicated" = "#009E73", "Not significant elsewhere" = "#0072B2",
         "Directionally discordant" = "#D55E00", "Absent from other feature universe" = "#777777",
         "Other-cohort status unavailable" = "#999999")
theme_set(theme_bw(base_size = 12) + theme(strip.text = element_text(face = "bold"),
  legend.position = "bottom", plot.title.position = "plot"))
dir.create(file.path(outdir, "figures"), recursive = TRUE)
dir.create(file.path(outdir, "figure_source"), recursive = TRUE)

p1 <- ggplot(d, aes(perturbation_persistence, effect_fidelity, colour = Replication)) +
  geom_point(alpha = 0.72, size = 2) +
  facet_grid(Population ~ interaction(Cohort, Profiler, sep = "\n")) +
  scale_colour_manual(values = pal, drop = FALSE) +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(title = "Controlled perturbations map CRC-biomarker reliability",
       subtitle = "Implanted taxa excluded; external replication is an independent annotation",
       x = "Persistence across spike stress tests", y = "Disease-effect fidelity",
       colour = "External evidence")

p2 <- ggplot(d, aes(Replication, perturbation_reliability_score, fill = Replication)) +
  geom_boxplot(width = 0.42, outlier.shape = NA, alpha = 0.60) +
  geom_jitter(width = 0.08, alpha = 0.35, size = 1) +
  facet_grid(Population ~ Profiler, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(title = "Perturbation reliability is reported as a continuum",
       subtitle = "Equal-weight component score; group differences are descriptive, not truth validation",
       x = NULL, y = "Perturbation reliability score", fill = "External evidence") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

for (spec in list(list(p = p1, name = "reliability_landscape", w = 14, h = 7.5),
                  list(p = p2, name = "reliability_by_replication", w = 12, h = 7.5))) {
  ggsave(file.path(outdir, "figures", paste0(spec$name, ".pdf")), spec$p,
         width = spec$w, height = spec$h, device = cairo_pdf)
  ggsave(file.path(outdir, "figures", paste0(spec$name, ".png")), spec$p,
         width = spec$w, height = spec$h, dpi = 300)
}
write.table(d, file.path(outdir, "figure_source", "perturbation_reliability_scores.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
writeLines(c(
  "status=DEVELOPMENT_ONLY", "automatic_feature_removal=NO",
  "score_is_false_positive_probability=NO"), file.path(outdir, "DEVELOPMENT_ONLY.txt"))
writeLines(c(
  "# Perturbation reliability figures", "",
  "The landscape jointly shows call persistence and effect-size fidelity. External replication",
  "is encoded independently and is not used to calculate the reliability score.", "",
  "The distribution panel asks whether directionally replicated biomarkers are more robust.",
  "These plots annotate reliability; they do not validate automatic feature removal."
), file.path(outdir, "FIGURE_GUIDE.md"))
files <- c(normalizePath(input), list.files(outdir, recursive = TRUE, full.names = TRUE))
files <- files[!grepl("SUCCESS$|reliability_figures.sha256$", files)]
status <- system2("sha256sum", files,
                  stdout = file.path(outdir, "reliability_figures.sha256"))
if (!identical(status, 0L)) stop("Could not seal reliability figures")
writeLines("status=PASS", file.path(outdir, "SUCCESS"))
message("[PASS] Perturbation reliability figures: ", normalizePath(outdir))
