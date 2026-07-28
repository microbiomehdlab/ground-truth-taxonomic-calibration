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
metadata_file <- get_arg("--metadata", "metadata_w_study.tsv")
qc_list <- get_arg("--qc-list", file.path(run_root, "qc_depth_original_samples", "qc_files.txt"))
recovery_file <- get_arg("--recovery-file", file.path(run_root, "spike_metrics", "target_member_errors_with_condition.csv"))
outdir <- get_arg("--outdir", file.path(run_root, "revised_manuscript_figures", "source_panels", "FigC8_community_depth_qc_recovery"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "models"), recursive = TRUE, showWarnings = FALSE)

message("[INFO] metadata: ", metadata_file)
message("[INFO] qc list: ", qc_list)
message("[INFO] recovery file: ", recovery_file)
message("[INFO] outdir: ", outdir)

meta <- read_tsv(metadata_file, show_col_types = FALSE) %>%
  mutate(
    sample_id = as.character(sample_id),
    Study = as.character(Study),
    Target_Condition = as.character(Target_Condition)
  )

sample_ids <- meta$sample_id

extract_sample_id <- function(x, ids) {
  x <- as.character(x)
  hits <- ids[vapply(ids, function(id) grepl(id, x, fixed = TRUE), logical(1))]
  if (length(hits) == 0) return(NA_character_)
  hits[which.max(nchar(hits))]
}

read_fastqc_total <- function(zipfile) {
  tryCatch({
    if (!file.exists(zipfile)) return(NA_real_)
    if (file.info(zipfile)$size <= 0) return(NA_real_)

    files <- unzip(zipfile, list = TRUE)$Name
    fqdata <- files[grepl("fastqc_data.txt$", files)]
    if (length(fqdata) < 1) return(NA_real_)

    con <- unz(zipfile, fqdata[1])
    on.exit(close(con), add = TRUE)

    txt <- readLines(con, warn = FALSE)
    line <- txt[grepl("^Total Sequences\\t", txt)]
    if (length(line) < 1) return(NA_real_)

    as.numeric(strsplit(line[1], "\t")[[1]][2])
  }, error = function(e) {
    message("[WARN] Skipping unreadable FastQC zip: ", zipfile)
    return(NA_real_)
  })
}

qc_files <- readLines(qc_list, warn = FALSE)
qc_files <- qc_files[grepl("_fastqc\\.zip$", qc_files)]

message("[INFO] Reading FastQC zip files: ", length(qc_files))

qc <- tibble(path = qc_files) %>%
  mutate(
    sample_id = vapply(path, extract_sample_id, character(1), ids = sample_ids),
    qc_stage = case_when(
      grepl("/qc_before/", path) ~ "qc_before",
      grepl("/qc_after/", path) ~ "qc_after",
      TRUE ~ NA_character_
    ),
    read_mate = case_when(
      grepl("_1_fastqc\\.zip$", path) ~ "R1",
      grepl("_2_fastqc\\.zip$", path) ~ "R2",
      TRUE ~ NA_character_
    ),
    total_sequences = vapply(path, read_fastqc_total, numeric(1))
  ) %>%
  filter(!is.na(sample_id), !is.na(qc_stage), !is.na(read_mate), !is.na(total_sequences))

write_tsv(qc, file.path(outdir, "fastqc_read_counts_long.tsv"))

qc_sample <- qc %>%
  group_by(sample_id, qc_stage) %>%
  summarise(
    read_pairs = min(total_sequences, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = qc_stage, values_from = read_pairs) %>%
  mutate(
    raw_read_pairs = qc_before,
    postqc_read_pairs = qc_after,
    retention_pct = 100 * postqc_read_pairs / raw_read_pairs,
    raw_read_pairs_million = raw_read_pairs / 1e6,
    postqc_read_pairs_million = postqc_read_pairs / 1e6
  ) %>%
  left_join(meta, by = "sample_id")

write_tsv(qc_sample, file.path(outdir, "fastqc_depth_by_sample.tsv"))

message("[INFO] Samples with depth table: ", nrow(qc_sample))
message("[INFO] Samples with post-QC depth: ", sum(!is.na(qc_sample$postqc_read_pairs_million)))
message("[INFO] Samples with retention: ", sum(!is.na(qc_sample$retention_pct)))

recovery <- read_csv(recovery_file, show_col_types = FALSE)

write_tsv(
  tibble(column = names(recovery)),
  file.path(outdir, "recovery_table_columns.tsv")
)

char_cols <- names(recovery)[vapply(recovery, function(x) is.character(x) || is.factor(x), logical(1))]

scores <- vapply(char_cols, function(cc) {
  ids <- vapply(recovery[[cc]], extract_sample_id, character(1), ids = sample_ids)
  sum(!is.na(ids))
}, numeric(1))

best_sample_col <- char_cols[which.max(scores)]
message("[INFO] Recovery sample column inferred as: ", best_sample_col)

# Choose the best target/taxon column.
candidate_taxon_cols <- names(recovery)[grepl("target|taxon|species|member|spike_label", names(recovery), ignore.case = TRUE)]

score_taxon_col <- function(cc) {
  x <- as.character(recovery[[cc]])
  x <- x[!is.na(x)]
  x <- x[!x %in% c("CRCpanel", "community", "independent")]
  length(unique(x))
}

taxon_scores <- vapply(candidate_taxon_cols, score_taxon_col, numeric(1))
best_taxon_col <- candidate_taxon_cols[which.max(taxon_scores)]

if (length(best_taxon_col) == 0 || is.na(best_taxon_col)) {
  stop("Could not infer target/taxon column. Check recovery_table_columns.tsv.")
}

message("[INFO] Target/taxon column inferred as: ", best_taxon_col)

tool_labs <- c(
  kraken2_bracken = "Kraken2 + Bracken",
  metaphlan4 = "MetaPhlAn 4"
)

# Good recovery definition:
# Prefer explicit recovery-class column if present; otherwise use abs(relative_error) <= 0.10.
class_col <- names(recovery)[grepl("recovery.*class|class", names(recovery), ignore.case = TRUE)]
class_col <- if (length(class_col) > 0) class_col[[1]] else NA_character_

message("[INFO] Recovery class column: ", ifelse(is.na(class_col), "none; using abs(relative_error) <= 0.10", class_col))

recovery2 <- recovery %>%
  mutate(
    sample_id = vapply(.data[[best_sample_col]], extract_sample_id, character(1), ids = sample_ids),
    target_taxon = as.character(.data[[best_taxon_col]]),
    tool_label = recode(as.character(tool), !!!tool_labs, .default = as.character(tool)),
    spike_fraction_total = as.numeric(spike_fraction_total),
    relative_error = suppressWarnings(as.numeric(relative_error)),
    observed_over_expected = suppressWarnings(as.numeric(observed_over_expected))
  )

if (!is.na(class_col)) {
  recovery2 <- recovery2 %>%
    mutate(
      good_recovery = tolower(as.character(.data[[class_col]])) %in% c("good", "well_recovered", "well recovered")
    )
} else {
  recovery2 <- recovery2 %>%
    mutate(
      good_recovery = !is.na(observed_over_expected) &
        !is.na(relative_error) &
        abs(relative_error) <= 0.10
    )
}

qc_join <- qc_sample %>%
  select(
    sample_id,
    raw_read_pairs,
    postqc_read_pairs,
    retention_pct,
    raw_read_pairs_million,
    postqc_read_pairs_million,
    Study,
    Target_Condition
  )

community <- recovery2 %>%
  select(-any_of(c(
    "Study",
    "Target_Condition",
    "raw_read_pairs",
    "postqc_read_pairs",
    "retention_pct",
    "raw_read_pairs_million",
    "postqc_read_pairs_million"
  ))) %>%
  filter(!is.na(sample_id)) %>%
  filter(spike_mode == "community") %>%
  filter(!is.na(spike_fraction_total)) %>%
  filter(!is.na(target_taxon)) %>%
  filter(!target_taxon %in% c("CRCpanel", "community", "Community", "NA", "")) %>%
  mutate(
    effective_fraction_per_species = spike_fraction_total / 10,
    effective_fraction_pct = 100 * effective_fraction_per_species
  ) %>%
  left_join(qc_join, by = "sample_id") %>%
  filter(!is.na(postqc_read_pairs), postqc_read_pairs > 0) %>%
  mutate(
    implanted_read_pairs_target = postqc_read_pairs * effective_fraction_per_species,
    log10_implanted = log10(pmax(implanted_read_pairs_target, 1)),
    Study = factor(Study),
    Target_Condition = factor(Target_Condition),
    tool_label = factor(tool_label, levels = c("Kraken2 + Bracken", "MetaPhlAn 4")),
    target_taxon = factor(target_taxon)
  )

write_tsv(community, file.path(outdir, "community_target_level_depth_recovery.tsv"))

message("[INFO] Community target-level rows: ", nrow(community))
message("[INFO] Community unique samples: ", n_distinct(community$sample_id))
message("[INFO] Community unique taxa: ", n_distinct(community$target_taxon))
message("[INFO] Community total fractions: ", paste(sort(unique(community$spike_fraction_total)), collapse = ", "))
message("[INFO] Community effective per-species %: ", paste(signif(sort(unique(community$effective_fraction_pct)), 4), collapse = ", "))

# Fractions for D-F:
# effective per-species 0.001%, 0.01%, 0.1%
wanted_eff <- c(0.00001, 0.0001, 0.001)

community_fig <- community %>%
  mutate(
    fraction_group = case_when(
      abs(effective_fraction_per_species - wanted_eff[1]) < 1e-12 ~ "0.001% effective per species",
      abs(effective_fraction_per_species - wanted_eff[2]) < 1e-12 ~ "0.01% effective per species",
      abs(effective_fraction_per_species - wanted_eff[3]) < 1e-12 ~ "0.1% effective per species",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(fraction_group)) %>%
  mutate(
    fraction_group = factor(
      fraction_group,
      levels = c(
        "0.001% effective per species",
        "0.01% effective per species",
        "0.1% effective per species"
      )
    )
  )

if (nrow(community_fig) == 0) {
  stop("No community rows found for the requested effective per-species fractions. Check community_target_level_depth_recovery.tsv.")
}

# Collapse duplicate target rows if any.
community_target <- community_fig %>%
  group_by(
    sample_id, Study, Target_Condition, tool_label,
    fraction_group, spike_fraction_total, effective_fraction_per_species,
    target_taxon
  ) %>%
  summarise(
    good_recovery = any(good_recovery, na.rm = TRUE),
    postqc_read_pairs = first(postqc_read_pairs),
    postqc_read_pairs_million = first(postqc_read_pairs_million),
    raw_read_pairs_million = first(raw_read_pairs_million),
    retention_pct = first(retention_pct),
    implanted_read_pairs_target = first(implanted_read_pairs_target),
    log10_implanted = first(log10_implanted),
    .groups = "drop"
  )

community_agg <- community_target %>%
  group_by(
    sample_id, Study, Target_Condition, tool_label,
    fraction_group, spike_fraction_total, effective_fraction_per_species,
    postqc_read_pairs, postqc_read_pairs_million,
    implanted_read_pairs_target, log10_implanted
  ) %>%
  summarise(
    n_targets = n_distinct(target_taxon),
    good_targets = sum(good_recovery, na.rm = TRUE),
    good_recovery_prop = good_targets / n_targets,
    good_recovery_pct = 100 * good_recovery_prop,
    .groups = "drop"
  )

write_tsv(community_target, file.path(outdir, "community_target_level_for_models.tsv"))
write_tsv(community_agg, file.path(outdir, "community_sample_level_good_recovery.tsv"))

message("[INFO] Community sample-level rows for figure: ", nrow(community_agg))

theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

qc_plot <- qc_sample %>%
  filter(!is.na(Study), !is.na(Target_Condition), !is.na(postqc_read_pairs_million)) %>%
  mutate(
    Study = factor(Study),
    Target_Condition = factor(Target_Condition)
  )

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
  scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
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

make_recovery_panel <- function(dat, frac_label, panel_letter) {
  d <- dat %>% filter(fraction_group == frac_label)

  ggplot(
    d,
    aes(
      x = implanted_read_pairs_target,
      y = good_recovery_prop,
      colour = tool_label,
      shape = Study,
      weight = n_targets
    )
  ) +
    geom_point(alpha = 0.65, size = 1.9) +
    geom_smooth(
      method = "glm",
      method.args = list(family = binomial),
      se = TRUE,
      linewidth = 0.7
    ) +
    scale_x_log10(labels = label_number(accuracy = 1)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(
      title = paste0(panel_letter, ". Community recovery at ", frac_label),
      x = "Implanted read pairs per target",
      y = "Good targets / 10",
      colour = "Profiler",
      shape = "Cohort"
    ) +
    theme_pub()
}

pD <- make_recovery_panel(community_agg, "0.001% effective per species", "D")
pE <- make_recovery_panel(community_agg, "0.01% effective per species", "E")
pF <- make_recovery_panel(community_agg, "0.1% effective per species", "F")

fig <- (pA | pB) / (pC | pD) / (pE | pF) +
  plot_annotation(
    title = "Supplementary Fig. C8. Sequencing depth, preprocessing QC, and community spike-in recovery",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))
  )

ggsave(file.path(outdir, "FigC8_community_depth_qc_recovery.pdf"), fig, width = 13.5, height = 14.5, bg = "white")
ggsave(file.path(outdir, "FigC8_community_depth_qc_recovery.png"), fig, width = 13.5, height = 14.5, dpi = 320, bg = "white")

# ----------------------------
# Models
# ----------------------------

sink(file.path(outdir, "models", "model_notes.txt"))
cat("Model notes\n")
cat("===========\n\n")
cat("N_implanted,target = post-QC read pairs x effective per-species spike fraction.\n")
cat("Effective per-species spike fraction = total community spike fraction / 10.\n")
cat("Good recovery uses the explicit recovery class if detected; otherwise abs(relative_error) <= 0.10.\n\n")
cat("Target-level GLMM is fitted at 0.01% effective per-species fraction.\n")
cat("Aggregated binomial GLM is fitted across the three plotted community fractions.\n")
sink()

# Primary target-level GLMM at 0.01% effective per species.
target_model_data <- community_target %>%
  filter(fraction_group == "0.01% effective per species") %>%
  mutate(
    good_recovery = as.integer(good_recovery),
    log10_implanted = as.numeric(log10_implanted),
    sample_id = factor(sample_id),
    tool_label = droplevels(factor(tool_label)),
    Study = droplevels(factor(Study)),
    Target_Condition = droplevels(factor(Target_Condition)),
    target_taxon = droplevels(factor(target_taxon))
  ) %>%
  filter(
    !is.na(good_recovery),
    !is.na(log10_implanted),
    !is.na(tool_label),
    !is.na(Study),
    !is.na(Target_Condition),
    !is.na(target_taxon)
  )

write_tsv(target_model_data, file.path(outdir, "models", "target_level_model_data_0p01_effective.tsv"))

if (requireNamespace("lme4", quietly = TRUE)) {
  message("[INFO] Fitting target-level GLMM with lme4.")

  m_glmm <- lme4::glmer(
    good_recovery ~ log10_implanted * tool_label + Study + Target_Condition + target_taxon + (1 | sample_id),
    data = target_model_data,
    family = binomial(),
    control = lme4::glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  )

  sink(file.path(outdir, "models", "target_level_glmm_0p01_effective_summary.txt"))
  print(summary(m_glmm))
  sink()

  glmm_coef <- as.data.frame(coef(summary(m_glmm)))
  glmm_coef$term <- rownames(glmm_coef)
  rownames(glmm_coef) <- NULL
  write_tsv(glmm_coef, file.path(outdir, "models", "target_level_glmm_0p01_effective_coefficients.tsv"))

} else {
  message("[WARN] lme4 not available. Fitting fixed-effect logistic GLM instead.")

  m_glm_target <- glm(
    good_recovery ~ log10_implanted * tool_label + Study + Target_Condition + target_taxon,
    data = target_model_data,
    family = binomial()
  )

  sink(file.path(outdir, "models", "target_level_glm_0p01_effective_summary_NO_RANDOM_EFFECT.txt"))
  print(summary(m_glm_target))
  sink()

  glm_coef <- as.data.frame(coef(summary(m_glm_target)))
  glm_coef$term <- rownames(glm_coef)
  rownames(glm_coef) <- NULL
  write_tsv(glm_coef, file.path(outdir, "models", "target_level_glm_0p01_effective_coefficients_NO_RANDOM_EFFECT.tsv"))
}

# Alternative aggregated binomial model across plotted fractions.
agg_model_data <- community_agg %>%
  mutate(
    failures = n_targets - good_targets,
    fraction_group = droplevels(factor(fraction_group)),
    tool_label = droplevels(factor(tool_label)),
    Study = droplevels(factor(Study)),
    Target_Condition = droplevels(factor(Target_Condition)),
    log10_implanted = as.numeric(log10_implanted)
  ) %>%
  filter(
    !is.na(good_targets),
    !is.na(failures),
    !is.na(log10_implanted),
    !is.na(tool_label),
    !is.na(Study),
    !is.na(Target_Condition),
    !is.na(fraction_group)
  )

write_tsv(agg_model_data, file.path(outdir, "models", "aggregated_binomial_model_data_all_plotted_fractions.tsv"))

m_binom <- glm(
  cbind(good_targets, failures) ~ log10_implanted * tool_label + fraction_group + Study + Target_Condition,
  data = agg_model_data,
  family = binomial()
)

sink(file.path(outdir, "models", "aggregated_binomial_glm_all_plotted_fractions_summary.txt"))
print(summary(m_binom))
sink()

binom_coef <- as.data.frame(coef(summary(m_binom)))
binom_coef$term <- rownames(binom_coef)
rownames(binom_coef) <- NULL
write_tsv(binom_coef, file.path(outdir, "models", "aggregated_binomial_glm_all_plotted_fractions_coefficients.tsv"))

message("[OK] Wrote figure and model outputs to: ", outdir)
