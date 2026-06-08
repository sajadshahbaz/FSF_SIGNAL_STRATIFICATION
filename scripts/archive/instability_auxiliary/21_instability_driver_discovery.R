#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 21: Instability Driver Discovery
#
# Purpose:
#   Identify condition-specific and shared instable genes,
#   compute instability burden, and rank instability drivers.
#
# Input:
#   results/real_data/condition_fsf/condition_wise_fsf_metrics.tsv
#
# Outputs:
#   results/real_data/instability_drivers/*_instable_genes.tsv
#   results/real_data/instability_drivers/instability_overlap_matrix.tsv
#   results/real_data/instability_drivers/gene_instability_burden.tsv
#   results/real_data/instability_drivers/top_instability_genes.tsv
#   results/real_data/figures/condition_instability_burden_barplot.png
#   results/logs/21_instability_driver_discovery.log
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

INPUT_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/condition_fsf/condition_wise_fsf_metrics.tsv"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/instability_drivers"
)

FIG_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/figures"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/21_instability_driver_discovery.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 21: Instability Driver Discovery\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

fsf <- read_tsv(INPUT_FILE, show_col_types = FALSE)

required_cols <- c(
  "condition",
  "feature_id",
  "SSI",
  "stability_deviation",
  "fsf_signal_class"
)

missing_cols <- setdiff(required_cols, colnames(fsf))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

fsf <- fsf %>%
  mutate(
    condition = as.character(condition),
    feature_id = as.character(feature_id),
    SSI = as.numeric(SSI),
    stability_deviation = as.numeric(stability_deviation),
    fsf_signal_class = as.character(fsf_signal_class)
  )

conditions <- sort(unique(fsf$condition))

cat("Conditions:", paste(conditions, collapse = ", "), "\n")
cat("Total rows:", nrow(fsf), "\n")
cat("Unique features:", n_distinct(fsf$feature_id), "\n\n")

# ------------------------------------------------------------
# A. Extract instable genes per condition
# ------------------------------------------------------------

instable <- fsf %>%
  filter(fsf_signal_class == "instable") %>%
  arrange(condition, desc(stability_deviation), feature_id)

for (cond in conditions) {
  cond_tbl <- instable %>%
    filter(condition == cond)

  out_file <- file.path(
    OUT_DIR,
    paste0(cond, "_instable_genes.tsv")
  )

  write_tsv(cond_tbl, out_file)

  cat("Saved", cond, "instable genes:",
      nrow(cond_tbl),
      "->",
      out_file,
      "\n")
}

cat("\n")

# ------------------------------------------------------------
# B. Overlap matrix between conditions
# ------------------------------------------------------------

instable_sets <- lapply(
  conditions,
  function(cond) {
    instable %>%
      filter(condition == cond) %>%
      pull(feature_id) %>%
      unique()
  }
)

names(instable_sets) <- conditions

overlap_matrix <- expand.grid(
  condition_a = conditions,
  condition_b = conditions,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    n_overlap = length(
      intersect(
        instable_sets[[condition_a]],
        instable_sets[[condition_b]]
      )
    ),
    n_a = length(instable_sets[[condition_a]]),
    n_b = length(instable_sets[[condition_b]]),
    jaccard = ifelse(
      length(union(instable_sets[[condition_a]], instable_sets[[condition_b]])) > 0,
      n_overlap / length(union(instable_sets[[condition_a]], instable_sets[[condition_b]])),
      NA_real_
    )
  ) %>%
  ungroup()

write_tsv(
  overlap_matrix,
  file.path(OUT_DIR, "instability_overlap_matrix.tsv")
)

cat("Saved instability overlap matrix.\n")

# ------------------------------------------------------------
# C. Instability burden per gene
# ------------------------------------------------------------

burden <- fsf %>%
  mutate(
    is_instable = fsf_signal_class == "instable"
  ) %>%
  group_by(feature_id) %>%
  summarise(
    instability_count = sum(is_instable),
    instability_conditions = paste(sort(condition[is_instable]), collapse = ";"),
    mean_SSI = mean(SSI),
    mean_stability_deviation = mean(stability_deviation),
    max_stability_deviation = max(stability_deviation),
    .groups = "drop"
  ) %>%
  arrange(desc(instability_count), desc(mean_stability_deviation), feature_id)

write_tsv(
  burden,
  file.path(OUT_DIR, "gene_instability_burden.tsv")
)

cat("Saved gene instability burden table.\n")

# ------------------------------------------------------------
# D. Top instability drivers
# ------------------------------------------------------------

top_instability <- burden %>%
  filter(instability_count > 0) %>%
  arrange(desc(instability_count), desc(mean_stability_deviation)) %>%
  mutate(rank = row_number()) %>%
  select(rank, everything())

write_tsv(
  top_instability,
  file.path(OUT_DIR, "top_instability_genes.tsv")
)

write_tsv(
  top_instability %>% slice_head(n = 50),
  file.path(OUT_DIR, "top50_instability_genes.tsv")
)

write_tsv(
  top_instability %>% slice_head(n = 100),
  file.path(OUT_DIR, "top100_instability_genes.tsv")
)

write_tsv(
  top_instability %>% slice_head(n = 200),
  file.path(OUT_DIR, "top200_instability_genes.tsv")
)

cat("Saved top instability driver tables.\n")

# ------------------------------------------------------------
# E. Condition instability burden summary
# ------------------------------------------------------------

condition_burden <- fsf %>%
  group_by(condition) %>%
  summarise(
    n_features = n(),
    n_instable = sum(fsf_signal_class == "instable"),
    instable_proportion = n_instable / n_features,
    mean_SSI = mean(SSI),
    mean_stability_deviation = mean(stability_deviation),
    median_stability_deviation = median(stability_deviation),
    .groups = "drop"
  ) %>%
  arrange(desc(instable_proportion))

write_tsv(
  condition_burden,
  file.path(OUT_DIR, "condition_instability_burden.tsv")
)

cat("Condition instability burden:\n")
print(condition_burden)
cat("\n")

# ------------------------------------------------------------
# F. Plot condition instability burden
# ------------------------------------------------------------

condition_order <- condition_burden %>%
  arrange(instable_proportion) %>%
  pull(condition)

p <- ggplot(
  condition_burden,
  aes(
    x = factor(condition, levels = condition_order),
    y = instable_proportion
  )
) +
  geom_col() +
  geom_text(
    aes(label = sprintf("%.1f%%", instable_proportion * 100)),
    vjust = -0.3,
    size = 3.5
  ) +
  labs(
    title = "Condition-wise instability burden",
    x = "Condition",
    y = "Instable feature proportion"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(FIG_DIR, "condition_instability_burden_barplot.png"),
  p,
  width = 7,
  height = 5,
  dpi = 300
)

cat("Saved condition instability burden plot.\n\n")

cat("Main outputs saved to:\n")
cat(OUT_DIR, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

