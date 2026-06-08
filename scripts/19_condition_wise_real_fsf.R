#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 19: Condition-wise Real FSF Stratification
#
# Purpose:
#   Compute SSI, Stability Deviation, and FSF signal classes
#   separately for each biological condition.
#
# Input:
#   results/real_data/effect_estimates/real_effect_estimates.tsv
#
# Outputs:
#   results/real_data/condition_fsf/condition_wise_fsf_metrics.tsv
#   results/real_data/condition_fsf/condition_wise_fsf_class_summary.tsv
#   results/real_data/condition_fsf/condition_wise_signal_state_summary.tsv
#   results/logs/19_condition_wise_real_fsf.log
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

INPUT_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/effect_estimates/real_effect_estimates.tsv"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/condition_fsf"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/19_condition_wise_real_fsf.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 19: Condition-wise Real FSF Stratification\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

tau <- 0.5

highly_stable_cutoff <- 0.90
stable_cutoff <- 0.75
weakly_stable_cutoff <- 0.60

cat("Input:", INPUT_FILE, "\n")
cat("Tau:", tau, "\n\n")

effects <- read_tsv(INPUT_FILE, show_col_types = FALSE)

required_cols <- c(
  "feature_id",
  "perturbation_id",
  "effect_estimate",
  "condition"
)

missing_cols <- setdiff(required_cols, colnames(effects))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

effects <- effects %>%
  mutate(
    feature_id = as.character(feature_id),
    perturbation_id = as.character(perturbation_id),
    condition = as.character(condition),
    effect_estimate = as.numeric(effect_estimate)
  )

cat("Rows:", nrow(effects), "\n")
cat("Features:", n_distinct(effects$feature_id), "\n")
cat("Conditions:", paste(sort(unique(effects$condition)), collapse = ", "), "\n\n")

# ------------------------------------------------------------
# Assign signal states
# ------------------------------------------------------------

states <- effects %>%
  mutate(
    signal_state = case_when(
      effect_estimate > tau  ~ "up",
      effect_estimate < -tau ~ "down",
      TRUE                   ~ "const"
    ),
    signal_numeric = case_when(
      signal_state == "up"    ~ 1L,
      signal_state == "const" ~ 0L,
      signal_state == "down"  ~ -1L,
      TRUE ~ NA_integer_
    )
  )

state_summary <- states %>%
  count(condition, signal_state, name = "n_assignments") %>%
  group_by(condition) %>%
  mutate(proportion = n_assignments / sum(n_assignments)) %>%
  ungroup() %>%
  arrange(condition, desc(n_assignments))

write_tsv(
  state_summary,
  file.path(OUT_DIR, "condition_wise_signal_state_summary.tsv")
)

cat("Condition-wise signal-state summary:\n")
print(state_summary, n = 100)
cat("\n")

# ------------------------------------------------------------
# Compute FSF per condition
# ------------------------------------------------------------

probs <- states %>%
  count(condition, feature_id, signal_state, name = "n_state") %>%
  group_by(condition, feature_id) %>%
  mutate(
    n_perturbations = sum(n_state),
    probability = n_state / n_perturbations
  ) %>%
  ungroup() %>%
  select(condition, feature_id, n_perturbations, signal_state, probability) %>%
  pivot_wider(
    names_from = signal_state,
    values_from = probability,
    values_fill = 0,
    names_prefix = "p_"
  )

for (col in c("p_up", "p_down", "p_const")) {
  if (!col %in% colnames(probs)) {
    probs[[col]] <- 0
  }
}

fsf <- probs %>%
  rowwise() %>%
  mutate(
    SSI = max(c(p_up, p_down, p_const), na.rm = TRUE),
    stability_deviation = 1 - SSI,

    dominant_state = case_when(
      p_up == SSI & p_up > p_down & p_up > p_const ~ "up",
      p_down == SSI & p_down > p_up & p_down > p_const ~ "down",
      p_const == SSI & p_const > p_up & p_const > p_down ~ "constant",
      TRUE ~ "mixed"
    ),

    stability_level = case_when(
      SSI >= highly_stable_cutoff ~ "highly_stable",
      SSI >= stable_cutoff ~ "stable",
      SSI >= weakly_stable_cutoff ~ "weakly_stable_transitional",
      TRUE ~ "instable"
    ),

    dominant_signal_identity = case_when(
      stability_level == "instable" ~ "instable_no_stable_identity",
      dominant_state == "up" ~ "stable_up_identity",
      dominant_state == "down" ~ "stable_down_identity",
      dominant_state == "constant" ~ "stable_constant_identity",
      TRUE ~ "mixed_identity"
    ),

    fsf_signal_class = case_when(
      stability_level == "instable" ~ "instable",
      dominant_state == "up" ~ paste0(stability_level, "_up"),
      dominant_state == "down" ~ paste0(stability_level, "_down"),
      dominant_state == "constant" ~ paste0(stability_level, "_constant"),
      TRUE ~ "mixed"
    )
  ) %>%
  ungroup() %>%
  select(
    condition,
    feature_id,
    n_perturbations,
    p_up,
    p_down,
    p_const,
    dominant_state,
    SSI,
    stability_deviation,
    stability_level,
    dominant_signal_identity,
    fsf_signal_class
  ) %>%
  arrange(condition, desc(SSI), stability_deviation, feature_id)

summary_tbl <- fsf %>%
  count(
    condition,
    stability_level,
    dominant_signal_identity,
    fsf_signal_class,
    name = "n_features"
  ) %>%
  group_by(condition) %>%
  mutate(proportion = n_features / sum(n_features)) %>%
  ungroup() %>%
  arrange(condition, stability_level, desc(n_features))

write_tsv(
  fsf,
  file.path(OUT_DIR, "condition_wise_fsf_metrics.tsv")
)

write_tsv(
  summary_tbl,
  file.path(OUT_DIR, "condition_wise_fsf_class_summary.tsv")
)

cat("Condition-wise FSF class summary:\n")
print(summary_tbl, n = 200)
cat("\n")

# ------------------------------------------------------------
# Save per-condition files
# ------------------------------------------------------------

conditions <- sort(unique(fsf$condition))

for (cond in conditions) {
  cond_dir <- file.path(OUT_DIR, cond)
  dir.create(cond_dir, recursive = TRUE, showWarnings = FALSE)

  write_tsv(
    fsf %>% filter(condition == cond),
    file.path(cond_dir, paste0("fsf_metrics_", cond, ".tsv"))
  )

  write_tsv(
    summary_tbl %>% filter(condition == cond),
    file.path(cond_dir, paste0("fsf_summary_", cond, ".tsv"))
  )
}

cat("Saved condition-wise FSF outputs to:\n")
cat(OUT_DIR, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

