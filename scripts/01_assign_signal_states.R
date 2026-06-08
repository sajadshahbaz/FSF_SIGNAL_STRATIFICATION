#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 01: Assign Signal States
#
# Input:
#   data/processed/effect_estimates.tsv
#
# Required columns:
#   feature_id
#   perturbation_id
#   effect_estimate
#
# Output:
#   results/signal_states/feature_signal_states.tsv
#   results/signal_states/signal_state_summary.tsv
#   results/logs/01_assign_signal_states.log
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# ----------------------------
# Locked project path
# ----------------------------

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

INPUT_FILE <- file.path(ROOT_DIR, "data", "processed", "effect_estimates.tsv")

OUT_DIR <- file.path(ROOT_DIR, "results", "signal_states")
LOG_DIR <- file.path(ROOT_DIR, "results", "logs")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

OUTPUT_FILE <- file.path(OUT_DIR, "feature_signal_states.tsv")
SUMMARY_FILE <- file.path(OUT_DIR, "signal_state_summary.tsv")
LOG_FILE <- file.path(LOG_DIR, "01_assign_signal_states.log")

# ----------------------------
# Start log
# ----------------------------

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 01: Assign Signal States\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

# ----------------------------
# Parameters
# ----------------------------

tau <- 0.5

cat("Project root:", ROOT_DIR, "\n")
cat("Input file:", INPUT_FILE, "\n")
cat("Effect threshold tau:", tau, "\n\n")

# ----------------------------
# Validate input
# ----------------------------

if (!file.exists(INPUT_FILE)) {
  stop(
    "Input file not found: ", INPUT_FILE, "\n",
    "Expected file with columns: feature_id, perturbation_id, effect_estimate"
  )
}

effect_df <- read_tsv(INPUT_FILE, show_col_types = FALSE)

required_cols <- c("feature_id", "perturbation_id", "effect_estimate")
missing_cols <- setdiff(required_cols, colnames(effect_df))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

effect_df <- effect_df %>%
  mutate(
    feature_id = as.character(feature_id),
    perturbation_id = as.character(perturbation_id),
    effect_estimate = as.numeric(effect_estimate)
  )

if (any(is.na(effect_df$effect_estimate))) {
  n_na <- sum(is.na(effect_df$effect_estimate))
  stop("effect_estimate contains NA/non-numeric values: ", n_na)
}

cat("Input rows:", nrow(effect_df), "\n")
cat("Unique features:", n_distinct(effect_df$feature_id), "\n")
cat("Unique perturbations:", n_distinct(effect_df$perturbation_id), "\n\n")

# ----------------------------
# Assign states
# ----------------------------

signal_states <- effect_df %>%
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
  ) %>%
  arrange(feature_id, perturbation_id)

# ----------------------------
# Summary
# ----------------------------

summary_table <- signal_states %>%
  count(signal_state, name = "n_assignments") %>%
  mutate(
    proportion = n_assignments / sum(n_assignments)
  ) %>%
  arrange(desc(n_assignments))

cat("Signal state summary:\n")
print(summary_table)
cat("\n")

# ----------------------------
# Save
# ----------------------------

write_tsv(signal_states, OUTPUT_FILE)
write_tsv(summary_table, SUMMARY_FILE)

cat("Saved signal states to:\n", OUTPUT_FILE, "\n\n")
cat("Saved summary to:\n", SUMMARY_FILE, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

