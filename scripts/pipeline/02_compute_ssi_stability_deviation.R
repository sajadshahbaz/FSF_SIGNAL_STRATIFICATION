#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 02: Compute SSI + Stability Deviation
#
# Input:
#   results/signal_states/feature_signal_states.tsv
#
# Output:
#   results/fsf_metrics/fsf_ssi_stability_deviation.tsv
#   results/fsf_metrics/fsf_signal_class_summary.tsv
#   results/logs/02_compute_ssi_stability_deviation.log
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

INPUT_FILE <- file.path(ROOT_DIR, "results", "signal_states", "feature_signal_states.tsv")

OUT_DIR <- file.path(ROOT_DIR, "results", "fsf_metrics")
LOG_DIR <- file.path(ROOT_DIR, "results", "logs")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

OUTPUT_FILE <- file.path(OUT_DIR, "fsf_ssi_stability_deviation.tsv")
SUMMARY_FILE <- file.path(OUT_DIR, "fsf_signal_class_summary.tsv")
LOG_FILE <- file.path(LOG_DIR, "02_compute_ssi_stability_deviation.log")

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 02: Compute SSI + Stability Deviation\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

highly_stable_cutoff <- 0.90
stable_cutoff <- 0.75
weakly_stable_cutoff <- 0.60

cat("Project root:", ROOT_DIR, "\n")
cat("Input file:", INPUT_FILE, "\n")
cat("Highly stable cutoff:", highly_stable_cutoff, "\n")
cat("Stable cutoff:", stable_cutoff, "\n")
cat("Weakly stable / transitional cutoff:", weakly_stable_cutoff, "\n\n")

if (!file.exists(INPUT_FILE)) {
  stop("Input file not found: ", INPUT_FILE)
}

states <- read_tsv(INPUT_FILE, show_col_types = FALSE)

required_cols <- c("feature_id", "perturbation_id", "signal_state")
missing_cols <- setdiff(required_cols, colnames(states))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

states <- states %>%
  mutate(
    feature_id = as.character(feature_id),
    perturbation_id = as.character(perturbation_id),
    signal_state = as.character(signal_state)
  )

allowed_states <- c("up", "down", "const")
bad_states <- setdiff(unique(states$signal_state), allowed_states)

if (length(bad_states) > 0) {
  stop("Invalid signal_state values: ", paste(bad_states, collapse = ", "))
}

cat("Input rows:", nrow(states), "\n")
cat("Unique features:", n_distinct(states$feature_id), "\n")
cat("Unique perturbations:", n_distinct(states$perturbation_id), "\n\n")

fsf_probs <- states %>%
  count(feature_id, signal_state, name = "n_state") %>%
  group_by(feature_id) %>%
  mutate(
    n_perturbations = sum(n_state),
    probability = n_state / n_perturbations
  ) %>%
  ungroup() %>%
  select(feature_id, n_perturbations, signal_state, probability) %>%
  pivot_wider(
    names_from = signal_state,
    values_from = probability,
    values_fill = 0,
    names_prefix = "p_"
  )

for (col in c("p_up", "p_down", "p_const")) {
  if (!col %in% colnames(fsf_probs)) {
    fsf_probs[[col]] <- 0
  }
}

fsf_results <- fsf_probs %>%
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
  arrange(desc(SSI), stability_deviation, feature_id)

summary_table <- fsf_results %>%
  count(stability_level, dominant_signal_identity, fsf_signal_class, name = "n_features") %>%
  arrange(stability_level, desc(n_features))

cat("FSF class summary:\n")
print(summary_table)
cat("\n")

write_tsv(fsf_results, OUTPUT_FILE)
write_tsv(summary_table, SUMMARY_FILE)

cat("Saved FSF metrics to:\n", OUTPUT_FILE, "\n\n")
cat("Saved class summary to:\n", SUMMARY_FILE, "\n\n")
cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

