#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 18: Run FSF on Real Effect Estimates
#
# Input:
#   data/processed/effect_estimates_real.tsv
#
# Outputs:
#   results/real_data/signal_states/real_feature_signal_states.tsv
#   results/real_data/signal_states/real_signal_state_summary.tsv
#   results/real_data/fsf_metrics/real_fsf_ssi_stability_deviation.tsv
#   results/real_data/fsf_metrics/real_fsf_signal_class_summary.tsv
#   results/logs/18_run_fsf_on_real_effect_estimates.log
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

INPUT_FILE <- file.path(
  ROOT_DIR,
  "data/processed/effect_estimates_real.tsv"
)

OUT_STATES_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/signal_states"
)

OUT_METRICS_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/fsf_metrics"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/18_run_fsf_on_real_effect_estimates.log"
)

dir.create(OUT_STATES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_METRICS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 18: Run FSF on Real Effect Estimates\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

tau <- 0.5

highly_stable_cutoff <- 0.90
stable_cutoff <- 0.75
weakly_stable_cutoff <- 0.60

cat("Input file:", INPUT_FILE, "\n")
cat("Tau:", tau, "\n")
cat("Highly stable cutoff:", highly_stable_cutoff, "\n")
cat("Stable cutoff:", stable_cutoff, "\n")
cat("Weakly stable/transitional cutoff:", weakly_stable_cutoff, "\n\n")

if (!file.exists(INPUT_FILE)) {
  stop("Missing real effect estimate input file: ", INPUT_FILE)
}

effects <- read_tsv(INPUT_FILE, show_col_types = FALSE)

required_cols <- c(
  "feature_id",
  "perturbation_id",
  "effect_estimate"
)

missing_cols <- setdiff(required_cols, colnames(effects))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

effects <- effects %>%
  mutate(
    feature_id = as.character(feature_id),
    perturbation_id = as.character(perturbation_id),
    effect_estimate = as.numeric(effect_estimate)
  )

if (any(is.na(effects$effect_estimate))) {
  stop("NA values detected in effect_estimate.")
}

cat("Rows:", nrow(effects), "\n")
cat("Unique features:", n_distinct(effects$feature_id), "\n")
cat("Unique perturbations:", n_distinct(effects$perturbation_id), "\n\n")

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
  ) %>%
  arrange(feature_id, perturbation_id)

state_summary <- states %>%
  count(signal_state, name = "n_assignments") %>%
  mutate(proportion = n_assignments / sum(n_assignments)) %>%
  arrange(desc(n_assignments))

write_tsv(
  states,
  file.path(OUT_STATES_DIR, "real_feature_signal_states.tsv")
)

write_tsv(
  state_summary,
  file.path(OUT_STATES_DIR, "real_signal_state_summary.tsv")
)

cat("Signal state summary:\n")
print(state_summary)
cat("\n")

# ------------------------------------------------------------
# Compute probabilities
# ------------------------------------------------------------

probs <- states %>%
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
  if (!col %in% colnames(probs)) {
    probs[[col]] <- 0
  }
}

# ------------------------------------------------------------
# Compute SSI + Stability Deviation
# ------------------------------------------------------------

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

fsf_summary <- fsf %>%
  count(
    stability_level,
    dominant_signal_identity,
    fsf_signal_class,
    name = "n_features"
  ) %>%
  mutate(
    proportion = n_features / sum(n_features)
  ) %>%
  arrange(stability_level, desc(n_features))

write_tsv(
  fsf,
  file.path(OUT_METRICS_DIR, "real_fsf_ssi_stability_deviation.tsv")
)

write_tsv(
  fsf_summary,
  file.path(OUT_METRICS_DIR, "real_fsf_signal_class_summary.tsv")
)

cat("FSF class summary:\n")
print(fsf_summary, n = 100)
cat("\n")

cat("Saved outputs:\n")
cat(file.path(OUT_STATES_DIR, "real_feature_signal_states.tsv"), "\n")
cat(file.path(OUT_STATES_DIR, "real_signal_state_summary.tsv"), "\n")
cat(file.path(OUT_METRICS_DIR, "real_fsf_ssi_stability_deviation.tsv"), "\n")
cat(file.path(OUT_METRICS_DIR, "real_fsf_signal_class_summary.tsv"), "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

