#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 08: Tau Sensitivity Analysis
#
# Purpose:
#   Evaluate how FSF classifications change across directional
#   effect thresholds tau.
#
# Input:
#   data/processed/effect_estimates.tsv
#
# Outputs:
#   results/synthetic/tau_sensitivity/tau_sensitivity_feature_metrics.tsv
#   results/synthetic/tau_sensitivity/tau_sensitivity_summary.tsv
#   results/synthetic/tau_sensitivity/tau_sensitivity_class_counts.tsv
#   results/figures/tau_sensitivity_mean_ssi.png
#   results/figures/tau_sensitivity_class_counts.png
#   results/logs/08_tau_sensitivity_analysis.log
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
  "data",
  "processed",
  "effect_estimates.tsv"
)

TRUTH_FILE <- file.path(
  ROOT_DIR,
  "results",
  "synthetic",
  "synthetic_truth_table.tsv"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results",
  "synthetic",
  "tau_sensitivity"
)

FIG_DIR <- file.path(ROOT_DIR, "results", "figures")
LOG_DIR <- file.path(ROOT_DIR, "results", "logs")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_FEATURE_METRICS <- file.path(OUT_DIR, "tau_sensitivity_feature_metrics.tsv")
OUT_SUMMARY <- file.path(OUT_DIR, "tau_sensitivity_summary.tsv")
OUT_CLASS_COUNTS <- file.path(OUT_DIR, "tau_sensitivity_class_counts.tsv")

FIG_MEAN_SSI <- file.path(FIG_DIR, "tau_sensitivity_mean_ssi.png")
FIG_CLASS_COUNTS <- file.path(FIG_DIR, "tau_sensitivity_class_counts.png")

LOG_FILE <- file.path(LOG_DIR, "08_tau_sensitivity_analysis.log")

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 08: Tau Sensitivity Analysis\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

tau_values <- c(0.25, 0.50, 0.75, 1.00)

highly_stable_cutoff <- 0.90
stable_cutoff <- 0.75
weakly_stable_cutoff <- 0.60

cat("Input file:", INPUT_FILE, "\n")
cat("Truth file:", TRUTH_FILE, "\n")
cat("Tau values:", paste(tau_values, collapse = ", "), "\n\n")

if (!file.exists(INPUT_FILE)) {
  stop("Missing input file: ", INPUT_FILE)
}

effects <- read_tsv(INPUT_FILE, show_col_types = FALSE)

required_cols <- c("feature_id", "perturbation_id", "effect_estimate")
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

cat("Rows:", nrow(effects), "\n")
cat("Features:", n_distinct(effects$feature_id), "\n")
cat("Perturbations:", n_distinct(effects$perturbation_id), "\n\n")

has_truth <- file.exists(TRUTH_FILE)

if (has_truth) {
  truth <- read_tsv(TRUTH_FILE, show_col_types = FALSE) %>%
    mutate(feature_id = as.character(feature_id))
  cat("Truth table loaded:", TRUTH_FILE, "\n\n")
} else {
  truth <- NULL
  cat("Truth table not found. Proceeding without truth-class summaries.\n\n")
}

assign_state <- function(x, tau) {
  case_when(
    x > tau ~ "up",
    x < -tau ~ "down",
    TRUE ~ "const"
  )
}

compute_fsf_for_tau <- function(effect_df, tau) {
  states <- effect_df %>%
    mutate(
      tau = tau,
      signal_state = assign_state(effect_estimate, tau)
    )

  probs <- states %>%
    count(tau, feature_id, signal_state, name = "n_state") %>%
    group_by(tau, feature_id) %>%
    mutate(
      n_perturbations = sum(n_state),
      probability = n_state / n_perturbations
    ) %>%
    ungroup() %>%
    select(tau, feature_id, n_perturbations, signal_state, probability) %>%
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

  probs %>%
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
    ungroup()
}

feature_metrics <- bind_rows(
  lapply(
    tau_values,
    function(tau) compute_fsf_for_tau(effects, tau)
  )
)

if (has_truth) {
  feature_metrics <- feature_metrics %>%
    left_join(truth, by = "feature_id")
} else {
  feature_metrics <- feature_metrics %>%
    mutate(truth_class = NA_character_)
}

write_tsv(feature_metrics, OUT_FEATURE_METRICS)

cat("Saved feature metrics:", OUT_FEATURE_METRICS, "\n")

summary_tbl <- feature_metrics %>%
  group_by(tau, truth_class) %>%
  summarise(
    n_features = n(),
    mean_ssi = mean(SSI),
    median_ssi = median(SSI),
    mean_sd = mean(stability_deviation),
    median_sd = median(stability_deviation),
    .groups = "drop"
  ) %>%
  arrange(tau, truth_class)

write_tsv(summary_tbl, OUT_SUMMARY)

cat("Saved summary:", OUT_SUMMARY, "\n")

class_counts <- feature_metrics %>%
  count(tau, truth_class, fsf_signal_class, name = "n_features") %>%
  group_by(tau, truth_class) %>%
  mutate(proportion = n_features / sum(n_features)) %>%
  ungroup() %>%
  arrange(tau, truth_class, desc(n_features))

write_tsv(class_counts, OUT_CLASS_COUNTS)

cat("Saved class counts:", OUT_CLASS_COUNTS, "\n")

# ------------------------------------------------------------
# Plot 1: Mean SSI across tau
# ------------------------------------------------------------

p1 <- ggplot(
  summary_tbl,
  aes(
    x = tau,
    y = mean_ssi,
    group = truth_class,
    linetype = truth_class
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(
    title = "Tau sensitivity: mean SSI",
    x = "Directional effect threshold (tau)",
    y = "Mean SSI",
    linetype = "Truth class"
  ) +
  theme_bw(base_size = 12)

ggsave(
  FIG_MEAN_SSI,
  p1,
  width = 8,
  height = 5,
  dpi = 300
)

cat("Saved figure:", FIG_MEAN_SSI, "\n")

# ------------------------------------------------------------
# Plot 2: Class composition across tau
# ------------------------------------------------------------

p2 <- ggplot(
  class_counts,
  aes(
    x = factor(tau),
    y = proportion,
    fill = fsf_signal_class
  )
) +
  geom_col() +
  facet_wrap(~ truth_class, nrow = 2) +
  labs(
    title = "Tau sensitivity: FSF class composition",
    x = "Directional effect threshold (tau)",
    y = "Proportion of features",
    fill = "FSF class"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  FIG_CLASS_COUNTS,
  p2,
  width = 12,
  height = 7,
  dpi = 300
)

cat("Saved figure:", FIG_CLASS_COUNTS, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

