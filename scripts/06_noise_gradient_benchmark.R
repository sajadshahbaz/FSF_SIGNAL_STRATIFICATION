#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 06: Noise-Gradient Synthetic Benchmark
#
# Purpose:
#   Evaluate SSI and Stability Deviation behavior under
#   increasing effect-estimate noise.
#
# Outputs:
#   results/synthetic/noise_gradient/noise_gradient_feature_metrics.tsv
#   results/synthetic/noise_gradient/noise_gradient_summary.tsv
#   results/synthetic/noise_gradient/noise_gradient_confusion.tsv
#   results/figures/noise_gradient_mean_ssi.png
#   results/figures/noise_gradient_mean_sd.png
#   results/figures/noise_gradient_class_recovery.png
#   results/logs/06_noise_gradient_benchmark.log
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

OUT_DIR <- file.path(ROOT_DIR, "results", "synthetic", "noise_gradient")
FIG_DIR <- file.path(ROOT_DIR, "results", "figures")
LOG_DIR <- file.path(ROOT_DIR, "results", "logs")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_FEATURE_METRICS <- file.path(OUT_DIR, "noise_gradient_feature_metrics.tsv")
OUT_SUMMARY <- file.path(OUT_DIR, "noise_gradient_summary.tsv")
OUT_CONFUSION <- file.path(OUT_DIR, "noise_gradient_confusion.tsv")

FIG_MEAN_SSI <- file.path(FIG_DIR, "noise_gradient_mean_ssi.png")
FIG_MEAN_SD <- file.path(FIG_DIR, "noise_gradient_mean_sd.png")
FIG_RECOVERY <- file.path(FIG_DIR, "noise_gradient_class_recovery.png")

LOG_FILE <- file.path(LOG_DIR, "06_noise_gradient_benchmark.log")

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 06: Noise-Gradient Synthetic Benchmark\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

set.seed(12345)

# ----------------------------
# Parameters
# ----------------------------

n_per_class <- 500
n_perturb <- 100
tau <- 0.5

noise_levels <- c(0.05, 0.10, 0.20, 0.35, 0.50, 0.75, 1.00)

highly_stable_cutoff <- 0.90
stable_cutoff <- 0.75
weakly_stable_cutoff <- 0.60

cat("Features per class:", n_per_class, "\n")
cat("Perturbations:", n_perturb, "\n")
cat("Tau:", tau, "\n")
cat("Noise levels:", paste(noise_levels, collapse = ", "), "\n\n")

# ----------------------------
# Helper: assign states
# ----------------------------

assign_state <- function(x, tau = 0.5) {
  dplyr::case_when(
    x > tau ~ "up",
    x < -tau ~ "down",
    TRUE ~ "const"
  )
}

# ----------------------------
# Helper: compute FSF metrics
# ----------------------------

compute_fsf <- function(effect_df) {
  states <- effect_df %>%
    mutate(signal_state = assign_state(effect_estimate, tau = tau))

  probs <- states %>%
    count(noise_sd, truth_class, feature_id, signal_state, name = "n_state") %>%
    group_by(noise_sd, truth_class, feature_id) %>%
    mutate(
      n_perturbations = sum(n_state),
      probability = n_state / n_perturbations
    ) %>%
    ungroup() %>%
    select(noise_sd, truth_class, feature_id, n_perturbations, signal_state, probability) %>%
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

# ----------------------------
# Generate synthetic data
# ----------------------------

all_effects <- list()

for (noise_sd in noise_levels) {

  stable_up <- expand.grid(
    feature_index = seq_len(n_per_class),
    perturbation = seq_len(n_perturb)
  ) %>%
    mutate(
      feature_id = paste0("SU_", feature_index),
      truth_class = "stable_up",
      noise_sd = noise_sd,
      effect_estimate = rnorm(n(), mean = 1.5, sd = noise_sd)
    )

  stable_down <- expand.grid(
    feature_index = seq_len(n_per_class),
    perturbation = seq_len(n_perturb)
  ) %>%
    mutate(
      feature_id = paste0("SD_", feature_index),
      truth_class = "stable_down",
      noise_sd = noise_sd,
      effect_estimate = rnorm(n(), mean = -1.5, sd = noise_sd)
    )

  stable_constant <- expand.grid(
    feature_index = seq_len(n_per_class),
    perturbation = seq_len(n_perturb)
  ) %>%
    mutate(
      feature_id = paste0("SC_", feature_index),
      truth_class = "stable_constant",
      noise_sd = noise_sd,
      effect_estimate = rnorm(n(), mean = 0, sd = noise_sd)
    )

  weakly_stable <- expand.grid(
    feature_index = seq_len(n_per_class),
    perturbation = seq_len(n_perturb)
  ) %>%
    mutate(
      feature_id = paste0("WT_", feature_index),
      truth_class = "weakly_stable",
      noise_sd = noise_sd,
      effect_estimate = rnorm(n(), mean = 0.8, sd = noise_sd)
    )

  instable <- expand.grid(
    feature_index = seq_len(n_per_class),
    perturbation = seq_len(n_perturb)
  ) %>%
    mutate(
      feature_id = paste0("IN_", feature_index),
      truth_class = "instable",
      noise_sd = noise_sd,
      hidden_state = sample(c("up", "down", "const"), size = n(), replace = TRUE),
      effect_estimate = case_when(
        hidden_state == "up" ~ rnorm(n(), mean = 1.5, sd = noise_sd),
        hidden_state == "down" ~ rnorm(n(), mean = -1.5, sd = noise_sd),
        hidden_state == "const" ~ rnorm(n(), mean = 0, sd = noise_sd)
      )
    ) %>%
    select(-hidden_state)

  all_effects[[as.character(noise_sd)]] <- bind_rows(
    stable_up,
    stable_down,
    stable_constant,
    weakly_stable,
    instable
  ) %>%
    transmute(
      noise_sd,
      truth_class,
      feature_id = paste0(truth_class, "_", feature_id),
      perturbation_id = paste0("P", perturbation),
      effect_estimate
    )
}

effects <- bind_rows(all_effects)

cat("Generated rows:", nrow(effects), "\n")
cat("Generated features per noise level:", n_distinct(effects$feature_id), "\n\n")

# ----------------------------
# Compute FSF
# ----------------------------

feature_metrics <- compute_fsf(effects)

write_tsv(feature_metrics, OUT_FEATURE_METRICS)

cat("Saved feature metrics:", OUT_FEATURE_METRICS, "\n")

# ----------------------------
# Summary by truth class and noise
# ----------------------------

summary_tbl <- feature_metrics %>%
  group_by(noise_sd, truth_class) %>%
  summarise(
    n_features = n(),
    mean_ssi = mean(SSI),
    median_ssi = median(SSI),
    mean_sd = mean(stability_deviation),
    median_sd = median(stability_deviation),
    .groups = "drop"
  )

write_tsv(summary_tbl, OUT_SUMMARY)

cat("Saved summary:", OUT_SUMMARY, "\n")

# ----------------------------
# Confusion by noise level
# ----------------------------

confusion <- feature_metrics %>%
  count(noise_sd, truth_class, fsf_signal_class, name = "n_features") %>%
  group_by(noise_sd, truth_class) %>%
  mutate(proportion = n_features / sum(n_features)) %>%
  ungroup()

write_tsv(confusion, OUT_CONFUSION)

cat("Saved confusion:", OUT_CONFUSION, "\n")

# ----------------------------
# Plot: mean SSI
# ----------------------------

p1 <- ggplot(summary_tbl, aes(x = noise_sd, y = mean_ssi, group = truth_class, linetype = truth_class)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(
    title = "Noise-gradient benchmark: mean SSI",
    x = "Noise standard deviation",
    y = "Mean SSI",
    linetype = "Truth class"
  ) +
  theme_bw(base_size = 12)

ggsave(FIG_MEAN_SSI, p1, width = 8, height = 5, dpi = 300)

cat("Saved figure:", FIG_MEAN_SSI, "\n")

# ----------------------------
# Plot: mean Stability Deviation
# ----------------------------

p2 <- ggplot(summary_tbl, aes(x = noise_sd, y = mean_sd, group = truth_class, linetype = truth_class)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(
    title = "Noise-gradient benchmark: mean Stability Deviation",
    x = "Noise standard deviation",
    y = "Mean Stability Deviation",
    linetype = "Truth class"
  ) +
  theme_bw(base_size = 12)

ggsave(FIG_MEAN_SD, p2, width = 8, height = 5, dpi = 300)

cat("Saved figure:", FIG_MEAN_SD, "\n")

# ----------------------------
# Plot: class recovery heatmap
# ----------------------------

p3 <- ggplot(confusion, aes(x = fsf_signal_class, y = truth_class, fill = proportion)) +
  geom_tile() +
  facet_wrap(~ noise_sd, nrow = 2) +
  labs(
    title = "Noise-gradient benchmark: FSF class recovery",
    x = "FSF class",
    y = "Truth class",
    fill = "Proportion"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(FIG_RECOVERY, p3, width = 13, height = 7, dpi = 300)

cat("Saved figure:", FIG_RECOVERY, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

