#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 07: Publication Benchmark Summary
#
# Inputs:
#   results/synthetic/noise_gradient/noise_gradient_summary.tsv
#   results/synthetic/noise_gradient/noise_gradient_confusion.tsv
#
# Outputs:
#   results/tables/publication_noise_gradient_summary.tsv
#   results/tables/publication_noise_gradient_dominant_recovery.tsv
#   results/logs/07_make_publication_benchmark_summary.log
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

SUMMARY_FILE <- file.path(
  ROOT_DIR,
  "results",
  "synthetic",
  "noise_gradient",
  "noise_gradient_summary.tsv"
)

CONFUSION_FILE <- file.path(
  ROOT_DIR,
  "results",
  "synthetic",
  "noise_gradient",
  "noise_gradient_confusion.tsv"
)

OUT_TABLE_1 <- file.path(
  ROOT_DIR,
  "results",
  "tables",
  "publication_noise_gradient_summary.tsv"
)

OUT_TABLE_2 <- file.path(
  ROOT_DIR,
  "results",
  "tables",
  "publication_noise_gradient_dominant_recovery.tsv"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results",
  "logs",
  "07_make_publication_benchmark_summary.log"
)

dir.create(dirname(OUT_TABLE_1), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 07: Publication Benchmark Summary\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

if (!file.exists(SUMMARY_FILE)) {
  stop("Missing summary file: ", SUMMARY_FILE)
}

if (!file.exists(CONFUSION_FILE)) {
  stop("Missing confusion file: ", CONFUSION_FILE)
}

summary_tbl <- read_tsv(SUMMARY_FILE, show_col_types = FALSE)
confusion_tbl <- read_tsv(CONFUSION_FILE, show_col_types = FALSE)

# ------------------------------------------------------------
# Table 1: Rounded benchmark summary
# ------------------------------------------------------------

publication_summary <- summary_tbl %>%
  mutate(
    noise_sd = round(noise_sd, 2),
    mean_ssi = round(mean_ssi, 3),
    median_ssi = round(median_ssi, 3),
    mean_stability_deviation = round(mean_sd, 3),
    median_stability_deviation = round(median_sd, 3)
  ) %>%
  select(
    noise_sd,
    truth_class,
    n_features,
    mean_ssi,
    median_ssi,
    mean_stability_deviation,
    median_stability_deviation
  ) %>%
  arrange(noise_sd, truth_class)

write_tsv(publication_summary, OUT_TABLE_1)

# ------------------------------------------------------------
# Table 2: Dominant recovered FSF class per truth class/noise
# ------------------------------------------------------------

dominant_recovery <- confusion_tbl %>%
  group_by(noise_sd, truth_class) %>%
  arrange(desc(proportion), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    noise_sd = round(noise_sd, 2),
    recovery_proportion = round(proportion, 3)
  ) %>%
  select(
    noise_sd,
    truth_class,
    dominant_recovered_fsf_class = fsf_signal_class,
    recovery_proportion,
    n_features
  ) %>%
  arrange(noise_sd, truth_class)

write_tsv(dominant_recovery, OUT_TABLE_2)

cat("Publication summary table:\n")
print(publication_summary)
cat("\n")

cat("Dominant recovery table:\n")
print(dominant_recovery)
cat("\n")

cat("Saved:\n")
cat(OUT_TABLE_1, "\n")
cat(OUT_TABLE_2, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

