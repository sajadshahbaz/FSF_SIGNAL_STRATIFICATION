#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 09: Final Benchmark Report Table
#
# Inputs:
#   results/synthetic/synthetic_recovery_confusion_matrix.tsv
#   results/synthetic/synthetic_recovery_summary.tsv
#   results/tables/publication_noise_gradient_summary.tsv
#   results/tables/publication_noise_gradient_dominant_recovery.tsv
#   results/synthetic/tau_sensitivity/tau_sensitivity_summary.tsv
#   results/synthetic/tau_sensitivity/tau_sensitivity_class_counts.tsv
#
# Outputs:
#   results/tables/FSF_v1_final_benchmark_report.tsv
#   results/tables/FSF_v1_key_results_for_manuscript.tsv
#   results/logs/09_make_final_benchmark_report_table.log
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

RECOVERY_CONF <- file.path(
  ROOT_DIR,
  "results",
  "synthetic",
  "synthetic_recovery_confusion_matrix.tsv"
)

RECOVERY_SUMMARY <- file.path(
  ROOT_DIR,
  "results",
  "synthetic",
  "synthetic_recovery_summary.tsv"
)

NOISE_SUMMARY <- file.path(
  ROOT_DIR,
  "results",
  "tables",
  "publication_noise_gradient_summary.tsv"
)

NOISE_RECOVERY <- file.path(
  ROOT_DIR,
  "results",
  "tables",
  "publication_noise_gradient_dominant_recovery.tsv"
)

TAU_SUMMARY <- file.path(
  ROOT_DIR,
  "results",
  "synthetic",
  "tau_sensitivity",
  "tau_sensitivity_summary.tsv"
)

TAU_COUNTS <- file.path(
  ROOT_DIR,
  "results",
  "synthetic",
  "tau_sensitivity",
  "tau_sensitivity_class_counts.tsv"
)

OUT_FINAL <- file.path(
  ROOT_DIR,
  "results",
  "tables",
  "FSF_v1_final_benchmark_report.tsv"
)

OUT_KEY <- file.path(
  ROOT_DIR,
  "results",
  "tables",
  "FSF_v1_key_results_for_manuscript.tsv"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results",
  "logs",
  "09_make_final_benchmark_report_table.log"
)

dir.create(dirname(OUT_FINAL), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 09: Final Benchmark Report Table\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

required_files <- c(
  RECOVERY_CONF,
  RECOVERY_SUMMARY,
  NOISE_SUMMARY,
  NOISE_RECOVERY,
  TAU_SUMMARY,
  TAU_COUNTS
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required files:\n",
    paste(missing_files, collapse = "\n")
  )
}

recovery_conf <- read_tsv(RECOVERY_CONF, show_col_types = FALSE)
recovery_summary <- read_tsv(RECOVERY_SUMMARY, show_col_types = FALSE)
noise_summary <- read_tsv(NOISE_SUMMARY, show_col_types = FALSE)
noise_recovery <- read_tsv(NOISE_RECOVERY, show_col_types = FALSE)
tau_summary <- read_tsv(TAU_SUMMARY, show_col_types = FALSE)
tau_counts <- read_tsv(TAU_COUNTS, show_col_types = FALSE)

cat("Loaded benchmark files successfully.\n\n")

# ------------------------------------------------------------
# Section 1: Baseline synthetic recovery
# ------------------------------------------------------------

baseline_recovery <- recovery_conf %>%
  group_by(truth_class) %>%
  mutate(
    total_truth_features = sum(n_features),
    recovery_proportion = n_features / total_truth_features
  ) %>%
  arrange(truth_class, desc(recovery_proportion)) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    benchmark = "baseline_synthetic_recovery",
    setting = "default_tau_0.5",
    truth_class,
    key_metric = "dominant_recovered_class",
    value = fsf_signal_class,
    numeric_value = round(recovery_proportion, 3),
    interpretation = paste0(
      truth_class,
      " dominantly recovered as ",
      fsf_signal_class,
      " with proportion ",
      round(recovery_proportion, 3)
    )
  )

# ------------------------------------------------------------
# Section 2: Baseline SSI / SD
# ------------------------------------------------------------

baseline_ssi <- recovery_summary %>%
  transmute(
    benchmark = "baseline_synthetic_recovery",
    setting = "default_tau_0.5",
    truth_class,
    key_metric = "mean_ssi",
    value = as.character(round(mean_ssi, 3)),
    numeric_value = round(mean_ssi, 3),
    interpretation = paste0(
      truth_class,
      " mean SSI = ",
      round(mean_ssi, 3),
      "; mean Stability Deviation = ",
      round(mean_sd, 3)
    )
  )

# ------------------------------------------------------------
# Section 3: Noise gradient low vs high noise
# ------------------------------------------------------------

noise_low_high <- noise_summary %>%
  filter(noise_sd %in% c(0.05, 1.00)) %>%
  transmute(
    benchmark = "noise_gradient",
    setting = paste0("noise_sd_", noise_sd),
    truth_class,
    key_metric = "mean_ssi_and_sd",
    value = paste0(
      "SSI=", mean_ssi,
      "; SD=", mean_stability_deviation
    ),
    numeric_value = mean_ssi,
    interpretation = paste0(
      truth_class,
      " at noise_sd=",
      noise_sd,
      " had mean SSI=",
      mean_ssi,
      " and mean Stability Deviation=",
      mean_stability_deviation
    )
  )

noise_recovery_low_high <- noise_recovery %>%
  filter(noise_sd %in% c(0.05, 1.00)) %>%
  transmute(
    benchmark = "noise_gradient",
    setting = paste0("noise_sd_", noise_sd),
    truth_class,
    key_metric = "dominant_recovered_class",
    value = dominant_recovered_fsf_class,
    numeric_value = recovery_proportion,
    interpretation = paste0(
      truth_class,
      " at noise_sd=",
      noise_sd,
      " dominantly recovered as ",
      dominant_recovered_fsf_class,
      " with proportion ",
      recovery_proportion
    )
  )

# ------------------------------------------------------------
# Section 4: Tau sensitivity default tau result
# ------------------------------------------------------------

tau_default_summary <- tau_summary %>%
  filter(abs(tau - 0.5) < 1e-9) %>%
  transmute(
    benchmark = "tau_sensitivity",
    setting = "tau_0.5",
    truth_class,
    key_metric = "mean_ssi_and_sd",
    value = paste0(
      "SSI=",
      round(mean_ssi, 3),
      "; SD=",
      round(mean_sd, 3)
    ),
    numeric_value = round(mean_ssi, 3),
    interpretation = paste0(
      truth_class,
      " at tau=0.5 had mean SSI=",
      round(mean_ssi, 3),
      " and mean Stability Deviation=",
      round(mean_sd, 3)
    )
  )

tau_default_recovery <- tau_counts %>%
  filter(abs(tau - 0.5) < 1e-9) %>%
  group_by(truth_class) %>%
  arrange(desc(proportion), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    benchmark = "tau_sensitivity",
    setting = "tau_0.5",
    truth_class,
    key_metric = "dominant_recovered_class",
    value = fsf_signal_class,
    numeric_value = round(proportion, 3),
    interpretation = paste0(
      truth_class,
      " at tau=0.5 dominantly recovered as ",
      fsf_signal_class,
      " with proportion ",
      round(proportion, 3)
    )
  )

# ------------------------------------------------------------
# Merge final report
# ------------------------------------------------------------

final_report <- bind_rows(
  baseline_recovery,
  baseline_ssi,
  noise_low_high,
  noise_recovery_low_high,
  tau_default_summary,
  tau_default_recovery
) %>%
  arrange(benchmark, setting, truth_class, key_metric)

write_tsv(final_report, OUT_FINAL)

# ------------------------------------------------------------
# Key manuscript statements table
# ------------------------------------------------------------

key_results <- tibble::tibble(
  result_id = c(
    "R1",
    "R2",
    "R3",
    "R4",
    "R5",
    "R6"
  ),
  manuscript_use = c(
    "Baseline synthetic recovery",
    "Baseline synthetic recovery",
    "Noise-gradient robustness",
    "Noise-gradient robustness",
    "Tau sensitivity",
    "Interpretation of Stability Deviation"
  ),
  result_statement = c(
    "FSF recovered all highly stable directional and constant synthetic classes under default settings.",
    "The deliberately weakly stable synthetic class was primarily assigned to the weakly stable/transitional stratum, with minor redistribution into adjacent stable or instable classes.",
    "Stable Up and Stable Down signals retained high SSI values under increasing noise.",
    "Instable synthetic signals remained consistently low-SSI and were recovered as instable across all tested noise levels.",
    "Tau sensitivity analysis supported tau = 0.5 as a balanced default threshold for separating stable, transitional, and instable synthetic classes.",
    "Stability Deviation increased as signal identity weakened, supporting its use as a continuous measure of departure from ideal stability."
  ),
  supporting_output = c(
    "synthetic_recovery_confusion_matrix.tsv",
    "synthetic_recovery_confusion_matrix.tsv",
    "publication_noise_gradient_summary.tsv",
    "publication_noise_gradient_dominant_recovery.tsv",
    "tau_sensitivity_summary.tsv; tau_sensitivity_class_counts.tsv",
    "synthetic_recovery_summary.tsv; noise_gradient_summary.tsv"
  )
)

write_tsv(key_results, OUT_KEY)

cat("Saved final benchmark report:\n")
cat(OUT_FINAL, "\n\n")

cat("Saved key manuscript results:\n")
cat(OUT_KEY, "\n\n")

cat("Preview final report:\n")
print(final_report, n = 40)

cat("\nPreview key manuscript results:\n")
print(key_results)

cat("\nFinished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

