#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 04: Generate Synthetic Effect Estimates
#
# Purpose:
#   Generate controlled synthetic benchmark dataset
#
# Output:
#   data/processed/effect_estimates.tsv
#   results/synthetic/synthetic_truth_table.tsv
#   results/logs/04_generate_synthetic_effect_estimates.log
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

OUT_DATA <- file.path(
  ROOT_DIR,
  "data",
  "processed",
  "effect_estimates.tsv"
)

OUT_TRUTH <- file.path(
  ROOT_DIR,
  "results",
  "synthetic",
  "synthetic_truth_table.tsv"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results",
  "logs",
  "04_generate_synthetic_effect_estimates.log"
)

dir.create(dirname(OUT_DATA), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_TRUTH), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("====================================================\n")
cat("FSF Synthetic Benchmark Generator\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("====================================================\n\n")

set.seed(12345)

n_per_class <- 1000
n_perturb <- 100

cat("Features per class:", n_per_class, "\n")
cat("Perturbations:", n_perturb, "\n\n")

# ----------------------------------------------------
# Stable Up
# ----------------------------------------------------

stable_up <- expand.grid(
  gene = paste0("SU_", seq_len(n_per_class)),
  perturbation = seq_len(n_perturb)
) %>%
  mutate(
    effect_estimate = rnorm(
      n(),
      mean = 1.5,
      sd = 0.15
    ),
    truth_class = "stable_up"
  )

# ----------------------------------------------------
# Stable Down
# ----------------------------------------------------

stable_down <- expand.grid(
  gene = paste0("SD_", seq_len(n_per_class)),
  perturbation = seq_len(n_perturb)
) %>%
  mutate(
    effect_estimate = rnorm(
      n(),
      mean = -1.5,
      sd = 0.15
    ),
    truth_class = "stable_down"
  )

# ----------------------------------------------------
# Stable Constant
# ----------------------------------------------------

stable_const <- expand.grid(
  gene = paste0("SC_", seq_len(n_per_class)),
  perturbation = seq_len(n_perturb)
) %>%
  mutate(
    effect_estimate = rnorm(
      n(),
      mean = 0,
      sd = 0.15
    ),
    truth_class = "stable_constant"
  )

# ----------------------------------------------------
# Weakly Stable Transitional
# ----------------------------------------------------

weakly_stable <- expand.grid(
  gene = paste0("WT_", seq_len(n_per_class)),
  perturbation = seq_len(n_perturb)
) %>%
  mutate(
    effect_estimate = rnorm(
      n(),
      mean = 0.8,
      sd = 0.7
    ),
    truth_class = "weakly_stable"
  )

# ----------------------------------------------------
# Instable
# ----------------------------------------------------

instable <- expand.grid(
  gene = paste0("IN_", seq_len(n_per_class)),
  perturbation = seq_len(n_perturb)
) %>%
  mutate(
    effect_estimate = sample(
      c(
        rnorm(n(), 1.5, 0.2),
        rnorm(n(), -1.5, 0.2),
        rnorm(n(), 0, 0.2)
      ),
      size = n(),
      replace = TRUE
    ),
    truth_class = "instable"
  )

# ----------------------------------------------------
# Merge
# ----------------------------------------------------

synthetic <- bind_rows(
  stable_up,
  stable_down,
  stable_const,
  weakly_stable,
  instable
)

effect_table <- synthetic %>%
  transmute(
    feature_id = gene,
    perturbation_id = paste0("P", perturbation),
    effect_estimate
  )

truth_table <- synthetic %>%
  distinct(
    gene,
    truth_class
  ) %>%
  rename(
    feature_id = gene
  )

write_tsv(
  effect_table,
  OUT_DATA
)

write_tsv(
  truth_table,
  OUT_TRUTH
)

cat("Synthetic rows:", nrow(effect_table), "\n")
cat("Features:", n_distinct(effect_table$feature_id), "\n\n")

cat("Saved:\n")
cat(OUT_DATA, "\n")
cat(OUT_TRUTH, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")

sink()

