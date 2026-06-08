#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 17: Generate Real Effect Estimates
#
# Effect estimate:
#   E(i,p) = mean(VST_treatment) - mean(VST_control)
#
# Outputs:
#   results/real_data/effect_estimates/real_effect_estimates.tsv
#   data/processed/effect_estimates_real.tsv
#   results/logs/17_generate_real_effect_estimates.log
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

VST_FILE <- "/media/saji/5E06441D0643F5152/phase_2/results/results/vst_counts_filtered.tsv"

CONTRAST_MAP_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/design/real_contrast_map.tsv"
)

SAMPLE_MAP_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/design/real_contrast_sample_map.tsv"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/effect_estimates"
)

OUT_EFFECTS <- file.path(
  OUT_DIR,
  "real_effect_estimates.tsv"
)

OUT_EFFECTS_PROCESSED <- file.path(
  ROOT_DIR,
  "data/processed/effect_estimates_real.tsv"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/17_generate_real_effect_estimates.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_EFFECTS_PROCESSED), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 17: Generate Real Effect Estimates\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

vst <- read_tsv(VST_FILE, show_col_types = FALSE)
contrast_map <- read_tsv(CONTRAST_MAP_FILE, show_col_types = FALSE)
sample_map <- read_tsv(SAMPLE_MAP_FILE, show_col_types = FALSE)

if (!"GeneID" %in% colnames(vst)) {
  stop("VST file must contain GeneID column.")
}

required_sample_cols <- c("contrast_id", "sample_id", "role")
missing_sample_cols <- setdiff(required_sample_cols, colnames(sample_map))

if (length(missing_sample_cols) > 0) {
  stop("Missing sample map columns: ", paste(missing_sample_cols, collapse = ", "))
}

vst_samples <- setdiff(colnames(vst), "GeneID")

bad_samples <- setdiff(unique(sample_map$sample_id), vst_samples)

if (length(bad_samples) > 0) {
  stop("Samples in contrast sample map missing from VST: ", paste(bad_samples, collapse = ", "))
}

cat("Genes:", nrow(vst), "\n")
cat("Samples:", length(vst_samples), "\n")
cat("Contrasts:", nrow(contrast_map), "\n\n")

effect_list <- list()

for (cid in contrast_map$contrast_id) {

  sm <- sample_map %>%
    filter(contrast_id == cid)

  ctrl_samples <- sm %>%
    filter(role == "control") %>%
    pull(sample_id)

  trt_samples <- sm %>%
    filter(role == "treatment") %>%
    pull(sample_id)

  if (length(ctrl_samples) == 0 || length(trt_samples) == 0) {
    stop("Missing control or treatment samples for contrast: ", cid)
  }

  ctrl_mat <- as.matrix(vst[, ctrl_samples, drop = FALSE])
  trt_mat  <- as.matrix(vst[, trt_samples, drop = FALSE])

  ctrl_mean <- rowMeans(ctrl_mat, na.rm = TRUE)
  trt_mean  <- rowMeans(trt_mat, na.rm = TRUE)

  effect <- trt_mean - ctrl_mean

  meta_row <- contrast_map %>%
    filter(contrast_id == cid)

  effect_list[[cid]] <- tibble(
    feature_id = vst$GeneID,
    perturbation_id = cid,
    effect_estimate = effect,
    condition = meta_row$condition[1],
    baseline_block = meta_row$baseline_block[1],
    perturbation_type = meta_row$perturbation_type[1],
    timepoint = ifelse("timepoint" %in% colnames(meta_row), meta_row$timepoint[1], NA),
    n_control = length(ctrl_samples),
    n_treatment = length(trt_samples)
  )

  cat("Processed:", cid,
      "| control:", length(ctrl_samples),
      "| treatment:", length(trt_samples),
      "\n")
}

effects <- bind_rows(effect_list)

write_tsv(effects, OUT_EFFECTS)

# Also save FSF-compatible minimal input
effects_minimal <- effects %>%
  select(
    feature_id,
    perturbation_id,
    effect_estimate
  )

write_tsv(effects_minimal, OUT_EFFECTS_PROCESSED)

cat("\nSaved full real effect estimates:\n")
cat(OUT_EFFECTS, "\n\n")

cat("Saved FSF-compatible minimal real effect estimates:\n")
cat(OUT_EFFECTS_PROCESSED, "\n\n")

cat("Rows:", nrow(effects), "\n")
cat("Unique features:", dplyr::n_distinct(effects$feature_id), "\n")
cat("Unique perturbations:", dplyr::n_distinct(effects$perturbation_id), "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

