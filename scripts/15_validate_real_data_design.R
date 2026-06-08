#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 15: Validate Real-Data Design
#
# Purpose:
#   Validate sample matching between metadata and VST matrix,
#   and summarize biological design for FSF perturbation setup.
#
# Outputs:
#   results/real_data/input_validation/sample_matching_summary.tsv
#   results/real_data/input_validation/metadata_group_summary.tsv
#   results/real_data/input_validation/valid_real_metadata.tsv
#   results/logs/15_validate_real_data_design.log
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

METADATA_FILE <- "/media/saji/5E06441D0643F5152/phase_1/parse_meta/revised_nature_method/Simplified_Metadata_Table.csv"
VST_FILE <- "/media/saji/5E06441D0643F5152/phase_2/results/results/vst_counts_filtered.tsv"

OUT_DIR <- file.path(ROOT_DIR, "results", "real_data", "input_validation")
LOG_FILE <- file.path(ROOT_DIR, "results", "logs", "15_validate_real_data_design.log")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 15: Validate Real-Data Design\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

meta <- read_csv(METADATA_FILE, show_col_types = FALSE)
vst <- read_tsv(VST_FILE, show_col_types = FALSE)

required_meta_cols <- c(
  "sample_id",
  "condition",
  "group_type",
  "timepoint",
  "baseline_block"
)

missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))

if (length(missing_meta_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta_cols, collapse = ", "))
}

if (!"GeneID" %in% colnames(vst)) {
  stop("VST matrix must contain GeneID column.")
}

meta <- meta %>%
  mutate(
    sample_id = as.character(sample_id),
    condition = as.character(condition),
    group_type = as.character(group_type),
    timepoint = as.character(timepoint),
    baseline_block = as.character(baseline_block)
  )

vst_samples <- setdiff(colnames(vst), "GeneID")
meta_samples <- meta$sample_id

samples_in_both <- intersect(meta_samples, vst_samples)
metadata_not_in_vst <- setdiff(meta_samples, vst_samples)
vst_not_in_metadata <- setdiff(vst_samples, meta_samples)

sample_summary <- tibble::tibble(
  item = c(
    "metadata_samples",
    "vst_samples",
    "samples_in_both",
    "metadata_not_in_vst",
    "vst_not_in_metadata"
  ),
  n = c(
    length(meta_samples),
    length(vst_samples),
    length(samples_in_both),
    length(metadata_not_in_vst),
    length(vst_not_in_metadata)
  ),
  values = c(
    paste(head(meta_samples, 20), collapse = ";"),
    paste(head(vst_samples, 20), collapse = ";"),
    paste(head(samples_in_both, 20), collapse = ";"),
    paste(metadata_not_in_vst, collapse = ";"),
    paste(vst_not_in_metadata, collapse = ";")
  )
)

write_tsv(
  sample_summary,
  file.path(OUT_DIR, "sample_matching_summary.tsv")
)

cat("Sample matching summary:\n")
print(sample_summary)
cat("\n")

if (length(samples_in_both) == 0) {
  stop("No overlapping samples between metadata and VST matrix.")
}

if (length(metadata_not_in_vst) > 0) {
  cat("WARNING: metadata samples missing from VST:\n")
  print(metadata_not_in_vst)
}

if (length(vst_not_in_metadata) > 0) {
  cat("WARNING: VST samples missing from metadata:\n")
  print(vst_not_in_metadata)
}

valid_meta <- meta %>%
  filter(sample_id %in% samples_in_both) %>%
  arrange(condition, group_type, timepoint, baseline_block, sample_id)

write_tsv(
  valid_meta,
  file.path(OUT_DIR, "valid_real_metadata.tsv")
)

group_summary <- valid_meta %>%
  count(
    condition,
    group_type,
    timepoint,
    baseline_block,
    name = "n_samples"
  ) %>%
  arrange(condition, group_type, timepoint, baseline_block)

write_tsv(
  group_summary,
  file.path(OUT_DIR, "metadata_group_summary.tsv")
)

cat("Group summary:\n")
print(group_summary, n = 200)
cat("\n")

condition_summary <- valid_meta %>%
  count(condition, group_type, name = "n_samples") %>%
  arrange(condition, group_type)

write_tsv(
  condition_summary,
  file.path(OUT_DIR, "condition_group_summary.tsv")
)

cat("Condition/group summary:\n")
print(condition_summary, n = 100)
cat("\n")

cat("Saved outputs to:\n")
cat(OUT_DIR, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

