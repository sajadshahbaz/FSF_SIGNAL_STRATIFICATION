#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 24: Define Clean Annotation Universe
#
# Purpose:
#   Define the exact feature universe used by FSF real-data
#   analysis and prepare it for clean de novo annotation.
#
# This replaces the previous noisy annotation-resource search.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

VST_FILE <- "/media/saji/5E06441D0643F5152/phase_2/results/results/vst_counts_filtered.tsv"

FSF_METRICS_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/condition_fsf/condition_wise_fsf_metrics.tsv"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/clean_annotation_universe"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/24_define_clean_annotation_universe.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 24: Define Clean Annotation Universe\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

vst <- read_tsv(VST_FILE, show_col_types = FALSE)
fsf <- read_tsv(FSF_METRICS_FILE, show_col_types = FALSE)

if (!"GeneID" %in% colnames(vst)) {
  stop("VST file must contain GeneID column.")
}

if (!"feature_id" %in% colnames(fsf)) {
  stop("FSF metrics file must contain feature_id column.")
}

vst_ids <- unique(as.character(vst$GeneID))
fsf_ids <- unique(as.character(fsf$feature_id))

shared_ids <- intersect(vst_ids, fsf_ids)
vst_not_fsf <- setdiff(vst_ids, fsf_ids)
fsf_not_vst <- setdiff(fsf_ids, vst_ids)

summary_tbl <- tibble::tibble(
  item = c(
    "vst_feature_ids",
    "fsf_feature_ids",
    "shared_feature_ids",
    "vst_not_in_fsf",
    "fsf_not_in_vst"
  ),
  n = c(
    length(vst_ids),
    length(fsf_ids),
    length(shared_ids),
    length(vst_not_fsf),
    length(fsf_not_vst)
  )
)

write_tsv(
  summary_tbl,
  file.path(OUT_DIR, "clean_annotation_universe_summary.tsv")
)

annotation_universe <- tibble::tibble(
  feature_id = sort(shared_ids),
  annotation_universe = "FSF_v1_real_data_universe",
  source = "vst_counts_filtered_and_condition_wise_fsf",
  include_for_annotation = TRUE
)

write_tsv(
  annotation_universe,
  file.path(OUT_DIR, "FSF_v1_clean_annotation_universe.tsv")
)

write_lines(
  annotation_universe$feature_id,
  file.path(OUT_DIR, "FSF_v1_clean_annotation_universe.ids.txt")
)

cat("Clean annotation universe summary:\n")
print(summary_tbl)
cat("\n")

cat("Saved:\n")
cat(file.path(OUT_DIR, "FSF_v1_clean_annotation_universe.tsv"), "\n")
cat(file.path(OUT_DIR, "FSF_v1_clean_annotation_universe.ids.txt"), "\n\n")

if (length(vst_not_fsf) > 0 || length(fsf_not_vst) > 0) {
  cat("WARNING: VST and FSF feature IDs are not perfectly matched.\n")
  cat("Inspect mismatch files before continuing.\n")

  write_lines(
    vst_not_fsf,
    file.path(OUT_DIR, "vst_not_in_fsf.ids.txt")
  )

  write_lines(
    fsf_not_vst,
    file.path(OUT_DIR, "fsf_not_in_vst.ids.txt")
  )
} else {
  cat("VST and FSF feature IDs match perfectly.\n")
}

cat("\nNext step:\n")
cat("Step 25 must extract transcript/protein sequences matching these 15,187 IDs.\n")
cat("Do not use Phase 4 shortlisted BRAKER annotations for FSF enrichment.\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

