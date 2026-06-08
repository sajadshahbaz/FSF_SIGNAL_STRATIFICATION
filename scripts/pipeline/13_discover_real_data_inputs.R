#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 13: Register Real-Data Input Paths
#
# Purpose:
#   Register the locked real-data files used for biological
#   FSF validation.
#
# This script replaces the earlier broad disk-scanning version.
# No recursive search is performed.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

METADATA_FILE <- "/media/saji/5E06441D0643F5152/phase_1/parse_meta/revised_nature_method/Simplified_Metadata_Table.csv"
VST_FILE <- "/media/saji/5E06441D0643F5152/phase_2/results/results/vst_counts_filtered.tsv"
COUNT_FILE <- "/media/saji/5E06441D0643F5152/phase_1/DE_Phase1B/x_count_filt.csv"

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/input_validation"
)

CONFIG_DIR <- file.path(
  ROOT_DIR,
  "config"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/13_discover_real_data_inputs.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CONFIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 13: Register Real-Data Input Paths\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

input_registry <- tibble(
  role = c(
    "metadata",
    "vst_matrix",
    "count_matrix"
  ),
  path = c(
    METADATA_FILE,
    VST_FILE,
    COUNT_FILE
  )
) %>%
  mutate(
    exists = file.exists(path),
    file_size_mb = ifelse(
      exists,
      round(file.info(path)$size / 1024^2, 3),
      NA_real_
    ),
    modified_time = ifelse(
      exists,
      as.character(file.info(path)$mtime),
      NA_character_
    )
  )

print(input_registry)

if (any(!input_registry$exists)) {
  stop("One or more locked real-data input files do not exist.")
}

write_tsv(
  input_registry,
  file.path(CONFIG_DIR, "real_data_inputs.tsv")
)

write_tsv(
  input_registry,
  file.path(OUT_DIR, "real_data_input_registry.tsv")
)

cat("\nSaved:\n")
cat(file.path(CONFIG_DIR, "real_data_inputs.tsv"), "\n")
cat(file.path(OUT_DIR, "real_data_input_registry.tsv"), "\n\n")

cat("Note:\n")
cat("This script intentionally does not search the whole disk.\n")
cat("Real-data inputs are locked manually to avoid accidental use of unrelated files.\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

