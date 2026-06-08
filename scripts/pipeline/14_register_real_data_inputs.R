#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 14: Register and Validate Real-Data Inputs
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

METADATA_FILE <- "/media/saji/5E06441D0643F5152/phase_1/parse_meta/revised_nature_method/Simplified_Metadata_Table.csv"
VST_FILE <- "/media/saji/5E06441D0643F5152/phase_2/results/results/vst_counts_filtered.tsv"
COUNT_FILE <- "/media/saji/5E06441D0643F5152/phase_1/DE_Phase1B/x_count_filt.csv"

OUT_DIR <- file.path(ROOT_DIR, "results", "real_data", "input_validation")
CONFIG_FILE <- file.path(ROOT_DIR, "config", "real_data_inputs.tsv")
LOG_FILE <- file.path(ROOT_DIR, "results", "logs", "14_register_real_data_inputs.log")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(CONFIG_FILE), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 14: Register and Validate Real-Data Inputs\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

inputs <- tibble::tibble(
  role = c("metadata", "vst_matrix", "count_matrix"),
  path = c(METADATA_FILE, VST_FILE, COUNT_FILE),
  exists = file.exists(c(METADATA_FILE, VST_FILE, COUNT_FILE)),
  file_size_mb = ifelse(
    file.exists(c(METADATA_FILE, VST_FILE, COUNT_FILE)),
    round(file.info(c(METADATA_FILE, VST_FILE, COUNT_FILE))$size / 1024^2, 3),
    NA_real_
  )
)

print(inputs)

if (any(!inputs$exists)) {
  stop("One or more required input files do not exist. Fix paths before continuing.")
}

write_tsv(inputs, CONFIG_FILE)

cat("\nSaved input registry:\n")
cat(CONFIG_FILE, "\n\n")

# ----------------------------
# Read metadata
# ----------------------------

meta <- read_csv(METADATA_FILE, show_col_types = FALSE)

cat("Metadata dimensions:", nrow(meta), "rows x", ncol(meta), "columns\n")
cat("Metadata columns:\n")
print(colnames(meta))
cat("\n")

write_tsv(
  tibble::tibble(metadata_column = colnames(meta)),
  file.path(OUT_DIR, "metadata_columns.tsv")
)

# ----------------------------
# Read VST matrix
# ----------------------------

vst <- read_tsv(VST_FILE, show_col_types = FALSE)

cat("VST dimensions:", nrow(vst), "rows x", ncol(vst), "columns\n")
cat("First VST columns:\n")
print(head(colnames(vst), 20))
cat("\n")

write_tsv(
  tibble::tibble(vst_column = colnames(vst)),
  file.path(OUT_DIR, "vst_columns.tsv")
)

# ----------------------------
# Read count matrix header only
# ----------------------------

count_header <- read_csv(COUNT_FILE, n_max = 5, show_col_types = FALSE)

cat("Count matrix preview dimensions:", nrow(count_header), "rows x", ncol(count_header), "columns\n")
cat("First count columns:\n")
print(head(colnames(count_header), 20))
cat("\n")

write_tsv(
  tibble::tibble(count_column = colnames(count_header)),
  file.path(OUT_DIR, "count_columns.tsv")
)

# ----------------------------
# Guess sample ID column in metadata
# ----------------------------

possible_sample_cols <- c("sample_id", "sample", "Sample", "Sample_ID", "Run", "run", "sra", "SRA", "accession")

detected_sample_cols <- intersect(possible_sample_cols, colnames(meta))

if (length(detected_sample_cols) == 0) {
  cat("WARNING: No obvious sample ID column detected in metadata.\n")
} else {
  cat("Detected possible metadata sample ID columns:\n")
  print(detected_sample_cols)
}

# ----------------------------
# Save validation summary
# ----------------------------

summary <- tibble::tibble(
  item = c(
    "metadata_rows",
    "metadata_columns",
    "vst_rows",
    "vst_columns",
    "count_preview_columns",
    "detected_sample_id_columns"
  ),
  value = c(
    as.character(nrow(meta)),
    as.character(ncol(meta)),
    as.character(nrow(vst)),
    as.character(ncol(vst)),
    as.character(ncol(count_header)),
    paste(detected_sample_cols, collapse = ";")
  )
)

write_tsv(
  summary,
  file.path(OUT_DIR, "real_data_input_validation_summary.tsv")
)

cat("\nSaved validation summary:\n")
cat(file.path(OUT_DIR, "real_data_input_validation_summary.tsv"), "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

