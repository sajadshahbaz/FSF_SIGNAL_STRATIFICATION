#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 24: Validate NCBI Annotation Universe
#
# Purpose:
#   Compare the FSF analyzed feature universe against the full
#   NCBI raw CDS/transcript FASTA universe.
#
# This step does NOT perform annotation or enrichment.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

RAW_DIR <- "/media/saji/5E06441D0643F5152/nature_method2/NM2_instability_signal/01_data/raw"

FSF_ID_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/clean_annotation_universe/FSF_v1_clean_annotation_universe.ids.txt"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/annotation_universe"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/24_validate_ncbi_annotation_universe.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 24: Validate NCBI Annotation Universe\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

if (!dir.exists(RAW_DIR)) {
  stop("Raw NCBI directory not found: ", RAW_DIR)
}

if (!file.exists(FSF_ID_FILE)) {
  stop("FSF ID file not found: ", FSF_ID_FILE)
}

fsf_ids <- read_lines(FSF_ID_FILE)
fsf_ids <- sort(unique(fsf_ids[fsf_ids != ""]))

cat("FSF IDs:", length(fsf_ids), "\n\n")

# ------------------------------------------------------------
# Locate FASTA files in raw directory
# ------------------------------------------------------------

fasta_files <- list.files(
  RAW_DIR,
  pattern = "\\.(fa|fasta|fna|ffn|faa|fas)(\\.gz)?$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(fasta_files) == 0) {
  stop("No FASTA files found in raw directory: ", RAW_DIR)
}

cat("FASTA files found:", length(fasta_files), "\n")
print(fasta_files)
cat("\n")

# ------------------------------------------------------------
# Function to parse gene IDs from FASTA headers
# ------------------------------------------------------------

parse_gene_ids_from_fasta <- function(fasta_file) {

  con <- if (grepl("\\.gz$", fasta_file, ignore.case = TRUE)) {
    gzfile(fasta_file, open = "rt")
  } else {
    file(fasta_file, open = "rt")
  }

  on.exit(close(con), add = TRUE)

  lines <- readLines(con, warn = FALSE)
  headers <- lines[startsWith(lines, ">")]

  gene_ids <- str_match(headers, "\\[gene=([^\\]]+)\\]")[, 2]
  locus_tags <- str_match(headers, "\\[locus_tag=([^\\]]+)\\]")[, 2]
  protein_ids <- str_match(headers, "\\[protein_id=([^\\]]+)\\]")[, 2]

  tibble(
    fasta_file = fasta_file,
    header = headers,
    gene_id = gene_ids,
    locus_tag = locus_tags,
    protein_id = protein_ids
  )
}

# ------------------------------------------------------------
# Parse all FASTA headers
# ------------------------------------------------------------

all_header_tbl <- bind_rows(
  lapply(fasta_files, parse_gene_ids_from_fasta)
)

write_tsv(
  all_header_tbl,
  file.path(OUT_DIR, "ncbi_raw_fasta_header_parsing.tsv")
)

fasta_summary <- all_header_tbl %>%
  group_by(fasta_file) %>%
  summarise(
    n_headers = n(),
    n_gene_ids = sum(!is.na(gene_id)),
    n_unique_gene_ids = n_distinct(gene_id[!is.na(gene_id)]),
    n_locus_tags = sum(!is.na(locus_tag)),
    n_unique_locus_tags = n_distinct(locus_tag[!is.na(locus_tag)]),
    n_protein_ids = sum(!is.na(protein_id)),
    n_unique_protein_ids = n_distinct(protein_id[!is.na(protein_id)]),
    .groups = "drop"
  ) %>%
  arrange(desc(n_unique_gene_ids), desc(n_headers))

write_tsv(
  fasta_summary,
  file.path(OUT_DIR, "ncbi_raw_fasta_file_summary.tsv")
)

cat("FASTA summary:\n")
print(fasta_summary, n = 100)
cat("\n")

# ------------------------------------------------------------
# Select best FASTA by gene_id coverage
# ------------------------------------------------------------

gene_tbl <- all_header_tbl %>%
  filter(!is.na(gene_id)) %>%
  distinct(fasta_file, gene_id, .keep_all = TRUE)

coverage_tbl <- gene_tbl %>%
  group_by(fasta_file) %>%
  summarise(
    ncbi_gene_ids = n_distinct(gene_id),
    shared_with_fsf = length(intersect(unique(gene_id), fsf_ids)),
    fsf_missing_from_this_fasta = length(setdiff(fsf_ids, unique(gene_id))),
    coverage_fraction = shared_with_fsf / length(fsf_ids),
    .groups = "drop"
  ) %>%
  arrange(desc(shared_with_fsf), desc(coverage_fraction), desc(ncbi_gene_ids))

write_tsv(
  coverage_tbl,
  file.path(OUT_DIR, "fasta_fsf_coverage_ranking.tsv")
)

cat("FASTA coverage ranking:\n")
print(coverage_tbl, n = 100)
cat("\n")

best_fasta <- coverage_tbl$fasta_file[1]

cat("Best FASTA selected:\n")
cat(best_fasta, "\n\n")

best_gene_ids <- gene_tbl %>%
  filter(fasta_file == best_fasta) %>%
  pull(gene_id) %>%
  unique() %>%
  sort()

shared_ids <- intersect(fsf_ids, best_gene_ids)
ncbi_only_ids <- setdiff(best_gene_ids, fsf_ids)
fsf_only_ids <- setdiff(fsf_ids, best_gene_ids)

comparison_tbl <- tibble(
  item = c(
    "ncbi_gene_ids_best_fasta",
    "fsf_gene_ids",
    "shared_gene_ids",
    "ncbi_only_gene_ids",
    "fsf_only_gene_ids"
  ),
  n = c(
    length(best_gene_ids),
    length(fsf_ids),
    length(shared_ids),
    length(ncbi_only_ids),
    length(fsf_only_ids)
  )
)

write_tsv(
  comparison_tbl,
  file.path(OUT_DIR, "annotation_universe_summary.tsv")
)

write_lines(fsf_ids, file.path(OUT_DIR, "fsf_gene_ids.txt"))
write_lines(best_gene_ids, file.path(OUT_DIR, "ncbi_gene_ids.txt"))
write_lines(shared_ids, file.path(OUT_DIR, "shared_gene_ids.txt"))
write_lines(ncbi_only_ids, file.path(OUT_DIR, "ncbi_only_gene_ids.txt"))
write_lines(fsf_only_ids, file.path(OUT_DIR, "fsf_only_gene_ids.txt"))

best_fasta_record <- tibble(
  selected_fasta = best_fasta,
  ncbi_gene_ids = length(best_gene_ids),
  fsf_gene_ids = length(fsf_ids),
  shared_gene_ids = length(shared_ids),
  ncbi_only_gene_ids = length(ncbi_only_ids),
  fsf_only_gene_ids = length(fsf_only_ids),
  fsf_coverage_fraction = length(shared_ids) / length(fsf_ids)
)

write_tsv(
  best_fasta_record,
  file.path(OUT_DIR, "selected_ncbi_fasta_for_fsf_annotation.tsv")
)

cat("Annotation universe comparison:\n")
print(comparison_tbl)
cat("\n")

cat("Selected FASTA record:\n")
print(best_fasta_record)
cat("\n")

if (length(fsf_only_ids) == 0) {
  cat("SUCCESS: All FSF IDs are present in selected NCBI FASTA.\n")
} else {
  cat("WARNING: Some FSF IDs are missing from selected NCBI FASTA.\n")
  cat("Inspect fsf_only_gene_ids.txt before sequence extraction.\n")
}

cat("\nNext step:\n")
cat("Step 25 should extract sequences for shared FSF IDs from selected NCBI FASTA.\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

