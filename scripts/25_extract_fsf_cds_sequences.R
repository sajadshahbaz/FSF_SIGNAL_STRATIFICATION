#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 25: Extract FSF CDS Sequences
#
# Purpose:
#   Extract CDS sequences for the 15,187 FSF-analyzed genes
#   from the validated NCBI CDS FASTA.
#
# Matching key:
#   FASTA header field [gene=...]
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

SELECTED_FASTA_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/annotation_universe/selected_ncbi_fasta_for_fsf_annotation.tsv"
)

FSF_ID_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/annotation_universe/shared_gene_ids.txt"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/annotation_sequences"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/25_extract_fsf_cds_sequences.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 25: Extract FSF CDS Sequences\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

if (!file.exists(SELECTED_FASTA_FILE)) {
  stop("Selected FASTA record not found: ", SELECTED_FASTA_FILE)
}

if (!file.exists(FSF_ID_FILE)) {
  stop("FSF shared ID file not found: ", FSF_ID_FILE)
}

selected <- read_tsv(SELECTED_FASTA_FILE, show_col_types = FALSE)

if (!"selected_fasta" %in% colnames(selected)) {
  stop("selected_ncbi_fasta_for_fsf_annotation.tsv must contain selected_fasta column.")
}

FASTA_FILE <- selected$selected_fasta[1]

if (!file.exists(FASTA_FILE)) {
  stop("Selected FASTA does not exist: ", FASTA_FILE)
}

fsf_ids <- read_lines(FSF_ID_FILE)
fsf_ids <- sort(unique(fsf_ids[fsf_ids != ""]))

cat("Selected FASTA:\n")
cat(FASTA_FILE, "\n\n")
cat("FSF IDs to extract:", length(fsf_ids), "\n\n")

# ------------------------------------------------------------
# Read FASTA records
# ------------------------------------------------------------

read_fasta_records <- function(fasta_file) {

  con <- if (grepl("\\.gz$", fasta_file, ignore.case = TRUE)) {
    gzfile(fasta_file, open = "rt")
  } else {
    file(fasta_file, open = "rt")
  }

  on.exit(close(con), add = TRUE)

  lines <- readLines(con, warn = FALSE)

  header_idx <- which(startsWith(lines, ">"))

  if (length(header_idx) == 0) {
    stop("No FASTA headers found in: ", fasta_file)
  }

  end_idx <- c(header_idx[-1] - 1, length(lines))

  records <- vector("list", length(header_idx))

  for (i in seq_along(header_idx)) {
    header <- lines[header_idx[i]]
    seq_lines <- lines[(header_idx[i] + 1):end_idx[i]]
    seq_lines <- seq_lines[!startsWith(seq_lines, ">")]
    sequence <- paste(seq_lines, collapse = "")

    gene_id <- str_match(header, "\\[gene=([^\\]]+)\\]")[, 2]
    locus_tag <- str_match(header, "\\[locus_tag=([^\\]]+)\\]")[, 2]
    protein_id <- str_match(header, "\\[protein_id=([^\\]]+)\\]")[, 2]

    records[[i]] <- tibble(
      header = header,
      gene_id = gene_id,
      locus_tag = locus_tag,
      protein_id = protein_id,
      sequence = sequence,
      sequence_length = nchar(sequence)
    )
  }

  bind_rows(records)
}

fasta_tbl <- read_fasta_records(FASTA_FILE)

cat("FASTA records read:", nrow(fasta_tbl), "\n")
cat("Unique gene IDs in FASTA:", n_distinct(fasta_tbl$gene_id[!is.na(fasta_tbl$gene_id)]), "\n\n")

write_tsv(
  fasta_tbl %>%
    select(header, gene_id, locus_tag, protein_id, sequence_length),
  file.path(OUT_DIR, "ncbi_cds_fasta_record_index.tsv")
)

# ------------------------------------------------------------
# Extract FSF sequences
# ------------------------------------------------------------

matched_tbl <- fasta_tbl %>%
  filter(gene_id %in% fsf_ids) %>%
  arrange(gene_id, desc(sequence_length))

# Keep one representative CDS per gene_id:
# longest CDS if multiple records exist for the same gene.
representative_tbl <- matched_tbl %>%
  group_by(gene_id) %>%
  arrange(desc(sequence_length), header, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  rename(feature_id = gene_id) %>%
  mutate(
    representative_rule = "longest_cds_per_gene_id"
  ) %>%
  select(
    feature_id,
    locus_tag,
    protein_id,
    sequence_length,
    representative_rule,
    header,
    sequence
  ) %>%
  arrange(feature_id)

matched_ids <- unique(representative_tbl$feature_id)
unmatched_ids <- setdiff(fsf_ids, matched_ids)

summary_tbl <- tibble(
  item = c(
    "fsf_ids_requested",
    "matched_fsf_gene_ids",
    "unmatched_fsf_gene_ids",
    "raw_matched_cds_records",
    "representative_cds_records"
  ),
  n = c(
    length(fsf_ids),
    length(matched_ids),
    length(unmatched_ids),
    nrow(matched_tbl),
    nrow(representative_tbl)
  )
)

write_tsv(
  summary_tbl,
  file.path(OUT_DIR, "fsf_cds_sequence_extraction_summary.tsv")
)

write_tsv(
  representative_tbl,
  file.path(OUT_DIR, "FSF_v1_clean_annotation_universe_cds_sequences.tsv")
)

write_lines(
  unmatched_ids,
  file.path(OUT_DIR, "unmatched_fsf_gene_ids.txt")
)

# ------------------------------------------------------------
# Write FASTA
# ------------------------------------------------------------

out_fasta <- file.path(
  OUT_DIR,
  "FSF_v1_clean_annotation_universe_cds_sequences.fna"
)

fasta_out_lines <- unlist(
  lapply(seq_len(nrow(representative_tbl)), function(i) {

    header <- paste0(
      ">",
      representative_tbl$feature_id[i],
      " locus_tag=",
      ifelse(is.na(representative_tbl$locus_tag[i]), "NA", representative_tbl$locus_tag[i]),
      " protein_id=",
      ifelse(is.na(representative_tbl$protein_id[i]), "NA", representative_tbl$protein_id[i]),
      " length=",
      representative_tbl$sequence_length[i]
    )

    seq <- representative_tbl$sequence[i]

    wrapped <- strwrap(seq, width = 80)

    c(header, wrapped)
  })
)

write_lines(fasta_out_lines, out_fasta)

cat("Extraction summary:\n")
print(summary_tbl)
cat("\n")

cat("Saved:\n")
cat(out_fasta, "\n")
cat(file.path(OUT_DIR, "FSF_v1_clean_annotation_universe_cds_sequences.tsv"), "\n")
cat(file.path(OUT_DIR, "fsf_cds_sequence_extraction_summary.tsv"), "\n\n")

if (length(unmatched_ids) == 0) {
  cat("SUCCESS: All FSF IDs were extracted from NCBI CDS FASTA.\n")
} else {
  cat("WARNING: Some FSF IDs were not extracted.\n")
  cat("Inspect unmatched_fsf_gene_ids.txt before annotation.\n")
}

cat("\nNext step:\n")
cat("Step 26 can translate CDS to protein sequences or run CDS-level annotation.\n")
cat("For eggNOG/InterPro, protein FASTA is preferred.\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

