#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 26: Translate FSF CDS to Protein FASTA
#
# Purpose:
#   Translate the 15,187 representative CDS sequences from
#   the FSF annotation universe into protein sequences for
#   downstream eggNOG / InterProScan / Pfam annotation.
#
# Input:
#   results/real_data/annotation_sequences/
#   FSF_v1_clean_annotation_universe_cds_sequences.tsv
#
# Outputs:
#   results/real_data/annotation_proteins/
#   FSF_v1_clean_annotation_universe_proteins.faa
#   FSF_v1_translation_summary.tsv
#   FSF_v1_translation_quality_flags.tsv
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

CDS_TSV <- file.path(
  ROOT_DIR,
  "results/real_data/annotation_sequences/FSF_v1_clean_annotation_universe_cds_sequences.tsv"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/annotation_proteins"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/26_translate_fsf_cds_to_proteins.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 26: Translate FSF CDS to Protein FASTA\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

if (!file.exists(CDS_TSV)) {
  stop("CDS TSV not found: ", CDS_TSV)
}

cds <- read_tsv(CDS_TSV, show_col_types = FALSE)

required_cols <- c(
  "feature_id",
  "locus_tag",
  "protein_id",
  "sequence_length",
  "header",
  "sequence"
)

missing_cols <- setdiff(required_cols, colnames(cds))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

cds <- cds %>%
  mutate(
    feature_id = as.character(feature_id),
    locus_tag = as.character(locus_tag),
    protein_id = as.character(protein_id),
    sequence = toupper(as.character(sequence)),
    sequence_length = nchar(sequence)
  )

cat("CDS records loaded:", nrow(cds), "\n")
cat("Unique feature IDs:", n_distinct(cds$feature_id), "\n\n")

# ------------------------------------------------------------
# Genetic code table
# Standard nuclear code.
# ------------------------------------------------------------

genetic_code <- c(
  TTT="F", TTC="F", TTA="L", TTG="L",
  TCT="S", TCC="S", TCA="S", TCG="S",
  TAT="Y", TAC="Y", TAA="*", TAG="*",
  TGT="C", TGC="C", TGA="*", TGG="W",

  CTT="L", CTC="L", CTA="L", CTG="L",
  CCT="P", CCC="P", CCA="P", CCG="P",
  CAT="H", CAC="H", CAA="Q", CAG="Q",
  CGT="R", CGC="R", CGA="R", CGG="R",

  ATT="I", ATC="I", ATA="I", ATG="M",
  ACT="T", ACC="T", ACA="T", ACG="T",
  AAT="N", AAC="N", AAA="K", AAG="K",
  AGT="S", AGC="S", AGA="R", AGG="R",

  GTT="V", GTC="V", GTA="V", GTG="V",
  GCT="A", GCC="A", GCA="A", GCG="A",
  GAT="D", GAC="D", GAA="E", GAG="E",
  GGT="G", GGC="G", GGA="G", GGG="G"
)

translate_cds <- function(seq) {
  seq <- toupper(seq)
  seq <- gsub("[^ACGTN]", "N", seq)

  n <- nchar(seq)
  usable_len <- n - (n %% 3)

  if (usable_len < 3) {
    return("")
  }

  seq_trim <- substr(seq, 1, usable_len)

  codons <- substring(
    seq_trim,
    seq(1, usable_len, by = 3),
    seq(3, usable_len, by = 3)
  )

  aa <- genetic_code[codons]
  aa[is.na(aa)] <- "X"

  protein <- paste0(aa, collapse = "")

  # Remove terminal stop only, keep internal stops for QC flagging.
  protein <- sub("\\*$", "", protein)

  return(protein)
}

# ------------------------------------------------------------
# Translate
# ------------------------------------------------------------

protein_tbl <- cds %>%
  rowwise() %>%
  mutate(
    cds_length = nchar(sequence),
    cds_mod3 = cds_length %% 3,
    starts_with_atg = str_starts(sequence, "ATG"),
    raw_stop_count = str_count(translate_cds(sequence), "\\*"),
    protein_sequence = translate_cds(sequence),
    protein_length = nchar(protein_sequence),
    contains_internal_stop = str_detect(protein_sequence, "\\*"),
    contains_unknown_aa = str_detect(protein_sequence, "X"),
    translation_status = case_when(
      protein_length == 0 ~ "failed_empty_translation",
      contains_internal_stop ~ "translated_with_internal_stop",
      contains_unknown_aa ~ "translated_with_unknown_aa",
      cds_mod3 != 0 ~ "translated_cds_length_not_multiple_of_3",
      !starts_with_atg ~ "translated_non_atg_start",
      TRUE ~ "translated_clean"
    )
  ) %>%
  ungroup()

# ------------------------------------------------------------
# Save protein table and QC flags
# ------------------------------------------------------------

write_tsv(
  protein_tbl %>%
    select(
      feature_id,
      locus_tag,
      protein_id,
      cds_length,
      cds_mod3,
      starts_with_atg,
      protein_length,
      contains_internal_stop,
      contains_unknown_aa,
      translation_status,
      header,
      protein_sequence
    ),
  file.path(OUT_DIR, "FSF_v1_clean_annotation_universe_proteins.tsv")
)

quality_flags <- protein_tbl %>%
  filter(translation_status != "translated_clean") %>%
  select(
    feature_id,
    locus_tag,
    protein_id,
    cds_length,
    cds_mod3,
    starts_with_atg,
    protein_length,
    contains_internal_stop,
    contains_unknown_aa,
    translation_status,
    header
  ) %>%
  arrange(translation_status, feature_id)

write_tsv(
  quality_flags,
  file.path(OUT_DIR, "FSF_v1_translation_quality_flags.tsv")
)

summary_tbl <- protein_tbl %>%
  count(translation_status, name = "n_features") %>%
  mutate(proportion = n_features / sum(n_features)) %>%
  arrange(desc(n_features))

overall_summary <- tibble(
  item = c(
    "input_cds_records",
    "unique_feature_ids",
    "translated_protein_records",
    "clean_translations",
    "flagged_translations",
    "min_protein_length",
    "median_protein_length",
    "mean_protein_length",
    "max_protein_length"
  ),
  value = c(
    nrow(cds),
    n_distinct(cds$feature_id),
    nrow(protein_tbl),
    sum(protein_tbl$translation_status == "translated_clean"),
    sum(protein_tbl$translation_status != "translated_clean"),
    min(protein_tbl$protein_length, na.rm = TRUE),
    median(protein_tbl$protein_length, na.rm = TRUE),
    round(mean(protein_tbl$protein_length, na.rm = TRUE), 2),
    max(protein_tbl$protein_length, na.rm = TRUE)
  )
)

write_tsv(
  summary_tbl,
  file.path(OUT_DIR, "FSF_v1_translation_status_summary.tsv")
)

write_tsv(
  overall_summary,
  file.path(OUT_DIR, "FSF_v1_translation_summary.tsv")
)

# ------------------------------------------------------------
# Write protein FASTA
# ------------------------------------------------------------

protein_fasta_tbl <- protein_tbl %>%
  filter(protein_length > 0) %>%
  arrange(feature_id)

out_fasta <- file.path(
  OUT_DIR,
  "FSF_v1_clean_annotation_universe_proteins.faa"
)

fasta_lines <- unlist(
  lapply(seq_len(nrow(protein_fasta_tbl)), function(i) {
    header <- paste0(
      ">",
      protein_fasta_tbl$feature_id[i],
      " locus_tag=",
      ifelse(is.na(protein_fasta_tbl$locus_tag[i]), "NA", protein_fasta_tbl$locus_tag[i]),
      " protein_id=",
      ifelse(is.na(protein_fasta_tbl$protein_id[i]), "NA", protein_fasta_tbl$protein_id[i]),
      " aa_length=",
      protein_fasta_tbl$protein_length[i],
      " status=",
      protein_fasta_tbl$translation_status[i]
    )

    seq <- protein_fasta_tbl$protein_sequence[i]
    wrapped <- strwrap(seq, width = 80)

    c(header, wrapped)
  })
)

write_lines(fasta_lines, out_fasta)

cat("Translation status summary:\n")
print(summary_tbl, n = 100)
cat("\n")

cat("Overall translation summary:\n")
print(overall_summary)
cat("\n")

cat("Saved:\n")
cat(out_fasta, "\n")
cat(file.path(OUT_DIR, "FSF_v1_clean_annotation_universe_proteins.tsv"), "\n")
cat(file.path(OUT_DIR, "FSF_v1_translation_quality_flags.tsv"), "\n")
cat(file.path(OUT_DIR, "FSF_v1_translation_summary.tsv"), "\n\n")

cat("Next step:\n")
cat("Inspect translation flags. If acceptable, Step 27 can prepare annotation commands for eggNOG/InterProScan/Pfam.\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

