#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

PROTEIN_TSV <- file.path(ROOT_DIR, "results/real_data/annotation_proteins/FSF_v1_clean_annotation_universe_proteins.tsv")
EGGNOG_FILE <- file.path(ROOT_DIR, "results/real_data/functional_annotation/eggnog/FSF_v1_eggnog.emapper.annotations")
INTERPRO_FILE <- file.path(ROOT_DIR, "results/real_data/functional_annotation/interproscan/FSF_v1_clean_annotation_universe_proteins_noSTOP.faa.tsv")
PFAM_FILE <- file.path(ROOT_DIR, "results/real_data/functional_annotation/pfam/FSF_v1_pfam.domtblout")

OUT_DIR <- file.path(ROOT_DIR, "results/real_data/annotation_master")
LOG_FILE <- file.path(ROOT_DIR, "results/logs/28_build_integrated_annotation_master.log")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 28: Build Integrated Annotation Master\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

clean_vec <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  x <- x[x != ""]
  x <- x[x != "-"]
  x <- unique(x)
  sort(x)
}

collapse_clean <- function(x) {
  y <- clean_vec(x)
  if (length(y) == 0) return(NA_character_)
  paste(y, collapse = ";")
}

protein_tbl <- read_tsv(PROTEIN_TSV, show_col_types = FALSE) %>%
  transmute(
    feature_id = as.character(feature_id),
    locus_tag = as.character(locus_tag),
    protein_id = as.character(protein_id),
    protein_length = as.integer(protein_length),
    translation_status = as.character(translation_status)
  ) %>%
  distinct(feature_id, .keep_all = TRUE)

cat("FSF protein universe:", nrow(protein_tbl), "\n")

eggnog <- read_tsv(EGGNOG_FILE, comment = "##", show_col_types = FALSE)
colnames(eggnog)[1] <- str_replace(colnames(eggnog)[1], "^#", "")

eggnog_clean <- eggnog %>%
  transmute(
    feature_id = as.character(query),
    eggnog_annotated = TRUE,
    eggnog_seed_ortholog = as.character(seed_ortholog),
    eggnog_evalue = as.character(evalue),
    eggnog_score = as.character(score),
    eggnog_ogs = as.character(eggNOG_OGs),
    eggnog_max_annot_lvl = as.character(max_annot_lvl),
    eggnog_cog_category = as.character(COG_category),
    eggnog_description = as.character(Description),
    eggnog_preferred_name = as.character(Preferred_name),
    eggnog_go = as.character(GOs),
    eggnog_ec = as.character(EC),
    eggnog_kegg_ko = as.character(KEGG_ko),
    eggnog_kegg_pathway = as.character(KEGG_Pathway),
    eggnog_kegg_module = as.character(KEGG_Module),
    eggnog_pfam = as.character(PFAMs)
  ) %>%
  distinct(feature_id, .keep_all = TRUE)

cat("eggNOG annotated features:", nrow(eggnog_clean), "\n")

interpro_raw <- read_tsv(INTERPRO_FILE, col_names = FALSE, show_col_types = FALSE)

if (ncol(interpro_raw) < 11) {
  stop("InterProScan TSV has fewer than expected columns.")
}

ipr_names <- c(
  "feature_id",
  "sequence_md5",
  "sequence_length",
  "analysis",
  "signature_accession",
  "signature_description",
  "start",
  "stop",
  "score",
  "status",
  "date",
  "interpro_accession",
  "interpro_description",
  "go_terms",
  "pathways"
)

colnames(interpro_raw)[seq_len(min(length(ipr_names), ncol(interpro_raw)))] <-
  ipr_names[seq_len(min(length(ipr_names), ncol(interpro_raw)))]

for (nm in c("interpro_accession", "interpro_description", "go_terms", "pathways")) {
  if (!nm %in% colnames(interpro_raw)) {
    interpro_raw[[nm]] <- NA_character_
  }
}

interpro_clean <- interpro_raw %>%
  mutate(
    feature_id = as.character(feature_id),
    analysis = as.character(analysis),
    signature_accession = as.character(signature_accession),
    signature_description = as.character(signature_description),
    interpro_accession = as.character(interpro_accession),
    interpro_description = as.character(interpro_description),
    go_terms = as.character(go_terms),
    pathways = as.character(pathways)
  )

interpro_summary <- interpro_clean %>%
  group_by(feature_id) %>%
  summarise(
    interpro_annotated = TRUE,
    interpro_analyses = collapse_clean(analysis),
    interpro_signature_accessions = collapse_clean(signature_accession),
    interpro_signature_descriptions = collapse_clean(signature_description),
    interpro_ids = collapse_clean(interpro_accession),
    interpro_descriptions = collapse_clean(interpro_description),
    interpro_go = collapse_clean(unlist(str_split(go_terms, "\\|"))),
    interpro_pathways = collapse_clean(unlist(str_split(pathways, "\\|"))),
    n_interpro_rows = n(),
    .groups = "drop"
  )

cat("InterProScan annotated features:", nrow(interpro_summary), "\n")
cat("InterProScan total rows:", nrow(interpro_raw), "\n")

# QC specifically for the previous suspected corruption
interpro_qc <- interpro_summary %>%
  summarise(
    n_features = n(),
    n_unique_interpro_id_strings = n_distinct(interpro_ids, na.rm = TRUE),
    n_unique_interpro_description_strings = n_distinct(interpro_descriptions, na.rm = TRUE),
    n_features_with_ipr035921 = sum(str_detect(coalesce(interpro_ids, ""), "IPR035921")),
    n_features_with_fv_atpase_description = sum(str_detect(coalesce(interpro_descriptions, ""), fixed("F/V-ATP synthase subunit C superfamily")))
  )

write_tsv(interpro_qc, file.path(OUT_DIR, "FSF_v1_interpro_aggregation_QC.tsv"))

cat("\nInterPro aggregation QC:\n")
print(interpro_qc)
cat("\n")

pfam_lines <- read_lines(PFAM_FILE)
pfam_data_lines <- pfam_lines[!str_starts(pfam_lines, "#")]
pfam_data_lines <- pfam_data_lines[pfam_data_lines != ""]

parse_pfam_line <- function(x) {
  parts <- str_split(x, "\\s+", simplify = TRUE)
  if (ncol(parts) < 22) return(NULL)

  tibble(
    target_name = parts[1],
    target_accession = parts[2],
    target_len = parts[3],
    feature_id = parts[4],
    query_accession = parts[5],
    query_len = parts[6],
    full_evalue = parts[7],
    full_score = parts[8],
    full_bias = parts[9],
    domain_number = parts[10],
    domain_total = parts[11],
    c_evalue = parts[12],
    i_evalue = parts[13],
    domain_score = parts[14],
    domain_bias = parts[15],
    hmm_from = parts[16],
    hmm_to = parts[17],
    ali_from = parts[18],
    ali_to = parts[19],
    env_from = parts[20],
    env_to = parts[21],
    acc = parts[22],
    description = ifelse(ncol(parts) > 22, paste(parts[23:ncol(parts)], collapse = " "), NA_character_)
  )
}

pfam_raw <- bind_rows(lapply(pfam_data_lines, parse_pfam_line))

pfam_summary <- pfam_raw %>%
  mutate(
    feature_id = as.character(feature_id),
    target_name = as.character(target_name),
    target_accession = as.character(target_accession),
    description = as.character(description)
  ) %>%
  group_by(feature_id) %>%
  summarise(
    pfam_annotated = TRUE,
    pfam_domains = collapse_clean(target_name),
    pfam_accessions = collapse_clean(target_accession),
    pfam_descriptions = collapse_clean(description),
    n_pfam_hits = n(),
    .groups = "drop"
  )

cat("Pfam annotated features:", nrow(pfam_summary), "\n")
cat("Pfam total domain hits:", nrow(pfam_raw), "\n\n")

annotation_master <- protein_tbl %>%
  left_join(eggnog_clean, by = "feature_id") %>%
  left_join(interpro_summary, by = "feature_id") %>%
  left_join(pfam_summary, by = "feature_id") %>%
  mutate(
    eggnog_annotated = ifelse(is.na(eggnog_annotated), FALSE, eggnog_annotated),
    interpro_annotated = ifelse(is.na(interpro_annotated), FALSE, interpro_annotated),
    pfam_annotated = ifelse(is.na(pfam_annotated), FALSE, pfam_annotated),
    annotation_source_count =
      as.integer(eggnog_annotated) +
      as.integer(interpro_annotated) +
      as.integer(pfam_annotated),
    annotation_status = case_when(
      annotation_source_count == 0 ~ "unannotated",
      annotation_source_count == 1 ~ "single_source_annotation",
      annotation_source_count == 2 ~ "two_source_annotation",
      annotation_source_count == 3 ~ "three_source_annotation",
      TRUE ~ "unknown"
    )
  ) %>%
  arrange(feature_id)

write_tsv(annotation_master, file.path(OUT_DIR, "FSF_v1_annotation_master.tsv"))

coverage_summary <- tibble(
  metric = c(
    "total_fsf_features",
    "eggnog_annotated_features",
    "interpro_annotated_features",
    "pfam_annotated_features",
    "features_with_at_least_one_annotation",
    "features_with_at_least_two_annotation_layers",
    "features_with_all_three_annotation_layers",
    "unannotated_features"
  ),
  n_features = c(
    nrow(annotation_master),
    sum(annotation_master$eggnog_annotated),
    sum(annotation_master$interpro_annotated),
    sum(annotation_master$pfam_annotated),
    sum(annotation_master$annotation_source_count >= 1),
    sum(annotation_master$annotation_source_count >= 2),
    sum(annotation_master$annotation_source_count == 3),
    sum(annotation_master$annotation_source_count == 0)
  )
) %>%
  mutate(proportion = n_features / nrow(annotation_master))

write_tsv(coverage_summary, file.path(OUT_DIR, "FSF_v1_annotation_coverage_summary.tsv"))

layer_statistics <- annotation_master %>%
  count(annotation_status, name = "n_features") %>%
  mutate(proportion = n_features / sum(n_features)) %>%
  arrange(desc(n_features))

write_tsv(layer_statistics, file.path(OUT_DIR, "FSF_v1_annotation_layer_statistics.tsv"))

unannotated_genes <- annotation_master %>%
  filter(annotation_source_count == 0) %>%
  select(feature_id, locus_tag, protein_id, protein_length, translation_status)

write_tsv(unannotated_genes, file.path(OUT_DIR, "FSF_v1_unannotated_genes.tsv"))

cat("Annotation coverage summary:\n")
print(coverage_summary, n = 100)
cat("\n")

cat("Annotation layer statistics:\n")
print(layer_statistics, n = 100)
cat("\n")

cat("Saved outputs:\n")
cat(file.path(OUT_DIR, "FSF_v1_annotation_master.tsv"), "\n")
cat(file.path(OUT_DIR, "FSF_v1_interpro_aggregation_QC.tsv"), "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()
