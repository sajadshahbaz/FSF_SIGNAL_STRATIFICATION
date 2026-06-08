#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

FSF_FILE <- file.path(ROOT_DIR, "results/real_data/condition_fsf/condition_wise_fsf_metrics.tsv")
ANNOTATION_FILE <- file.path(ROOT_DIR, "results/real_data/annotation_master/FSF_v1_annotation_master.tsv")
OUT_DIR <- file.path(ROOT_DIR, "results/real_data/class_annotation")
LOG_FILE <- file.path(ROOT_DIR, "results/logs/29_class_annotation_summary.log")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 29: Class Annotation Summary\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

fsf <- read_tsv(FSF_FILE, show_col_types = FALSE)
anno <- read_tsv(ANNOTATION_FILE, show_col_types = FALSE)

dat <- fsf %>%
  mutate(
    feature_id = as.character(feature_id),
    condition = as.character(condition),
    fsf_signal_class = as.character(fsf_signal_class),
    stability_level = as.character(stability_level),
    dominant_state = as.character(dominant_state)
  ) %>%
  left_join(
    anno %>% mutate(feature_id = as.character(feature_id)),
    by = "feature_id"
  ) %>%
  mutate(
    has_annotation = annotation_source_count > 0,
    condition_class = paste(condition, fsf_signal_class, sep = "__")
  )

if (any(is.na(dat$annotation_source_count))) {
  write_tsv(
    dat %>% filter(is.na(annotation_source_count)),
    file.path(OUT_DIR, "ERROR_fsf_rows_missing_annotation_join.tsv")
  )
  stop("Some FSF rows failed annotation join.")
}

write_tsv(dat, file.path(OUT_DIR, "FSF_condition_class_annotation_master.tsv"))

class_summary <- dat %>%
  group_by(condition, stability_level, dominant_state, fsf_signal_class) %>%
  summarise(
    n_features = n(),
    mean_SSI = mean(SSI, na.rm = TRUE),
    median_SSI = median(SSI, na.rm = TRUE),
    mean_stability_deviation = mean(stability_deviation, na.rm = TRUE),
    median_stability_deviation = median(stability_deviation, na.rm = TRUE),
    eggnog_annotated = sum(eggnog_annotated, na.rm = TRUE),
    interpro_annotated = sum(interpro_annotated, na.rm = TRUE),
    pfam_annotated = sum(pfam_annotated, na.rm = TRUE),
    at_least_one_annotation = sum(annotation_source_count >= 1, na.rm = TRUE),
    at_least_two_annotation_layers = sum(annotation_source_count >= 2, na.rm = TRUE),
    all_three_annotation_layers = sum(annotation_source_count == 3, na.rm = TRUE),
    unannotated = sum(annotation_source_count == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    eggnog_fraction = eggnog_annotated / n_features,
    interpro_fraction = interpro_annotated / n_features,
    pfam_fraction = pfam_annotated / n_features,
    annotated_fraction = at_least_one_annotation / n_features,
    unannotated_fraction = unannotated / n_features
  ) %>%
  arrange(condition, stability_level, dominant_state, fsf_signal_class)

write_tsv(class_summary, file.path(OUT_DIR, "FSF_condition_class_annotation_summary.tsv"))

split_terms <- function(tbl, term_col, sep_pattern) {
  tbl %>%
    filter(!is.na(.data[[term_col]]), .data[[term_col]] != "", .data[[term_col]] != "-") %>%
    select(condition, stability_level, dominant_state, fsf_signal_class, feature_id, term = all_of(term_col)) %>%
    separate_rows(term, sep = sep_pattern) %>%
    mutate(term = str_trim(term)) %>%
    filter(!is.na(term), term != "", term != "-")
}

make_term_summary <- function(long_tbl, term_name = "term") {
  long_tbl %>%
    count(condition, stability_level, dominant_state, fsf_signal_class, .data[[term_name]], name = "n_genes") %>%
    group_by(condition, stability_level, dominant_state, fsf_signal_class) %>%
    mutate(rank = dense_rank(desc(n_genes))) %>%
    ungroup() %>%
    arrange(condition, stability_level, dominant_state, fsf_signal_class, rank)
}

eggnog_go_summary <- split_terms(dat, "eggnog_go", ",") %>%
  make_term_summary("term")
write_tsv(eggnog_go_summary, file.path(OUT_DIR, "FSF_condition_class_eggNOG_GO_summary.tsv"))

interpro_go_summary <- split_terms(dat, "interpro_go", ";") %>%
  make_term_summary("term")
write_tsv(interpro_go_summary, file.path(OUT_DIR, "FSF_condition_class_InterPro_GO_summary.tsv"))

interpro_id_summary <- dat %>%
  filter(!is.na(interpro_ids), interpro_ids != "", interpro_ids != "-") %>%
  select(condition, stability_level, dominant_state, fsf_signal_class, feature_id, interpro_ids) %>%
  separate_rows(interpro_ids, sep = ";") %>%
  mutate(interpro_ids = str_trim(interpro_ids)) %>%
  filter(interpro_ids != "", interpro_ids != "-") %>%
  count(condition, stability_level, dominant_state, fsf_signal_class, interpro_ids, name = "n_genes") %>%
  group_by(condition, stability_level, dominant_state, fsf_signal_class) %>%
  mutate(rank = dense_rank(desc(n_genes))) %>%
  ungroup() %>%
  arrange(condition, stability_level, dominant_state, fsf_signal_class, rank, interpro_ids)

interpro_desc_lookup <- dat %>%
  filter(!is.na(interpro_ids), !is.na(interpro_descriptions)) %>%
  select(interpro_ids, interpro_descriptions) %>%
  separate_rows(interpro_ids, sep = ";") %>%
  mutate(interpro_ids = str_trim(interpro_ids)) %>%
  group_by(interpro_ids) %>%
  summarise(
    interpro_descriptions = first(na.omit(interpro_descriptions)),
    .groups = "drop"
  )

interpro_summary <- interpro_id_summary %>%
  left_join(interpro_desc_lookup, by = "interpro_ids")

write_tsv(interpro_summary, file.path(OUT_DIR, "FSF_condition_class_InterPro_summary.tsv"))

pfam_summary <- dat %>%
  filter(!is.na(pfam_domains), pfam_domains != "", pfam_domains != "-") %>%
  select(condition, stability_level, dominant_state, fsf_signal_class, feature_id, pfam_domains) %>%
  separate_rows(pfam_domains, sep = ";") %>%
  mutate(pfam_domains = str_trim(pfam_domains)) %>%
  filter(pfam_domains != "", pfam_domains != "-") %>%
  count(condition, stability_level, dominant_state, fsf_signal_class, pfam_domains, name = "n_genes") %>%
  group_by(condition, stability_level, dominant_state, fsf_signal_class) %>%
  mutate(rank = dense_rank(desc(n_genes))) %>%
  ungroup() %>%
  arrange(condition, stability_level, dominant_state, fsf_signal_class, rank, pfam_domains)

pfam_desc_lookup <- dat %>%
  filter(!is.na(pfam_domains), !is.na(pfam_descriptions)) %>%
  select(pfam_domains, pfam_descriptions) %>%
  separate_rows(pfam_domains, sep = ";") %>%
  mutate(pfam_domains = str_trim(pfam_domains)) %>%
  group_by(pfam_domains) %>%
  summarise(
    pfam_descriptions = first(na.omit(pfam_descriptions)),
    .groups = "drop"
  )

pfam_summary <- pfam_summary %>%
  left_join(pfam_desc_lookup, by = "pfam_domains")

write_tsv(pfam_summary, file.path(OUT_DIR, "FSF_condition_class_Pfam_summary.tsv"))

kegg_ko_summary <- split_terms(dat, "eggnog_kegg_ko", ",") %>%
  make_term_summary("term")
write_tsv(kegg_ko_summary, file.path(OUT_DIR, "FSF_condition_class_eggNOG_KEGG_KO_summary.tsv"))

kegg_pathway_summary <- split_terms(dat, "eggnog_kegg_pathway", ",") %>%
  make_term_summary("term")
write_tsv(kegg_pathway_summary, file.path(OUT_DIR, "FSF_condition_class_eggNOG_KEGG_pathway_summary.tsv"))

top_interpro <- interpro_summary %>%
  group_by(condition, stability_level, dominant_state, fsf_signal_class) %>%
  slice_head(n = 20) %>%
  ungroup()
write_tsv(top_interpro, file.path(OUT_DIR, "TOP20_FSF_condition_class_InterPro_summary.tsv"))

top_pfam <- pfam_summary %>%
  group_by(condition, stability_level, dominant_state, fsf_signal_class) %>%
  slice_head(n = 20) %>%
  ungroup()
write_tsv(top_pfam, file.path(OUT_DIR, "TOP20_FSF_condition_class_Pfam_summary.tsv"))

top_go <- interpro_go_summary %>%
  group_by(condition, stability_level, dominant_state, fsf_signal_class) %>%
  slice_head(n = 20) %>%
  ungroup()
write_tsv(top_go, file.path(OUT_DIR, "TOP20_FSF_condition_class_InterPro_GO_summary.tsv"))

directional_stable <- dat %>%
  filter(
    dominant_state %in% c("up", "down"),
    stability_level %in% c("highly_stable", "stable"),
    annotation_source_count >= 1
  ) %>%
  arrange(
    condition,
    dominant_state,
    desc(stability_level == "highly_stable"),
    desc(annotation_source_count),
    desc(SSI),
    stability_deviation,
    feature_id
  )

write_tsv(directional_stable, file.path(OUT_DIR, "FSF_directional_stable_annotated_genes.tsv"))

top_directional_stable <- directional_stable %>%
  group_by(condition, dominant_state) %>%
  slice_head(n = 100) %>%
  ungroup()

write_tsv(top_directional_stable, file.path(OUT_DIR, "TOP100_FSF_directional_stable_annotated_genes.tsv"))

cat("Class annotation summary:\n")
print(class_summary, n = 100)
cat("\n")

cat("Rows written:\n")
cat("Full joined master:", nrow(dat), "\n")
cat("Class summary:", nrow(class_summary), "\n")
cat("eggNOG GO summary:", nrow(eggnog_go_summary), "\n")
cat("InterPro GO summary:", nrow(interpro_go_summary), "\n")
cat("InterPro summary:", nrow(interpro_summary), "\n")
cat("Pfam summary:", nrow(pfam_summary), "\n")
cat("KEGG KO summary:", nrow(kegg_ko_summary), "\n")
cat("KEGG pathway summary:", nrow(kegg_pathway_summary), "\n")
cat("Directional stable annotated genes:", nrow(directional_stable), "\n\n")

cat("Saved outputs to:\n")
cat(OUT_DIR, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()
