#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

KEGG_SIG <- file.path(ROOT_DIR, "results/real_data/enrichment_KEGG/FSF_significant_KEGG_enrichment_FDR005.tsv")
KEGG_SUM <- file.path(ROOT_DIR, "results/real_data/enrichment_KEGG/FSF_KEGG_enrichment_summary_by_gene_set.tsv")
GO_BIO <- file.path(ROOT_DIR, "results/real_data/enrichment_GO_interpretation/FSF_GO_biological_theme_summary.tsv")
MANIFEST <- file.path(ROOT_DIR, "results/real_data/enrichment_gene_sets/FSF_enrichment_gene_set_manifest.tsv")

OUT_DIR <- file.path(ROOT_DIR, "results/real_data/enrichment_KEGG_interpretation")
FINAL_DIR <- file.path(ROOT_DIR, "results/real_data/final_manuscript_tables")
LOG_FILE <- file.path(ROOT_DIR, "results/logs/34_build_biological_profiles.log")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FINAL_DIR, recursive = TRUE, showWarnings = FALSE)
sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 34: KEGG Theme Summary + Final Biological Profiles\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

for (f in c(KEGG_SIG, KEGG_SUM, GO_BIO, MANIFEST)) {
  if (!file.exists(f)) stop("Missing file: ", f)
}

kegg_sig <- read_tsv(KEGG_SIG, show_col_types = FALSE)
kegg_sum <- read_tsv(KEGG_SUM, show_col_types = FALSE)
go_bio <- read_tsv(GO_BIO, show_col_types = FALSE)
manifest <- read_tsv(MANIFEST, show_col_types = FALSE)

assign_kegg_theme <- function(term) {
  case_when(
    str_detect(term, "ko03010|map03010|Ribosome|K029|K028") ~ "Translation / ribosome",
    str_detect(term, "ko03030|map03030|ko034|map034|K107|K108|K109|DNA|repair|replication|recombination") ~ "Genome maintenance / DNA repair",
    str_detect(term, "ko04110|map04110|ko04111|map04111|ko04113|map04113|cell cycle") ~ "Cell cycle regulation",
    str_detect(term, "ko03050|map03050|ko04120|map04120|proteasome|ubiquitin|K056") ~ "Protein turnover / proteostasis",
    str_detect(term, "ko00190|map00190|ko01100|map01100|ko01110|map01110|ko01120|map01120|ko01200|map01200|oxidative|metabolism") ~ "Energy / core metabolism",
    str_detect(term, "ko04010|map04010|ko04013|map04013|ko04022|map04022|ko04151|map04151|ko046|map046|signaling|MAPK") ~ "Stress signaling / signal transduction",
    str_detect(term, "ko045|map045|adhesion|junction|ECM") ~ "Cell adhesion / extracellular interaction",
    str_detect(term, "ko048|map048|cytoskeleton|actin") ~ "Cytoskeleton / cellular remodeling",
    str_detect(term, "ko047|map047|synapse|neuro|neuron") ~ "Neural signaling",
    str_detect(term, "ko049|map049|endocrine|hormone|insulin") ~ "Endocrine / physiological regulation",
    str_detect(term, "ko052|map052|cancer") ~ "Broad proliferation / regulatory pathway",
    TRUE ~ "Other KEGG / manual review"
  )
}

top_kegg <- kegg_sig %>%
  mutate(KEGG_theme = assign_kegg_theme(kegg_term)) %>%
  group_by(enrichment_type, condition, stability_level, dominant_state, fsf_signal_class, gene_set_id) %>%
  arrange(p_adjust_BH, p_value, desc(enrichment_ratio), .by_group = TRUE) %>%
  slice_head(n = 30) %>%
  ungroup()

kegg_theme_summary <- top_kegg %>%
  group_by(condition, stability_level, dominant_state, fsf_signal_class, gene_set_id, KEGG_theme) %>%
  summarise(
    n_top_terms = n(),
    best_KEGG_FDR = min(p_adjust_BH, na.rm = TRUE),
    median_KEGG_enrichment_ratio = median(enrichment_ratio, na.rm = TRUE),
    representative_KEGG_terms = paste(head(unique(kegg_term), 10), collapse = ";"),
    KEGG_types = paste(sort(unique(enrichment_type)), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(condition, gene_set_id, best_KEGG_FDR)

kegg_bio <- kegg_theme_summary %>%
  group_by(condition, stability_level, dominant_state, fsf_signal_class, gene_set_id) %>%
  arrange(best_KEGG_FDR, desc(n_top_terms), .by_group = TRUE) %>%
  summarise(
    KEGG_theme = first(KEGG_theme),
    KEGG_theme_secondary = paste(head(unique(KEGG_theme), 5), collapse = ";"),
    n_KEGG_themes = n_distinct(KEGG_theme),
    best_KEGG_FDR = min(best_KEGG_FDR, na.rm = TRUE),
    representative_KEGG_terms = paste(head(unique(representative_KEGG_terms), 5), collapse = "|"),
    .groups = "drop"
  ) %>%
  right_join(
    manifest %>% select(condition, stability_level, dominant_state, fsf_signal_class, gene_set_id, n_genes, n_annotated, annotated_fraction),
    by = c("condition","stability_level","dominant_state","fsf_signal_class","gene_set_id")
  ) %>%
  left_join(
    kegg_sum %>%
      group_by(gene_set_id) %>%
      summarise(
        significant_KEGG_terms_FDR005 = sum(significant_terms_FDR005, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "gene_set_id"
  ) %>%
  mutate(
    KEGG_theme = ifelse(is.na(KEGG_theme), "No significant KEGG enrichment", KEGG_theme),
    KEGG_theme_secondary = ifelse(is.na(KEGG_theme_secondary), NA_character_, KEGG_theme_secondary),
    n_KEGG_themes = ifelse(is.na(n_KEGG_themes), 0, n_KEGG_themes),
    significant_KEGG_terms_FDR005 = ifelse(is.na(significant_KEGG_terms_FDR005), 0, significant_KEGG_terms_FDR005)
  ) %>%
  arrange(condition, stability_level, dominant_state, fsf_signal_class)

make_interpretation <- function(condition, stability_level, dominant_state, fsf_signal_class, go_theme, kegg_theme) {
  cls <- paste(stability_level, dominant_state, fsf_signal_class, sep = " / ")
  paste0(
    condition, " class ", cls,
    " shows a signal-class-specific biological profile. ",
    "GO theme: ", go_theme, ". ",
    "KEGG theme: ", kegg_theme, "."
  )
}

final_table <- go_bio %>%
  select(
    condition, stability_level, dominant_state, fsf_signal_class, gene_set_id,
    n_genes, n_annotated, annotated_fraction,
    GO_theme, GO_theme_secondary, significant_GO_terms_FDR005,
    best_GO_FDR, representative_GO_terms
  ) %>%
  left_join(
    kegg_bio %>%
      select(
        gene_set_id,
        KEGG_theme, KEGG_theme_secondary,
        significant_KEGG_terms_FDR005,
        best_KEGG_FDR, representative_KEGG_terms
      ),
    by = "gene_set_id"
  ) %>%
  mutate(
    KEGG_theme = ifelse(is.na(KEGG_theme), "No significant KEGG enrichment", KEGG_theme),
    significant_KEGG_terms_FDR005 = ifelse(is.na(significant_KEGG_terms_FDR005), 0, significant_KEGG_terms_FDR005),
    biological_support = case_when(
      significant_GO_terms_FDR005 > 0 & significant_KEGG_terms_FDR005 > 0 ~ "GO_and_KEGG_supported",
      significant_GO_terms_FDR005 > 0 & significant_KEGG_terms_FDR005 == 0 ~ "GO_supported",
      significant_GO_terms_FDR005 == 0 & significant_KEGG_terms_FDR005 > 0 ~ "KEGG_supported",
      TRUE ~ "no_significant_enrichment"
    ),
    manuscript_interpretation = mapply(
      make_interpretation,
      condition, stability_level, dominant_state, fsf_signal_class,
      GO_theme, KEGG_theme,
      USE.NAMES = FALSE
    )
  ) %>%
  arrange(condition, stability_level, dominant_state, fsf_signal_class)

condition_architecture <- final_table %>%
  group_by(condition) %>%
  summarise(
    n_signal_classes = n(),
    total_genes = sum(n_genes, na.rm = TRUE),
    classes_with_GO_support = sum(significant_GO_terms_FDR005 > 0, na.rm = TRUE),
    classes_with_KEGG_support = sum(significant_KEGG_terms_FDR005 > 0, na.rm = TRUE),
    dominant_GO_themes = paste(head(unique(GO_theme[GO_theme != "No significant GO enrichment"]), 8), collapse = ";"),
    dominant_KEGG_themes = paste(head(unique(KEGG_theme[KEGG_theme != "No significant KEGG enrichment"]), 8), collapse = ";"),
    example_gene_sets = paste(head(unique(gene_set_id), 12), collapse = ";"),
    .groups = "drop"
  )

write_tsv(top_kegg, file.path(OUT_DIR, "TOP30_FSF_KEGG_terms_by_class.tsv"))
write_tsv(kegg_theme_summary, file.path(OUT_DIR, "FSF_KEGG_theme_summary_by_class.tsv"))
write_tsv(kegg_bio, file.path(OUT_DIR, "FSF_KEGG_biological_theme_summary.tsv"))

write_tsv(final_table, file.path(FINAL_DIR, "Table_SignalClass_Biology.tsv"))
write_tsv(condition_architecture, file.path(FINAL_DIR, "Table_Condition_SignalArchitecture.tsv"))

cat("KEGG biological theme rows:", nrow(kegg_bio), "\n")
cat("Final manuscript biology rows:", nrow(final_table), "\n")
cat("Condition architecture rows:", nrow(condition_architecture), "\n\n")

cat("Saved:\n")
cat(file.path(OUT_DIR, "FSF_KEGG_biological_theme_summary.tsv"), "\n")
cat(file.path(FINAL_DIR, "Table_SignalClass_Biology.tsv"), "\n")
cat(file.path(FINAL_DIR, "Table_Condition_SignalArchitecture.tsv"), "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")
sink()
