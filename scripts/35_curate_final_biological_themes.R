#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

INPUT_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/final_manuscript_tables/Table_SignalClass_Biology.tsv"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/final_biological_theme_curation"
)

FINAL_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/final_manuscript_tables"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/35_curate_final_biological_themes.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FINAL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 35: Curate Final Biological Themes\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

if (!file.exists(INPUT_FILE)) {
  stop("Missing input file: ", INPUT_FILE)
}

bio <- read_tsv(INPUT_FILE, show_col_types = FALSE)

collapse_theme <- function(go_theme, kegg_theme, go_terms, kegg_terms) {
  
  txt <- paste(
    go_theme,
    kegg_theme,
    go_terms,
    kegg_terms,
    sep = " "
  ) %>%
    str_to_lower()
  
  case_when(
    str_detect(txt, "translation|ribosome|ribosomal|protein targeting|secretion") ~
      "Translation / protein targeting",
    
    str_detect(txt, "dna repair|genome maintenance|replication|recombination|chromosome|nuclear organization|nucleotide excision|homologous recombination|mismatch repair") ~
      "Genome maintenance / nuclear regulation",
    
    str_detect(txt, "cell cycle|cell division|mitotic|meiosis|spindle|chromosome segregation") ~
      "Cell cycle / division",
    
    str_detect(txt, "proteostasis|protein turnover|proteasome|ubiquitin|folding|chaperone") ~
      "Proteostasis / protein turnover",
    
    str_detect(txt, "energy|metabolism|oxidation|oxidative|oxidoreduction|mitochond|carbon metabolism|oxidative phosphorylation") ~
      "Energy metabolism / redox regulation",
    
    str_detect(txt, "stress signaling|signal transduction|mapk|calcium|phosphatidylinositol|signaling") ~
      "Stress signaling / signal transduction",
    
    str_detect(txt, "adhesion|junction|extracellular|ecm|cell adhesion") ~
      "Cell adhesion / extracellular interaction",
    
    str_detect(txt, "development|morphogenesis|tissue|structural remodeling|organ development") ~
      "Development / morphogenesis",
    
    str_detect(txt, "neural|neuron|synapse|axon|projection|behavior|movement|taxis|locomotion") ~
      "Neural / behavioral organization",
    
    str_detect(txt, "no significant go enrichment no significant kegg enrichment") ~
      "No clear enrichment",
    
    str_detect(txt, "no significant") & str_detect(txt, "broad regulatory|other kegg|manual review|unclassified") ~
      "Weak / nonspecific enrichment",
    
    str_detect(txt, "broad regulatory|other kegg|manual review|unclassified") ~
      "Broad cellular regulation",
    
    TRUE ~
      "Mixed functional program"
  )
}

make_architecture_type <- function(condition, stability_level, dominant_state, fsf_signal_class) {
  case_when(
    condition %in% c("lt", "osm", "ht") ~
      "Highly stable-dominated architecture",
    
    condition == "des" ~
      "Mixed stable-instabile architecture",
    
    condition %in% c("uv", "gam") ~
      "Multi-layer signal architecture",
    
    TRUE ~
      "Unassigned architecture"
  )
}

make_result_priority <- function(n_genes, n_annotated, go_n, kegg_n, final_theme) {
  case_when(
    final_theme == "No clear enrichment" ~ "low",
    go_n > 100 & kegg_n > 10 ~ "high",
    go_n > 20 | kegg_n > 5 ~ "moderate",
    go_n > 0 | kegg_n > 0 ~ "limited",
    TRUE ~ "low"
  )
}

make_interpretation <- function(condition, stability_level, dominant_state, fsf_signal_class, theme, support) {
  
  class_label <- paste(stability_level, dominant_state, fsf_signal_class, sep = " / ")
  
  paste0(
    toupper(condition),
    " ",
    class_label,
    " is assigned to ",
    theme,
    " based on combined GO and KEGG evidence. ",
    "Support level: ",
    support,
    "."
  )
}

curated <- bio %>%
  mutate(
    GO_theme = ifelse(is.na(GO_theme), "No significant GO enrichment", GO_theme),
    KEGG_theme = ifelse(is.na(KEGG_theme), "No significant KEGG enrichment", KEGG_theme),
    representative_GO_terms = ifelse(is.na(representative_GO_terms), "", representative_GO_terms),
    representative_KEGG_terms = ifelse(is.na(representative_KEGG_terms), "", representative_KEGG_terms),
    significant_GO_terms_FDR005 = ifelse(is.na(significant_GO_terms_FDR005), 0, significant_GO_terms_FDR005),
    significant_KEGG_terms_FDR005 = ifelse(is.na(significant_KEGG_terms_FDR005), 0, significant_KEGG_terms_FDR005),
    
    final_biological_theme = mapply(
      collapse_theme,
      GO_theme,
      KEGG_theme,
      representative_GO_terms,
      representative_KEGG_terms,
      USE.NAMES = FALSE
    ),
    
    signal_architecture_type = mapply(
      make_architecture_type,
      condition,
      stability_level,
      dominant_state,
      fsf_signal_class,
      USE.NAMES = FALSE
    ),
    
    evidence_priority = mapply(
      make_result_priority,
      n_genes,
      n_annotated,
      significant_GO_terms_FDR005,
      significant_KEGG_terms_FDR005,
      final_biological_theme,
      USE.NAMES = FALSE
    ),
    
    final_interpretation = mapply(
      make_interpretation,
      condition,
      stability_level,
      dominant_state,
      fsf_signal_class,
      final_biological_theme,
      evidence_priority,
      USE.NAMES = FALSE
    )
  ) %>%
  select(
    condition,
    signal_architecture_type,
    stability_level,
    dominant_state,
    fsf_signal_class,
    gene_set_id,
    n_genes,
    n_annotated,
    annotated_fraction,
    final_biological_theme,
    evidence_priority,
    GO_theme,
    KEGG_theme,
    significant_GO_terms_FDR005,
    significant_KEGG_terms_FDR005,
    best_GO_FDR,
    best_KEGG_FDR,
    representative_GO_terms,
    representative_KEGG_terms,
    final_interpretation,
    everything()
  )

condition_summary <- curated %>%
  group_by(condition, signal_architecture_type) %>%
  summarise(
    n_signal_classes = n(),
    total_genes_across_classes = sum(n_genes, na.rm = TRUE),
    classes_with_GO_support = sum(significant_GO_terms_FDR005 > 0, na.rm = TRUE),
    classes_with_KEGG_support = sum(significant_KEGG_terms_FDR005 > 0, na.rm = TRUE),
    high_priority_classes = sum(evidence_priority == "high", na.rm = TRUE),
    moderate_priority_classes = sum(evidence_priority == "moderate", na.rm = TRUE),
    dominant_biological_themes = paste(
      head(unique(final_biological_theme[final_biological_theme != "No clear enrichment"]), 10),
      collapse = ";"
    ),
    representative_gene_sets = paste(head(unique(gene_set_id), 20), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(condition)

theme_frequency <- curated %>%
  count(final_biological_theme, evidence_priority, sort = TRUE)

paper_focus_table <- curated %>%
  filter(
    evidence_priority %in% c("high", "moderate") |
      condition %in% c("des", "uv", "gam")
  ) %>%
  arrange(
    condition,
    desc(evidence_priority == "high"),
    desc(significant_GO_terms_FDR005 + significant_KEGG_terms_FDR005),
    stability_level,
    dominant_state
  )

write_tsv(
  curated,
  file.path(OUT_DIR, "FSF_curated_signal_class_biology.tsv")
)

write_tsv(
  condition_summary,
  file.path(OUT_DIR, "FSF_curated_condition_architecture_summary.tsv")
)

write_tsv(
  theme_frequency,
  file.path(OUT_DIR, "FSF_curated_theme_frequency.tsv")
)

write_tsv(
  paper_focus_table,
  file.path(OUT_DIR, "FSF_paper_focus_signal_classes.tsv")
)

write_tsv(
  curated,
  file.path(FINAL_DIR, "Table_SignalClass_Biology_CURATED.tsv")
)

write_tsv(
  condition_summary,
  file.path(FINAL_DIR, "Table_Condition_SignalArchitecture_CURATED.tsv")
)

write_tsv(
  paper_focus_table,
  file.path(FINAL_DIR, "Table_PaperFocus_SignalClasses_CURATED.tsv")
)

cat("Curated signal-class rows:", nrow(curated), "\n")
cat("Condition summary rows:", nrow(condition_summary), "\n")
cat("Paper-focus rows:", nrow(paper_focus_table), "\n\n")

cat("Theme frequency:\n")
print(theme_frequency, n = 100)

cat("\nSaved outputs:\n")
cat(file.path(FINAL_DIR, "Table_SignalClass_Biology_CURATED.tsv"), "\n")
cat(file.path(FINAL_DIR, "Table_Condition_SignalArchitecture_CURATED.tsv"), "\n")
cat(file.path(FINAL_DIR, "Table_PaperFocus_SignalClasses_CURATED.tsv"), "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()
