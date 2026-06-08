#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

GO_SIG <- file.path(ROOT_DIR, "results/real_data/enrichment_GO/FSF_significant_GO_enrichment_FDR005.tsv")
GO_SUM <- file.path(ROOT_DIR, "results/real_data/enrichment_GO/FSF_GO_enrichment_summary_by_gene_set.tsv")
MANIFEST <- file.path(ROOT_DIR, "results/real_data/enrichment_gene_sets/FSF_enrichment_gene_set_manifest.tsv")
OUT_DIR <- file.path(ROOT_DIR, "results/real_data/enrichment_GO_interpretation")
LOG_FILE <- file.path(ROOT_DIR, "results/logs/32_GO_enrichment_interpretation_tables.log")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 32: GO Biological Theme Summary\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

for (f in c(GO_SIG, GO_SUM, MANIFEST)) {
  if (!file.exists(f)) stop("Missing file: ", f)
}

go_sig <- read_tsv(GO_SIG, show_col_types = FALSE)
go_sum <- read_tsv(GO_SUM, show_col_types = FALSE)
manifest <- read_tsv(MANIFEST, show_col_types = FALSE)

assign_theme <- function(go_term) {
  case_when(
    go_term %in% c("GO:0022626","GO:0022625","GO:0022627","GO:0005840","GO:0003735","GO:0002181","GO:0006412","GO:0043043") ~ "Translation / ribosome",
    go_term %in% c("GO:0006612","GO:0006613","GO:0006614","GO:0045047","GO:0072599","GO:0070972") ~ "Protein targeting / secretion",
    go_term %in% c("GO:0006259","GO:0006260","GO:0006271","GO:0006273","GO:0006281","GO:0006974","GO:0000723","GO:0000724","GO:0000725","GO:0000726") ~ "Genome maintenance / DNA repair",
    go_term %in% c("GO:0005694","GO:0051276","GO:0000228","GO:0044427","GO:0044454","GO:0098813","GO:0005634","GO:0005654") ~ "Chromosome / nuclear organization",
    go_term %in% c("GO:0007049","GO:0022402","GO:0000278","GO:0000280","GO:0048285","GO:0140013","GO:0051321","GO:1903046","GO:1903047","GO:0000070","GO:0000075","GO:0000077") ~ "Cell cycle / cell division",
    go_term %in% c("GO:0035239","GO:0048729","GO:0060562","GO:0009887","GO:0048598","GO:0060541","GO:0001501","GO:0002009") ~ "Development / morphogenesis",
    go_term %in% c("GO:0007411","GO:0097485","GO:0031175","GO:0048812","GO:0061564") ~ "Neural / projection organization",
    go_term %in% c("GO:0040011","GO:0006935","GO:0050920","GO:0007610","GO:0008045","GO:0008038") ~ "Movement / behavior / taxis",
    go_term %in% c("GO:0007155","GO:0007156","GO:0098609","GO:0005912","GO:0098742") ~ "Cell adhesion / junction",
    go_term %in% c("GO:0016491","GO:0055114","GO:0044281","GO:0006082","GO:0019752","GO:0043436","GO:0071704","GO:0044238") ~ "Metabolism / oxidation-reduction",
    str_detect(go_term, "^GO:00") ~ "Broad regulatory / cellular process",
    TRUE ~ "Unclassified / manual review"
  )
}

top_go <- go_sig %>%
  mutate(biological_theme = assign_theme(go_term)) %>%
  group_by(condition, stability_level, dominant_state, fsf_signal_class, gene_set_id) %>%
  arrange(p_adjust_BH, p_value, desc(enrichment_ratio), .by_group = TRUE) %>%
  slice_head(n = 50) %>%
  ungroup()

theme_summary <- top_go %>%
  group_by(condition, stability_level, dominant_state, fsf_signal_class, gene_set_id, biological_theme) %>%
  summarise(
    n_top_terms = n(),
    best_FDR = min(p_adjust_BH, na.rm = TRUE),
    median_enrichment_ratio = median(enrichment_ratio, na.rm = TRUE),
    representative_GO_terms = paste(head(unique(go_term), 10), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(condition, gene_set_id, best_FDR)

go_bio <- theme_summary %>%
  group_by(condition, stability_level, dominant_state, fsf_signal_class, gene_set_id) %>%
  arrange(best_FDR, desc(n_top_terms), .by_group = TRUE) %>%
  summarise(
    GO_theme = first(biological_theme),
    GO_theme_secondary = paste(head(unique(biological_theme), 5), collapse = ";"),
    n_GO_themes = n_distinct(biological_theme),
    best_GO_FDR = min(best_FDR, na.rm = TRUE),
    representative_GO_terms = paste(head(unique(representative_GO_terms), 5), collapse = "|"),
    .groups = "drop"
  ) %>%
  right_join(
    manifest %>% select(condition, stability_level, dominant_state, fsf_signal_class, gene_set_id, n_genes, n_annotated, annotated_fraction),
    by = c("condition","stability_level","dominant_state","fsf_signal_class","gene_set_id")
  ) %>%
  left_join(
    go_sum %>% select(gene_set_id, significant_GO_terms_FDR005),
    by = "gene_set_id"
  ) %>%
  mutate(
    GO_theme = ifelse(is.na(GO_theme), "No significant GO enrichment", GO_theme),
    GO_theme_secondary = ifelse(is.na(GO_theme_secondary), NA_character_, GO_theme_secondary),
    n_GO_themes = ifelse(is.na(n_GO_themes), 0, n_GO_themes),
    significant_GO_terms_FDR005 = ifelse(is.na(significant_GO_terms_FDR005), 0, significant_GO_terms_FDR005)
  ) %>%
  arrange(condition, stability_level, dominant_state, fsf_signal_class)

write_tsv(top_go, file.path(OUT_DIR, "TOP50_FSF_GO_terms_by_class.tsv"))
write_tsv(theme_summary, file.path(OUT_DIR, "FSF_GO_theme_summary_by_class.tsv"))
write_tsv(go_bio, file.path(OUT_DIR, "FSF_GO_biological_theme_summary.tsv"))

cat("GO biological theme rows:", nrow(go_bio), "\n")
cat("Saved:", file.path(OUT_DIR, "FSF_GO_biological_theme_summary.tsv"), "\n")
cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")
sink()
