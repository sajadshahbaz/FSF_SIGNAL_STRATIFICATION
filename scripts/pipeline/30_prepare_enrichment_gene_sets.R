#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

INPUT_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/class_annotation/FSF_condition_class_annotation_master.tsv"
)

OUT_DIR <- file.path(ROOT_DIR, "results/real_data/enrichment_gene_sets")

MAIN_DIR <- file.path(OUT_DIR, "gene_sets/main_class_sets")
COLLAPSED_DIR <- file.path(OUT_DIR, "gene_sets/collapsed_stability_sets")
DIRECTIONAL_DIR <- file.path(OUT_DIR, "gene_sets/directional_stable_sets")

LOG_FILE <- file.path(ROOT_DIR, "results/logs/30_prepare_enrichment_gene_sets.log")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

if (dir.exists(file.path(OUT_DIR, "gene_sets"))) {
  unlink(file.path(OUT_DIR, "gene_sets"), recursive = TRUE, force = TRUE)
}

dir.create(MAIN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(COLLAPSED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DIRECTIONAL_DIR, recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 30: Prepare Enrichment-Ready Gene Sets\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

dat <- read_tsv(INPUT_FILE, show_col_types = FALSE)

dat <- dat %>%
  mutate(
    condition = as.character(condition),
    feature_id = as.character(feature_id),
    stability_level = as.character(stability_level),
    dominant_state = as.character(dominant_state),
    fsf_signal_class = as.character(fsf_signal_class),
    gene_set_id = paste(
      condition,
      stability_level,
      dominant_state,
      fsf_signal_class,
      sep = "__"
    ),
    gene_set_id = str_replace_all(gene_set_id, "[^A-Za-z0-9_]+", "_")
  )

cat("Input rows:", nrow(dat), "\n")
cat("Unique features:", n_distinct(dat$feature_id), "\n")
cat("Main class gene sets:", n_distinct(dat$gene_set_id), "\n\n")

universe_all <- dat %>% distinct(feature_id) %>% arrange(feature_id)

universe_annotated <- dat %>%
  filter(annotation_source_count >= 1) %>%
  distinct(feature_id) %>%
  arrange(feature_id)

write_lines(universe_all$feature_id, file.path(OUT_DIR, "FSF_universe_all_15187.genes.txt"))
write_lines(universe_annotated$feature_id, file.path(OUT_DIR, "FSF_universe_annotated.genes.txt"))

main_summary <- dat %>%
  group_by(condition, stability_level, dominant_state, fsf_signal_class, gene_set_id) %>%
  summarise(
    n_genes = n_distinct(feature_id),
    n_annotated = n_distinct(feature_id[annotation_source_count >= 1]),
    annotated_fraction = n_annotated / n_genes,
    .groups = "drop"
  ) %>%
  arrange(condition, stability_level, dominant_state, fsf_signal_class)

write_tsv(main_summary, file.path(OUT_DIR, "FSF_enrichment_gene_set_summary.tsv"))

for (gs in main_summary$gene_set_id) {
  genes_all <- dat %>%
    filter(gene_set_id == gs) %>%
    distinct(feature_id) %>%
    arrange(feature_id)

  genes_annotated <- dat %>%
    filter(gene_set_id == gs, annotation_source_count >= 1) %>%
    distinct(feature_id) %>%
    arrange(feature_id)

  write_lines(genes_all$feature_id, file.path(MAIN_DIR, paste0(gs, ".genes.txt")))
  write_lines(genes_annotated$feature_id, file.path(MAIN_DIR, paste0(gs, ".annotated.genes.txt")))
}

collision_check <- main_summary %>%
  count(gene_set_id, name = "n_rows") %>%
  filter(n_rows > 1)

write_tsv(collision_check, file.path(OUT_DIR, "FSF_gene_set_id_collision_check.tsv"))

if (nrow(collision_check) > 0) {
  stop("gene_set_id collision detected.")
}

collapsed <- dat %>%
  mutate(
    collapsed_gene_set_id = paste(condition, stability_level, sep = "__"),
    collapsed_gene_set_id = str_replace_all(collapsed_gene_set_id, "[^A-Za-z0-9_]+", "_")
  )

collapsed_summary <- collapsed %>%
  group_by(condition, stability_level, collapsed_gene_set_id) %>%
  summarise(
    n_genes = n_distinct(feature_id),
    n_annotated = n_distinct(feature_id[annotation_source_count >= 1]),
    annotated_fraction = n_annotated / n_genes,
    .groups = "drop"
  ) %>%
  arrange(condition, stability_level)

write_tsv(collapsed_summary, file.path(OUT_DIR, "FSF_collapsed_stability_gene_set_summary.tsv"))

for (gs in collapsed_summary$collapsed_gene_set_id) {
  genes_all <- collapsed %>%
    filter(collapsed_gene_set_id == gs) %>%
    distinct(feature_id) %>%
    arrange(feature_id)

  genes_annotated <- collapsed %>%
    filter(collapsed_gene_set_id == gs, annotation_source_count >= 1) %>%
    distinct(feature_id) %>%
    arrange(feature_id)

  write_lines(genes_all$feature_id, file.path(COLLAPSED_DIR, paste0(gs, ".genes.txt")))
  write_lines(genes_annotated$feature_id, file.path(COLLAPSED_DIR, paste0(gs, ".annotated.genes.txt")))
}

directional <- dat %>%
  filter(
    stability_level %in% c("highly_stable", "stable"),
    dominant_state %in% c("up", "down")
  ) %>%
  mutate(
    directional_gene_set_id = paste(
      condition,
      stability_level,
      dominant_state,
      fsf_signal_class,
      sep = "__"
    ),
    directional_gene_set_id = str_replace_all(directional_gene_set_id, "[^A-Za-z0-9_]+", "_")
  )

directional_summary <- directional %>%
  group_by(condition, stability_level, dominant_state, fsf_signal_class, directional_gene_set_id) %>%
  summarise(
    n_genes = n_distinct(feature_id),
    n_annotated = n_distinct(feature_id[annotation_source_count >= 1]),
    annotated_fraction = n_annotated / n_genes,
    .groups = "drop"
  ) %>%
  arrange(condition, stability_level, dominant_state, fsf_signal_class)

write_tsv(directional_summary, file.path(OUT_DIR, "FSF_directional_stable_gene_set_summary.tsv"))

for (gs in directional_summary$directional_gene_set_id) {
  genes_all <- directional %>%
    filter(directional_gene_set_id == gs) %>%
    distinct(feature_id) %>%
    arrange(feature_id)

  genes_annotated <- directional %>%
    filter(directional_gene_set_id == gs, annotation_source_count >= 1) %>%
    distinct(feature_id) %>%
    arrange(feature_id)

  write_lines(genes_all$feature_id, file.path(DIRECTIONAL_DIR, paste0(gs, ".genes.txt")))
  write_lines(genes_annotated$feature_id, file.path(DIRECTIONAL_DIR, paste0(gs, ".annotated.genes.txt")))
}

manifest <- main_summary %>%
  mutate(
    gene_file = file.path(MAIN_DIR, paste0(gene_set_id, ".genes.txt")),
    annotated_gene_file = file.path(MAIN_DIR, paste0(gene_set_id, ".annotated.genes.txt"))
  )

write_tsv(manifest, file.path(OUT_DIR, "FSF_enrichment_gene_set_manifest.tsv"))

cat("Main class gene sets:", nrow(main_summary), "\n")
cat("Collapsed stability gene sets:", nrow(collapsed_summary), "\n")
cat("Directional stable gene sets:", nrow(directional_summary), "\n")
cat("Universe all:", nrow(universe_all), "\n")
cat("Universe annotated:", nrow(universe_annotated), "\n\n")

cat("Main gene-set directory:\n", MAIN_DIR, "\n\n")
cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()
