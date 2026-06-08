#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 23: Prepare Stable Signal Annotation Inputs
#
# Purpose:
#   Export clean condition-wise stable signal strata for
#   downstream functional annotation and enrichment.
#
# Input:
#   results/real_data/stable_signal_catalog/condition_stable_signal_catalog.tsv
#
# Outputs:
#   results/real_data/stable_signal_annotation/stable_signal_annotation_master.tsv
#   results/real_data/stable_signal_annotation/stable_signal_gene_sets/
#   results/real_data/stable_signal_annotation/stable_signal_gene_set_summary.tsv
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

INPUT_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/stable_signal_catalog/condition_stable_signal_catalog.tsv"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/stable_signal_annotation"
)

GENESET_DIR <- file.path(
  OUT_DIR,
  "stable_signal_gene_sets"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/23_prepare_stable_signal_annotation_inputs.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GENESET_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 23: Prepare Stable Signal Annotation Inputs\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

stable <- read_tsv(INPUT_FILE, show_col_types = FALSE)

required_cols <- c(
  "condition",
  "feature_id",
  "dominant_state",
  "SSI",
  "stability_deviation",
  "stability_level",
  "fsf_signal_class"
)

missing_cols <- setdiff(required_cols, colnames(stable))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

stable <- stable %>%
  mutate(
    condition = as.character(condition),
    feature_id = as.character(feature_id),
    dominant_state = as.character(dominant_state),
    stability_level = as.character(stability_level),
    fsf_signal_class = as.character(fsf_signal_class),
    SSI = as.numeric(SSI),
    stability_deviation = as.numeric(stability_deviation),
    gene_set_id = paste(condition, fsf_signal_class, sep = "__")
  ) %>%
  arrange(condition, fsf_signal_class, desc(SSI), stability_deviation, feature_id)

write_tsv(
  stable,
  file.path(OUT_DIR, "stable_signal_annotation_master.tsv")
)

summary_tbl <- stable %>%
  count(
    condition,
    stability_level,
    dominant_state,
    fsf_signal_class,
    gene_set_id,
    name = "n_genes"
  ) %>%
  arrange(condition, stability_level, dominant_state)

write_tsv(
  summary_tbl,
  file.path(OUT_DIR, "stable_signal_gene_set_summary.tsv")
)

cat("Stable signal gene-set summary:\n")
print(summary_tbl, n = 100)
cat("\n")

gene_sets <- unique(stable$gene_set_id)

for (gs in gene_sets) {

  gs_tbl <- stable %>%
    filter(gene_set_id == gs) %>%
    arrange(desc(SSI), stability_deviation, feature_id)

  safe_name <- str_replace_all(gs, "[^A-Za-z0-9_\\-]+", "_")

  write_tsv(
    gs_tbl,
    file.path(GENESET_DIR, paste0(safe_name, ".tsv"))
  )

  write_lines(
    gs_tbl$feature_id,
    file.path(GENESET_DIR, paste0(safe_name, ".genes.txt"))
  )

  cat("Saved gene set:", gs, "| n =", nrow(gs_tbl), "\n")
}

cat("\nSaved outputs to:\n")
cat(OUT_DIR, "\n")
cat(GENESET_DIR, "\n\n")

cat("Next required input for enrichment:\n")
cat("- gene-to-GO / gene-to-KEGG / eggNOG / InterPro annotation table\n")
cat("Without annotation mapping, enrichment cannot be done honestly.\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

