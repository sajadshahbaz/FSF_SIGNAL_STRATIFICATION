#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(fs)
})

ROOT <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

MANUSCRIPT_DIR <- file.path(ROOT, "results/manuscript")

DIRS <- c(
  "figures/main",
  "figures/supplementary",
  "tables/main",
  "tables/supplementary",
  "logs",
  "source_index"
)

for (d in DIRS) {
  dir_create(file.path(MANUSCRIPT_DIR, d))
}

log_file <- file.path(MANUSCRIPT_DIR, "logs", "step36_collect_manuscript_materials.log")
sink(log_file, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 36: Collect Manuscript Materials\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

copy_if_exists <- function(src, dest_dir, new_name = NULL) {
  if (!file.exists(src)) {
    cat("[MISSING]", src, "\n")
    return(tibble(source = src, destination = NA_character_, status = "missing"))
  }
  dest <- file.path(dest_dir, ifelse(is.null(new_name), basename(src), new_name))
  file.copy(src, dest, overwrite = TRUE)
  cat("[COPIED]", src, "->", dest, "\n")
  tibble(source = src, destination = dest, status = "copied")
}

copy_pattern <- function(pattern, dest_dir, tag = NULL) {
  files <- Sys.glob(pattern)
  if (length(files) == 0) {
    cat("[NO MATCH]", pattern, "\n")
    return(tibble(source = pattern, destination = NA_character_, status = "no_match"))
  }
  bind_rows(lapply(files, function(f) {
    new_name <- if (is.null(tag)) basename(f) else paste0(tag, "__", basename(f))
    copy_if_exists(f, dest_dir, new_name)
  }))
}

main_tables_dir <- file.path(MANUSCRIPT_DIR, "tables/main")
supp_tables_dir <- file.path(MANUSCRIPT_DIR, "tables/supplementary")
main_fig_dir <- file.path(MANUSCRIPT_DIR, "figures/main")
supp_fig_dir <- file.path(MANUSCRIPT_DIR, "figures/supplementary")

records <- list()

# =========================
# MAIN TABLES
# =========================

main_tables <- c(
  "results/real_data/final_manuscript_tables/Table_SignalClass_Biology_CURATED.tsv",
  "results/real_data/final_manuscript_tables/Table_Condition_SignalArchitecture_CURATED.tsv",
  "results/real_data/final_manuscript_tables/Table_PaperFocus_SignalClasses_CURATED.tsv",
  "results/real_data/enrichment_gene_sets/FSF_enrichment_gene_set_summary.tsv",
  "results/real_data/enrichment_GO_interpretation/FSF_GO_biological_theme_summary.tsv",
  "results/real_data/enrichment_KEGG_interpretation/FSF_KEGG_biological_theme_summary.tsv"
)

for (f in main_tables) {
  records[[length(records) + 1]] <- copy_if_exists(file.path(ROOT, f), main_tables_dir)
}

# =========================
# SUPPLEMENTARY TABLES
# =========================

supp_tables <- c(
  "results/real_data/class_annotation/FSF_condition_class_annotation_master.tsv",
  "results/real_data/class_annotation/FSF_condition_class_annotation_summary.tsv",
  "results/real_data/annotation_master/FSF_v1_annotation_master.tsv",
  "results/real_data/annotation_master/FSF_v1_annotation_coverage_summary.tsv",
  "results/real_data/annotation_master/FSF_v1_annotation_layer_statistics.tsv",
  "results/real_data/enrichment_GO/FSF_GO_annotation_QC.tsv",
  "results/real_data/enrichment_GO/FSF_GO_enrichment_summary_by_gene_set.tsv",
  "results/real_data/enrichment_GO/FSF_significant_GO_enrichment_FDR005.tsv",
  "results/real_data/enrichment_KEGG/FSF_KEGG_annotation_QC.tsv",
  "results/real_data/enrichment_KEGG/FSF_KEGG_enrichment_summary_by_gene_set.tsv",
  "results/real_data/enrichment_KEGG/FSF_significant_KEGG_enrichment_FDR005.tsv",
  "results/real_data/final_biological_theme_curation/FSF_curated_signal_class_biology.tsv",
  "results/real_data/final_biological_theme_curation/FSF_curated_condition_architecture_summary.tsv",
  "results/real_data/final_biological_theme_curation/FSF_curated_theme_frequency.tsv"
)

for (f in supp_tables) {
  records[[length(records) + 1]] <- copy_if_exists(file.path(ROOT, f), supp_tables_dir)
}

# Large raw enrichment tables stay supplementary, because apparently reviewers enjoy megatables.
records[[length(records) + 1]] <- copy_if_exists(
  file.path(ROOT, "results/real_data/enrichment_GO/FSF_all_GO_enrichment_results.tsv"),
  supp_tables_dir
)

records[[length(records) + 1]] <- copy_if_exists(
  file.path(ROOT, "results/real_data/enrichment_KEGG/FSF_all_KEGG_enrichment_results.tsv"),
  supp_tables_dir
)

# =========================
# MAIN FIGURES
# Existing figures only. Generation comes after collection.
# =========================

main_figure_patterns <- c(
  "results/real_data/stability_landscape/*.{pdf,png}",
  "results/real_data/*figure*.{pdf,png}",
  "results/figures/*.{pdf,png}",
  "results/tables/*probability*triangle*.{pdf,png}"
)

for (p in main_figure_patterns) {
  records[[length(records) + 1]] <- copy_pattern(file.path(ROOT, p), main_fig_dir)
}

# =========================
# SUPPLEMENTARY FIGURES
# =========================

supp_figure_patterns <- c(
  "results/real_data/**/plots/*.{pdf,png}",
  "results/real_data/**/*plot*.{pdf,png}",
  "results/real_data/**/*heatmap*.{pdf,png}",
  "results/real_data/**/*barplot*.{pdf,png}",
  "results/synthetic/**/*.{pdf,png}",
  "results/benchmark*/**/*.{pdf,png}"
)

for (p in supp_figure_patterns) {
  records[[length(records) + 1]] <- copy_pattern(file.path(ROOT, p), supp_fig_dir)
}

# =========================
# SYNTHETIC / BENCHMARK TABLES
# =========================

synthetic_patterns <- c(
  "results/synthetic/**/*.tsv",
  "results/benchmark*/**/*.tsv",
  "results/tables/FSF_v1_final_benchmark_report.tsv",
  "results/tables/FSF_v1_key_results_for_manuscript.tsv",
  "results/tables/fsf_probability_triangle_ready_table.tsv"
)

for (p in synthetic_patterns) {
  records[[length(records) + 1]] <- copy_pattern(file.path(ROOT, p), supp_tables_dir, tag = "synthetic_or_benchmark")
}

# =========================
# Make manuscript index
# =========================

index_tbl <- bind_rows(records) %>%
  mutate(
    exists_after_copy = ifelse(!is.na(destination), file.exists(destination), FALSE),
    file_size_bytes = ifelse(exists_after_copy, file.info(destination)$size, NA_real_)
  )

write_tsv(
  index_tbl,
  file.path(MANUSCRIPT_DIR, "source_index", "manuscript_materials_index.tsv")
)

summary_tbl <- index_tbl %>%
  count(status, exists_after_copy, name = "n_files")

write_tsv(
  summary_tbl,
  file.path(MANUSCRIPT_DIR, "source_index", "manuscript_materials_collection_summary.tsv")
)

cat("\nCollection summary:\n")
print(summary_tbl, n = 100)

cat("\nManuscript directory:\n")
cat(MANUSCRIPT_DIR, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()
