#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 29: Join Stable Signal Strata with Annotation Master
#
# Purpose:
#   Combine condition-wise stable signal strata with the clean
#   integrated annotation master.
#
# Inputs:
#   - condition_stable_signal_catalog.tsv
#   - FSF_v1_annotation_master.tsv
#
# Outputs:
#   - condition_stable_signal_annotated_catalog.tsv
#   - condition_stable_signal_annotation_coverage.tsv
#   - top_annotated_stable_directional_signals.tsv
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

STABLE_SIGNAL_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/stable_signal_catalog/condition_stable_signal_catalog.tsv"
)

ANNOTATION_MASTER_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/annotation_master/FSF_v1_annotation_master.tsv"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/stable_signal_annotation_integrated"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/29_join_stable_signals_with_annotation.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 29: Join Stable Signals with Annotation\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

if (!file.exists(STABLE_SIGNAL_FILE)) {
  stop("Stable signal file not found: ", STABLE_SIGNAL_FILE)
}

if (!file.exists(ANNOTATION_MASTER_FILE)) {
  stop("Annotation master file not found: ", ANNOTATION_MASTER_FILE)
}

stable <- read_tsv(STABLE_SIGNAL_FILE, show_col_types = FALSE)
anno <- read_tsv(ANNOTATION_MASTER_FILE, show_col_types = FALSE)

required_stable_cols <- c(
  "condition",
  "feature_id",
  "dominant_state",
  "SSI",
  "stability_deviation",
  "stability_level",
  "fsf_signal_class"
)

required_anno_cols <- c(
  "feature_id",
  "locus_tag",
  "protein_id",
  "protein_length",
  "translation_status",
  "eggnog_annotated",
  "interpro_annotated",
  "pfam_annotated",
  "annotation_source_count",
  "annotation_status"
)

missing_stable <- setdiff(required_stable_cols, colnames(stable))
missing_anno <- setdiff(required_anno_cols, colnames(anno))

if (length(missing_stable) > 0) {
  stop("Missing stable signal columns: ", paste(missing_stable, collapse = ", "))
}

if (length(missing_anno) > 0) {
  stop("Missing annotation master columns: ", paste(missing_anno, collapse = ", "))
}

cat("Stable signal rows:", nrow(stable), "\n")
cat("Unique stable signal features:", n_distinct(stable$feature_id), "\n")
cat("Annotation master rows:", nrow(anno), "\n")
cat("Unique annotation features:", n_distinct(anno$feature_id), "\n\n")

# ------------------------------------------------------------
# 1. Join
# ------------------------------------------------------------

annotated_stable <- stable %>%
  mutate(
    feature_id = as.character(feature_id),
    condition = as.character(condition),
    dominant_state = as.character(dominant_state),
    stability_level = as.character(stability_level),
    fsf_signal_class = as.character(fsf_signal_class)
  ) %>%
  left_join(
    anno %>% mutate(feature_id = as.character(feature_id)),
    by = "feature_id"
  ) %>%
  mutate(
    has_annotation = annotation_source_count > 0,
    stable_signal_group = paste(condition, fsf_signal_class, sep = "__")
  ) %>%
  arrange(
    condition,
    fsf_signal_class,
    desc(annotation_source_count),
    desc(SSI),
    stability_deviation,
    feature_id
  )

missing_annotation_rows <- annotated_stable %>%
  filter(is.na(annotation_source_count))

if (nrow(missing_annotation_rows) > 0) {
  write_tsv(
    missing_annotation_rows,
    file.path(OUT_DIR, "stable_signals_missing_from_annotation_master.tsv")
  )
  stop("Some stable signal rows did not match annotation master. Inspect stable_signals_missing_from_annotation_master.tsv")
}

write_tsv(
  annotated_stable,
  file.path(OUT_DIR, "condition_stable_signal_annotated_catalog.tsv")
)

# ------------------------------------------------------------
# 2. Annotation coverage per stable signal stratum
# ------------------------------------------------------------

coverage <- annotated_stable %>%
  group_by(
    condition,
    stability_level,
    dominant_state,
    fsf_signal_class,
    stable_signal_group
  ) %>%
  summarise(
    n_genes = n(),
    eggnog_annotated = sum(eggnog_annotated, na.rm = TRUE),
    interpro_annotated = sum(interpro_annotated, na.rm = TRUE),
    pfam_annotated = sum(pfam_annotated, na.rm = TRUE),
    at_least_one_annotation = sum(annotation_source_count >= 1, na.rm = TRUE),
    at_least_two_annotation_layers = sum(annotation_source_count >= 2, na.rm = TRUE),
    all_three_annotation_layers = sum(annotation_source_count == 3, na.rm = TRUE),
    unannotated = sum(annotation_source_count == 0, na.rm = TRUE),
    mean_annotation_source_count = mean(annotation_source_count, na.rm = TRUE),
    mean_SSI = mean(SSI, na.rm = TRUE),
    median_SSI = median(SSI, na.rm = TRUE),
    mean_stability_deviation = mean(stability_deviation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    eggnog_fraction = eggnog_annotated / n_genes,
    interpro_fraction = interpro_annotated / n_genes,
    pfam_fraction = pfam_annotated / n_genes,
    annotated_fraction = at_least_one_annotation / n_genes,
    unannotated_fraction = unannotated / n_genes
  ) %>%
  arrange(condition, stability_level, dominant_state)

write_tsv(
  coverage,
  file.path(OUT_DIR, "condition_stable_signal_annotation_coverage.tsv")
)

# ------------------------------------------------------------
# 3. Directional stable signals only
# ------------------------------------------------------------

directional_annotated <- annotated_stable %>%
  filter(
    dominant_state %in% c("up", "down"),
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

write_tsv(
  directional_annotated,
  file.path(OUT_DIR, "condition_stable_directional_annotated_signals.tsv")
)

top_directional <- directional_annotated %>%
  group_by(condition, dominant_state) %>%
  slice_head(n = 100) %>%
  ungroup()

write_tsv(
  top_directional,
  file.path(OUT_DIR, "top100_annotated_stable_directional_signals.tsv")
)

# ------------------------------------------------------------
# 4. Compact biological interpretation table
# ------------------------------------------------------------

compact_cols <- c(
  "condition",
  "dominant_state",
  "feature_id",
  "SSI",
  "stability_deviation",
  "stability_level",
  "fsf_signal_class",
  "protein_id",
  "protein_length",
  "annotation_source_count",
  "annotation_status",
  "eggnog_description",
  "eggnog_preferred_name",
  "eggnog_go",
  "eggnog_kegg_ko",
  "eggnog_kegg_pathway",
  "interpro_ids",
  "interpro_descriptions",
  "interpro_go",
  "pfam_domains",
  "pfam_descriptions"
)

compact <- annotated_stable %>%
  select(any_of(compact_cols)) %>%
  arrange(
    condition,
    dominant_state,
    desc(annotation_source_count),
    desc(SSI),
    stability_deviation,
    feature_id
  )

write_tsv(
  compact,
  file.path(OUT_DIR, "condition_stable_signal_annotated_compact.tsv")
)

top_compact_directional <- compact %>%
  filter(dominant_state %in% c("up", "down"), annotation_source_count >= 1) %>%
  group_by(condition, dominant_state) %>%
  slice_head(n = 50) %>%
  ungroup()

write_tsv(
  top_compact_directional,
  file.path(OUT_DIR, "top50_annotated_stable_directional_signals_compact.tsv")
)

# ------------------------------------------------------------
# 5. Print report
# ------------------------------------------------------------

cat("Stable signal annotation coverage:\n")
print(coverage, n = 100)
cat("\n")

cat("Directional annotated signal rows:", nrow(directional_annotated), "\n")
cat("Top directional rows:", nrow(top_directional), "\n\n")

cat("Saved outputs:\n")
cat(file.path(OUT_DIR, "condition_stable_signal_annotated_catalog.tsv"), "\n")
cat(file.path(OUT_DIR, "condition_stable_signal_annotation_coverage.tsv"), "\n")
cat(file.path(OUT_DIR, "condition_stable_directional_annotated_signals.tsv"), "\n")
cat(file.path(OUT_DIR, "top100_annotated_stable_directional_signals.tsv"), "\n")
cat(file.path(OUT_DIR, "condition_stable_signal_annotated_compact.tsv"), "\n")
cat(file.path(OUT_DIR, "top50_annotated_stable_directional_signals_compact.tsv"), "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

