#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

ANNOTATION_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/annotation_master/FSF_v1_annotation_master.tsv"
)

MANIFEST_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/enrichment_gene_sets/FSF_enrichment_gene_set_manifest.tsv"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/enrichment_KEGG"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/33_run_KEGG_enrichment.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 33: KEGG / KO Enrichment Analysis\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

for (f in c(ANNOTATION_FILE, MANIFEST_FILE)) {
  if (!file.exists(f)) stop("Missing required file: ", f)
}

annot <- read_tsv(ANNOTATION_FILE, show_col_types = FALSE)
manifest <- read_tsv(MANIFEST_FILE, show_col_types = FALSE)

required_manifest_cols <- c(
  "condition",
  "stability_level",
  "dominant_state",
  "fsf_signal_class",
  "gene_set_id",
  "n_genes",
  "n_annotated",
  "annotated_gene_file"
)

missing_manifest_cols <- setdiff(required_manifest_cols, colnames(manifest))
if (length(missing_manifest_cols) > 0) {
  stop("Manifest missing required columns: ", paste(missing_manifest_cols, collapse = ", "))
}

missing_gene_files <- manifest$annotated_gene_file[!file.exists(manifest$annotated_gene_file)]
if (length(missing_gene_files) > 0) {
  stop("Missing annotated gene-set files:\n", paste(missing_gene_files, collapse = "\n"))
}

if (n_distinct(manifest$gene_set_id) != nrow(manifest)) {
  stop("Manifest gene_set_id values are not unique.")
}

cat("Manifest gene sets:", nrow(manifest), "\n")
cat("Unique manifest gene_set_id:", n_distinct(manifest$gene_set_id), "\n\n")

# -----------------------------
# Build KO mapping
# -----------------------------

gene2ko <- annot %>%
  select(feature_id, eggnog_kegg_ko) %>%
  filter(
    !is.na(eggnog_kegg_ko),
    eggnog_kegg_ko != "",
    eggnog_kegg_ko != "-"
  ) %>%
  separate_rows(eggnog_kegg_ko, sep = "[,;]") %>%
  mutate(
    kegg_term = str_trim(eggnog_kegg_ko),
    kegg_term = gsub("^ko:", "", kegg_term)
  ) %>%
  filter(kegg_term != "", kegg_term != "-") %>%
  select(feature_id, kegg_term) %>%
  distinct()

# -----------------------------
# Build KEGG pathway mapping
# -----------------------------

gene2path <- annot %>%
  select(feature_id, eggnog_kegg_pathway) %>%
  filter(
    !is.na(eggnog_kegg_pathway),
    eggnog_kegg_pathway != "",
    eggnog_kegg_pathway != "-"
  ) %>%
  separate_rows(eggnog_kegg_pathway, sep = "[,;]") %>%
  mutate(
    kegg_term = str_trim(eggnog_kegg_pathway)
  ) %>%
  filter(kegg_term != "", kegg_term != "-") %>%
  select(feature_id, kegg_term) %>%
  distinct()

write_tsv(gene2ko, file.path(OUT_DIR, "FSF_gene2KO.tsv"))
write_tsv(gene2path, file.path(OUT_DIR, "FSF_gene2Pathway.tsv"))

qc <- tibble(
  metric = c(
    "total_FSF_features",
    "KO_annotated_genes",
    "unique_KOs",
    "pathway_annotated_genes",
    "unique_pathways",
    "main_gene_sets_tested"
  ),
  value = c(
    n_distinct(annot$feature_id),
    n_distinct(gene2ko$feature_id),
    n_distinct(gene2ko$kegg_term),
    n_distinct(gene2path$feature_id),
    n_distinct(gene2path$kegg_term),
    nrow(manifest)
  )
)

write_tsv(qc, file.path(OUT_DIR, "FSF_KEGG_annotation_QC.tsv"))

cat("KEGG annotation QC:\n")
print(qc, n = 100)
cat("\n")

run_enrichment <- function(query_genes, gene2term, gene_set_id, enrichment_type) {

  universe <- unique(gene2term$feature_id)
  query_genes <- intersect(unique(query_genes), universe)

  if (length(query_genes) < 5) return(tibble())

  term_universe <- gene2term %>%
    filter(feature_id %in% universe) %>%
    distinct(feature_id, kegg_term)

  N <- length(universe)
  n <- length(query_genes)

  term_counts <- term_universe %>%
    group_by(kegg_term) %>%
    summarise(
      background_term_genes = n_distinct(feature_id),
      .groups = "drop"
    )

  query_counts <- term_universe %>%
    filter(feature_id %in% query_genes) %>%
    group_by(kegg_term) %>%
    summarise(
      overlap_genes = n_distinct(feature_id),
      query_gene_ids = paste(sort(unique(feature_id)), collapse = ";"),
      .groups = "drop"
    )

  if (nrow(query_counts) == 0) return(tibble())

  query_counts %>%
    inner_join(term_counts, by = "kegg_term") %>%
    mutate(
      enrichment_type = enrichment_type,
      gene_set_id = gene_set_id,
      query_size = n,
      universe_size = N,
      enrichment_ratio =
        (overlap_genes / query_size) /
        (background_term_genes / universe_size),
      p_value = phyper(
        q = overlap_genes - 1,
        m = background_term_genes,
        n = universe_size - background_term_genes,
        k = query_size,
        lower.tail = FALSE
      ),
      p_adjust_BH = p.adjust(p_value, method = "BH")
    ) %>%
    arrange(p_adjust_BH, p_value, desc(enrichment_ratio))
}

all_results <- list()

for (i in seq_len(nrow(manifest))) {

  gs <- manifest$gene_set_id[i]
  query <- read_lines(manifest$annotated_gene_file[i])

  ko_res <- run_enrichment(
    query_genes = query,
    gene2term = gene2ko,
    gene_set_id = gs,
    enrichment_type = "KEGG_KO"
  )

  path_res <- run_enrichment(
    query_genes = query,
    gene2term = gene2path,
    gene_set_id = gs,
    enrichment_type = "KEGG_pathway"
  )

  res <- bind_rows(ko_res, path_res)

  if (nrow(res) == 0) next

  res <- res %>%
    mutate(
      condition = manifest$condition[i],
      stability_level = manifest$stability_level[i],
      dominant_state = manifest$dominant_state[i],
      fsf_signal_class = manifest$fsf_signal_class[i],
      n_genes_in_set = manifest$n_genes[i],
      n_annotated_in_set = manifest$n_annotated[i]
    ) %>%
    select(
      enrichment_type,
      condition,
      stability_level,
      dominant_state,
      fsf_signal_class,
      gene_set_id,
      kegg_term,
      query_size,
      universe_size,
      overlap_genes,
      background_term_genes,
      enrichment_ratio,
      p_value,
      p_adjust_BH,
      query_gene_ids,
      n_genes_in_set,
      n_annotated_in_set
    )

  all_results[[gs]] <- res

  write_tsv(
    res,
    file.path(OUT_DIR, paste0(gs, "__KEGG_enrichment.tsv"))
  )
}

if (length(all_results) == 0) {
  stop("No KEGG enrichment results were generated.")
}

kegg_all <- bind_rows(all_results)

kegg_sig <- kegg_all %>%
  filter(p_adjust_BH <= 0.05) %>%
  arrange(enrichment_type, condition, stability_level, dominant_state, p_adjust_BH, p_value)

top20 <- kegg_all %>%
  group_by(enrichment_type, condition, stability_level, dominant_state, fsf_signal_class, gene_set_id) %>%
  arrange(p_adjust_BH, p_value, desc(enrichment_ratio), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()

summary_tbl <- kegg_all %>%
  group_by(enrichment_type, condition, stability_level, dominant_state, fsf_signal_class, gene_set_id) %>%
  summarise(
    tested_terms = n_distinct(kegg_term),
    significant_terms_FDR005 = sum(p_adjust_BH <= 0.05, na.rm = TRUE),
    min_p_adjust_BH = min(p_adjust_BH, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(enrichment_type, condition, stability_level, dominant_state, fsf_signal_class)

write_tsv(kegg_all, file.path(OUT_DIR, "FSF_all_KEGG_enrichment_results.tsv"))
write_tsv(kegg_sig, file.path(OUT_DIR, "FSF_significant_KEGG_enrichment_FDR005.tsv"))
write_tsv(top20, file.path(OUT_DIR, "TOP20_FSF_KEGG_enrichment_by_gene_set.tsv"))
write_tsv(summary_tbl, file.path(OUT_DIR, "FSF_KEGG_enrichment_summary_by_gene_set.tsv"))

cat("KEGG enrichment completed.\n")
cat("Total enrichment rows:", nrow(kegg_all), "\n")
cat("Significant FDR<=0.05 rows:", nrow(kegg_sig), "\n")
cat("Gene sets with results:", n_distinct(kegg_all$gene_set_id), "\n")
cat("Gene sets with significant results:", n_distinct(kegg_sig$gene_set_id), "\n\n")

cat("Saved outputs to:\n")
cat(OUT_DIR, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()
