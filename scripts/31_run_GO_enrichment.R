#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 31: GO Enrichment Analysis
#
# Clean final version:
#   - Uses only FSF_enrichment_gene_set_manifest.tsv
#   - Uses only manifest$annotated_gene_file
#   - Uses only 40 main class gene sets
#   - Canonicalizes GO IDs
#   - Replaces old output in-place
# ============================================================

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

UNIVERSE_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/enrichment_gene_sets/FSF_universe_annotated.genes.txt"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/enrichment_GO"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/31_run_GO_enrichment.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 31: GO Enrichment Analysis\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

for (f in c(ANNOTATION_FILE, MANIFEST_FILE, UNIVERSE_FILE)) {
  if (!file.exists(f)) {
    stop("Missing required file: ", f)
  }
}

anno <- read_tsv(ANNOTATION_FILE, show_col_types = FALSE)
manifest <- read_tsv(MANIFEST_FILE, show_col_types = FALSE)
universe <- read_lines(UNIVERSE_FILE)

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
  stop(
    "Manifest missing required columns: ",
    paste(missing_manifest_cols, collapse = ", ")
  )
}

missing_gene_files <- manifest$annotated_gene_file[!file.exists(manifest$annotated_gene_file)]

if (length(missing_gene_files) > 0) {
  stop(
    "Some annotated gene-set files are missing:\n",
    paste(missing_gene_files, collapse = "\n")
  )
}

cat("Universe genes:", length(universe), "\n")
cat("Manifest gene sets:", nrow(manifest), "\n")
cat("Unique manifest gene_set_id:", n_distinct(manifest$gene_set_id), "\n\n")

if (n_distinct(manifest$gene_set_id) != nrow(manifest)) {
  stop("Manifest gene_set_id values are not unique.")
}

canonicalize_go <- function(x) {
  str_extract(as.character(x), "GO:[0-9]{7}")
}

extract_go_terms <- function(tbl, source_col, sep_pattern, source_name) {
  tbl %>%
    filter(feature_id %in% universe) %>%
    filter(
      !is.na(.data[[source_col]]),
      .data[[source_col]] != "",
      .data[[source_col]] != "-"
    ) %>%
    select(feature_id, raw_go_term = all_of(source_col)) %>%
    separate_rows(raw_go_term, sep = sep_pattern) %>%
    mutate(
      raw_go_term = str_trim(raw_go_term),
      go_term = canonicalize_go(raw_go_term),
      source = source_name
    ) %>%
    filter(!is.na(go_term), go_term != "") %>%
    distinct(feature_id, go_term, source)
}

eggnog_go <- extract_go_terms(
  anno,
  source_col = "eggnog_go",
  sep_pattern = ",",
  source_name = "eggNOG"
)

interpro_go <- extract_go_terms(
  anno,
  source_col = "interpro_go",
  sep_pattern = ";",
  source_name = "InterPro"
)

gene2go_with_source <- bind_rows(eggnog_go, interpro_go) %>%
  distinct(feature_id, go_term, source)

gene2go <- gene2go_with_source %>%
  distinct(feature_id, go_term)

write_tsv(
  gene2go_with_source,
  file.path(OUT_DIR, "FSF_gene2GO_with_source.tsv")
)

write_tsv(
  gene2go,
  file.path(OUT_DIR, "FSF_gene2GO_combined.tsv")
)

go_qc <- tibble(
  metric = c(
    "universe_genes",
    "GO_annotated_genes",
    "unique_GO_terms_combined",
    "unique_GO_terms_eggNOG",
    "unique_GO_terms_InterPro",
    "gene_GO_pairs_combined",
    "gene_GO_pairs_with_source",
    "main_gene_sets_tested"
  ),
  value = c(
    length(universe),
    n_distinct(gene2go$feature_id),
    n_distinct(gene2go$go_term),
    n_distinct(eggnog_go$go_term),
    n_distinct(interpro_go$go_term),
    nrow(gene2go),
    nrow(gene2go_with_source),
    nrow(manifest)
  )
)

write_tsv(
  go_qc,
  file.path(OUT_DIR, "FSF_GO_annotation_QC.tsv")
)

cat("GO annotation QC:\n")
print(go_qc, n = 100)
cat("\n")

run_enrichment <- function(query_genes, universe_genes, gene2term, gene_set_id) {

  query_genes <- intersect(unique(query_genes), universe_genes)

  if (length(query_genes) < 5) {
    return(tibble())
  }

  term_universe <- gene2term %>%
    filter(feature_id %in% universe_genes) %>%
    distinct(feature_id, go_term)

  N <- length(unique(universe_genes))
  n <- length(unique(query_genes))

  term_counts <- term_universe %>%
    group_by(go_term) %>%
    summarise(
      background_term_genes = n_distinct(feature_id),
      .groups = "drop"
    )

  query_counts <- term_universe %>%
    filter(feature_id %in% query_genes) %>%
    group_by(go_term) %>%
    summarise(
      overlap_genes = n_distinct(feature_id),
      query_gene_ids = paste(sort(unique(feature_id)), collapse = ";"),
      .groups = "drop"
    )

  if (nrow(query_counts) == 0) {
    return(tibble())
  }

  query_counts %>%
    inner_join(term_counts, by = "go_term") %>%
    mutate(
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
  gene_file <- manifest$annotated_gene_file[i]

  query <- read_lines(gene_file)

  res <- run_enrichment(
    query_genes = query,
    universe_genes = universe,
    gene2term = gene2go,
    gene_set_id = gs
  )

  if (nrow(res) == 0) {
    next
  }

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
      condition,
      stability_level,
      dominant_state,
      fsf_signal_class,
      gene_set_id,
      go_term,
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
    file.path(OUT_DIR, paste0(gs, "__GO_enrichment.tsv"))
  )
}

if (length(all_results) == 0) {
  stop("No GO enrichment results were generated. Check gene2GO and gene-set files.")
}

go_all <- bind_rows(all_results)

go_sig <- go_all %>%
  filter(p_adjust_BH <= 0.05) %>%
  arrange(condition, stability_level, dominant_state, p_adjust_BH, p_value)

top20 <- go_all %>%
  group_by(condition, stability_level, dominant_state, fsf_signal_class, gene_set_id) %>%
  arrange(p_adjust_BH, p_value, desc(enrichment_ratio), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()

summary_tbl <- go_all %>%
  group_by(condition, stability_level, dominant_state, fsf_signal_class, gene_set_id) %>%
  summarise(
    tested_GO_terms = n_distinct(go_term),
    significant_GO_terms_FDR005 = sum(p_adjust_BH <= 0.05, na.rm = TRUE),
    min_p_adjust_BH = min(p_adjust_BH, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(condition, stability_level, dominant_state, fsf_signal_class)

write_tsv(
  go_all,
  file.path(OUT_DIR, "FSF_all_GO_enrichment_results.tsv")
)

write_tsv(
  go_sig,
  file.path(OUT_DIR, "FSF_significant_GO_enrichment_FDR005.tsv")
)

write_tsv(
  top20,
  file.path(OUT_DIR, "TOP20_FSF_GO_enrichment_by_gene_set.tsv")
)

write_tsv(
  summary_tbl,
  file.path(OUT_DIR, "FSF_GO_enrichment_summary_by_gene_set.tsv")
)

cat("GO enrichment completed.\n")
cat("Total enrichment rows:", nrow(go_all), "\n")
cat("Significant FDR<=0.05 rows:", nrow(go_sig), "\n")
cat("Gene sets with results:", n_distinct(go_all$gene_set_id), "\n")
cat("Gene sets with significant results:", n_distinct(go_sig$gene_set_id), "\n\n")

cat("Saved outputs to:\n")
cat(OUT_DIR, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()
