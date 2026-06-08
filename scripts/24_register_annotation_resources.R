#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 24: Register Annotation Resources
#
# Purpose:
#   Locate and register available functional annotation files
#   for downstream stable-signal interpretation.
#
# This step does NOT perform enrichment.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/annotation_registry"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/24_register_annotation_resources.log"
)

CONFIG_FILE <- file.path(
  ROOT_DIR,
  "config/annotation_resources.tsv"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(CONFIG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 24: Register Annotation Resources\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

search_roots <- c(
  "/media/saji/5E06441D0643F5152/phase_1",
  "/media/saji/5E06441D0643F5152/phase_2",
  "/media/saji/5E06441D0643F5152/phase_3",
  "/media/saji/5E06441D0643F5152/phase_4",
  "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"
)

search_roots <- search_roots[dir.exists(search_roots)]

cat("Search roots:\n")
print(search_roots)
cat("\n")

patterns <- c(
  "eggnog",
  "emapper",
  "interpro",
  "iprscan",
  "pfam",
  "go",
  "kegg",
  "annotation",
  "functional"
)

all_files <- unlist(
  lapply(search_roots, function(root) {
    list.files(
      root,
      recursive = TRUE,
      full.names = TRUE,
      all.files = FALSE
    )
  }),
  use.names = FALSE
)

candidate_files <- all_files[
  str_detect(
    basename(all_files),
    regex(paste(patterns, collapse = "|"), ignore_case = TRUE)
  )
]

candidate_files <- candidate_files[
  file.info(candidate_files)$isdir == FALSE
]

registry <- tibble(
  file_path = candidate_files,
  file_name = basename(candidate_files),
  file_size_mb = round(file.info(candidate_files)$size / 1024^2, 3),
  modified_time = as.character(file.info(candidate_files)$mtime),
  guessed_type = case_when(
    str_detect(file_name, regex("eggnog|emapper", ignore_case = TRUE)) ~ "eggNOG/emapper",
    str_detect(file_name, regex("interpro|iprscan", ignore_case = TRUE)) ~ "InterProScan",
    str_detect(file_name, regex("pfam", ignore_case = TRUE)) ~ "Pfam",
    str_detect(file_name, regex("kegg", ignore_case = TRUE)) ~ "KEGG",
    str_detect(file_name, regex("go", ignore_case = TRUE)) ~ "GO",
    str_detect(file_name, regex("annotation|functional", ignore_case = TRUE)) ~ "general_annotation",
    TRUE ~ "unknown"
  )
) %>%
  arrange(guessed_type, file_name)

write_tsv(
  registry,
  file.path(OUT_DIR, "annotation_resource_candidates.tsv")
)

write_tsv(
  registry,
  CONFIG_FILE
)

cat("Candidate annotation files found:", nrow(registry), "\n\n")
print(registry, n = 200)

cat("\nSaved:\n")
cat(file.path(OUT_DIR, "annotation_resource_candidates.tsv"), "\n")
cat(CONFIG_FILE, "\n\n")

cat("Next step:\n")
cat("Manually inspect the registry and select the correct gene-to-function mapping file.\n")
cat("Do not run enrichment until feature IDs can be matched to annotation IDs.\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

