#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 27: Prepare Annotation Commands
#
# Purpose:
#   Prepare reproducible eggNOG / InterProScan / Pfam annotation
#   commands for the clean FSF protein universe.
#
# This step does NOT run heavy annotation.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

PROTEIN_FASTA <- file.path(
  ROOT_DIR,
  "results/real_data/annotation_proteins/FSF_v1_clean_annotation_universe_proteins.faa"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/annotation_commands"
)

ANNOT_OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/functional_annotation"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/27_prepare_annotation_commands.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ANNOT_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 27: Prepare Annotation Commands\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

if (!file.exists(PROTEIN_FASTA)) {
  stop("Protein FASTA not found: ", PROTEIN_FASTA)
}

n_proteins <- as.integer(system(
  paste("grep -c '^>'", shQuote(PROTEIN_FASTA)),
  intern = TRUE
))

cat("Protein FASTA:\n")
cat(PROTEIN_FASTA, "\n")
cat("Protein records:", n_proteins, "\n\n")

# ------------------------------------------------------------
# Tool availability check
# ------------------------------------------------------------

check_tool <- function(cmd) {
  path <- Sys.which(cmd)
  tibble(
    tool = cmd,
    available = path != "",
    path = ifelse(path == "", NA_character_, path)
  )
}

tool_check <- bind_rows(
  check_tool("emapper.py"),
  check_tool("interproscan.sh"),
  check_tool("hmmscan"),
  check_tool("diamond")
)

write_tsv(
  tool_check,
  file.path(OUT_DIR, "annotation_tool_availability.tsv")
)

cat("Tool availability:\n")
print(tool_check)
cat("\n")

# ------------------------------------------------------------
# Output locations
# ------------------------------------------------------------

EGGNOG_DIR <- file.path(ANNOT_OUT_DIR, "eggnog")
INTERPRO_DIR <- file.path(ANNOT_OUT_DIR, "interproscan")
PFAM_DIR <- file.path(ANNOT_OUT_DIR, "pfam")

dir.create(EGGNOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(INTERPRO_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PFAM_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Command templates
# ------------------------------------------------------------

eggnog_cmd <- paste(
  "emapper.py",
  "-i", shQuote(PROTEIN_FASTA),
  "--itype proteins",
  "--cpu 16",
  "-m diamond",
  "--output FSF_v1_eggnog",
  "--output_dir", shQuote(EGGNOG_DIR),
  "--override"
)

interpro_cmd <- paste(
  "interproscan.sh",
  "-i", shQuote(PROTEIN_FASTA),
  "-d", shQuote(INTERPRO_DIR),
  "-f TSV,GFF3",
  "-goterms",
  "-pa",
  "-cpu 16"
)

pfam_cmd <- paste(
  "hmmscan",
  "--cpu 16",
  "--domtblout", shQuote(file.path(PFAM_DIR, "FSF_v1_pfam.domtblout")),
  "PFAM_A_HMM_PATH_TO_REPLACE",
  shQuote(PROTEIN_FASTA),
  ">", shQuote(file.path(PFAM_DIR, "FSF_v1_pfam.hmmscan.log")),
  "2>&1"
)

commands_tbl <- tibble(
  annotation_layer = c("eggNOG", "InterProScan", "Pfam"),
  command = c(eggnog_cmd, interpro_cmd, pfam_cmd),
  output_dir = c(EGGNOG_DIR, INTERPRO_DIR, PFAM_DIR),
  status = c(
    "ready_if_emapper_database_configured",
    "ready_if_interproscan_installed",
    "requires_PFAM_A_HMM_PATH_TO_REPLACE"
  )
)

write_tsv(
  commands_tbl,
  file.path(OUT_DIR, "FSF_v1_annotation_command_templates.tsv")
)

# ------------------------------------------------------------
# Write runnable shell template
# ------------------------------------------------------------

shell_file <- file.path(OUT_DIR, "run_FSF_v1_functional_annotation_TEMPLATE.sh")

shell_lines <- c(
  "#!/usr/bin/env bash",
  "set -euo pipefail",
  "",
  "# ============================================================",
  "# FSF SIGNAL STRATIFICATION v1.0",
  "# Functional annotation command template",
  "#",
  "# Edit PFAM_A_HMM path before running Pfam.",
  "# Run only after confirming tool/database installation.",
  "# ============================================================",
  "",
  paste0("ROOT_DIR=", shQuote(ROOT_DIR)),
  paste0("PROTEIN_FASTA=", shQuote(PROTEIN_FASTA)),
  paste0("ANNOT_OUT_DIR=", shQuote(ANNOT_OUT_DIR)),
  "",
  "mkdir -p \"$ANNOT_OUT_DIR\"",
  "",
  "echo 'Protein FASTA:' \"$PROTEIN_FASTA\"",
  "grep -c '^>' \"$PROTEIN_FASTA\"",
  "",
  "echo 'Checking tools...'",
  "command -v emapper.py || true",
  "command -v interproscan.sh || true",
  "command -v hmmscan || true",
  "command -v diamond || true",
  "",
  "echo '============================================================'",
  "echo 'eggNOG-mapper command'",
  "echo '============================================================'",
  eggnog_cmd,
  "",
  "echo '============================================================'",
  "echo 'InterProScan command'",
  "echo '============================================================'",
  interpro_cmd,
  "",
  "echo '============================================================'",
  "echo 'Pfam hmmscan command'",
  "echo '============================================================'",
  "# Replace PFAM_A_HMM_PATH_TO_REPLACE with your Pfam-A.hmm path before running.",
  paste0("# ", pfam_cmd),
  "",
  "echo 'Annotation template finished.'"
)

write_lines(shell_lines, shell_file)
Sys.chmod(shell_file, mode = "0755")

# ------------------------------------------------------------
# Write README
# ------------------------------------------------------------

readme_file <- file.path(OUT_DIR, "README_annotation_commands.md")

readme_lines <- c(
  "# FSF v1.0 Functional Annotation Command Preparation",
  "",
  "## Input protein FASTA",
  paste0("`", PROTEIN_FASTA, "`"),
  "",
  paste0("Protein records: ", n_proteins),
  "",
  "## Purpose",
  "",
  "This directory contains reproducible command templates for annotating the clean FSF v1.0 protein universe.",
  "",
  "The protein universe was derived from the 15,187 FSF-analyzed features and should be used instead of old Phase 4/DTC-FSF shortlisted annotation outputs.",
  "",
  "## Main command template",
  "",
  paste0("`", shell_file, "`"),
  "",
  "## Annotation layers",
  "",
  "- eggNOG-mapper: broad orthology, GO, KEGG, COG-style annotation",
  "- InterProScan: domain/family signatures, GO terms, pathways where available",
  "- Pfam/HMMER: Pfam domain validation",
  "",
  "## Important warning",
  "",
  "Do not use the old Phase 4 shortlisted BRAKER protein annotation for this FSF paper. It belongs to a different feature-selection universe and would bias the biological interpretation.",
  "",
  "## Next step",
  "",
  "Run or adapt the command template after confirming database paths and installed tools."
)

write_lines(readme_lines, readme_file)

cat("Annotation commands prepared:\n")
print(commands_tbl)
cat("\n")

cat("Saved:\n")
cat(file.path(OUT_DIR, "annotation_tool_availability.tsv"), "\n")
cat(file.path(OUT_DIR, "FSF_v1_annotation_command_templates.tsv"), "\n")
cat(shell_file, "\n")
cat(readme_file, "\n\n")

cat("Next step:\n")
cat("Inspect tool availability and confirm database paths before running annotation.\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

