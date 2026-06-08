#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Install it with: install.packages('yaml')")
  }
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required. Install it with: install.packages('digest')")
  }
})

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

config_file <- file.path("config", "paths.yml")

if (!file.exists(config_file)) {
  stop("Missing config file: ", config_file)
}

cfg <- yaml::read_yaml(config_file)

project_root <- normalizePath(cfg$project$root, mustWork = TRUE)

if (!dir.exists(project_root)) {
  stop("Project root does not exist: ", project_root)
}

setwd(project_root)

required_dirs <- unlist(cfg$directories, use.names = TRUE)

for (d in required_dirs) {
  dir.create(file.path(project_root, d), recursive = TRUE, showWarnings = FALSE)
}

log_dir <- file.path(project_root, cfg$directories$logs)
script_dir <- file.path(project_root, cfg$directories$scripts)

log_file <- file.path(log_dir, paste0("00_setup_project_", timestamp, ".log"))
session_file <- file.path(log_dir, paste0("sessionInfo_00_setup_project_", timestamp, ".txt"))
checksum_file <- file.path(log_dir, paste0("script_checksum_00_setup_project_", timestamp, ".txt"))

script_path <- file.path(script_dir, "00_setup_project.R")

cat("FSF Signal Stratification Project setup\n", file = log_file)
cat("Timestamp: ", timestamp, "\n", file = log_file, append = TRUE)
cat("Project root: ", project_root, "\n", file = log_file, append = TRUE)
cat("Config file: ", normalizePath(config_file), "\n\n", file = log_file, append = TRUE)

cat("Directory check:\n", file = log_file, append = TRUE)

for (nm in names(required_dirs)) {
  d <- file.path(project_root, required_dirs[[nm]])
  status <- ifelse(dir.exists(d), "OK", "MISSING")
  cat(nm, "\t", d, "\t", status, "\n", file = log_file, append = TRUE)
}

cat("\nRules:\n", file = log_file, append = TRUE)
cat("Hard-coded paths allowed: ", cfg$rules$hardcoded_paths_allowed, "\n", file = log_file, append = TRUE)
cat("Old project import requires review: ", cfg$rules$old_project_import_requires_review, "\n", file = log_file, append = TRUE)
cat("Framework role: ", cfg$rules$framework_role, "\n", file = log_file, append = TRUE)

capture.output(sessionInfo(), file = session_file)

if (file.exists(script_path)) {
  checksum <- digest::digest(file = script_path, algo = "sha256")
  cat("Script: ", script_path, "\n", file = checksum_file)
  cat("SHA256: ", checksum, "\n", file = checksum_file, append = TRUE)
}

cat("Setup completed successfully.\n")
cat("Log saved to: ", log_file, "\n")
cat("Session info saved to: ", session_file, "\n")
cat("Checksum saved to: ", checksum_file, "\n")
