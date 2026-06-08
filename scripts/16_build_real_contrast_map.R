#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 16: Build Real Contrast Map
#
# Purpose:
#   Construct biologically valid perturbation definitions
#   from the tardigrade metadata.
#
# Outputs:
#   real_contrast_map.tsv
#   real_contrast_sample_map.tsv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

META_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/input_validation/valid_real_metadata.tsv"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/design"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/16_build_real_contrast_map.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 16: Build Real Contrast Map\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

meta <- read_tsv(META_FILE, show_col_types = FALSE)

cat("Samples loaded:", nrow(meta), "\n\n")

contrast_map <- list()
sample_map <- list()

contrast_counter <- 1

# ----------------------------------------------------------
# DES
# ----------------------------------------------------------

for(block in c("fast","slow")) {

  ctrl <- meta %>%
    filter(
      condition == "des",
      group_type == "control",
      baseline_block == block
    )

  trt <- meta %>%
    filter(
      condition == "des",
      group_type == "treatment",
      baseline_block == block
    )

  cid <- sprintf("DES_%s", toupper(block))

  contrast_map[[length(contrast_map)+1]] <- tibble(
    contrast_id = cid,
    condition = "des",
    baseline_block = block,
    perturbation_type = "treatment_vs_control",
    n_control = nrow(ctrl),
    n_treatment = nrow(trt)
  )

  sample_map[[length(sample_map)+1]] <- bind_rows(
    ctrl %>% mutate(contrast_id = cid,
                    role = "control"),
    trt %>% mutate(contrast_id = cid,
                   role = "treatment")
  )
}

# ----------------------------------------------------------
# UV
# ----------------------------------------------------------

for(block in c("short","long")) {

  ctrl <- meta %>%
    filter(
      condition == "uv",
      group_type == "control",
      baseline_block == block
    )

  trt <- meta %>%
    filter(
      condition == "uv",
      group_type == "treatment",
      baseline_block == block
    )

  times <- sort(unique(trt$timepoint))

  for(tp in times) {

    trt_tp <- trt %>%
      filter(timepoint == tp)

    cid <- paste(
      "UV",
      toupper(block),
      tp,
      sep = "_"
    )

    contrast_map[[length(contrast_map)+1]] <- tibble(
      contrast_id = cid,
      condition = "uv",
      baseline_block = block,
      perturbation_type = "timepoint_vs_control",
      timepoint = tp,
      n_control = nrow(ctrl),
      n_treatment = nrow(trt_tp)
    )

    sample_map[[length(sample_map)+1]] <- bind_rows(
      ctrl %>% mutate(
        contrast_id = cid,
        role = "control"
      ),
      trt_tp %>% mutate(
        contrast_id = cid,
        role = "treatment"
      )
    )
  }
}

# ----------------------------------------------------------
# GAM
# ----------------------------------------------------------

ctrl <- meta %>%
  filter(
    condition == "gam",
    group_type == "control"
  )

trt_times <- meta %>%
  filter(
    condition == "gam",
    group_type == "treatment"
  ) %>%
  pull(timepoint) %>%
  unique() %>%
  sort()

for(tp in trt_times) {

  trt_tp <- meta %>%
    filter(
      condition == "gam",
      group_type == "treatment",
      timepoint == tp
    )

  cid <- paste0("GAM_", tp)

  contrast_map[[length(contrast_map)+1]] <- tibble(
    contrast_id = cid,
    condition = "gam",
    baseline_block = "pool",
    perturbation_type = "timepoint_vs_00h",
    timepoint = tp,
    n_control = nrow(ctrl),
    n_treatment = nrow(trt_tp)
  )

  sample_map[[length(sample_map)+1]] <- bind_rows(
    ctrl %>% mutate(
      contrast_id = cid,
      role = "control"
    ),
    trt_tp %>% mutate(
      contrast_id = cid,
      role = "treatment"
    )
  )
}

# ----------------------------------------------------------
# LT
# ----------------------------------------------------------

for(tp in c("02h","24h")) {

  ctrl <- meta %>%
    filter(
      condition == "lt",
      group_type == "control",
      timepoint == tp
    )

  trt <- meta %>%
    filter(
      condition == "lt",
      group_type == "treatment",
      timepoint == tp
    )

  cid <- paste0("LT_", tp)

  contrast_map[[length(contrast_map)+1]] <- tibble(
    contrast_id = cid,
    condition = "lt",
    baseline_block = "pool",
    perturbation_type = "timepoint_vs_control",
    timepoint = tp,
    n_control = nrow(ctrl),
    n_treatment = nrow(trt)
  )

  sample_map[[length(sample_map)+1]] <- bind_rows(
    ctrl %>% mutate(
      contrast_id = cid,
      role = "control"
    ),
    trt %>% mutate(
      contrast_id = cid,
      role = "treatment"
    )
  )
}

# ----------------------------------------------------------
# HT
# ----------------------------------------------------------

for(cond in c("ht","osm")) {

  ctrl <- meta %>%
    filter(
      condition == cond,
      group_type == "control"
    )

  trt <- meta %>%
    filter(
      condition == cond,
      group_type == "treatment"
    )

  cid <- toupper(cond)

  contrast_map[[length(contrast_map)+1]] <- tibble(
    contrast_id = cid,
    condition = cond,
    baseline_block = "pool",
    perturbation_type = "treatment_vs_control",
    n_control = nrow(ctrl),
    n_treatment = nrow(trt)
  )

  sample_map[[length(sample_map)+1]] <- bind_rows(
    ctrl %>% mutate(
      contrast_id = cid,
      role = "control"
    ),
    trt %>% mutate(
      contrast_id = cid,
      role = "treatment"
    )
  )
}

contrast_tbl <- bind_rows(contrast_map)
sample_tbl <- bind_rows(sample_map)

write_tsv(
  contrast_tbl,
  file.path(
    OUT_DIR,
    "real_contrast_map.tsv"
  )
)

write_tsv(
  sample_tbl,
  file.path(
    OUT_DIR,
    "real_contrast_sample_map.tsv"
  )
)

cat("Number of contrasts:",
    nrow(contrast_tbl),
    "\n\n")

print(contrast_tbl, n = 200)

cat("\nSaved:\n")
cat(file.path(
  OUT_DIR,
  "real_contrast_map.tsv"
), "\n")

cat(file.path(
  OUT_DIR,
  "real_contrast_sample_map.tsv"
), "\n\n")

cat("Finished:",
    as.character(Sys.time()),
    "\n")

cat("============================================================\n")

sink()

