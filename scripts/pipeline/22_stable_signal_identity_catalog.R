#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 22: Stable Signal Identity Catalog
#
# Purpose:
#   Extract condition-wise stable signal identities from FSF
#   strata, focusing on stable up, stable down, and stable
#   constant classes.
#
# Defensive version:
#   Reconstructs dominant_signal_identity and fsf_signal_class
#   instead of trusting malformed upstream headers.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

INPUT_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/condition_fsf/condition_wise_fsf_metrics.tsv"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/stable_signal_catalog"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/22_stable_signal_identity_catalog.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 22: Stable Signal Identity Catalog\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

fsf_raw <- read_tsv(INPUT_FILE, show_col_types = FALSE)

cat("Input columns:\n")
print(colnames(fsf_raw))
cat("\n")

required_core_cols <- c(
  "condition",
  "feature_id",
  "n_perturbations",
  "p_up",
  "p_down",
  "p_const",
  "dominant_state",
  "SSI",
  "stability_deviation",
  "stability_level"
)

missing_core_cols <- setdiff(required_core_cols, colnames(fsf_raw))

if (length(missing_core_cols) > 0) {
  stop("Missing required core columns: ", paste(missing_core_cols, collapse = ", "))
}

# ------------------------------------------------------------
# Reconstruct clean identity columns
# ------------------------------------------------------------

clean_fsf <- fsf_raw %>%
  transmute(
    condition = as.character(condition),
    feature_id = as.character(feature_id),
    n_perturbations = as.integer(n_perturbations),
    p_up = as.numeric(p_up),
    p_down = as.numeric(p_down),
    p_const = as.numeric(p_const),
    dominant_state = as.character(dominant_state),
    SSI = as.numeric(SSI),
    stability_deviation = as.numeric(stability_deviation),
    stability_level = as.character(stability_level),

    dominant_signal_identity = case_when(
      stability_level == "instable" ~ "instable_no_stable_identity",
      dominant_state == "up" ~ "stable_up_identity",
      dominant_state == "down" ~ "stable_down_identity",
      dominant_state == "constant" ~ "stable_constant_identity",
      TRUE ~ "mixed_identity"
    ),

    fsf_signal_class = case_when(
      stability_level == "instable" ~ "instable",
      dominant_state == "up" ~ paste0(stability_level, "_up"),
      dominant_state == "down" ~ paste0(stability_level, "_down"),
      dominant_state == "constant" ~ paste0(stability_level, "_constant"),
      TRUE ~ "mixed"
    )
  )

# ------------------------------------------------------------
# Stable catalog
# ------------------------------------------------------------

stable_catalog <- clean_fsf %>%
  filter(stability_level %in% c("highly_stable", "stable")) %>%
  arrange(
    condition,
    dominant_state,
    desc(SSI),
    stability_deviation,
    feature_id
  )

write_tsv(
  stable_catalog,
  file.path(OUT_DIR, "condition_stable_signal_catalog.tsv")
)

# ------------------------------------------------------------
# Stable directional signals
# ------------------------------------------------------------

directional_stable <- stable_catalog %>%
  filter(dominant_state %in% c("up", "down")) %>%
  arrange(
    condition,
    dominant_state,
    desc(SSI),
    stability_deviation,
    feature_id
  )

write_tsv(
  directional_stable,
  file.path(OUT_DIR, "condition_stable_directional_signals.tsv")
)

# ------------------------------------------------------------
# Stable constant signals
# ------------------------------------------------------------

constant_stable <- stable_catalog %>%
  filter(dominant_state == "constant") %>%
  arrange(
    condition,
    desc(SSI),
    stability_deviation,
    feature_id
  )

write_tsv(
  constant_stable,
  file.path(OUT_DIR, "condition_stable_constant_signals.tsv")
)

# ------------------------------------------------------------
# Summary table
# ------------------------------------------------------------

summary_tbl <- stable_catalog %>%
  count(
    condition,
    stability_level,
    dominant_state,
    fsf_signal_class,
    name = "n_features"
  ) %>%
  group_by(condition) %>%
  mutate(
    proportion = n_features / sum(n_features)
  ) %>%
  ungroup() %>%
  arrange(condition, stability_level, dominant_state)

write_tsv(
  summary_tbl,
  file.path(OUT_DIR, "condition_stable_signal_summary.tsv")
)

# ------------------------------------------------------------
# Top directional stable signals
# ------------------------------------------------------------

top_directional <- directional_stable %>%
  group_by(condition, dominant_state) %>%
  arrange(
    desc(SSI),
    stability_deviation,
    feature_id,
    .by_group = TRUE
  ) %>%
  slice_head(n = 100) %>%
  ungroup() %>%
  select(
    condition,
    dominant_state,
    feature_id,
    n_perturbations,
    p_up,
    p_down,
    p_const,
    SSI,
    stability_deviation,
    stability_level,
    dominant_signal_identity,
    fsf_signal_class
  )

write_tsv(
  top_directional,
  file.path(OUT_DIR, "top100_condition_stable_directional_signals.tsv")
)

cat("Stable signal summary:\n")
print(summary_tbl, n = 100)

cat("\nOutput column check:\n")
print(colnames(top_directional))

cat("\nSaved outputs to:\n")
cat(OUT_DIR, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

