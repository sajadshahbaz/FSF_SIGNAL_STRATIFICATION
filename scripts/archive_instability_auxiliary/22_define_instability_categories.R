#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 22: Define Shared and Condition-Specific Instability Categories
#
# Purpose:
#   Classify instable genes into condition-specific and shared
#   instability categories.
#
# Input:
#   results/real_data/instability_drivers/gene_instability_burden.tsv
#
# Outputs:
#   results/real_data/instability_categories/instability_category_table.tsv
#   results/real_data/instability_categories/instability_category_summary.tsv
#   results/real_data/instability_categories/*_genes.tsv
#   results/real_data/figures/instability_category_barplot.png
#   results/logs/22_define_instability_categories.log
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

INPUT_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/instability_drivers/gene_instability_burden.tsv"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/instability_categories"
)

FIG_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/figures"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/22_define_instability_categories.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 22: Define Instability Categories\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

burden <- read_tsv(INPUT_FILE, show_col_types = FALSE)

required_cols <- c(
  "feature_id",
  "instability_count",
  "instability_conditions",
  "mean_SSI",
  "mean_stability_deviation",
  "max_stability_deviation"
)

missing_cols <- setdiff(required_cols, colnames(burden))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

burden <- burden %>%
  mutate(
    feature_id = as.character(feature_id),
    instability_conditions = as.character(instability_conditions),
    instability_count = as.integer(instability_count),
    mean_SSI = as.numeric(mean_SSI),
    mean_stability_deviation = as.numeric(mean_stability_deviation),
    max_stability_deviation = as.numeric(max_stability_deviation)
  )

# ------------------------------------------------------------
# Add condition flags
# ------------------------------------------------------------

cat_tbl <- burden %>%
  mutate(
    is_des = str_detect(instability_conditions, "(^|;)des(;|$)"),
    is_uv  = str_detect(instability_conditions, "(^|;)uv(;|$)"),
    is_gam = str_detect(instability_conditions, "(^|;)gam(;|$)"),
    is_lt  = str_detect(instability_conditions, "(^|;)lt(;|$)"),
    is_ht  = str_detect(instability_conditions, "(^|;)ht(;|$)"),
    is_osm = str_detect(instability_conditions, "(^|;)osm(;|$)")
  ) %>%
  mutate(
    instability_category = case_when(
      instability_count == 0 ~ "not_instable",

      is_uv & !is_des & !is_gam & !is_lt & !is_ht & !is_osm ~ "uv_specific_instable",
      is_des & !is_uv & !is_gam & !is_lt & !is_ht & !is_osm ~ "des_specific_instable",
      is_gam & !is_uv & !is_des & !is_lt & !is_ht & !is_osm ~ "gam_specific_instable",
      is_lt & !is_uv & !is_des & !is_gam & !is_ht & !is_osm ~ "lt_specific_instable",

      is_uv & is_des & !is_gam & !is_lt & !is_ht & !is_osm ~ "uv_des_shared_instable",

      is_uv & is_gam & !is_des & !is_lt & !is_ht & !is_osm ~ "uv_gam_shared_instable",
      is_des & is_gam & !is_uv & !is_lt & !is_ht & !is_osm ~ "des_gam_shared_instable",
      is_uv & is_lt & !is_des & !is_gam & !is_ht & !is_osm ~ "uv_lt_shared_instable",
      is_des & is_lt & !is_uv & !is_gam & !is_ht & !is_osm ~ "des_lt_shared_instable",
      is_gam & is_lt & !is_uv & !is_des & !is_ht & !is_osm ~ "gam_lt_shared_instable",

      instability_count >= 3 ~ "multi_condition_instable",

      TRUE ~ "other_shared_instable"
    )
  ) %>%
  arrange(
    desc(instability_count),
    desc(mean_stability_deviation),
    feature_id
  )

write_tsv(
  cat_tbl,
  file.path(OUT_DIR, "instability_category_table.tsv")
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary_tbl <- cat_tbl %>%
  count(instability_category, name = "n_genes") %>%
  mutate(proportion = n_genes / sum(n_genes)) %>%
  arrange(desc(n_genes))

write_tsv(
  summary_tbl,
  file.path(OUT_DIR, "instability_category_summary.tsv")
)

cat("Instability category summary:\n")
print(summary_tbl, n = 100)
cat("\n")

# ------------------------------------------------------------
# Export category-specific gene files
# ------------------------------------------------------------

categories <- unique(cat_tbl$instability_category)

for (cat_name in categories) {
  out_file <- file.path(
    OUT_DIR,
    paste0(cat_name, "_genes.tsv")
  )

  write_tsv(
    cat_tbl %>%
      filter(instability_category == cat_name),
    out_file
  )

  cat("Saved:", cat_name, "->", out_file, "\n")
}

cat("\n")

# ------------------------------------------------------------
# Focused manuscript categories
# ------------------------------------------------------------

focused_categories <- c(
  "uv_specific_instable",
  "des_specific_instable",
  "uv_des_shared_instable",
  "gam_specific_instable",
  "lt_specific_instable",
  "multi_condition_instable"
)

focused_tbl <- cat_tbl %>%
  filter(instability_category %in% focused_categories)

write_tsv(
  focused_tbl,
  file.path(OUT_DIR, "focused_instability_categories.tsv")
)

# ------------------------------------------------------------
# Barplot
# ------------------------------------------------------------

plot_tbl <- summary_tbl %>%
  filter(instability_category != "not_instable") %>%
  mutate(
    instability_category = factor(
      instability_category,
      levels = rev(instability_category)
    )
  )

p <- ggplot(
  plot_tbl,
  aes(
    x = instability_category,
    y = n_genes
  )
) +
  geom_col() +
  coord_flip() +
  geom_text(
    aes(label = n_genes),
    hjust = -0.15,
    size = 3.5
  ) +
  labs(
    title = "Shared and condition-specific instability categories",
    x = "Instability category",
    y = "Number of genes"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(FIG_DIR, "instability_category_barplot.png"),
  p,
  width = 8,
  height = 6,
  dpi = 300
)

cat("Saved instability category barplot.\n\n")

cat("Saved main outputs to:\n")
cat(OUT_DIR, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

