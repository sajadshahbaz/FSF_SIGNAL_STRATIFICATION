#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

ROOT <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

INPUT <- file.path(
  ROOT,
  "results/real_data/class_annotation/FSF_condition_class_annotation_master.tsv"
)

OUTDIR <- file.path(
  ROOT,
  "results/archive/pre_repair_manuscript/figures/main"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

df <- read_tsv(INPUT, show_col_types = FALSE)

needed <- c("condition", "feature_id", "SSI")
missing <- setdiff(needed, colnames(df))

if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

df_clean <- df %>%
  distinct(condition, feature_id, SSI) %>%
  mutate(
    condition = factor(condition, levels = c("ht", "lt", "osm", "des", "gam", "uv"))
  )

p <- ggplot(
  df_clean,
  aes(
    x = SSI,
    colour = condition
  )
) +
  stat_ecdf(linewidth = 0.9) +
  scale_x_continuous(
    limits = c(0.3, 1.02),
    breaks = c(0.4, 0.6, 0.8, 1.0)
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = "Figure 4. Empirical cumulative distribution of SSI by condition",
    x = "Signal Stratification Index (SSI)",
    y = "Cumulative proportion of features",
    colour = "Condition"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 17, hjust = 0),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.major.y = element_line(color = "grey94", linewidth = 0.25)
  )

ggsave(
  file.path(OUTDIR, "Figure4_ECDF_SSI_by_condition.pdf"),
  p,
  width = 8.4,
  height = 5.4,
  device = "pdf"
)

if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(
    file.path(OUTDIR, "Figure4_ECDF_SSI_by_condition.png"),
    p,
    width = 8.4,
    height = 5.4,
    dpi = 300,
    device = ragg::agg_png
  )
} else {
  ggsave(
    file.path(OUTDIR, "Figure4_ECDF_SSI_by_condition.png"),
    p,
    width = 8.4,
    height = 5.4,
    dpi = 300,
    device = "png"
  )
}

write_tsv(
  df_clean,
  file.path(OUTDIR, "Figure4_source_data.tsv")
)

cat("Generated final Figure 4:\n")
cat(file.path(OUTDIR, "Figure4_ECDF_SSI_by_condition.pdf"), "\n")
cat(file.path(OUTDIR, "Figure4_ECDF_SSI_by_condition.png"), "\n")
