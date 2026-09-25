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

ARCHIVE <- file.path(
  ROOT,
  "results/archive/pre_repair_manuscript/figures/archive/Figure4_replaced_versions"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ARCHIVE, recursive = TRUE, showWarnings = FALSE)

old_fig4 <- list.files(OUTDIR, pattern = "^Figure4", full.names = TRUE)
if (length(old_fig4) > 0) {
  file.copy(old_fig4, ARCHIVE, overwrite = TRUE)
  file.remove(old_fig4)
}

df <- read_tsv(INPUT, show_col_types = FALSE)

needed <- c("condition", "feature_id", "SSI", "stability_level")
missing <- setdiff(needed, colnames(df))

if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

df_clean <- df %>%
  distinct(condition, feature_id, SSI, stability_level) %>%
  mutate(
    condition = factor(condition, levels = c("lt", "ht", "osm", "des", "uv", "gam")),
    stability_level_clean = case_when(
      stability_level == "highly_stable" ~ "Highly stable",
      stability_level == "stable" ~ "Stable",
      stability_level == "weakly_stable_transitional" ~ "Transitional",
      stability_level == "instable" ~ "Instability",
      TRUE ~ stability_level
    ),
    stability_level_clean = factor(
      stability_level_clean,
      levels = c("Highly stable", "Stable", "Transitional", "Instability")
    )
  )

cols <- c(
  "Highly stable" = "#4D9221",
  "Stable" = "#A1D99B",
  "Transitional" = "#7B3294",
  "Instability" = "#B2182B"
)

p <- ggplot(
  df_clean,
  aes(
    x = SSI,
    y = condition,
    colour = stability_level_clean
  )
) +
  geom_jitter(
    height = 0.18,
    width = 0,
    alpha = 0.24,
    size = 0.45
  ) +
  stat_summary(
    aes(group = condition),
    fun = median,
    geom = "point",
    shape = 21,
    size = 3.1,
    fill = "white",
    colour = "black",
    stroke = 0.8
  ) +
  scale_colour_manual(values = cols, drop = FALSE, name = "Stability level") +
  scale_x_continuous(
    limits = c(0.3, 1.02),
    breaks = c(0.4, 0.6, 0.8, 1.0)
  ) +
  labs(
    title = "Figure 4. Condition-specific SSI distributions",
    x = "Signal Stratification Index (SSI)",
    y = "Condition"
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
  file.path(OUTDIR, "Figure4_condition_specific_SSI_distributions.pdf"),
  p,
  width = 8.8,
  height = 5.4,
  device = cairo_pdf
)

ggsave(
  file.path(OUTDIR, "Figure4_condition_specific_SSI_distributions.png"),
  p,
  width = 8.8,
  height = 5.4,
  dpi = 300
)

write_tsv(
  df_clean,
  file.path(OUTDIR, "Figure4_source_data.tsv")
)

cat("Generated final Figure 4:\n")
cat(file.path(OUTDIR, "Figure4_condition_specific_SSI_distributions.pdf"), "\n")
cat(file.path(OUTDIR, "Figure4_condition_specific_SSI_distributions.png"), "\n")
cat("Old Figure4 files archived in:\n")
cat(ARCHIVE, "\n")
