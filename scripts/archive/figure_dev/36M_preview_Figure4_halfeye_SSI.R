#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

ROOT <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"
INFILE <- file.path(ROOT, "results/real_data/class_annotation/FSF_condition_class_annotation_master.tsv")

OUT <- file.path(ROOT, "results/archive/pre_repair_manuscript/figures/supplementary")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Clean only this preview output
file.remove(list.files(OUT, pattern = "^FigureS4_alt_halfeye_SSI", full.names = TRUE))

df <- read_tsv(INFILE, show_col_types = FALSE)

df_clean <- df %>%
  distinct(condition, feature_id, SSI, stability_deviation,
           stability_level, dominant_state, fsf_signal_class) %>%
  mutate(
    condition = factor(condition, levels = rev(c("lt", "ht", "osm", "des", "uv", "gam"))),
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

level_cols <- c(
  "Highly stable" = "#4D9221",
  "Stable" = "#A1D99B",
  "Transitional" = "#7B3294",
  "Instability" = "#B2182B"
)

theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 4, hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor.x = element_line(color = "grey95", linewidth = 0.15)
    )
}

# If ggdist is installed, use real half-eye.
# If not, fall back to clean horizontal density ridges using base ggplot2.
has_ggdist <- requireNamespace("ggdist", quietly = TRUE)

if (has_ggdist) {

  suppressPackageStartupMessages(library(ggdist))

  p <- ggplot(df_clean, aes(x = SSI, y = condition, fill = stability_level_clean)) +
    ggdist::stat_halfeye(
      adjust = 0.45,
      width = 0.75,
      justification = -0.15,
      .width = c(0.5, 0.8, 0.95),
      point_interval = median_qi,
      alpha = 0.75,
      slab_colour = "grey25",
      slab_linewidth = 0.25
    ) +
    scale_fill_manual(values = level_cols, drop = FALSE, name = "Stability level") +
    scale_x_continuous(limits = c(0.3, 1.02), breaks = c(0.4, 0.6, 0.8, 1.0)) +
    labs(
      title = "Figure S4 alternative. Half-eye SSI landscape by condition",
      x = "Signal Stratification Index (SSI)",
      y = "Condition"
    ) +
    theme_pub(12)

} else {

  message("Package 'ggdist' not installed. Using fallback jitter + boxplot preview.")

  p <- ggplot(df_clean, aes(x = SSI, y = condition, colour = stability_level_clean)) +
    geom_jitter(
      height = 0.16,
      width = 0,
      alpha = 0.14,
      size = 0.45
    ) +
    geom_boxplot(
      aes(fill = stability_level_clean),
      width = 0.28,
      outlier.shape = NA,
      alpha = 0.55,
      colour = "black",
      linewidth = 0.25
    ) +
    scale_colour_manual(values = level_cols, drop = FALSE, name = "Stability level") +
    scale_fill_manual(values = level_cols, drop = FALSE, name = "Stability level") +
    scale_x_continuous(limits = c(0.3, 1.02), breaks = c(0.4, 0.6, 0.8, 1.0)) +
    labs(
      title = "Figure S4 alternative. SSI landscape by condition",
      x = "Signal Stratification Index (SSI)",
      y = "Condition"
    ) +
    theme_pub(12)
}

ggsave(
  file.path(OUT, "FigureS4_alt_halfeye_SSI_landscape.pdf"),
  p,
  width = 8.8,
  height = 5.4,
  device = cairo_pdf
)

ggsave(
  file.path(OUT, "FigureS4_alt_halfeye_SSI_landscape.png"),
  p,
  width = 8.8,
  height = 5.4,
  dpi = 300
)

write_tsv(df_clean, file.path(OUT, "FigureS4_alt_halfeye_SSI_source_data.tsv"))

cat("Generated half-eye SSI preview:\n")
cat(file.path(OUT, "FigureS4_alt_halfeye_SSI_landscape.pdf"), "\n")
cat(file.path(OUT, "FigureS4_alt_halfeye_SSI_landscape.png"), "\n")
