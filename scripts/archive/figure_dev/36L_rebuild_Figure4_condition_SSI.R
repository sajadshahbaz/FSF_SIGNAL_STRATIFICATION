#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

ROOT <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"
INFILE <- file.path(ROOT, "results/real_data/class_annotation/FSF_condition_class_annotation_master.tsv")

MAIN_OUT <- file.path(ROOT, "results/archive/pre_repair_manuscript/figures/main")
SUPP_OUT <- file.path(ROOT, "results/archive/pre_repair_manuscript/figures/supplementary")
ARCHIVE_OUT <- file.path(ROOT, "results/archive/pre_repair_manuscript/figures/archive/Figure4_old_versions")

dir.create(MAIN_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(ARCHIVE_OUT, recursive = TRUE, showWarnings = FALSE)

# Archive old Figure 4/S4 files instead of deleting them, because apparently civilization requires guardrails.
old_main <- list.files(MAIN_OUT, pattern = "^Figure4", full.names = TRUE)
old_supp <- list.files(SUPP_OUT, pattern = "^FigureS4", full.names = TRUE)

if (length(old_main) > 0) {
  file.copy(old_main, ARCHIVE_OUT, overwrite = TRUE)
  file.remove(old_main)
}

if (length(old_supp) > 0) {
  file.copy(old_supp, ARCHIVE_OUT, overwrite = TRUE)
  file.remove(old_supp)
}

df <- read_tsv(INFILE, show_col_types = FALSE)

df_clean <- df %>%
  distinct(
    condition,
    feature_id,
    SSI,
    stability_deviation,
    stability_level,
    dominant_state,
    fsf_signal_class
  ) %>%
  mutate(
    condition = factor(condition, levels = c("ht", "lt", "osm", "des", "gam", "uv")),
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
    ),
    fsf_class_clean = case_when(
      stability_level == "highly_stable" & dominant_state == "up" ~ "Highly stable up",
      stability_level == "highly_stable" & dominant_state == "constant" ~ "Highly stable constant",
      stability_level == "highly_stable" & dominant_state == "down" ~ "Highly stable down",
      stability_level == "stable" & dominant_state == "up" ~ "Stable up",
      stability_level == "stable" & dominant_state == "constant" ~ "Stable constant",
      stability_level == "stable" & dominant_state == "down" ~ "Stable down",
      stability_level == "weakly_stable_transitional" & dominant_state == "up" ~ "Transitional up",
      stability_level == "weakly_stable_transitional" & dominant_state == "constant" ~ "Transitional constant",
      stability_level == "weakly_stable_transitional" & dominant_state == "down" ~ "Transitional down",
      stability_level == "instable" ~ "Instability",
      TRUE ~ fsf_signal_class
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
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank()
    )
}

save_fig <- function(p, outdir, name, w, h) {
  ggsave(file.path(outdir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(outdir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

# ------------------------------------------------------------
# Figure 4A: SSI strip landscape
# ------------------------------------------------------------

set.seed(20260608)

p4a <- ggplot(
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
    alpha = 0.22,
    size = 0.45
  ) +
  stat_summary(
    aes(group = condition),
    fun = median,
    geom = "point",
    shape = 21,
    size = 3.2,
    fill = "white",
    colour = "black",
    stroke = 0.8
  ) +
  scale_colour_manual(values = level_cols, drop = FALSE, name = "Stability level") +
  scale_x_continuous(
    limits = c(0.3, 1.02),
    breaks = c(0.4, 0.6, 0.8, 1.0)
  ) +
  labs(
    title = "Figure 4A. Condition-specific SSI landscape",
    x = "Signal Stratification Index (SSI)",
    y = "Condition"
  ) +
  theme_pub(12)

save_fig(p4a, MAIN_OUT, "Figure4A_condition_specific_SSI_landscape", 8.8, 5.4)

# ------------------------------------------------------------
# Figure 4B: stability architecture composition
# ------------------------------------------------------------

p4b_data <- df_clean %>%
  count(condition, stability_level_clean, name = "n_features") %>%
  group_by(condition) %>%
  mutate(proportion = n_features / sum(n_features)) %>%
  ungroup()

p4b <- ggplot(
  p4b_data,
  aes(
    x = condition,
    y = proportion,
    fill = stability_level_clean
  )
) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = level_cols, drop = FALSE, name = "Stability level") +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = "Figure 4B. Condition-level stability architecture",
    x = "Condition",
    y = "Proportion of features"
  ) +
  theme_pub(12)

save_fig(p4b, MAIN_OUT, "Figure4B_condition_level_stability_architecture", 8.2, 5.2)

# ------------------------------------------------------------
# Figure S4: SSI ECDF
# ------------------------------------------------------------

pS4 <- ggplot(
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
    title = "Figure S4. Empirical cumulative distribution of SSI by condition",
    x = "Signal Stratification Index (SSI)",
    y = "Cumulative proportion of features",
    colour = "Condition"
  ) +
  theme_pub(12)

save_fig(pS4, SUPP_OUT, "FigureS4_ECDF_SSI_by_condition", 8.4, 5.4)

# ------------------------------------------------------------
# Source data
# ------------------------------------------------------------

write_tsv(df_clean, file.path(MAIN_OUT, "Figure4A_source_data.tsv"))
write_tsv(p4b_data, file.path(MAIN_OUT, "Figure4B_source_data.tsv"))

cat("Generated final Figure 4 set:\n")
cat(file.path(MAIN_OUT, "Figure4A_condition_specific_SSI_landscape.pdf"), "\n")
cat(file.path(MAIN_OUT, "Figure4B_condition_level_stability_architecture.pdf"), "\n")
cat(file.path(SUPP_OUT, "FigureS4_ECDF_SSI_by_condition.pdf"), "\n")
cat("Old Figure 4/S4 files archived in:\n")
cat(ARCHIVE_OUT, "\n")
