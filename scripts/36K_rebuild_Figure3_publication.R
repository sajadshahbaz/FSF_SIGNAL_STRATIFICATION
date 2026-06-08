#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
})

ROOT <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"
INFILE <- file.path(ROOT, "results/real_data/class_annotation/FSF_condition_class_annotation_master.tsv")
OUT <- file.path(ROOT, "results/manuscript/figures/main")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

df <- read_tsv(INFILE, show_col_types = FALSE)

df_clean <- df %>%
  distinct(condition, feature_id, stability_level, dominant_state, fsf_signal_class) %>%
  mutate(
    condition = factor(condition, levels = c("lt", "ht", "osm", "des", "uv", "gam")),
    fsf_signal_class = case_when(
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
    ),
    stability_level_clean = case_when(
      stability_level == "highly_stable" ~ "Highly stable",
      stability_level == "stable" ~ "Stable",
      stability_level == "weakly_stable_transitional" ~ "Transitional",
      stability_level == "instable" ~ "Instability",
      TRUE ~ stability_level
    )
  )

class_order <- c(
  "Highly stable up",
  "Highly stable constant",
  "Highly stable down",
  "Stable up",
  "Stable constant",
  "Stable down",
  "Transitional up",
  "Transitional constant",
  "Transitional down",
  "Instability"
)

class_cols <- c(
  "Highly stable up" = "#2166AC",
  "Highly stable constant" = "#1B7837",
  "Highly stable down" = "#E66101",
  "Stable up" = "#67A9CF",
  "Stable constant" = "#7FBF7B",
  "Stable down" = "#FDB863",
  "Transitional up" = "#8073AC",
  "Transitional constant" = "#B2ABD2",
  "Transitional down" = "#E08214",
  "Instability" = "#B2182B"
)

level_order <- c("Highly stable", "Stable", "Transitional", "Instability")

level_cols <- c(
  "Highly stable" = "#4D9221",
  "Stable" = "#A1D99B",
  "Transitional" = "#7B3294",
  "Instability" = "#B2182B"
)

theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 3, hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      strip.background = element_rect(fill = "grey92", color = "grey50"),
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25)
    )
}

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

# Figure 3A: full FSF class composition
fig3a_data <- df_clean %>%
  count(condition, fsf_signal_class, name = "n_genes") %>%
  group_by(condition) %>%
  mutate(proportion = n_genes / sum(n_genes)) %>%
  ungroup() %>%
  mutate(fsf_signal_class = factor(fsf_signal_class, levels = class_order))

p3a <- ggplot(fig3a_data, aes(x = condition, y = proportion, fill = fsf_signal_class)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = class_cols, drop = FALSE, name = "FSF class") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(
    title = "Figure 3A. Gene-weighted FSF signal-class composition",
    x = "Condition",
    y = "Proportion of genes"
  ) +
  theme_pub(12)

save_fig(p3a, "Figure3A_gene_weighted_signal_class_composition", 8.5, 5.2)

# Figure 3B: collapsed stability-level composition
fig3b_data <- df_clean %>%
  count(condition, stability_level_clean, name = "n_genes") %>%
  group_by(condition) %>%
  mutate(proportion = n_genes / sum(n_genes)) %>%
  ungroup() %>%
  mutate(stability_level_clean = factor(stability_level_clean, levels = level_order))

p3b <- ggplot(fig3b_data, aes(x = condition, y = proportion, fill = stability_level_clean)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = level_cols, drop = FALSE, name = "Stability level") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(
    title = "Figure 3B. Gene-weighted stability-level composition",
    x = "Condition",
    y = "Proportion of genes"
  ) +
  theme_pub(12)

save_fig(p3b, "Figure3B_gene_weighted_stability_level_composition", 7.8, 5.2)

# Figure 3C: condition architecture heatmap, replacing weak class-count plot
fig3c_data <- fig3a_data %>%
  complete(condition, fsf_signal_class, fill = list(n_genes = 0, proportion = 0)) %>%
  mutate(
    fsf_signal_class = factor(fsf_signal_class, levels = rev(class_order))
  )

p3c <- ggplot(fig3c_data, aes(x = condition, y = fsf_signal_class, fill = proportion)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = ifelse(proportion >= 0.02, sprintf("%.2f", proportion), "")),
    size = 3.2,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_gradient(
    low = "grey95",
    high = "#2166AC",
    name = "Gene\nproportion"
  ) +
  labs(
    title = "Figure 3C. Condition-level FSF signal architecture",
    x = "Condition",
    y = "FSF class"
  ) +
  theme_pub(12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid = element_blank()
  )

save_fig(p3c, "Figure3C_condition_level_signal_architecture_heatmap", 8.8, 5.8)

# Save source tables
write_tsv(fig3a_data, file.path(OUT, "Figure3A_source_data.tsv"))
write_tsv(fig3b_data, file.path(OUT, "Figure3B_source_data.tsv"))
write_tsv(fig3c_data, file.path(OUT, "Figure3C_source_data.tsv"))

cat("Figure 3 publication files generated:\n")
print(list.files(OUT, pattern = "^Figure3", full.names = FALSE))
