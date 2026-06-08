#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
})

ROOT <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"
OUT  <- file.path(ROOT, "results/manuscript/figures/main")
LOG  <- file.path(ROOT, "results/manuscript/logs/36E_fix_publication_figures.log")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

sink(LOG, split = TRUE)

cat("============================================================\n")
cat("Step 36E: Fix publication figures using gene-weighted data\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

bio_file <- file.path(
  ROOT,
  "results/real_data/final_manuscript_tables/Table_SignalClass_Biology_CURATED.tsv"
)

arch_file <- file.path(
  ROOT,
  "results/real_data/final_manuscript_tables/Table_Condition_SignalArchitecture_CURATED.tsv"
)

bio <- read_tsv(bio_file, show_col_types = FALSE)
arch <- read_tsv(arch_file, show_col_types = FALSE)

theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.title = element_text(face = "bold")
    )
}

save_main <- function(p, name, w = 8, h = 5) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

# ============================================================
# Figure 1: cleaner workflow, no insane overlap
# ============================================================

workflow <- tibble(
  step = factor(
    c("Input", "Perturbation", "State probabilities", "SSI", "Signal classes", "Biological profile"),
    levels = c("Input", "Perturbation", "State probabilities", "SSI", "Signal classes", "Biological profile")
  ),
  x = 1:6,
  y = 1,
  label = c(
    "Expression matrix\nor candidate genes",
    "Repeated\nperturbation",
    "P(up), P(down),\nP(const)",
    "Signal\nStratification\nIndex",
    "Stable,\ntransitional,\ninstable",
    "GO / KEGG\ninterpretation"
  )
)

p1 <- ggplot(workflow, aes(x = x, y = y)) +
  geom_segment(
    data = workflow |> filter(x < 6),
    aes(x = x + 0.28, xend = x + 0.72, y = y, yend = y),
    arrow = arrow(length = unit(0.18, "cm")),
    linewidth = 0.5
  ) +
  geom_label(aes(label = label), size = 3.3, label.size = 0.25, fill = "white") +
  scale_x_continuous(limits = c(0.6, 6.4), breaks = NULL) +
  scale_y_continuous(limits = c(0.9, 1.1), breaks = NULL) +
  labs(title = "Figure 1. FSF workflow") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0))

save_main(p1, "Figure1_Workflow_FIXED", 11, 2.8)

# ============================================================
# Figure 3A: corrected gene-weighted signal-class composition
# ============================================================

p3a <- bio %>%
  group_by(condition, fsf_signal_class) %>%
  summarise(n_genes = sum(n_genes), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(prop = n_genes / sum(n_genes)) %>%
  ungroup() %>%
  ggplot(aes(x = condition, y = prop, fill = fsf_signal_class)) +
  geom_col(width = 0.75) +
  labs(
    title = "Figure 3A. Gene-weighted FSF signal-class composition",
    x = "Condition",
    y = "Proportion of genes",
    fill = "FSF class"
  ) +
  theme_pub(12)

save_main(p3a, "Figure3A_gene_weighted_signal_class_composition_FIXED", 9, 5.5)

# ============================================================
# Figure 3B: corrected gene-weighted stability-level composition
# ============================================================

p3b <- bio %>%
  group_by(condition, stability_level) %>%
  summarise(n_genes = sum(n_genes), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(prop = n_genes / sum(n_genes)) %>%
  ungroup() %>%
  ggplot(aes(x = condition, y = prop, fill = stability_level)) +
  geom_col(width = 0.75) +
  labs(
    title = "Figure 3B. Gene-weighted stability-level composition",
    x = "Condition",
    y = "Proportion of genes",
    fill = "Stability level"
  ) +
  theme_pub(12)

save_main(p3b, "Figure3B_gene_weighted_stability_level_composition_FIXED", 8, 5.2)

# ============================================================
# Figure 3C: architecture, keep but cleaner
# ============================================================

p3c <- arch %>%
  ggplot(aes(x = condition, y = n_signal_classes, fill = signal_architecture_type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_signal_classes), vjust = -0.3, size = 4) +
  labs(
    title = "Figure 3C. Condition-level FSF architecture",
    x = "Condition",
    y = "Number of FSF classes",
    fill = "Architecture"
  ) +
  theme_pub(12)

save_main(p3c, "Figure3C_condition_architecture_FIXED", 8.5, 5)

# ============================================================
# Figure 5: corrected biological theme figures weighted by genes
# ============================================================

p5a <- bio %>%
  group_by(final_biological_theme) %>%
  summarise(n_genes = sum(n_genes), .groups = "drop") %>%
  mutate(final_biological_theme = fct_reorder(final_biological_theme, n_genes)) %>%
  ggplot(aes(x = final_biological_theme, y = n_genes)) +
  geom_col(fill = "grey35") +
  coord_flip() +
  labs(
    title = "Figure 5A. Gene-weighted curated biological themes",
    x = "Theme",
    y = "Number of genes"
  ) +
  theme_pub(12) +
  theme(axis.text.x = element_text(angle = 0))

save_main(p5a, "Figure5A_gene_weighted_theme_distribution_FIXED", 8.5, 5.5)

p5b <- bio %>%
  group_by(condition, final_biological_theme) %>%
  summarise(n_genes = sum(n_genes), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(prop = n_genes / sum(n_genes)) %>%
  ungroup() %>%
  ggplot(aes(x = final_biological_theme, y = condition, fill = prop)) +
  geom_tile(color = "white") +
  labs(
    title = "Figure 5B. Gene-weighted biological theme landscape",
    x = "Theme",
    y = "Condition",
    fill = "Gene proportion"
  ) +
  theme_pub(11)

save_main(p5b, "Figure5B_gene_weighted_theme_heatmap_FIXED", 11, 5.5)

p5c <- bio %>%
  filter(condition %in% c("des", "uv", "gam")) %>%
  ggplot(aes(x = fsf_signal_class, y = n_genes, fill = final_biological_theme)) +
  geom_col(width = 0.75) +
  facet_wrap(~condition, scales = "free_x") +
  labs(
    title = "Figure 5C. Gene-weighted DES, UV and GAM biological profiles",
    x = "FSF class",
    y = "Number of genes",
    fill = "Theme"
  ) +
  theme_pub(10)

save_main(p5c, "Figure5C_gene_weighted_DES_UV_GAM_profiles_FIXED", 12, 6)

# ============================================================
# Figure 7: simpler final model, no absurd overlap
# ============================================================

model <- tibble(
  x = c(1, 2, 3),
  y = 1,
  label = c(
    "Stable-dominated\nLT, HT, OSM\nMostly high-stability signal classes",
    "Mixed stable-instable\nDES\nStable and instable programs coexist",
    "Multi-layer architecture\nUV, GAM\nStable, transitional and instable layers coexist"
  )
)

p7 <- ggplot(model, aes(x = x, y = y)) +
  geom_segment(
    data = model |> filter(x < 3),
    aes(x = x + 0.25, xend = x + 0.75, y = y, yend = y),
    arrow = arrow(length = unit(0.18, "cm")),
    linewidth = 0.5
  ) +
  geom_label(aes(label = label), size = 3.4, label.size = 0.25, fill = "white") +
  scale_x_continuous(limits = c(0.6, 3.4), breaks = NULL) +
  scale_y_continuous(limits = c(0.9, 1.1), breaks = NULL) +
  labs(title = "Figure 7. Final FSF interpretation model") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0))

save_main(p7, "Figure7_Final_Model_FIXED", 10, 3)

cat("\nFixed figures generated:\n")
print(list.files(OUT, pattern = "FIXED", full.names = FALSE))

cat("\nFinished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()
