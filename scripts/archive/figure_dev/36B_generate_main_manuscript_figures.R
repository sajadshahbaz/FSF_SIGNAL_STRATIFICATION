#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(forcats)
})

ROOT <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"
FIG_MAIN <- file.path(ROOT, "results/manuscript/figures/main")
FIG_SUPP <- file.path(ROOT, "results/manuscript/figures/supplementary")
LOG <- file.path(ROOT, "results/manuscript/logs/36B_generate_main_manuscript_figures.log")

dir.create(FIG_MAIN, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_SUPP, recursive = TRUE, showWarnings = FALSE)

sink(LOG, split = TRUE)

cat("============================================================\n")
cat("Step 36B: Generate Main Manuscript Figures\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

bio <- read_tsv(
  file.path(ROOT, "results/manuscript/tables/main/Table_SignalClass_Biology_CURATED.tsv"),
  show_col_types = FALSE
)

arch <- read_tsv(
  file.path(ROOT, "results/manuscript/tables/main/Table_Condition_SignalArchitecture_CURATED.tsv"),
  show_col_types = FALSE
)

bench <- read_tsv(
  file.path(ROOT, "results/manuscript/tables/supplementary/synthetic_or_benchmark__FSF_v1_final_benchmark_report.tsv"),
  show_col_types = FALSE
)

# -------------------------
# Figure 3: Real-data architecture
# -------------------------

p1 <- bio %>%
  count(condition, fsf_signal_class) %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = condition, y = prop, fill = fsf_signal_class)) +
  geom_col(width = 0.75) +
  labs(
    title = "Figure 3A. FSF signal-class composition across real-data conditions",
    x = "Condition",
    y = "Proportion of signal classes",
    fill = "FSF class"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(FIG_MAIN, "Figure3A_signal_class_composition.pdf"), p1, width = 8, height = 5)
ggsave(file.path(FIG_MAIN, "Figure3A_signal_class_composition.png"), p1, width = 8, height = 5, dpi = 300)

p2 <- bio %>%
  count(condition, stability_level) %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = condition, y = prop, fill = stability_level)) +
  geom_col(width = 0.75) +
  labs(
    title = "Figure 3B. Stability-level composition across real-data conditions",
    x = "Condition",
    y = "Proportion",
    fill = "Stability level"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(FIG_MAIN, "Figure3B_stability_level_composition.pdf"), p2, width = 8, height = 5)
ggsave(file.path(FIG_MAIN, "Figure3B_stability_level_composition.png"), p2, width = 8, height = 5, dpi = 300)

p3 <- arch %>%
  ggplot(aes(x = condition, y = n_signal_classes, fill = signal_architecture_type)) +
  geom_col(width = 0.7) +
  labs(
    title = "Figure 3C. Condition-level signal architecture",
    x = "Condition",
    y = "Number of signal classes",
    fill = "Architecture type"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(FIG_MAIN, "Figure3C_condition_architecture_type.pdf"), p3, width = 8, height = 5)
ggsave(file.path(FIG_MAIN, "Figure3C_condition_architecture_type.png"), p3, width = 8, height = 5, dpi = 300)

# -------------------------
# Figure 4: Cross-condition landscape
# -------------------------

heat_data <- bio %>%
  count(condition, final_biological_theme) %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p4 <- ggplot(heat_data, aes(x = final_biological_theme, y = condition, fill = prop)) +
  geom_tile(color = "white") +
  labs(
    title = "Figure 4A. Cross-condition biological theme landscape",
    x = "Curated biological theme",
    y = "Condition",
    fill = "Class proportion"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(FIG_MAIN, "Figure4A_cross_condition_theme_heatmap.pdf"), p4, width = 11, height = 5.5)
ggsave(file.path(FIG_MAIN, "Figure4A_cross_condition_theme_heatmap.png"), p4, width = 11, height = 5.5, dpi = 300)

p5 <- bio %>%
  mutate(total_support = significant_GO_terms_FDR005 + significant_KEGG_terms_FDR005) %>%
  ggplot(aes(x = condition, y = total_support, fill = fsf_signal_class)) +
  geom_col(position = "stack") +
  labs(
    title = "Figure 4B. Functional support across signal classes",
    x = "Condition",
    y = "Significant GO + KEGG terms",
    fill = "FSF class"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(FIG_MAIN, "Figure4B_functional_support_by_class.pdf"), p5, width = 9, height = 5)
ggsave(file.path(FIG_MAIN, "Figure4B_functional_support_by_class.png"), p5, width = 9, height = 5, dpi = 300)

# -------------------------
# Figure 5: Biological interpretation
# -------------------------

p6 <- bio %>%
  count(final_biological_theme, sort = TRUE) %>%
  mutate(final_biological_theme = fct_reorder(final_biological_theme, n)) %>%
  ggplot(aes(x = final_biological_theme, y = n)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Figure 5A. Distribution of curated biological themes",
    x = "Curated biological theme",
    y = "Number of FSF classes"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(FIG_MAIN, "Figure5A_curated_theme_distribution.pdf"), p6, width = 8, height = 5.5)
ggsave(file.path(FIG_MAIN, "Figure5A_curated_theme_distribution.png"), p6, width = 8, height = 5.5, dpi = 300)

p7 <- bio %>%
  filter(condition %in% c("des", "uv", "gam")) %>%
  ggplot(aes(x = fsf_signal_class, y = n_genes, fill = final_biological_theme)) +
  geom_col(width = 0.75) +
  facet_wrap(~condition, scales = "free_x") +
  labs(
    title = "Figure 5B. Biological profiles of DES, UV and GAM signal classes",
    x = "FSF signal class",
    y = "Number of genes",
    fill = "Theme"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(FIG_MAIN, "Figure5B_DES_UV_GAM_biological_profiles.pdf"), p7, width = 12, height = 6)
ggsave(file.path(FIG_MAIN, "Figure5B_DES_UV_GAM_biological_profiles.png"), p7, width = 12, height = 6, dpi = 300)

# -------------------------
# Figure 6: Synthetic benchmark
# -------------------------

if (nrow(bench) > 0) {
  numeric_cols <- names(bench)[sapply(bench, is.numeric)]
  
  if (length(numeric_cols) >= 2) {
    bench_long <- bench %>%
      select(any_of(numeric_cols)) %>%
      pivot_longer(everything(), names_to = "metric", values_to = "value")
    
    p8 <- ggplot(bench_long, aes(x = metric, y = value)) +
      geom_boxplot() +
      coord_flip() +
      labs(
        title = "Figure 6A. Synthetic benchmark metric distribution",
        x = "Benchmark metric",
        y = "Value"
      ) +
      theme_bw(base_size = 12)
    
    ggsave(file.path(FIG_MAIN, "Figure6A_synthetic_benchmark_metrics.pdf"), p8, width = 8, height = 5)
    ggsave(file.path(FIG_MAIN, "Figure6A_synthetic_benchmark_metrics.png"), p8, width = 8, height = 5, dpi = 300)
  }
}

# -------------------------
# Supplementary figures
# -------------------------

p9 <- bio %>%
  ggplot(aes(x = annotated_fraction, y = n_genes, color = condition)) +
  geom_point(size = 3) +
  labs(
    title = "Supplementary Figure. Annotation fraction by FSF class size",
    x = "Annotated fraction",
    y = "Number of genes",
    color = "Condition"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(FIG_SUPP, "Supplementary_annotation_fraction_vs_class_size.pdf"), p9, width = 7, height = 5)
ggsave(file.path(FIG_SUPP, "Supplementary_annotation_fraction_vs_class_size.png"), p9, width = 7, height = 5, dpi = 300)

cat("Generated main figures:\n")
print(list.files(FIG_MAIN, pattern = "\\.(pdf|png)$", full.names = FALSE))

cat("\nGenerated supplementary figures:\n")
print(list.files(FIG_SUPP, pattern = "\\.(pdf|png)$", full.names = FALSE))

cat("\nFinished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()
