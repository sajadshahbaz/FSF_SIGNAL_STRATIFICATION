#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(grid)
})

ROOT <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"
FIG_MAIN <- file.path(ROOT, "results/archive/pre_repair_manuscript/figures/main")
FIG_SUPP <- file.path(ROOT, "results/archive/pre_repair_manuscript/figures/supplementary")
LOG <- file.path(ROOT, "results/archive/pre_repair_manuscript/logs/36C_rebuild_publication_figures.log")

dir.create(FIG_MAIN, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_SUPP, recursive = TRUE, showWarnings = FALSE)

sink(LOG, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 36C: Rebuild Publication-Ready Figures\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold")
    )
}

save_plot <- function(p, name, width = 8, height = 5) {
  ggsave(file.path(FIG_MAIN, paste0(name, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(FIG_MAIN, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
}

save_supp <- function(p, name, width = 8, height = 5) {
  ggsave(file.path(FIG_SUPP, paste0(name, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(FIG_SUPP, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
}

read_if_exists <- function(path) {
  if (!file.exists(path)) {
    cat("[MISSING]", path, "\n")
    return(NULL)
  }
  read_tsv(path, show_col_types = FALSE)
}

# ==============================
# INPUT TABLES
# ==============================

bio <- read_if_exists(file.path(ROOT, "results/real_data/final_manuscript_tables/Table_SignalClass_Biology_CURATED.tsv"))
arch <- read_if_exists(file.path(ROOT, "results/real_data/final_manuscript_tables/Table_Condition_SignalArchitecture_CURATED.tsv"))

fsf_metrics <- read_if_exists(file.path(ROOT, "results/fsf_metrics/fsf_ssi_stability_deviation.tsv"))
fsf_summary <- read_if_exists(file.path(ROOT, "results/fsf_metrics/fsf_signal_class_summary.tsv"))
signal_states <- read_if_exists(file.path(ROOT, "results/signal_states/feature_signal_states.tsv"))

truth <- read_if_exists(file.path(ROOT, "results/synthetic/synthetic_truth_table.tsv"))
conf <- read_if_exists(file.path(ROOT, "results/synthetic/synthetic_recovery_confusion_matrix.tsv"))
rec <- read_if_exists(file.path(ROOT, "results/synthetic/synthetic_recovery_summary.tsv"))

noise_summary <- read_if_exists(file.path(ROOT, "results/synthetic/noise_gradient/publication_noise_gradient_summary.tsv"))
noise_recovery <- read_if_exists(file.path(ROOT, "results/synthetic/noise_gradient/publication_noise_gradient_dominant_recovery.tsv"))

tau_files <- list.files(file.path(ROOT, "results/synthetic/tau_sensitivity"), pattern = "\\.tsv$", full.names = TRUE)
tau_data <- if (length(tau_files) > 0) {
  bind_rows(lapply(tau_files, function(x) {
    read_tsv(x, show_col_types = FALSE) %>% mutate(source_file = basename(x))
  }))
} else NULL

# ==============================
# FIGURE 1: WORKFLOW
# ==============================

workflow <- tibble(
  step = factor(
    c("Input matrix", "Perturbation", "State probabilities", "SSI", "Signal classes", "Biological interpretation"),
    levels = c("Input matrix", "Perturbation", "State probabilities", "SSI", "Signal classes", "Biological interpretation")
  ),
  x = 1:6,
  y = 1,
  label = c(
    "Expression / candidate features",
    "Repeated perturbation",
    "P(up), P(down), P(const)",
    "Signal Stratification Index",
    "Stable / transitional / instable",
    "GO / KEGG / theme profile"
  )
)

p1 <- ggplot(workflow, aes(x = x, y = y)) +
  geom_segment(aes(x = x, xend = x + 0.75, y = y, yend = y), data = workflow %>% filter(x < 6),
               arrow = arrow(length = unit(0.25, "cm")), linewidth = 0.5) +
  geom_label(aes(label = label), size = 3.5, label.size = 0.25, fill = "white") +
  scale_x_continuous(limits = c(0.5, 6.5), breaks = workflow$x, labels = workflow$step) +
  scale_y_continuous(limits = c(0.8, 1.2)) +
  labs(title = "Figure 1. FSF workflow", x = NULL, y = NULL) +
  theme_void(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0)
  )

save_plot(p1, "Figure1_Workflow", 11, 3)

# ==============================
# FIGURE 2: THEORY / PROBABILITY SIMPLEX
# ==============================

simplex <- expand.grid(
  p_up = seq(0, 1, by = 0.02),
  p_down = seq(0, 1, by = 0.02)
) %>%
  mutate(
    p_const = 1 - p_up - p_down
  ) %>%
  filter(p_const >= 0) %>%
  mutate(
    max_p = pmax(p_up, p_down, p_const),
    state = case_when(
      max_p == p_up ~ "Stable-up region",
      max_p == p_down ~ "Stable-down region",
      max_p == p_const ~ "Stable-constant region",
      TRUE ~ "Mixed"
    ),
    x = p_down + 0.5 * p_const,
    y = p_const * sqrt(3) / 2
  )

p2 <- ggplot(simplex, aes(x = x, y = y, fill = state)) +
  geom_tile(width = 0.015, height = 0.015, alpha = 0.85) +
  annotate("text", x = 0, y = 0, label = "P(up)", hjust = 1.1, size = 4) +
  annotate("text", x = 1, y = 0, label = "P(down)", hjust = -0.1, size = 4) +
  annotate("text", x = 0.5, y = sqrt(3)/2, label = "P(const)", vjust = -0.8, size = 4) +
  labs(title = "Figure 2. FSF probability space", x = NULL, y = NULL, fill = "Dominant region") +
  coord_equal() +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "right"
  )

save_plot(p2, "Figure2_Probability_Simplex", 7, 6)

# ==============================
# FIGURE 3: REAL-DATA ARCHITECTURE
# ==============================

if (!is.null(bio)) {
  p3a <- bio %>%
    count(condition, fsf_signal_class) %>%
    group_by(condition) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = condition, y = prop, fill = fsf_signal_class)) +
    geom_col(width = 0.75) +
    labs(
      title = "Figure 3A. Signal-class composition",
      x = "Condition", y = "Proportion", fill = "FSF class"
    ) +
    theme_pub(11)

  save_plot(p3a, "Figure3A_signal_class_composition", 8, 5)

  p3b <- bio %>%
    count(condition, stability_level) %>%
    group_by(condition) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = condition, y = prop, fill = stability_level)) +
    geom_col(width = 0.75) +
    labs(
      title = "Figure 3B. Stability-level composition",
      x = "Condition", y = "Proportion", fill = "Stability level"
    ) +
    theme_pub(11)

  save_plot(p3b, "Figure3B_stability_level_composition", 8, 5)
}

if (!is.null(arch)) {
  p3c <- arch %>%
    ggplot(aes(x = condition, y = n_signal_classes, fill = signal_architecture_type)) +
    geom_col(width = 0.7) +
    labs(
      title = "Figure 3C. Condition-level architecture",
      x = "Condition", y = "Number of signal classes", fill = "Architecture"
    ) +
    theme_pub(11)

  save_plot(p3c, "Figure3C_condition_architecture", 8, 5)
}

# ==============================
# FIGURE 4: SSI / DEVIATION / STATE SPACE
# ==============================

if (!is.null(fsf_metrics)) {

  names(fsf_metrics) <- make.names(names(fsf_metrics))

  ssi_col <- names(fsf_metrics)[str_detect(names(fsf_metrics), regex("^SSI$|ssi", ignore_case = TRUE))][1]
  dev_col <- names(fsf_metrics)[str_detect(names(fsf_metrics), regex("deviation", ignore_case = TRUE))][1]
  class_col <- names(fsf_metrics)[str_detect(names(fsf_metrics), regex("class", ignore_case = TRUE))][1]

  if (!is.na(ssi_col)) {
    p4a <- ggplot(fsf_metrics, aes(x = .data[[ssi_col]])) +
      geom_histogram(bins = 50, fill = "grey35", color = "white") +
      labs(title = "Figure 4A. SSI distribution", x = "SSI", y = "Number of features") +
      theme_pub(11)

    save_plot(p4a, "Figure4A_SSI_distribution", 7, 5)
  }

  if (!is.na(dev_col)) {
    p4b <- ggplot(fsf_metrics, aes(x = .data[[dev_col]])) +
      geom_histogram(bins = 50, fill = "grey35", color = "white") +
      labs(title = "Figure 4B. Stability-deviation distribution", x = "Stability deviation", y = "Number of features") +
      theme_pub(11)

    save_plot(p4b, "Figure4B_stability_deviation_distribution", 7, 5)
  }

  if (!is.na(ssi_col) && !is.na(dev_col) && !is.na(class_col)) {
    p4c <- ggplot(fsf_metrics, aes(x = .data[[ssi_col]], y = .data[[dev_col]], color = .data[[class_col]])) +
      geom_point(alpha = 0.35, size = 0.8) +
      labs(title = "Figure 4C. SSI versus stability deviation", x = "SSI", y = "Stability deviation", color = "Class") +
      theme_pub(10)

    save_plot(p4c, "Figure4C_SSI_vs_deviation", 8, 5.5)
  }
}

if (!is.null(signal_states)) {
  names(signal_states) <- make.names(names(signal_states))
  up_col <- names(signal_states)[str_detect(names(signal_states), regex("p.*up|up.*prob", ignore_case = TRUE))][1]
  down_col <- names(signal_states)[str_detect(names(signal_states), regex("p.*down|down.*prob", ignore_case = TRUE))][1]
  const_col <- names(signal_states)[str_detect(names(signal_states), regex("p.*const|const.*prob", ignore_case = TRUE))][1]
  class_col <- names(signal_states)[str_detect(names(signal_states), regex("class", ignore_case = TRUE))][1]

  if (!is.na(up_col) && !is.na(down_col)) {
    plot_data <- signal_states
    if (nrow(plot_data) > 20000) set.seed(1); plot_data <- plot_data %>% slice_sample(n = min(20000, n()))

    p4d <- ggplot(plot_data, aes(x = .data[[up_col]], y = .data[[down_col]], color = .data[[class_col]])) +
      geom_point(alpha = 0.35, size = 0.7) +
      labs(title = "Figure 4D. State-probability landscape", x = "P(up)", y = "P(down)", color = "FSF class") +
      theme_pub(10)

    save_plot(p4d, "Figure4D_state_probability_landscape", 8, 6)
  }
}

# ==============================
# FIGURE 5: BIOLOGICAL THEMES
# ==============================

if (!is.null(bio)) {
  p5a <- bio %>%
    count(final_biological_theme, sort = TRUE) %>%
    mutate(final_biological_theme = fct_reorder(final_biological_theme, n)) %>%
    ggplot(aes(x = final_biological_theme, y = n)) +
    geom_col(fill = "grey35") +
    coord_flip() +
    labs(title = "Figure 5A. Curated biological themes", x = "Theme", y = "Number of signal classes") +
    theme_pub(11) +
    theme(axis.text.x = element_text(angle = 0))

  save_plot(p5a, "Figure5A_curated_theme_distribution", 8, 5.5)

  p5b <- bio %>%
    count(condition, final_biological_theme) %>%
    group_by(condition) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = final_biological_theme, y = condition, fill = prop)) +
    geom_tile(color = "white") +
    labs(title = "Figure 5B. Biological theme landscape", x = "Theme", y = "Condition", fill = "Proportion") +
    theme_pub(10)

  save_plot(p5b, "Figure5B_theme_heatmap", 11, 5.5)

  p5c <- bio %>%
    filter(condition %in% c("des", "uv", "gam")) %>%
    ggplot(aes(x = fsf_signal_class, y = n_genes, fill = final_biological_theme)) +
    geom_col(width = 0.75) +
    facet_wrap(~condition, scales = "free_x") +
    labs(title = "Figure 5C. DES, UV and GAM biological profiles", x = "FSF class", y = "Number of genes", fill = "Theme") +
    theme_pub(10)

  save_plot(p5c, "Figure5C_DES_UV_GAM_profiles", 12, 6)
}

# ==============================
# FIGURE 6: SYNTHETIC BENCHMARKS
# ==============================

if (!is.null(conf)) {
  names(conf) <- make.names(names(conf))
  if (ncol(conf) >= 3) {
    p6a <- conf %>%
      ggplot(aes(x = .data[[names(conf)[1]]], y = .data[[names(conf)[2]]], fill = .data[[names(conf)[3]]])) +
      geom_tile(color = "white") +
      geom_text(aes(label = .data[[names(conf)[3]]]), size = 3) +
      labs(title = "Figure 6A. Synthetic truth versus FSF class", x = "Predicted FSF class", y = "Synthetic truth", fill = "Count") +
      theme_pub(10)

    save_plot(p6a, "Figure6A_synthetic_truth_confusion", 8, 6)
  }
}

if (!is.null(noise_summary)) {
  names(noise_summary) <- make.names(names(noise_summary))
  noise_col <- names(noise_summary)[str_detect(names(noise_summary), regex("noise", ignore_case = TRUE))][1]
  numeric_cols <- names(noise_summary)[sapply(noise_summary, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, noise_col)

  if (!is.na(noise_col) && length(numeric_cols) > 0) {
    p6b <- noise_summary %>%
      select(all_of(c(noise_col, numeric_cols))) %>%
      pivot_longer(cols = all_of(numeric_cols), names_to = "metric", values_to = "value") %>%
      ggplot(aes(x = .data[[noise_col]], y = value, linetype = metric)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      labs(title = "Figure 6B. Noise-gradient benchmark", x = "Noise level", y = "Metric value", linetype = "Metric") +
      theme_pub(10)

    save_plot(p6b, "Figure6B_noise_gradient_metrics", 8, 5)
  }
}

if (!is.null(noise_recovery)) {
  names(noise_recovery) <- make.names(names(noise_recovery))
  noise_col <- names(noise_recovery)[str_detect(names(noise_recovery), regex("noise", ignore_case = TRUE))][1]
  class_col <- names(noise_recovery)[str_detect(names(noise_recovery), regex("class|truth|state", ignore_case = TRUE))][1]
  value_col <- names(noise_recovery)[str_detect(names(noise_recovery), regex("recover|accuracy|rate|proportion|mean", ignore_case = TRUE))][1]

  if (!is.na(noise_col) && !is.na(value_col)) {
    p6c <- ggplot(noise_recovery, aes(x = .data[[noise_col]], y = .data[[value_col]], linetype = .data[[class_col]])) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      labs(title = "Figure 6C. Class recovery across noise", x = "Noise level", y = "Recovery", linetype = "Class") +
      theme_pub(10)

    save_plot(p6c, "Figure6C_noise_class_recovery", 8, 5)
  }
}

if (!is.null(tau_data)) {
  names(tau_data) <- make.names(names(tau_data))
  tau_col <- names(tau_data)[str_detect(names(tau_data), regex("tau|threshold", ignore_case = TRUE))][1]
  class_col <- names(tau_data)[str_detect(names(tau_data), regex("class", ignore_case = TRUE))][1]
  value_col <- names(tau_data)[str_detect(names(tau_data), regex("count|n$|mean|ssi", ignore_case = TRUE))][1]

  if (!is.na(tau_col) && !is.na(value_col)) {
    p6d <- ggplot(tau_data, aes(x = .data[[tau_col]], y = .data[[value_col]], linetype = .data[[class_col]])) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      labs(title = "Figure 6D. Tau sensitivity", x = "Tau / threshold", y = "Value", linetype = "Class") +
      theme_pub(10)

    save_plot(p6d, "Figure6D_tau_sensitivity", 8, 5)
  }
}

# ==============================
# FIGURE 7: FINAL MODEL
# ==============================

model <- tibble(
  architecture = c("Stable-dominated", "Mixed stable-instable", "Multi-layer"),
  conditions = c("LT, HT, OSM", "DES", "UV, GAM"),
  interpretation = c(
    "Most classes converge toward stable signal behavior",
    "Stable and instable classes coexist",
    "Multiple signal regimes coexist across stability levels"
  ),
  x = c(1, 2, 3),
  y = 1
)

p7 <- ggplot(model, aes(x = x, y = y)) +
  geom_label(aes(label = paste0(architecture, "\n", conditions, "\n", interpretation)),
             size = 3.4, label.size = 0.25, fill = "white") +
  geom_segment(aes(x = x, xend = x + 0.75, y = y, yend = y), data = model %>% filter(x < 3),
               arrow = arrow(length = unit(0.25, "cm")), linewidth = 0.5) +
  scale_x_continuous(limits = c(0.5, 3.5), breaks = NULL) +
  scale_y_continuous(limits = c(0.8, 1.2), breaks = NULL) +
  labs(title = "Figure 7. Final FSF interpretation model", x = NULL, y = NULL) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0))

save_plot(p7, "Figure7_Final_Model", 11, 3.5)

# ==============================
# SUPPLEMENTARY
# ==============================

if (!is.null(bio)) {
  pS1 <- bio %>%
    ggplot(aes(x = annotated_fraction, y = n_genes, color = condition)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = "Supplementary Figure S1. Annotation coverage by class size",
         x = "Annotated fraction", y = "Number of genes", color = "Condition") +
    theme_pub(11)

  save_supp(pS1, "Supplementary_FigureS1_annotation_coverage", 7, 5)
}

cat("\nGenerated main figures:\n")
print(list.files(FIG_MAIN, pattern = "\\.(pdf|png)$", full.names = FALSE))

cat("\nGenerated supplementary figures:\n")
print(list.files(FIG_SUPP, pattern = "\\.(pdf|png)$", full.names = FALSE))

cat("\nFinished:", as.character(Sys.time()), "\n")
cat("============================================================\n")
sink()
