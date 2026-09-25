#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
})

ROOT <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"
MAIN <- file.path(ROOT, "results/archive/pre_repair_manuscript/figures/main")
SUPP <- file.path(ROOT, "results/archive/pre_repair_manuscript/figures/supplementary")
LOG  <- file.path(ROOT, "results/archive/pre_repair_manuscript/logs/36F_cleanup_and_fix_figures.log")

dir.create(MAIN, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP, recursive = TRUE, showWarnings = FALSE)

sink(LOG, split = TRUE)

cat("============================================================\n")
cat("Step 36F: Clean duplicated/corrupted figures and repair Figure 1/2/6\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Remove corrupted / obsolete duplicate figures
# ------------------------------------------------------------

bad_patterns <- c(
  "Figure1_Workflow\\.",
  "Figure2_Probability_Simplex\\.",
  "Figure3A_signal_class_composition\\.",
  "Figure3B_stability_level_composition\\.",
  "Figure3C_condition_architecture\\.",
  "Figure5A_curated_theme_distribution\\.",
  "Figure5B_theme_heatmap\\.",
  "Figure5C_DES_UV_GAM_profiles\\.",
  "Figure7_Final_Model\\."
)

all_main <- list.files(MAIN, full.names = TRUE)

for (pat in bad_patterns) {
  f <- all_main[grepl(pat, basename(all_main))]
  if (length(f) > 0) {
    file.remove(f)
    cat("Removed obsolete:", basename(f), "\n")
  }
}

# Rename FIXED versions to canonical names
rename_pair <- function(old_base, new_base) {
  for (ext in c("pdf", "png")) {
    old <- file.path(MAIN, paste0(old_base, ".", ext))
    new <- file.path(MAIN, paste0(new_base, ".", ext))
    if (file.exists(old)) {
      file.copy(old, new, overwrite = TRUE)
      file.remove(old)
      cat("Renamed:", basename(old), "->", basename(new), "\n")
    }
  }
}

rename_pair("Figure3A_gene_weighted_signal_class_composition_FIXED", "Figure3A_gene_weighted_signal_class_composition")
rename_pair("Figure3B_gene_weighted_stability_level_composition_FIXED", "Figure3B_gene_weighted_stability_level_composition")
rename_pair("Figure3C_condition_architecture_FIXED", "Figure3C_condition_architecture")
rename_pair("Figure5A_gene_weighted_theme_distribution_FIXED", "Figure5A_gene_weighted_theme_distribution")
rename_pair("Figure5B_gene_weighted_theme_heatmap_FIXED", "Figure5B_gene_weighted_theme_heatmap")
rename_pair("Figure5C_gene_weighted_DES_UV_GAM_profiles_FIXED", "Figure5C_gene_weighted_DES_UV_GAM_profiles")
rename_pair("Figure7_Final_Model_FIXED", "Figure7_Final_Model")

# ------------------------------------------------------------
# Plot theme
# ------------------------------------------------------------

theme_clean <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.title = element_text(face = "bold")
    )
}

save_main <- function(p, name, w, h) {
  ggsave(file.path(MAIN, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(MAIN, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  cat("Saved:", name, "\n")
}

# ------------------------------------------------------------
# 2. Rebuild Figure 1 properly
# ------------------------------------------------------------

workflow <- tibble(
  x = 1:6,
  y = 1,
  label = c(
    "Input matrix\nor candidate set",
    "Repeated\nperturbation",
    "State probabilities\nP(up), P(down), P(const)",
    "Signal Stratification\nIndex (SSI)",
    "Signal class\nassignment",
    "Biological\ninterpretation"
  )
)

p1 <- ggplot(workflow, aes(x = x, y = y)) +
  geom_segment(
    data = workflow |> filter(x < 6),
    aes(x = x + 0.22, xend = x + 0.78, y = y, yend = y),
    arrow = arrow(length = unit(0.18, "cm")),
    linewidth = 0.45
  ) +
  geom_label(
    aes(label = label),
    size = 3.25,
    label.size = 0.25,
    fill = "white",
    label.padding = unit(0.25, "lines")
  ) +
  scale_x_continuous(limits = c(0.6, 6.4), breaks = NULL) +
  scale_y_continuous(limits = c(0.92, 1.08), breaks = NULL) +
  labs(title = "Figure 1. FSF workflow") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0))

save_main(p1, "Figure1_Workflow", 12, 2.6)

# ------------------------------------------------------------
# 3. Rebuild Figure 2 with labels inside boundaries
# ------------------------------------------------------------

simplex <- expand.grid(
  p_up = seq(0, 1, by = 0.025),
  p_down = seq(0, 1, by = 0.025)
) |>
  mutate(p_const = 1 - p_up - p_down) |>
  filter(p_const >= -1e-9) |>
  mutate(
    max_p = pmax(p_up, p_down, p_const),
    region = case_when(
      p_const == max_p ~ "Stable-constant region",
      p_down  == max_p ~ "Stable-down region",
      p_up    == max_p ~ "Stable-up region",
      TRUE ~ "Boundary"
    ),
    x = p_down + 0.5 * p_const,
    y = p_const * sqrt(3) / 2
  )

p2 <- ggplot(simplex, aes(x = x, y = y, fill = region)) +
  geom_point(shape = 22, size = 2.2, color = "white", stroke = 0.15) +
  annotate("segment", x = 0, xend = 1, y = 0, yend = 0, linewidth = 0.4) +
  annotate("segment", x = 0, xend = 0.5, y = 0, yend = sqrt(3)/2, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 0.5, y = 0, yend = sqrt(3)/2, linewidth = 0.4) +
  annotate("text", x = 0.08, y = 0.035, label = "P(up)", size = 4, hjust = 0) +
  annotate("text", x = 0.92, y = 0.035, label = "P(down)", size = 4, hjust = 1) +
  annotate("text", x = 0.5, y = 0.82, label = "P(const)", size = 4, vjust = 1) +
  scale_x_continuous(limits = c(-0.05, 1.05), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(-0.04, 0.92), expand = expansion(mult = 0)) +
  coord_equal() +
  labs(
    title = "Figure 2. FSF probability space",
    x = NULL,
    y = NULL,
    fill = "Dominant region"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "right"
  )

save_main(p2, "Figure2_Probability_Simplex", 7.5, 6)

# ------------------------------------------------------------
# 4. Try to rebuild missing Figure 6B and 6C
# ------------------------------------------------------------

noise_summary_file <- file.path(ROOT, "results/synthetic/noise_gradient/publication_noise_gradient_summary.tsv")
noise_recovery_file <- file.path(ROOT, "results/synthetic/noise_gradient/publication_noise_gradient_dominant_recovery.tsv")

cat("\nChecking synthetic noise files:\n")
cat("noise summary exists:", file.exists(noise_summary_file), "\n")
cat("noise recovery exists:", file.exists(noise_recovery_file), "\n")

if (file.exists(noise_summary_file)) {
  x <- read_tsv(noise_summary_file, show_col_types = FALSE)
  cat("Noise summary columns:\n")
  print(names(x))

  names(x) <- make.names(names(x))
  noise_col <- names(x)[grepl("noise", names(x), ignore.case = TRUE)][1]
  numeric_cols <- names(x)[sapply(x, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, noise_col)

  if (!is.na(noise_col) && length(numeric_cols) > 0) {
    long <- x |>
      select(all_of(c(noise_col, numeric_cols))) |>
      pivot_longer(cols = all_of(numeric_cols), names_to = "metric", values_to = "value")

    p6b <- ggplot(long, aes(x = .data[[noise_col]], y = value, linetype = metric)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      labs(
        title = "Figure 6B. Noise-gradient robustness",
        x = "Noise level",
        y = "Metric value",
        linetype = "Metric"
      ) +
      theme_clean(12)

    save_main(p6b, "Figure6B_noise_gradient_robustness", 8.5, 5)
  } else {
    cat("Cannot build Figure 6B: no usable noise/numeric columns.\n")
  }
}

if (file.exists(noise_recovery_file)) {
  x <- read_tsv(noise_recovery_file, show_col_types = FALSE)
  cat("Noise recovery columns:\n")
  print(names(x))

  names(x) <- make.names(names(x))
  noise_col <- names(x)[grepl("noise", names(x), ignore.case = TRUE)][1]
  class_col <- names(x)[grepl("class|truth|state", names(x), ignore.case = TRUE)][1]
  value_col <- names(x)[grepl("recover|accuracy|rate|proportion|mean", names(x), ignore.case = TRUE)][1]

  if (!is.na(noise_col) && !is.na(class_col) && !is.na(value_col)) {
    p6c <- ggplot(x, aes(x = .data[[noise_col]], y = .data[[value_col]], linetype = .data[[class_col]])) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      labs(
        title = "Figure 6C. Class-specific recovery across noise",
        x = "Noise level",
        y = "Recovery",
        linetype = "Signal class"
      ) +
      theme_clean(12)

    save_main(p6c, "Figure6C_class_specific_noise_recovery", 8.5, 5)
  } else {
    cat("Cannot build Figure 6C: required columns not detected.\n")
  }
}

# ------------------------------------------------------------
# 5. Final inventory
# ------------------------------------------------------------

cat("\nFinal main figure files:\n")
print(sort(list.files(MAIN, pattern = "\\.(pdf|png)$", full.names = FALSE)))

cat("\nFinished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()
