#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 20: Plot Condition-wise Real FSF Results
#
# Inputs:
#   results/real_data/condition_fsf/condition_wise_fsf_metrics.tsv
#   results/real_data/condition_fsf/condition_wise_fsf_class_summary.tsv
#
# Outputs:
#   results/real_data/figures/condition_fsf_stacked_barplot.png
#   results/real_data/figures/condition_instable_proportion.png
#   results/real_data/figures/condition_ssi_distribution.png
#   results/real_data/figures/condition_stability_deviation_distribution.png
#   results/logs/20_plot_condition_wise_fsf.log
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

METRICS_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/condition_fsf/condition_wise_fsf_metrics.tsv"
)

SUMMARY_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/condition_fsf/condition_wise_fsf_class_summary.tsv"
)

FIG_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/figures"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/20_plot_condition_wise_fsf.log"
)

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 20: Plot Condition-wise Real FSF Results\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

metrics <- read_tsv(METRICS_FILE, show_col_types = FALSE)
summary_tbl <- read_tsv(SUMMARY_FILE, show_col_types = FALSE)

condition_order <- c("osm", "lt", "ht", "gam", "uv", "des")

class_order <- c(
  "highly_stable_constant",
  "stable_constant",
  "weakly_stable_transitional_constant",
  "highly_stable_up",
  "stable_up",
  "weakly_stable_transitional_up",
  "highly_stable_down",
  "stable_down",
  "weakly_stable_transitional_down",
  "instable"
)

metrics <- metrics %>%
  mutate(
    condition = factor(condition, levels = condition_order),
    fsf_signal_class = factor(fsf_signal_class, levels = class_order),
    stability_level = factor(
      stability_level,
      levels = c(
        "highly_stable",
        "stable",
        "weakly_stable_transitional",
        "instable"
      )
    )
  )

summary_tbl <- summary_tbl %>%
  mutate(
    condition = factor(condition, levels = condition_order),
    fsf_signal_class = factor(fsf_signal_class, levels = class_order)
  )

cat("Metrics rows:", nrow(metrics), "\n")
cat("Summary rows:", nrow(summary_tbl), "\n\n")

# ------------------------------------------------------------
# Plot 1: Stacked class composition
# ------------------------------------------------------------

p1 <- ggplot(
  summary_tbl,
  aes(
    x = condition,
    y = proportion,
    fill = fsf_signal_class
  )
) +
  geom_col() +
  labs(
    title = "Condition-wise FSF signal class composition",
    x = "Condition",
    y = "Proportion of features",
    fill = "FSF signal class"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(FIG_DIR, "condition_fsf_stacked_barplot.png"),
  p1,
  width = 10,
  height = 6,
  dpi = 300
)

cat("Saved condition_fsf_stacked_barplot.png\n")

# ------------------------------------------------------------
# Plot 2: Instable proportion by condition
# ------------------------------------------------------------

instable_tbl <- summary_tbl %>%
  filter(fsf_signal_class == "instable") %>%
  select(condition, proportion, n_features)

p2 <- ggplot(
  instable_tbl,
  aes(
    x = condition,
    y = proportion
  )
) +
  geom_col() +
  geom_text(
    aes(label = sprintf("%.1f%%", proportion * 100)),
    vjust = -0.3,
    size = 3.5
  ) +
  ylim(0, max(instable_tbl$proportion) * 1.15) +
  labs(
    title = "Condition-wise instable signal proportion",
    x = "Condition",
    y = "Instable feature proportion"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(FIG_DIR, "condition_instable_proportion.png"),
  p2,
  width = 7,
  height = 5,
  dpi = 300
)

cat("Saved condition_instable_proportion.png\n")

# ------------------------------------------------------------
# Plot 3: SSI distribution
# ------------------------------------------------------------

p3 <- ggplot(
  metrics,
  aes(
    x = SSI
  )
) +
  geom_histogram(
    bins = 30,
    boundary = 0,
    closed = "left"
  ) +
  facet_wrap(~ condition, nrow = 2) +
  geom_vline(xintercept = 0.90, linetype = "dashed") +
  geom_vline(xintercept = 0.75, linetype = "dashed") +
  geom_vline(xintercept = 0.60, linetype = "dashed") +
  labs(
    title = "Condition-wise SSI distribution",
    x = "Signal Stratification Index (SSI)",
    y = "Number of features"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(FIG_DIR, "condition_ssi_distribution.png"),
  p3,
  width = 10,
  height = 7,
  dpi = 300
)

cat("Saved condition_ssi_distribution.png\n")

# ------------------------------------------------------------
# Plot 4: Stability Deviation distribution
# ------------------------------------------------------------

p4 <- ggplot(
  metrics,
  aes(
    x = stability_deviation
  )
) +
  geom_histogram(
    bins = 30,
    boundary = 0,
    closed = "left"
  ) +
  facet_wrap(~ condition, nrow = 2) +
  geom_vline(xintercept = 0.10, linetype = "dashed") +
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = 0.40, linetype = "dashed") +
  labs(
    title = "Condition-wise Stability Deviation distribution",
    x = "Stability Deviation (1 - SSI)",
    y = "Number of features"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(FIG_DIR, "condition_stability_deviation_distribution.png"),
  p4,
  width = 10,
  height = 7,
  dpi = 300
)

cat("Saved condition_stability_deviation_distribution.png\n\n")

cat("Saved all figures to:\n")
cat(FIG_DIR, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

