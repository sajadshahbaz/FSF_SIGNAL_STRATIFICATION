#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 03: Diagnostic Plots
#
# Input:
#   results/fsf_metrics/fsf_ssi_stability_deviation.tsv
#
# Outputs:
#   results/figures/fsf_ssi_distribution.png
#   results/figures/fsf_stability_deviation_distribution.png
#   results/figures/fsf_signal_class_barplot.png
#   results/figures/fsf_state_probability_scatter.png
#   results/tables/fsf_probability_triangle_ready_table.tsv
#   results/logs/03_fsf_diagnostic_plots.log
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

INPUT_FILE <- file.path(ROOT_DIR, "results", "fsf_metrics", "fsf_ssi_stability_deviation.tsv")

FIG_DIR <- file.path(ROOT_DIR, "results", "figures")
TABLE_DIR <- file.path(ROOT_DIR, "results", "tables")
LOG_DIR <- file.path(ROOT_DIR, "results", "logs")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

SSI_PLOT <- file.path(FIG_DIR, "fsf_ssi_distribution.png")
SD_PLOT <- file.path(FIG_DIR, "fsf_stability_deviation_distribution.png")
CLASS_BARPLOT <- file.path(FIG_DIR, "fsf_signal_class_barplot.png")
STATE_SCATTER <- file.path(FIG_DIR, "fsf_state_probability_scatter.png")
TRIANGLE_TABLE <- file.path(TABLE_DIR, "fsf_probability_triangle_ready_table.tsv")
LOG_FILE <- file.path(LOG_DIR, "03_fsf_diagnostic_plots.log")

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 03: Diagnostic Plots\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

cat("Project root:", ROOT_DIR, "\n")
cat("Input file:", INPUT_FILE, "\n\n")

if (!file.exists(INPUT_FILE)) {
  stop("Input file not found: ", INPUT_FILE)
}

fsf <- read_tsv(INPUT_FILE, show_col_types = FALSE)

required_cols <- c(
  "feature_id",
  "p_up",
  "p_down",
  "p_const",
  "SSI",
  "stability_deviation",
  "stability_level",
  "fsf_signal_class"
)

missing_cols <- setdiff(required_cols, colnames(fsf))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

fsf <- fsf %>%
  mutate(
    feature_id = as.character(feature_id),
    stability_level = factor(
      stability_level,
      levels = c(
        "highly_stable",
        "stable",
        "weakly_stable_transitional",
        "instable"
      )
    ),
    fsf_signal_class = as.factor(fsf_signal_class),
    SSI = as.numeric(SSI),
    stability_deviation = as.numeric(stability_deviation),
    p_up = as.numeric(p_up),
    p_down = as.numeric(p_down),
    p_const = as.numeric(p_const)
  )

cat("Loaded features:", nrow(fsf), "\n")
cat("FSF classes:\n")
print(table(fsf$fsf_signal_class))
cat("\n")

# ------------------------------------------------------------
# Plot 1: SSI distribution
# ------------------------------------------------------------

p1 <- ggplot(fsf, aes(x = SSI)) +
  geom_histogram(bins = 30, boundary = 0, closed = "left") +
  geom_vline(xintercept = 0.90, linetype = "dashed") +
  geom_vline(xintercept = 0.75, linetype = "dashed") +
  geom_vline(xintercept = 0.60, linetype = "dashed") +
  labs(
    title = "FSF Signal Stratification Index Distribution",
    x = "Signal Stratification Index (SSI)",
    y = "Number of features"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = SSI_PLOT,
  plot = p1,
  width = 7,
  height = 5,
  dpi = 300
)

cat("Saved:", SSI_PLOT, "\n")

# ------------------------------------------------------------
# Plot 2: Stability deviation distribution
# ------------------------------------------------------------

p2 <- ggplot(fsf, aes(x = stability_deviation)) +
  geom_histogram(bins = 30, boundary = 0, closed = "left") +
  geom_vline(xintercept = 0.10, linetype = "dashed") +
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = 0.40, linetype = "dashed") +
  labs(
    title = "FSF Stability Deviation Distribution",
    x = "Stability Deviation (1 - SSI)",
    y = "Number of features"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = SD_PLOT,
  plot = p2,
  width = 7,
  height = 5,
  dpi = 300
)

cat("Saved:", SD_PLOT, "\n")

# ------------------------------------------------------------
# Plot 3: Signal class barplot
# ------------------------------------------------------------

class_counts <- fsf %>%
  count(fsf_signal_class, name = "n_features") %>%
  arrange(desc(n_features))

p3 <- ggplot(class_counts, aes(x = reorder(fsf_signal_class, n_features), y = n_features)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "FSF Signal Class Composition",
    x = "FSF signal class",
    y = "Number of features"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = CLASS_BARPLOT,
  plot = p3,
  width = 8,
  height = 5,
  dpi = 300
)

cat("Saved:", CLASS_BARPLOT, "\n")

# ------------------------------------------------------------
# Plot 4: State probability scatter
# Simple 2D diagnostic: P(up) vs P(down), point size = P(const)
# ------------------------------------------------------------

p4 <- ggplot(fsf, aes(x = p_up, y = p_down, size = p_const, shape = stability_level)) +
  geom_point(alpha = 0.75) +
  labs(
    title = "FSF State Probability Diagnostic",
    x = "P(up)",
    y = "P(down)",
    size = "P(const)",
    shape = "Stability level"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = STATE_SCATTER,
  plot = p4,
  width = 7,
  height = 5,
  dpi = 300
)

cat("Saved:", STATE_SCATTER, "\n")

# ------------------------------------------------------------
# Triangle-ready table
# ------------------------------------------------------------

triangle_ready <- fsf %>%
  transmute(
    feature_id,
    p_up,
    p_down,
    p_const,
    SSI,
    stability_deviation,
    stability_level,
    fsf_signal_class
  )

write_tsv(triangle_ready, TRIANGLE_TABLE)

cat("Saved:", TRIANGLE_TABLE, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

