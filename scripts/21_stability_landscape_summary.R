#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 21: Stability Landscape Summary
#
# Purpose:
#   Summarize condition-level stability architecture from
#   FSF signal strata and generate publication-readable figures.
#
# This step replaces the archived instability-driver steps.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

INPUT_FILE <- file.path(
  ROOT_DIR,
  "results/real_data/condition_fsf/condition_wise_fsf_metrics.tsv"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/stability_landscape"
)

FIG_DIR <- file.path(
  ROOT_DIR,
  "results/real_data/figures"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results/logs/21_stability_landscape_summary.log"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 21: Stability Landscape Summary\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

fsf <- read_tsv(INPUT_FILE, show_col_types = FALSE)

required_cols <- c(
  "condition",
  "feature_id",
  "stability_level",
  "dominant_state",
  "fsf_signal_class",
  "SSI",
  "stability_deviation"
)

missing_cols <- setdiff(required_cols, colnames(fsf))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

condition_order <- c("osm", "lt", "ht", "gam", "des", "uv")

condition_labels <- c(
  osm = "OSM",
  lt  = "LT",
  ht  = "HT",
  gam = "GAM",
  des = "DES",
  uv  = "UV"
)

level_order <- c(
  "highly_stable",
  "stable",
  "weakly_stable_transitional",
  "instable"
)

level_labels <- c(
  highly_stable = "Highly stable",
  stable = "Stable",
  weakly_stable_transitional = "Weakly stable",
  instable = "Unstable"
)

class_order <- c(
  "highly_stable_constant",
  "highly_stable_up",
  "highly_stable_down",
  "stable_constant",
  "stable_up",
  "stable_down",
  "weakly_stable_transitional_constant",
  "weakly_stable_transitional_up",
  "weakly_stable_transitional_down",
  "instable"
)

class_labels <- c(
  highly_stable_constant = "Highly stable constant",
  highly_stable_up = "Highly stable up",
  highly_stable_down = "Highly stable down",
  stable_constant = "Stable constant",
  stable_up = "Stable up",
  stable_down = "Stable down",
  weakly_stable_transitional_constant = "Weakly stable constant",
  weakly_stable_transitional_up = "Weakly stable up",
  weakly_stable_transitional_down = "Weakly stable down",
  instable = "Unstable"
)

class_labels_short <- c(
  highly_stable_constant = "HSC",
  highly_stable_up = "HSU",
  highly_stable_down = "HSD",
  stable_constant = "SC",
  stable_up = "SU",
  stable_down = "SD",
  weakly_stable_transitional_constant = "WSC",
  weakly_stable_transitional_up = "WSU",
  weakly_stable_transitional_down = "WSD",
  instable = "UNS"
)

base_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

fsf <- fsf %>%
  mutate(
    condition = factor(condition, levels = condition_order),
    condition_label = factor(
      condition_labels[as.character(condition)],
      levels = condition_labels[condition_order]
    ),
    stability_level = factor(stability_level, levels = level_order),
    stability_level_label = factor(
      level_labels[as.character(stability_level)],
      levels = level_labels[level_order]
    ),
    fsf_signal_class = factor(fsf_signal_class, levels = class_order),
    fsf_class_label = factor(
      class_labels[as.character(fsf_signal_class)],
      levels = class_labels[class_order]
    ),
    fsf_class_short = factor(
      class_labels_short[as.character(fsf_signal_class)],
      levels = class_labels_short[class_order]
    ),
    SSI = as.numeric(SSI),
    stability_deviation = as.numeric(stability_deviation)
  )

cat("Rows:", nrow(fsf), "\n")
cat("Features:", n_distinct(fsf$feature_id), "\n")
cat("Conditions:", paste(levels(droplevels(fsf$condition)), collapse = ", "), "\n\n")

# ------------------------------------------------------------
# 1. Full FSF class landscape
# ------------------------------------------------------------

class_landscape <- fsf %>%
  count(
    condition,
    condition_label,
    fsf_signal_class,
    fsf_class_label,
    fsf_class_short,
    name = "n_features"
  ) %>%
  group_by(condition, condition_label) %>%
  mutate(
    proportion = n_features / sum(n_features),
    percent = proportion * 100
  ) %>%
  ungroup() %>%
  arrange(condition, fsf_signal_class)

write_tsv(
  class_landscape,
  file.path(OUT_DIR, "condition_fsf_class_landscape.tsv")
)

# ------------------------------------------------------------
# 2. Collapsed stability architecture
# ------------------------------------------------------------

stability_architecture <- fsf %>%
  count(
    condition,
    condition_label,
    stability_level,
    stability_level_label,
    name = "n_features"
  ) %>%
  group_by(condition, condition_label) %>%
  mutate(
    proportion = n_features / sum(n_features),
    percent = proportion * 100
  ) %>%
  ungroup() %>%
  arrange(condition, stability_level)

write_tsv(
  stability_architecture,
  file.path(OUT_DIR, "condition_stability_architecture.tsv")
)

# ------------------------------------------------------------
# 3. Directional stability architecture
# ------------------------------------------------------------

directional_architecture <- fsf %>%
  count(
    condition,
    condition_label,
    dominant_state,
    stability_level,
    stability_level_label,
    name = "n_features"
  ) %>%
  group_by(condition, condition_label) %>%
  mutate(
    proportion = n_features / sum(n_features),
    percent = proportion * 100
  ) %>%
  ungroup() %>%
  arrange(condition, dominant_state, stability_level)

write_tsv(
  directional_architecture,
  file.path(OUT_DIR, "condition_directional_stability_architecture.tsv")
)

# ------------------------------------------------------------
# 4. Condition-level summary statistics
# ------------------------------------------------------------

condition_summary <- fsf %>%
  group_by(condition, condition_label) %>%
  summarise(
    n_features = n(),
    mean_SSI = mean(SSI),
    median_SSI = median(SSI),
    mean_stability_deviation = mean(stability_deviation),
    median_stability_deviation = median(stability_deviation),
    highly_stable_fraction = mean(stability_level == "highly_stable"),
    stable_fraction = mean(stability_level == "stable"),
    weakly_stable_fraction = mean(stability_level == "weakly_stable_transitional"),
    unstable_fraction = mean(stability_level == "instable"),
    stable_architecture_fraction = mean(stability_level %in% c("highly_stable", "stable")),
    non_unstable_fraction = mean(stability_level != "instable"),
    .groups = "drop"
  ) %>%
  arrange(condition)

write_tsv(
  condition_summary,
  file.path(OUT_DIR, "condition_stability_summary.tsv")
)

cat("Condition stability summary:\n")
print(condition_summary, n = 100)
cat("\n")

# ------------------------------------------------------------
# 5. Abbreviation key
# ------------------------------------------------------------

abbr_key <- tibble::tibble(
  abbreviation = c("HSC", "HSU", "HSD", "SC", "SU", "SD", "WSC", "WSU", "WSD", "UNS"),
  meaning = c(
    "Highly stable constant",
    "Highly stable up",
    "Highly stable down",
    "Stable constant",
    "Stable up",
    "Stable down",
    "Weakly stable constant",
    "Weakly stable up",
    "Weakly stable down",
    "Unstable"
  )
)

write_tsv(
  abbr_key,
  file.path(OUT_DIR, "FSF_class_abbreviation_key.tsv")
)

# ------------------------------------------------------------
# 6. Heatmap of full FSF class landscape
# ------------------------------------------------------------

heatmap_tbl <- class_landscape %>%
  mutate(
    fsf_class_label = factor(
      fsf_class_label,
      levels = rev(class_labels[class_order])
    )
  )

p1 <- ggplot(
  heatmap_tbl,
  aes(
    x = condition_label,
    y = fsf_class_label,
    fill = percent
  )
) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = sprintf("%.1f", percent)),
    size = 3
  ) +
  labs(
    title = "FSF stability landscape across stress conditions",
    x = "Condition",
    y = "FSF stratum",
    fill = "Genes (%)"
  ) +
  base_theme +
  theme(
    axis.text.y = element_text(size = 9)
  )

ggsave(
  file.path(FIG_DIR, "stability_landscape_heatmap.png"),
  p1,
  width = 9.5,
  height = 7.2,
  dpi = 300
)

ggsave(
  file.path(FIG_DIR, "stability_landscape_heatmap.pdf"),
  p1,
  width = 9.5,
  height = 7.2
)

# ------------------------------------------------------------
# 7. Collapsed stability architecture stacked barplot
# ------------------------------------------------------------

bar_label_tbl <- stability_architecture %>%
  mutate(
    label = ifelse(
      percent >= 4,
      sprintf("%.1f%%", percent),
      ""
    )
  )

p2 <- ggplot(
  bar_label_tbl,
  aes(
    x = condition_label,
    y = proportion,
    fill = stability_level_label
  )
) +
  geom_col(color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  labs(
    title = "FSF stability architecture across conditions",
    x = "Condition",
    y = "Proportion of genes",
    fill = "Stability level"
  ) +
  base_theme

ggsave(
  file.path(FIG_DIR, "stability_architecture_stacked_barplot.png"),
  p2,
  width = 8.8,
  height = 5.4,
  dpi = 300
)

ggsave(
  file.path(FIG_DIR, "stability_architecture_stacked_barplot.pdf"),
  p2,
  width = 8.8,
  height = 5.4
)

# ------------------------------------------------------------
# 8. Mean Stability Deviation by condition
# ------------------------------------------------------------

p3 <- ggplot(
  condition_summary,
  aes(
    x = condition_label,
    y = mean_stability_deviation
  )
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = sprintf("%.3f", mean_stability_deviation)),
    vjust = -0.3,
    size = 3.5
  ) +
  labs(
    title = "Average Stability Deviation across conditions",
    x = "Condition",
    y = "Mean Stability Deviation"
  ) +
  expand_limits(y = max(condition_summary$mean_stability_deviation) * 1.15) +
  base_theme

ggsave(
  file.path(FIG_DIR, "mean_stability_deviation_by_condition.png"),
  p3,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(FIG_DIR, "mean_stability_deviation_by_condition.pdf"),
  p3,
  width = 7,
  height = 5
)

cat("Saved outputs to:\n")
cat(OUT_DIR, "\n")
cat(FIG_DIR, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

