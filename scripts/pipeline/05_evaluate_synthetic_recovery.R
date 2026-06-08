#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 05: Evaluate Synthetic Recovery
#
# Inputs:
#   results/synthetic/synthetic_truth_table.tsv
#   results/fsf_metrics/fsf_ssi_stability_deviation.tsv
#
# Outputs:
#   results/synthetic/synthetic_recovery_confusion_matrix.tsv
#   results/synthetic/synthetic_recovery_summary.tsv
#   results/figures/synthetic_truth_vs_fsf_class_heatmap.png
#   results/logs/05_evaluate_synthetic_recovery.log
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

TRUTH_FILE <- file.path(
  ROOT_DIR,
  "results",
  "synthetic",
  "synthetic_truth_table.tsv"
)

FSF_FILE <- file.path(
  ROOT_DIR,
  "results",
  "fsf_metrics",
  "fsf_ssi_stability_deviation.tsv"
)

OUT_CONF <- file.path(
  ROOT_DIR,
  "results",
  "synthetic",
  "synthetic_recovery_confusion_matrix.tsv"
)

OUT_SUMMARY <- file.path(
  ROOT_DIR,
  "results",
  "synthetic",
  "synthetic_recovery_summary.tsv"
)

OUT_FIG <- file.path(
  ROOT_DIR,
  "results",
  "figures",
  "synthetic_truth_vs_fsf_class_heatmap.png"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results",
  "logs",
  "05_evaluate_synthetic_recovery.log"
)

sink(LOG_FILE, split = TRUE)

cat("====================================================\n")
cat("FSF Synthetic Recovery Evaluation\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("====================================================\n\n")

truth <- read_tsv(TRUTH_FILE, show_col_types = FALSE)
fsf <- read_tsv(FSF_FILE, show_col_types = FALSE)

cat("Truth features:", nrow(truth), "\n")
cat("FSF features:", nrow(fsf), "\n\n")

merged <- truth %>%
  inner_join(
    fsf %>%
      select(
        feature_id,
        fsf_signal_class,
        stability_level,
        SSI,
        stability_deviation
      ),
    by = "feature_id"
  )

cat("Merged features:", nrow(merged), "\n\n")

# ----------------------------------------------------
# Confusion matrix
# ----------------------------------------------------

confusion <- merged %>%
  count(
    truth_class,
    fsf_signal_class,
    name = "n_features"
  ) %>%
  arrange(
    truth_class,
    desc(n_features)
  )

write_tsv(
  confusion,
  OUT_CONF
)

# ----------------------------------------------------
# Recovery summary
# ----------------------------------------------------

summary_tbl <- merged %>%
  group_by(truth_class) %>%
  summarise(
    n_features = n(),
    mean_ssi = mean(SSI),
    median_ssi = median(SSI),
    mean_sd = mean(stability_deviation),
    median_sd = median(stability_deviation),
    .groups = "drop"
  )

write_tsv(
  summary_tbl,
  OUT_SUMMARY
)

# ----------------------------------------------------
# Heatmap
# ----------------------------------------------------

heatmap_df <- confusion %>%
  group_by(truth_class) %>%
  mutate(
    proportion = n_features / sum(n_features)
  ) %>%
  ungroup()

p <- ggplot(
  heatmap_df,
  aes(
    x = fsf_signal_class,
    y = truth_class,
    fill = proportion
  )
) +
  geom_tile() +
  geom_text(
    aes(
      label = sprintf("%.2f", proportion)
    ),
    size = 3
  ) +
  labs(
    title = "Synthetic Truth vs FSF Classification",
    x = "FSF class",
    y = "Truth class",
    fill = "Proportion"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  OUT_FIG,
  p,
  width = 9,
  height = 6,
  dpi = 300
)

cat("Saved:\n")
cat(OUT_CONF, "\n")
cat(OUT_SUMMARY, "\n")
cat(OUT_FIG, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")

sink()

