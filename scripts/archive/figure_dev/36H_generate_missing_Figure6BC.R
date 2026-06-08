#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

ROOT <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"
OUT  <- file.path(ROOT, "results/manuscript/figures/main")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

theme_pub <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.title = element_text(face = "bold")
    )
}

savep <- function(p, name, w = 8, h = 5) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  cat("Saved:", name, "\n")
}

# ============================================================
# Figure 6B: mean SSI and stability deviation across noise
# ============================================================

summary_file <- file.path(
  ROOT,
  "results/synthetic/noise_gradient/noise_gradient_summary.tsv"
)

summary_tbl <- read_tsv(summary_file, show_col_types = FALSE)

p6b <- summary_tbl %>%
  select(noise_sd, truth_class, mean_ssi, mean_sd) %>%
  pivot_longer(
    cols = c(mean_ssi, mean_sd),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      mean_ssi = "Mean SSI",
      mean_sd  = "Mean stability deviation"
    )
  ) %>%
  ggplot(aes(x = noise_sd, y = value, linetype = truth_class)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~metric, nrow = 1) +
  labs(
    title = "Figure 6B. Noise-gradient effect on FSF stability metrics",
    x = "Noise standard deviation",
    y = "Metric value",
    linetype = "Synthetic truth class"
  ) +
  theme_pub()

savep(p6b, "Figure6B_noise_gradient_stability_metrics", 10, 5)

# ============================================================
# Figure 6C: class recovery across noise
# Correct recovery = mapped true class recovered as compatible FSF class
# ============================================================

conf_file <- file.path(
  ROOT,
  "results/synthetic/noise_gradient/noise_gradient_confusion.tsv"
)

conf <- read_tsv(conf_file, show_col_types = FALSE)

recovery <- conf %>%
  mutate(
    correct = case_when(
      truth_class == "instable" &
        fsf_signal_class == "instable" ~ TRUE,

      truth_class == "stable_constant" &
        fsf_signal_class %in% c(
          "stable_constant",
          "highly_stable_constant",
          "weakly_stable_transitional_constant"
        ) ~ TRUE,

      truth_class == "stable_down" &
        fsf_signal_class %in% c(
          "stable_down",
          "highly_stable_down",
          "weakly_stable_transitional_down"
        ) ~ TRUE,

      truth_class == "stable_up" &
        fsf_signal_class %in% c(
          "stable_up",
          "highly_stable_up",
          "weakly_stable_transitional_up"
        ) ~ TRUE,

      TRUE ~ FALSE
    )
  ) %>%
  group_by(noise_sd, truth_class) %>%
  summarise(
    recovered = sum(n_features[correct], na.rm = TRUE),
    total = sum(n_features, na.rm = TRUE),
    recovery_rate = recovered / total,
    .groups = "drop"
  )

write_tsv(
  recovery,
  file.path(ROOT, "results/synthetic/noise_gradient/noise_gradient_recovery_rate.tsv")
)

p6c <- ggplot(
  recovery,
  aes(x = noise_sd, y = recovery_rate, linetype = truth_class)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Figure 6C. Synthetic class recovery across noise",
    x = "Noise standard deviation",
    y = "Recovery rate",
    linetype = "Synthetic truth class"
  ) +
  theme_pub()

savep(p6c, "Figure6C_noise_gradient_class_recovery", 8.5, 5)

cat("Figure 6B and 6C completed.\n")
