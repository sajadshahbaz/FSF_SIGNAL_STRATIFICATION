#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
})

ROOT <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"
OUT  <- file.path(ROOT, "results/manuscript/figures/main")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

ternary_xy <- function(p_up, p_down, p_const) {
  tibble(
    x = p_down + 0.5 * p_const,
    y = p_const * sqrt(3) / 2
  )
}

grid <- expand.grid(
  p_up = seq(0, 1, by = 0.015),
  p_down = seq(0, 1, by = 0.015)
) %>%
  mutate(p_const = 1 - p_up - p_down) %>%
  filter(p_const >= 0) %>%
  mutate(
    SSI = pmax(p_up, p_down, p_const),
    stability_deviation = 1 - SSI,
    dominant_state = case_when(
      p_up >= p_down & p_up >= p_const ~ "Stable up",
      p_down >= p_up & p_down >= p_const ~ "Stable down",
      p_const >= p_up & p_const >= p_down ~ "Stable constant"
    ),
    fsf_region = case_when(
      SSI >= 0.80 & dominant_state == "Stable up" ~ "Highly stable up",
      SSI >= 0.80 & dominant_state == "Stable down" ~ "Highly stable down",
      SSI >= 0.80 & dominant_state == "Stable constant" ~ "Highly stable constant",
      SSI < 0.45 ~ "Instability region",
      SSI < 0.80 ~ "Transitional region",
      TRUE ~ "Other"
    )
  ) %>%
  bind_cols(ternary_xy(.$p_up, .$p_down, .$p_const))

boundary <- tibble(
  p_up = c(1, 0, 0, 1),
  p_down = c(0, 1, 0, 0),
  p_const = c(0, 0, 1, 0)
) %>%
  bind_cols(ternary_xy(.$p_up, .$p_down, .$p_const))

make_line <- function(type, value) {
  s <- seq(0, 1 - value, by = 0.01)

  if (type == "up") {
    df <- tibble(p_up = value, p_down = s, p_const = 1 - value - s)
  } else if (type == "down") {
    df <- tibble(p_down = value, p_up = s, p_const = 1 - value - s)
  } else {
    df <- tibble(p_const = value, p_up = s, p_down = 1 - value - s)
  }

  bind_cols(df, ternary_xy(df$p_up, df$p_down, df$p_const))
}

guides <- bind_rows(
  lapply(c(0.2, 0.4, 0.6, 0.8), \(v) make_line("up", v)),
  lapply(c(0.2, 0.4, 0.6, 0.8), \(v) make_line("down", v)),
  lapply(c(0.2, 0.4, 0.6, 0.8), \(v) make_line("const", v))
)

colors <- c(
  "Highly stable up" = "#2166AC",
  "Highly stable constant" = "#1B7837",
  "Highly stable down" = "#E66101",
  "Transitional region" = "#7B3294",
  "Instability region" = "#B2182B"
)

p <- ggplot() +
  geom_tile(
    data = grid %>% filter(fsf_region %in% names(colors)),
    aes(x = x, y = y, fill = fsf_region),
    width = 0.017,
    height = 0.017,
    alpha = 0.82
  ) +
  geom_path(
    data = guides,
    aes(x = x, y = y, group = interaction(round(p_up, 2), round(p_down, 2), round(p_const, 2))),
    color = "grey82",
    linewidth = 0.18,
    alpha = 0.55
  ) +
  geom_path(
    data = boundary,
    aes(x = x, y = y),
    linewidth = 0.75,
    color = "black"
  ) +
  annotate("text", x = 0.50, y = 0.91, label = "P(const)", size = 5, fontface = "bold") +
  annotate("text", x = -0.035, y = -0.035, label = "P(up)", size = 5, fontface = "bold", hjust = 0) +
  annotate("text", x = 1.035, y = -0.035, label = "P(down)", size = 5, fontface = "bold", hjust = 1) +

  annotate("text", x = 0.16, y = 0.10, label = "Highly\nstable up", size = 4.2, color = "white", fontface = "bold") +
  annotate("text", x = 0.84, y = 0.10, label = "Highly\nstable down", size = 4.2, color = "white", fontface = "bold") +
  annotate("text", x = 0.50, y = 0.685, label = "Highly\nstable constant", size = 4, color = "white", fontface = "bold") +
  annotate("text", x = 0.50, y = 0.43, label = "Transitional\nregion", size = 4.1, color = "white", fontface = "bold") +
  annotate("text", x = 0.50, y = 0.27, label = "Instability\nregion", size = 4.1, color = "white", fontface = "bold") +

  scale_fill_manual(values = colors, name = "FSF region") +
  coord_equal(
    xlim = c(-0.08, 1.18),
    ylim = c(-0.08, 0.96),
    clip = "off"
  ) +
  labs(
    title = "Figure 2. FSF probability simplex and signal-state regions",
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0, size = 17),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(20, 30, 20, 20)
  )

ggsave(
  file.path(OUT, "Figure2_FSF_probability_simplex_regions.pdf"),
  p,
  width = 9.5,
  height = 7,
  device = cairo_pdf
)

ggsave(
  file.path(OUT, "Figure2_FSF_probability_simplex_regions.png"),
  p,
  width = 9.5,
  height = 7,
  dpi = 300
)

cat("Saved Figure2_FSF_probability_simplex_regions.pdf/png\n")
