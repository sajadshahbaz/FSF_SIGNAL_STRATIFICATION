#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(dplyr)
  library(tibble)
})

ROOT <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"
OUT  <- file.path(ROOT, "results/archive/pre_repair_manuscript/figures/main")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

savep <- function(p, name, w = 12, h = 7) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

boxes <- tibble(
  id = c(
    "A", "B", "C", "D", "E", "F",
    "G", "H", "I"
  ),
  x = c(1, 3, 5, 7, 9, 11, 3, 7, 10),
  y = c(5, 5, 5, 5, 5, 5, 2.3, 2.3, 2.3),
  label = c(
    "Input data\nexpression matrix\nor candidate genes",
    "Perturbation layer\nresampling / noise / threshold variation",
    "State estimation\nP(up), P(down), P(const)",
    "Signal score\nSSI and stability deviation",
    "Signal stratification\nstable / transitional / instable",
    "Condition architecture\nclass composition per condition",
    "Functional annotation\nGO / KEGG / curated themes",
    "Synthetic validation\ntruth recovery and robustness",
    "Manuscript outputs\nfigures, tables, interpretation"
  )
)

arrows <- tibble(
  x = c(1.75, 3.75, 5.75, 7.75, 9.75, 5, 7.8),
  xend = c(2.25, 4.25, 6.25, 8.25, 10.25, 3.4, 7.4),
  y = c(5, 5, 5, 5, 5, 4.55, 4.55),
  yend = c(5, 5, 5, 5, 5, 2.8, 2.8)
)

p <- ggplot() +
  geom_segment(
    data = arrows,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.18, "cm")),
    linewidth = 0.45
  ) +
  geom_label(
    data = boxes,
    aes(x = x, y = y, label = label),
    size = 3.35,
    label.size = 0.35,
    label.padding = unit(0.35, "lines"),
    fill = "white"
  ) +
  annotate(
    "text",
    x = 1,
    y = 6.25,
    label = "A",
    fontface = "bold",
    size = 6
  ) +
  annotate(
    "text",
    x = 0.45,
    y = 5.7,
    label = "FSF computational workflow",
    fontface = "bold",
    hjust = 0,
    size = 5
  ) +
  annotate(
    "text",
    x = 0.45,
    y = 1.35,
    label = "The workflow converts feature-level perturbation behavior into interpretable signal classes and condition-level biological profiles.",
    hjust = 0,
    size = 3.6
  ) +
  scale_x_continuous(limits = c(0, 12), expand = c(0, 0)) +
  scale_y_continuous(limits = c(1, 6.6), expand = c(0, 0)) +
  theme_void(base_size = 12)

savep(p, "Figure1_Workflow_Schematic", 12, 7)

cat("Saved Figure1_Workflow_Schematic.pdf/png\n")
