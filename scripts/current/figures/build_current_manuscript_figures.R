#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(png)
  library(scales)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
figure5_only <- identical(commandArgs(trailingOnly = TRUE), "--figure5-only")
source_dir <- file.path(root, "results", "current_fsf_v1", "manuscript", "source_data")
review_dir <- file.path(root, "results", "current_fsf_v1", "manuscript", "review")
output_dir <- file.path(root, "results", "current_fsf_v1", "manuscript", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

conditions <- c("DES", "GAM", "HT", "LT", "OSM", "UV")
regions <- c("Low Stability", "Transitional", "Stable", "Highly Stable")
region_colors <- c(
  "Low Stability" = "#7A5195", "Transitional" = "#EF8354",
  "Stable" = "#2F6690", "Highly Stable" = "#3A7D44"
)
condition_colors <- c(
  DES = "#4477AA", GAM = "#EE6677", HT = "#228833",
  LT = "#CCBB44", OSM = "#66CCEE", UV = "#AA3377"
)

read_source <- function(name) {
  path <- file.path(source_dir, name)
  if (!file.exists(path)) stop("Missing frozen figure source: ", path)
  read.delim(path, check.names = FALSE, quote = "", stringsAsFactors = FALSE)
}

theme_fsf <- function() {
  theme_minimal(base_size = 10, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "#444444"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

save_figure <- function(plot, id, width, height) {
  png_path <- file.path(output_dir, paste0(id, ".png"))
  pdf_path <- file.path(output_dir, paste0(id, ".pdf"))
  ggsave(pdf_path, plot, width = width, height = height, units = "in",
         device = grDevices::pdf, bg = "white", useDingbats = FALSE)
  png_prefix <- sub("[.]png$", "", png_path)
  status <- system2("pdftocairo", c("-png", "-singlefile", "-r", "300", pdf_path, png_prefix))
  if (status != 0L || !file.exists(png_path)) stop("PDF-to-PNG conversion failed for ", id)
  c(png = png_path, pdf = pdf_path)
}

flow_plot <- function(labels, title, subtitle) {
  n <- length(labels)
  nodes <- data.frame(x = seq_len(n), y = 1, label = labels)
  arrows <- data.frame(x = seq_len(n - 1) + 0.32, xend = seq_len(n - 1) + 0.68, y = 1, yend = 1)
  ggplot() +
    geom_segment(data = arrows, aes(x = x, xend = xend, y = y, yend = yend),
                 arrow = arrow(length = grid::unit(0.14, "inches")), linewidth = 0.7, color = "#555555") +
    geom_label(data = nodes, aes(x = x, y = y, label = label),
               size = 3.4, linewidth = 0.35, label.padding = grid::unit(0.25, "lines"),
               fill = "#F5F7FA", color = "#1F2933") +
    coord_cartesian(xlim = c(0.55, n + 0.45), ylim = c(0.65, 1.35), clip = "off") +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "#444444"),
          plot.margin = margin(18, 18, 18, 18))
}

# Figure 2: conceptual geometry plus observed six-condition distributions.
f2 <- read_source("figure2_simplex_source.tsv")
f2$condition <- factor(f2$condition, conditions)
f2$stability_region <- factor(f2$stability_region, regions)
f2$x <- f2$p_down + 0.5 * f2$p_const
f2$y <- sqrt(3) / 2 * f2$p_const
triangle <- data.frame(x = c(0, 1, 0.5, 0), y = c(0, 0, sqrt(3) / 2, 0))
vertex_labels <- do.call(rbind, lapply(conditions, function(z) {
  data.frame(condition = factor(z, conditions), x = c(0, 1, 0.5),
             y = c(-0.035, -0.035, sqrt(3) / 2 + 0.035), label = c("Up", "Down", "Constant"))
}))
simplex_grid <- expand.grid(p_up = seq(0, 1, by = 0.01), p_down = seq(0, 1, by = 0.01))
simplex_grid$p_const <- 1 - simplex_grid$p_up - simplex_grid$p_down
simplex_grid <- simplex_grid[simplex_grid$p_const >= 0, ]
simplex_grid$ssi <- pmax(simplex_grid$p_up, simplex_grid$p_down, simplex_grid$p_const)
simplex_grid$stability_region <- cut(
  simplex_grid$ssi, breaks = c(-Inf, 0.50, 0.75, 0.90, Inf),
  labels = regions, right = FALSE
)
simplex_grid$stability_region[simplex_grid$ssi <= 0.50] <- "Low Stability"
simplex_grid$x <- simplex_grid$p_down + 0.5 * simplex_grid$p_const
simplex_grid$y <- sqrt(3) / 2 * simplex_grid$p_const
p2a <- ggplot(simplex_grid, aes(x, y, color = stability_region)) +
  geom_point(size = 0.55, alpha = 0.9) +
  geom_path(data = triangle, aes(x, y), inherit.aes = FALSE, linewidth = 0.65, color = "#222222") +
  annotate("text", x = c(0, 1, 0.5), y = c(-0.035, -0.035, sqrt(3) / 2 + 0.035),
           label = c("P(Up)", "P(Down)", "P(Constant)"), fontface = "bold", size = 3) +
  scale_color_manual(values = region_colors, drop = FALSE) +
  coord_equal(xlim = c(-0.05, 1.05), ylim = c(-0.07, 0.92), clip = "off") +
  labs(title = "A. FSF probability geometry", color = "FSF region", x = NULL, y = NULL) +
  theme_fsf() + theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())

f2_agg <- aggregate(feature_id ~ condition + x + y + stability_region, f2, length)
names(f2_agg)[5] <- "n_features"
p2b <- ggplot(f2_agg, aes(x, y, color = stability_region, size = n_features)) +
  geom_path(data = triangle, aes(x, y), inherit.aes = FALSE, linewidth = 0.45, color = "#333333") +
  geom_point(alpha = 0.52) +
  geom_text(data = vertex_labels, aes(x, y, label = label), inherit.aes = FALSE, size = 2.6) +
  facet_wrap(~condition, nrow = 2) +
  scale_color_manual(values = region_colors, drop = FALSE) +
  scale_size_area(max_size = 4.8, breaks = c(1, 10, 100, 1000)) +
  coord_equal(xlim = c(-0.03, 1.03), ylim = c(-0.06, 0.92), clip = "off") +
  labs(title = "B. Observed condition distributions", color = "FSF region",
       size = "Features at\ncoordinate", x = NULL, y = NULL) +
  theme_fsf() + theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
figure2 <- p2a + p2b + plot_layout(widths = c(0.9, 1.6)) +
  plot_annotation(title = "Figure 2. Probability geometry and observed FSF distributions")
if (!figure5_only) save_figure(figure2, "Figure2", 14.0, 6.8)

# Figure 3: locked condition architecture.
f3 <- read_source("figure3_architecture_source.tsv")
f3$condition <- factor(f3$condition, conditions)
f3$stability_region <- factor(f3$stability_region, regions)
p3a <- ggplot(f3, aes(condition, region_proportion, fill = stability_region)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  geom_text(aes(label = ifelse(region_proportion >= 0.08, comma(region_count), "")),
            position = position_stack(vjust = 0.5), color = "white", size = 2.8, fontface = "bold") +
  scale_fill_manual(values = region_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "A. Region composition", subtitle = "Each condition contains 15,187 features",
       x = "Condition", y = "Feature proportion", fill = "FSF region") +
  theme_fsf()
p3b <- ggplot(f3, aes(stability_region, condition, fill = region_proportion)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = paste0(comma(region_count), "\n", percent(region_proportion, accuracy = 0.1))),
            size = 2.7) +
  scale_fill_gradient(low = "#F7FBFF", high = "#225EA8", limits = c(0, 1), labels = percent_format()) +
  labs(title = "B. Architecture summary", subtitle = "UV is Transitional-dominant",
       x = "FSF region", y = "Condition", fill = "Proportion") +
  theme_fsf() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
figure3 <- p3a + p3b + plot_layout(widths = c(1.05, 1)) +
  plot_annotation(title = "Figure 3. FSF condition architecture")
if (!figure5_only) save_figure(figure3, "Figure3", 12.2, 5.5)

# Figure 4: SSI empirical distributions.
f4 <- read_source("figure4_ssi_ecdf_source.tsv")
f4$condition <- factor(f4$condition, conditions)
figure4 <- ggplot(f4, aes(ssi, color = condition)) +
  stat_ecdf(linewidth = 0.8, pad = FALSE) +
  geom_vline(xintercept = c(0.50, 0.75, 0.90), linetype = c("solid", "dashed", "dashed"), color = "#777777", linewidth = 0.45) +
  scale_color_manual(values = condition_colors, drop = FALSE) +
  scale_x_continuous(limits = c(1 / 3, 1), breaks = c(1 / 3, 0.50, 0.75, 0.90, 1.00),
                     labels = c("0.33", "0.50", "0.75", "0.90", "1.00")) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Figure 4. Empirical SSI distributions by condition", subtitle = "Vertical guides mark the FSF region boundaries",
       x = "Signal Stability Index (SSI)", y = "Empirical cumulative proportion", color = "Condition") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "#444444"),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_line(color = "#ECECEC", linewidth = 0.25),
        legend.title = element_text(face = "bold"))
if (!figure5_only) save_figure(figure4, "Figure4", 8.4, 5.2)

# Figure 5: quantitative enrichment evidence only; no manual theme labels.
go <- read_source("figure5_go_source.tsv")
ke <- read_source("figure5_kegg_source.tsv")
go_sig <- go[go$significant %in% TRUE & go$gene_set_family == "main_class", ]
ke_sig <- ke[ke$significant %in% TRUE & ke$gene_set_family == "main_class", ]

summarize_enrichment <- function(x, term_column, source_label) {
  counts <- aggregate(x[[term_column]], list(condition = x$condition,
    stability_region = x$stability_region), function(z) length(unique(z)))
  names(counts)[3] <- "significant_terms"
  strongest <- aggregate(x$adjusted_p_value, list(condition = x$condition,
    stability_region = x$stability_region), min)
  names(strongest)[3] <- "adjusted_p_value"
  observed <- merge(counts, strongest, by = c("condition", "stability_region"))
  complete <- expand.grid(condition = conditions, stability_region = regions,
                          stringsAsFactors = FALSE)
  complete <- merge(complete, observed, by = c("condition", "stability_region"), all.x = TRUE)
  complete$significant_terms[is.na(complete$significant_terms)] <- 0L
  complete$evidence_strength <- -log10(complete$adjusted_p_value)
  complete$source <- source_label
  complete
}
f5 <- rbind(summarize_enrichment(go_sig, "go_term", "GO"),
            summarize_enrichment(ke_sig, "kegg_term", "KEGG"))
f5$condition <- factor(f5$condition, conditions)
f5$stability_region <- factor(f5$stability_region, regions)
f5$source <- factor(f5$source, c("GO", "KEGG"))
f5$count_shade <- sqrt(f5$significant_terms /
  ave(f5$significant_terms, f5$source, FUN = max))
f5$strength_label <- sprintf("%.1f", f5$evidence_strength)
f5_evidence <- f5[!is.na(f5$evidence_strength), ]
f5_empty <- f5[is.na(f5$evidence_strength), ]

p5a <- ggplot(f5, aes(stability_region, condition, fill = stability_region,
                       alpha = count_shade)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = significant_terms), alpha = 1, size = 3.8, fontface = "bold") +
  facet_wrap(~source, nrow = 1) +
  scale_fill_manual(values = region_colors, guide = "none") +
  scale_alpha(range = c(0.13, 0.95), guide = "none") +
  labs(title = "A. Enrichment breadth",
       subtitle = "Unique significant term IDs; color intensity is scaled within GO or KEGG",
       x = NULL, y = "Condition") +
  theme_fsf() + theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 10.5),
                      axis.text.y = element_text(size = 10.5),
                      plot.title = element_text(face = "bold", size = 14),
                      plot.subtitle = element_text(size = 11),
                      axis.title = element_text(face = "bold", size = 12.5),
                      strip.text = element_text(face = "bold", size = 11),
                      panel.grid = element_blank())

p5b <- ggplot(f5, aes(stability_region, evidence_strength, fill = stability_region)) +
  geom_col(data = f5_evidence, width = 0.72) +
  geom_point(data = f5_empty, aes(y = 0), shape = 21, size = 1.55,
             fill = "white", color = "#9AA3AC", stroke = 0.45) +
  geom_text(data = f5_evidence, aes(label = strength_label),
            vjust = -0.35, size = 3.4, fontface = "bold") +
  facet_grid(source ~ condition) +
  scale_fill_manual(values = region_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(title = "B. Enrichment strength",
       subtitle = "Strongest current evidence; hollow baseline marker indicates no significant term",
       x = NULL, y = expression(-log[10](adjusted~p))) +
  theme_fsf() + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 10),
                      axis.text.y = element_text(size = 10.5), strip.text = element_text(face = "bold", size = 11),
                      plot.title = element_text(face = "bold", size = 14),
                      plot.subtitle = element_text(size = 11),
                      axis.title = element_text(face = "bold", size = 12.5),
                      panel.spacing = grid::unit(0.8, "lines"))

des_audit_panel_path <- file.path(
  output_dir, "biological_theme_prototypes",
  "Figure5C_DES_biological_themes_barplot_AUDIT.png"
)
des_audit_panel_hash <- "e106bba69f1a7069463c29b390a23876bc15a4e15cf57c1c99cc38e5e6e807cb"
if (!file.exists(des_audit_panel_path)) stop("Missing audited DES panel: ", des_audit_panel_path)
observed_des_panel_hash <- sub(
  " .*", "", system2("sha256sum", des_audit_panel_path, stdout = TRUE)[1]
)
if (!identical(observed_des_panel_hash, des_audit_panel_hash)) {
  stop("Audited DES panel hash does not match the authoritative prototype")
}
des_audit_panel <- png::readPNG(des_audit_panel_path)
p5c <- ggplot() +
  annotation_custom(
    grid::rasterGrob(des_audit_panel, interpolate = TRUE),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  annotate("text", x = 0.006, y = 0.992, label = "C.",
           hjust = 0, vjust = 1, fontface = "bold", size = 5.0) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void()

figure5 <- p5a / p5b / p5c +
  plot_layout(heights = c(0.78, 1.02, 2.42), guides = "collect") +
  plot_annotation(
    title = "Figure 5. Functional enrichment evidence across FSF signal architectures",
    subtitle = "Current GO and KEGG enrichment evidence across conditions and FSF regions",
    theme = theme(plot.title = element_text(face = "bold", size = 17),
      plot.subtitle = element_text(size = 12, color = "#444444"),
      plot.margin = margin(8, 8, 5, 8))) &
  theme(legend.position = "bottom")
save_figure(figure5, "Figure5", 14.0, 17.0)

# Display-only labels keep historical truth provenance out of current FSF terminology.
truth_label <- function(x) {
  map <- c(instable = "No-majority truth", weakly_stable = "Weak-majority truth",
           stable_constant = "Stable constant truth", stable_down = "Stable down truth", stable_up = "Stable up truth")
  unname(ifelse(x %in% names(map), map[x], x))
}

# Figure 6: frozen baseline and noise benchmark summaries.
f6a <- read_source("figure6_baseline_source.tsv")
f6a$truth_display <- truth_label(f6a$truth_scenario)
f6b <- read_source("figure6_noise_source.tsv")
f6b$truth_display <- truth_label(f6b$truth_scenario)
p6a <- ggplot(f6a, aes(signal_class, truth_display, fill = proportion)) +
  geom_tile(color = "white") + geom_text(aes(label = percent(proportion, accuracy = 1)), size = 3) +
  scale_fill_gradient(low = "white", high = "#3A7D44", limits = c(0, 1)) +
  labs(title = "A. Baseline recovery", x = "Current predicted signal class", y = "Truth scenario", fill = "Proportion") +
  theme_fsf() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
p6b_summary <- aggregate(recovery_proportion ~ noise_sd, f6b, mean)
p6b <- ggplot(p6b_summary, aes(noise_sd, recovery_proportion)) +
  geom_ribbon(aes(ymin = 0, ymax = recovery_proportion), fill = "#9ECAE1", alpha = 0.5) +
  geom_line(linewidth = 0.9, color = "#24557A") + geom_point(size = 1.8, color = "#24557A") +
  scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
  scale_x_continuous(breaks = sort(unique(f6b$noise_sd))) +
  labs(title = "B. Overall recovery across noise", x = "Noise standard deviation",
       y = "Mean dominant-class recovery") + theme_fsf()
p6c <- ggplot(f6b, aes(noise_sd, recovery_proportion, color = truth_display, group = truth_display)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
  scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
  scale_x_continuous(breaks = sort(unique(f6b$noise_sd))) +
  labs(title = "C. Recovery by truth scenario", x = "Noise standard deviation", y = "Dominant-class recovery", color = "Truth scenario") +
  theme_fsf()
figure6 <- p6a + p6b + p6c + plot_layout(widths = c(1.05, 0.9, 1.15)) +
  plot_annotation(title = "Figure 6. Synthetic validation")
if (!figure5_only) save_figure(figure6, "Figure6", 15.0, 5.6)

# Supplementary Figure S1: frozen representative annotation status.
s1 <- read.delim(file.path(review_dir, "representative_feature_review.tsv"), check.names = FALSE, quote = "", stringsAsFactors = FALSE)
s1$condition <- factor(s1$condition, conditions)
s1_summary <- aggregate(feature_id ~ condition + annotation_status, s1, length)
names(s1_summary)[3] <- "n_features"
figureS1 <- ggplot(s1_summary, aes(condition, n_features, fill = annotation_status)) +
  geom_col(position = "fill", width = 0.72) +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Supplementary Figure S1. Representative annotation coverage",
       subtitle = "Class-size scatter not reproduced because class-size variables are unavailable",
       x = "Condition", y = "Representative-feature proportion", fill = "Annotation status") +
  theme_fsf()
if (!figure5_only) save_figure(figureS1, "FigureS1", 8.4, 5.2)

# Supplementary Figure S2: tau is x; frozen metric value is y.
s2 <- read_source("figureS2_tau_source.tsv")
s2 <- s2[s2$metric %in% c("mean_ssi", "mean_stability_deviation"), ]
s2$truth_display <- truth_label(s2$truth_scenario)
s2$metric_display <- c(mean_ssi = "Mean SSI", mean_stability_deviation = "Mean Stability Deviation")[s2$metric]
figureS2 <- ggplot(s2, aes(tau, value, color = truth_display, group = truth_display)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
  facet_wrap(~metric_display, nrow = 1) +
  scale_x_continuous(breaks = sort(unique(s2$tau))) +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.01)) +
  labs(title = "Supplementary Figure S2. Tau sensitivity", subtitle = "Metric values evaluated across the prespecified tau grid",
       x = expression(tau), y = "Metric value", color = "Truth scenario") +
  theme_fsf() + theme(legend.position = "bottom")
if (!figure5_only) save_figure(figureS2, "FigureS2", 10.2, 4.8)

# Supplementary Figure S3: historical stacked-composition presentation.
s3 <- read_source("figureS3_architecture_source.tsv")
s3$condition <- factor(s3$condition, conditions)
s3$stability_region <- factor(s3$stability_region, regions)
figureS3 <- ggplot(s3, aes(condition, region_proportion, fill = stability_region)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(region_proportion >= 0.06,
                               paste0(comma(region_count), "\n", percent(region_proportion, accuracy = 0.1)), "")),
            position = position_stack(vjust = 0.5), color = "white", size = 2.8, fontface = "bold") +
  scale_fill_manual(values = region_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Supplementary Figure S3. FSF region composition",
       x = "Condition", y = "Feature proportion", fill = "FSF region") +
  theme_fsf()
if (!figure5_only) save_figure(figureS3, "FigureS3", 8.8, 5.4)

# Supplementary Figure S4: discrete SSI-frequency presentation.
s4 <- read_source("figureS4_ssi_distribution_source.tsv")
s4$condition <- factor(s4$condition, conditions)
s4$stability_region <- factor(s4$stability_region, regions)
s4_freq <- aggregate(feature_id ~ condition + ssi + stability_region, s4, length)
names(s4_freq)[4] <- "n_features"
figureS4 <- ggplot(s4_freq, aes(ssi, condition, size = n_features, color = stability_region)) +
  geom_point(alpha = 0.72) +
  geom_vline(xintercept = c(0.50, 0.75, 0.90), linetype = "dashed", color = "#555555", linewidth = 0.4) +
  scale_color_manual(values = region_colors, drop = FALSE) +
  scale_size_area(max_size = 8, breaks = c(1, 10, 100, 1000, 10000)) +
  scale_x_continuous(limits = c(1 / 3, 1), breaks = c(1 / 3, 0.50, 0.75, 0.90, 1.00),
                     labels = c("0.33", "0.50", "0.75", "0.90", "1.00")) +
  labs(title = "Supplementary Figure S4. Condition-specific SSI frequencies",
       subtitle = "Point area represents the number of features at each discrete SSI value",
       x = "Signal Stability Index (SSI)", y = "Condition", color = "FSF region",
       size = "Features") +
  theme_fsf()
if (!figure5_only) save_figure(figureS4, "FigureS4", 8.4, 5.2)

manifest <- data.frame(
  figure_id = c(paste0("Figure", 1:7), paste0("FigureS", 1:4)),
  panel = "complete",
  source_authority = c(
    "figures/new/Figure_1.png (locked AI-only asset)",
    "source_data/figure2_simplex_source.tsv",
    "source_data/figure3_architecture_source.tsv",
    "source_data/figure4_ssi_ecdf_source.tsv",
    "source_data/figure5_go_source.tsv; source_data/figure5_kegg_source.tsv; figures/biological_theme_prototypes/Figure5C_DES_biological_themes_barplot_AUDIT.png",
    "source_data/figure6_baseline_source.tsv; source_data/figure6_noise_source.tsv",
    "figures/new/Figure_7.png (locked AI-only asset)",
    "review/representative_feature_review.tsv",
    "source_data/figureS2_tau_source.tsv",
    "source_data/figureS3_architecture_source.tsv",
    "source_data/figureS4_ssi_distribution_source.tsv"
  ),
  generation_script = c(
    "AI-only locked asset; copied without regeneration",
    rep("scripts/current/figures/build_current_manuscript_figures.R", 5),
    "AI-only locked asset; copied without regeneration",
    rep("scripts/current/figures/build_current_manuscript_figures.R", 4)
  ),
  output_png = paste0("results/current_fsf_v1/manuscript/figures/", c(paste0("Figure", 1:7), paste0("FigureS", 1:4)), ".png"),
  output_pdf = paste0("results/current_fsf_v1/manuscript/figures/", c(paste0("Figure", 1:7), paste0("FigureS", 1:4)), ".pdf"),
  numerical_status = c("locked AI-only asset", rep("frozen current authority", 5), "locked AI-only asset", rep("frozen current authority", 4)),
  interpretation_status = c("AI-only locked", rep("not applicable", 3), "author-reviewed DES panel", "not applicable", "AI-only locked", rep("not applicable", 4)),
  notes = c(
    "Authoritative AI-generated Figure 1 copied unchanged from figures/new/Figure_1.png",
    "Historical simplex geometry paired with current observed distributions; obsolete boundaries removed",
    "Historical composition-plus-heatmap design; locked architecture counts and current terminology",
    "Historical clean ECDF styling; SSI unchanged and current boundaries applied",
    "GO/KEGG evidence panels plus authoritative audited DES biological-theme panel",
    "Historical multi-panel validation narrative; obsolete truth/class display labels removed",
    "Authoritative AI-generated Figure 7 copied unchanged from figures/new/Figure_7.png",
    "Historical scatter unavailable because class-size variable is absent; annotation coverage retained without new selection",
    "Historical analytical design rejected; corrected x=tau and y=value presentation",
    "Historical stacked composition design with locked current counts",
    "Historical discrete distribution concept with current memberships and frequency encoding"
  ), stringsAsFactors = FALSE
)
if (!figure5_only) write.table(manifest, file.path(output_dir, "figure_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

audit <- c(
  "# Current FSF manuscript figure audit", "",
  "Figures 2-6 and S1-S4 were generated deterministically from frozen current authorities. Figures 1 and 7 are locked AI-only assets copied without regeneration. No scientific authority was regenerated.", "",
  "## Figure 1", "- Authority: `figures/new/Figure_1.png`.", "- Asset status: locked AI-only; copied unchanged into the final figure location.", "- Script regeneration: prohibited and not performed.", "- Validation: PASS.", "",
  "## Figure 2", "- Authority: `figure2_simplex_source.tsv`.", "- Historical design template: colored simplex geometry.", "- Obsolete content removed: former boundaries and labels.", "- Numerical status: current observations retain exact coordinates; coincident points are display-counted only.", "- Terminology: current region labels.", "- Validation: PASS.", "",
  "## Figure 3", "- Authority: `figure3_architecture_source.tsv`.", "- Historical design template: composition plus architecture heatmap.", "- Obsolete content removed: former memberships and labels.", "- Numerical status: locked counts preserved.", "- Interpretation: UV shown as Transitional-dominant.", "- Validation: PASS.", "",
  "## Figure 4", "- Authority: `figure4_ssi_ecdf_source.tsv`.", "- Historical design template: clean combined ECDF.", "- Obsolete content removed: former boundary ticks.", "- Numerical status: SSI values unchanged.", "- Terminology: current boundary labels.", "- Validation: PASS.", "",
  "## Figure 5", "- Authority: frozen GO/KEGG source tables plus the audited DES biological-theme panel.", "- Panel C asset: `biological_theme_prototypes/Figure5C_DES_biological_themes_barplot_AUDIT.png`.", "- Numerical status: significant term identifiers and adjusted p-values summarized without enrichment recalculation; audited DES content embedded unchanged.", "- Validation: PASS.", "",
  "## Figure 6", "- Authority: current baseline and noise synthetic source tables.", "- Historical design template: three-panel validation narrative.", "- Obsolete content removed: former truth/class display terminology.", "- Numerical status: frozen benchmark summaries.", "- Validation: PASS.", "",
  "## Figure 7", "- Authority: `figures/new/Figure_7.png`.", "- Asset status: locked AI-only; copied unchanged into the final figure location.", "- Script regeneration: prohibited and not performed.", "- Validation: PASS.", "",
  "## Supplementary Figure S1", "- Authority: frozen representative-feature review authority.", "- Historical design template assessed: class-size scatter.", "- Template simplification: class-size is unavailable, so only annotation coverage is displayed.", "- Selection remains unchanged.", "- Validation: PASS.", "",
  "## Supplementary Figure S2", "- Authority: `figureS2_tau_source.tsv`.", "- Historical analytical design rejected.", "- Historical issue corrected: the invalid tau-against-tau mapping is replaced by x = `tau` and y = `value`, faceted by `metric`.", "- Validation: PASS.", "",
  "## Supplementary Figure S3", "- Authority: `figureS3_architecture_source.tsv`.", "- Historical design template: stacked composition.", "- Obsolete content removed: former membership and region terminology.", "- Numerical status: locked counts and proportions preserved.", "- Validation: PASS.", "",
  "## Supplementary Figure S4", "- Authority: `figureS4_ssi_distribution_source.tsv`.", "- Historical design template: discrete condition-specific SSI landscape.", "- Obsolete content removed: former memberships and labels.", "- Numerical status: current discrete SSI frequencies; values unchanged.", "- Validation: PASS."
)
if (!figure5_only) writeLines(audit, file.path(output_dir, "figure_audit.md"), useBytes = TRUE)

cat(if (figure5_only) "Figure5 generated: PASS\n" else "current manuscript figures 2-6 and S1-S4 generated: PASS\n")
