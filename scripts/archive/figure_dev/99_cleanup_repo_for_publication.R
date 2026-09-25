#!/usr/bin/env Rscript

ROOT <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

main_fig <- file.path(ROOT, "results/archive/pre_repair_manuscript/figures/main")
supp_fig <- file.path(ROOT, "results/archive/pre_repair_manuscript/figures/supplementary")
archive <- file.path(ROOT, "results/archive/pre_repair_manuscript/figures/archive/removed_before_publication_freeze")
dir.create(archive, recursive = TRUE, showWarnings = FALSE)

keep_main <- c(
  "Figure1_Signal-aware stratification in the FSF framework.png",
  "Figure2_FSF_probability_simplex_regions.pdf",
  "Figure2_FSF_probability_simplex_regions.png",
  "Figure3A_gene_weighted_signal_class_composition.pdf",
  "Figure3A_gene_weighted_signal_class_composition.png",
  "Figure3B_condition_level_signal_architecture_heatmap.pdf",
  "Figure3B_condition_level_signal_architecture_heatmap.png",
  "Figure4A_condition_level_stability_architecture.pdf",
  "Figure4A_condition_level_stability_architecture.png",
  "Figure5A_gene_weighted_theme_distribution.pdf",
  "Figure5A_gene_weighted_theme_distribution.png",
  "Figure5B_gene_weighted_theme_heatmap.pdf",
  "Figure5B_gene_weighted_theme_heatmap.png",
  "Figure5C_gene_weighted_DES_UV_GAM_profiles.pdf",
  "Figure5C_gene_weighted_DES_UV_GAM_profiles.png",
  "Figure6A_synthetic_truth_confusion.pdf",
  "Figure6A_synthetic_truth_confusion.png",
  "Figure6B_noise_gradient_stability_metrics.pdf",
  "Figure6B_noise_gradient_stability_metrics.png",
  "Figure6C_noise_gradient_class_recovery.pdf",
  "Figure6C_noise_gradient_class_recovery.png",
  "Figure7_Final_Model.pdf",
  "Figure7_Final_Model.png"
)

keep_supp <- c(
  "Supplementary_FigureS1_annotation_coverage.pdf",
  "Supplementary_FigureS1_annotation_coverage.png",
  "Supplementary_FigureS2_tau_sensitivity.pdf",
  "Supplementary_FigureS2_tau_sensitivity.png",
  "FigureS3_gene_weighted_stability_level_composition.pdf",
  "FigureS3_gene_weighted_stability_level_composition.png",
  "FigureS4_ECDF_SSI_by_condition.pdf",
  "FigureS4_ECDF_SSI_by_condition.png"
)

clean_dir <- function(path, keep) {
  files <- list.files(path, full.names = TRUE)
  remove <- files[!basename(files) %in% keep]
  if (length(remove) > 0) {
    file.copy(remove, archive, overwrite = TRUE)
    file.remove(remove)
  }
}

clean_dir(main_fig, keep_main)
clean_dir(supp_fig, keep_supp)

cat("Main figures after cleanup:\n")
print(list.files(main_fig))

cat("\nSupplementary figures after cleanup:\n")
print(list.files(supp_fig))

cat("\nRemoved files archived in:\n")
cat(archive, "\n")
