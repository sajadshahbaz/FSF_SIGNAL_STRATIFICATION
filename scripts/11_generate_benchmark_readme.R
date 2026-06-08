#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 11: Generate Benchmark README
# ============================================================

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

OUT_README <- file.path(
  ROOT_DIR,
  "results",
  "manuscript",
  "FSF_v1_benchmark_README.md"
)

OUT_DOC_COPY <- file.path(
  ROOT_DIR,
  "docs",
  "FSF_v1_benchmark_README.md"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results",
  "logs",
  "11_generate_benchmark_readme.log"
)

dir.create(dirname(OUT_README), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_DOC_COPY), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 11: Generate Benchmark README\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

readme <- c(
"# FSF v1.0 Benchmark README",
"",
"## Project",
"",
"Feature Stability Framework (FSF): Signal Stratification in Heterogeneous Transcriptomic Data",
"",
"## Locked Project Root",
"",
"`/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION`",
"",
"## Core Method",
"",
"FSF v1.0 uses:",
"",
"- directional effect estimates",
"- signal state assignment: up, down, constant",
"- signal-state probabilities: P(up), P(down), P(const)",
"- Signal Stratification Index (SSI)",
"- Stability Deviation: 1 - SSI",
"- final FSF signal classes",
"",
"## Scope",
"",
"FSF v1.0 is a signal stratification framework.",
"",
"FSF v1.0 does not perform:",
"",
"- differential expression analysis",
"- feature selection",
"- machine learning classification",
"- transition-point detection",
"- switching analysis",
"- shifting analysis",
"- regime detection",
"",
"## Default Parameters",
"",
"- directional threshold tau = 0.5",
"- highly stable cutoff: SSI >= 0.90",
"- stable cutoff: 0.75 <= SSI < 0.90",
"- weakly stable / transitional cutoff: 0.60 <= SSI < 0.75",
"- instable cutoff: SSI < 0.60",
"",
"## Script Order",
"",
"Run scripts in this order:",
"",
"```bash",
"Rscript scripts/04_generate_synthetic_effect_estimates.R",
"Rscript scripts/01_assign_signal_states.R",
"Rscript scripts/02_compute_ssi_stability_deviation.R",
"Rscript scripts/03_fsf_diagnostic_plots.R",
"Rscript scripts/05_evaluate_synthetic_recovery.R",
"Rscript scripts/06_noise_gradient_benchmark.R",
"Rscript scripts/07_make_publication_benchmark_summary.R",
"Rscript scripts/08_tau_sensitivity_analysis.R",
"Rscript scripts/09_make_final_benchmark_report_table.R",
"Rscript scripts/10_generate_results_text.R",
"Rscript scripts/11_generate_benchmark_readme.R",
"```",
"",
"## Main Inputs",
"",
"### Synthetic effect estimates",
"",
"`data/processed/effect_estimates.tsv`",
"",
"Required columns:",
"",
"- feature_id",
"- perturbation_id",
"- effect_estimate",
"",
"## Main Outputs",
"",
"### Signal states",
"",
"`results/signal_states/feature_signal_states.tsv`",
"",
"### FSF metrics",
"",
"`results/fsf_metrics/fsf_ssi_stability_deviation.tsv`",
"",
"Contains:",
"",
"- feature_id",
"- P(up)",
"- P(down)",
"- P(const)",
"- dominant state",
"- SSI",
"- Stability Deviation",
"- stability level",
"- FSF signal class",
"",
"### Synthetic benchmark outputs",
"",
"`results/synthetic/synthetic_recovery_confusion_matrix.tsv`",
"",
"`results/synthetic/synthetic_recovery_summary.tsv`",
"",
"### Noise-gradient outputs",
"",
"`results/synthetic/noise_gradient/noise_gradient_summary.tsv`",
"",
"`results/synthetic/noise_gradient/noise_gradient_confusion.tsv`",
"",
"### Tau sensitivity outputs",
"",
"`results/synthetic/tau_sensitivity/tau_sensitivity_summary.tsv`",
"",
"`results/synthetic/tau_sensitivity/tau_sensitivity_class_counts.tsv`",
"",
"### Manuscript-ready outputs",
"",
"`results/tables/FSF_v1_final_benchmark_report.tsv`",
"",
"`results/tables/FSF_v1_key_results_for_manuscript.tsv`",
"",
"`results/manuscript/benchmark_results_text.txt`",
"",
"## Main Figures",
"",
"- `results/figures/synthetic_truth_vs_fsf_class_heatmap.png`",
"- `results/figures/noise_gradient_mean_ssi.png`",
"- `results/figures/noise_gradient_mean_sd.png`",
"- `results/figures/noise_gradient_class_recovery.png`",
"- `results/figures/tau_sensitivity_mean_ssi.png`",
"- `results/figures/tau_sensitivity_class_counts.png`",
"",
"## Key Benchmark Findings",
"",
"1. FSF recovered highly stable up, down, and constant synthetic classes under default settings.",
"2. Instable synthetic signals remained low-SSI and were consistently recovered as instable.",
"3. Weakly stable synthetic signals were mainly assigned to the weakly stable/transitional stratum.",
"4. Stable Up and Stable Down signals retained high SSI under increasing noise.",
"5. Stable Constant signals degraded under high noise because constant signals sit near the zero-effect boundary.",
"6. Tau sensitivity supported tau = 0.5 as the default threshold for FSF v1.0.",
"7. Stability Deviation behaved as a continuous measure of departure from ideal stability.",
"",
"## Manuscript Position",
"",
"These benchmarks support FSF v1.0 as a reproducible signal stratification framework.",
"",
"They do not support claims that FSF is a feature-selection algorithm, a differential expression method, or a machine-learning classifier.",
"",
"## Generated",
"",
paste0("Generated on: ", Sys.time())
)

writeLines(readme, OUT_README)
writeLines(readme, OUT_DOC_COPY)

cat("Saved README to:\n")
cat(OUT_README, "\n")
cat(OUT_DOC_COPY, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

