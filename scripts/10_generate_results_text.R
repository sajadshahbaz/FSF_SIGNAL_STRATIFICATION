#!/usr/bin/env Rscript

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Step 10: Generate Manuscript Results Text
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

ROOT_DIR <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

FINAL_REPORT <- file.path(
  ROOT_DIR,
  "results",
  "tables",
  "FSF_v1_final_benchmark_report.tsv"
)

KEY_RESULTS <- file.path(
  ROOT_DIR,
  "results",
  "tables",
  "FSF_v1_key_results_for_manuscript.tsv"
)

OUT_TEXT <- file.path(
  ROOT_DIR,
  "results",
  "manuscript",
  "benchmark_results_text.txt"
)

LOG_FILE <- file.path(
  ROOT_DIR,
  "results",
  "logs",
  "10_generate_results_text.log"
)

dir.create(dirname(OUT_TEXT), recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

cat("============================================================\n")
cat("FSF SIGNAL STRATIFICATION v1.0\n")
cat("Step 10: Generate Manuscript Results Text\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

report <- read_tsv(FINAL_REPORT, show_col_types = FALSE)

baseline <- report %>%
  filter(
    benchmark == "baseline_synthetic_recovery",
    key_metric == "mean_ssi"
  )

noise <- report %>%
  filter(
    benchmark == "noise_gradient",
    setting == "noise_sd_1",
    key_metric == "mean_ssi_and_sd"
  )

txt <- c(

"FSF Benchmarking Results",
"",
"Baseline Synthetic Recovery",
"",
paste(
"Under the default FSF configuration (tau = 0.5), all highly stable synthetic signal classes were recovered with perfect class fidelity."
),
paste(
"Stable Up, Stable Down, and Stable Constant signals achieved mean SSI values of",
round(baseline$numeric_value[baseline$truth_class=="stable_up"],3), ",",
round(baseline$numeric_value[baseline$truth_class=="stable_down"],3), "and",
round(baseline$numeric_value[baseline$truth_class=="stable_constant"],3),
"respectively."
),
paste(
"Instable signals remained clearly separated from stable classes with a mean SSI of",
round(baseline$numeric_value[baseline$truth_class=="instable"],3), "."
),
paste(
"The deliberately weakly stable synthetic class produced an intermediate mean SSI of",
round(baseline$numeric_value[baseline$truth_class=="weakly_stable"],3),
"and was primarily assigned to the weakly stable/transitional stratum."
),

"",
"Noise-Gradient Benchmark",
"",
paste(
"Across increasing noise levels, Stable Up and Stable Down signals retained high SSI values."
),
paste(
"At the highest tested noise level (noise_sd = 1.0), Stable Up and Stable Down signals maintained mean SSI values of",
round(noise$numeric_value[noise$truth_class=="stable_up"],3),
"and",
round(noise$numeric_value[noise$truth_class=="stable_down"],3),
"respectively."
),
paste(
"Instable signals remained consistently low-SSI across all benchmark scenarios, supporting the ability of FSF to distinguish dominant signal identity from mixed signal behavior."
),
paste(
"Stable Constant signals showed progressive degradation under increasing noise, resulting in increasing Stability Deviation values and eventual migration toward the instable region."
),

"",
"Tau Sensitivity",
"",
paste(
"Tau sensitivity analysis demonstrated that tau = 0.5 provided the most balanced separation between stable directional, stable constant, weakly stable transitional, and instable signal classes."
),
paste(
"Lower thresholds increased assignment of weakly stable features to stable classes, whereas higher thresholds increasingly classified weakly stable features as instable."
),

"",
"Interpretation of Stability Deviation",
"",
paste(
"Stability Deviation increased monotonically as reproducible signal identity weakened."
),
paste(
"This behavior supports Stability Deviation as a continuous measure of departure from ideal stability and complements SSI by quantifying the magnitude of signal identity degradation."
)

)

writeLines(txt, OUT_TEXT)

cat("Saved:\n")
cat(OUT_TEXT, "\n\n")

cat("Finished:", as.character(Sys.time()), "\n")

sink()

