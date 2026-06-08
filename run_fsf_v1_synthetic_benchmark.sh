#!/usr/bin/env bash

# ============================================================
# FSF SIGNAL STRATIFICATION v1.0
# Master Runner: Synthetic Benchmark
# ============================================================

set -euo pipefail

ROOT_DIR="/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"
cd "$ROOT_DIR"

MASTER_LOG="$ROOT_DIR/results/logs/run_fsf_v1_synthetic_benchmark.log"

mkdir -p "$ROOT_DIR/results/logs"

{
echo "============================================================"
echo "FSF SIGNAL STRATIFICATION v1.0"
echo "Master Runner: Synthetic Benchmark"
echo "Started: $(date)"
echo "Root: $ROOT_DIR"
echo "============================================================"
echo ""

echo "[Step 04] Generate synthetic effect estimates"
Rscript scripts/04_generate_synthetic_effect_estimates.R

echo ""
echo "[Step 01] Assign signal states"
Rscript scripts/01_assign_signal_states.R

echo ""
echo "[Step 02] Compute SSI + Stability Deviation"
Rscript scripts/02_compute_ssi_stability_deviation.R

echo ""
echo "[Step 03] Generate diagnostic plots"
Rscript scripts/03_fsf_diagnostic_plots.R

echo ""
echo "[Step 05] Evaluate synthetic recovery"
Rscript scripts/05_evaluate_synthetic_recovery.R

echo ""
echo "[Step 06] Noise-gradient benchmark"
Rscript scripts/06_noise_gradient_benchmark.R

echo ""
echo "[Step 07] Publication benchmark summary"
Rscript scripts/07_make_publication_benchmark_summary.R

echo ""
echo "[Step 08] Tau sensitivity analysis"
Rscript scripts/08_tau_sensitivity_analysis.R

echo ""
echo "[Step 09] Final benchmark report table"
Rscript scripts/09_make_final_benchmark_report_table.R

echo ""
echo "[Step 10] Generate manuscript Results text"
Rscript scripts/10_generate_results_text.R

echo ""
echo "[Step 11] Generate benchmark README"
Rscript scripts/11_generate_benchmark_readme.R

echo ""
echo "============================================================"
echo "FSF v1.0 synthetic benchmark completed successfully."
echo "Finished: $(date)"
echo "============================================================"

} 2>&1 | tee "$MASTER_LOG"

