#!/usr/bin/env bash
set -euo pipefail

ROOT="/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"
MAIN="$ROOT/results/manuscript/figures/main"
SUPP="$ROOT/results/manuscript/figures/supplementary"
ARCH="$ROOT/results/manuscript/figures/archive/final_cleanup_removed"

mkdir -p "$ARCH"

# Archive redundant Figure 3 files produced by 36K
for f in \
  "$MAIN/Figure3B_gene_weighted_stability_level_composition.pdf" \
  "$MAIN/Figure3B_gene_weighted_stability_level_composition.png" \
  "$MAIN/Figure3C_condition_level_signal_architecture_heatmap.pdf" \
  "$MAIN/Figure3C_condition_level_signal_architecture_heatmap.png" \
  "$MAIN/Figure3C_source_data.tsv"
do
  [ -f "$f" ] && mv "$f" "$ARCH/"
done

# Keep Figure3B architecture heatmap as the official Figure3B.
# Keep Figure3B_source_data.tsv only if it belongs to the heatmap.
# Current 36K overwrites Figure3B_source_data.tsv with stability-level source data,
# so archive it and keep Figure3A source only in main.
[ -f "$MAIN/Figure3B_source_data.tsv" ] && mv "$MAIN/Figure3B_source_data.tsv" "$ARCH/"

# Archive duplicate old supplementary folder
if [ -d "$SUPP/old" ]; then
  mkdir -p "$ARCH/supplementary_old"
  mv "$SUPP/old/"* "$ARCH/supplementary_old/" 2>/dev/null || true
  rmdir "$SUPP/old" 2>/dev/null || true
fi

# Final expected main files
echo "FINAL MAIN FIGURES:"
find "$MAIN" -maxdepth 1 -type f | sort

echo
echo "FINAL SUPPLEMENTARY FIGURES:"
find "$SUPP" -maxdepth 1 -type f | sort

echo
echo "ARCHIVED REMOVED FILES:"
find "$ARCH" -type f | sort
