#!/usr/bin/env bash
set -euo pipefail

ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../.." && pwd )"
cd "$ROOT"

if ! command -v pdftoppm >/dev/null 2>&1; then
  echo "WARNING: pdftoppm not found."
  echo "PDF figures are preserved, but PNG conversion is skipped."
  echo "Install with: sudo apt install poppler-utils"
  exit 0
fi

FIG_DIRS=(
  "results/manuscript/figures/main"
  "results/manuscript/figures/supplementary"
)

for dir in "${FIG_DIRS[@]}"; do
  [ -d "$dir" ] || continue

  for pdf in "$dir"/*.pdf; do
    [ -e "$pdf" ] || continue

    base="${pdf%.pdf}"
    tmp="${base}_tmp"

    echo "Converting: $pdf"

    pdftoppm -png -r 300 "$pdf" "$tmp"

    if [ -f "${tmp}-1.png" ]; then
      mv "${tmp}-1.png" "${base}.png"
      echo "Saved: ${base}.png"
    else
      echo "WARNING: PNG not created for $pdf"
    fi
  done
done

echo "PDF-to-PNG conversion completed."
