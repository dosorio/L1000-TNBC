#!/usr/bin/env bash
# Runs the revision-3 analyses in order. Requires R >= 4.3 with ccdata, preprocessCore, fgsea, readxl,
# ggplot2, ggrepel and patchwork, and metaLINCS 0.9.0 (R CMD INSTALL metaLINCS after 00_download.sh).
set -euo pipefail
cd "$(dirname "$0")"
bash 00_download.sh
mkdir -p results figures
for s in 01_lincs_metadata 02_signatures 03_benchmark 04_stability 05_public_sensitivity 06_viability 07_threshold 08_deg_direction 09_enrichment 10_figures 11_public_sensitivity_figure; do
  echo "== $s"; Rscript "$s.R" > "results/$s.log" 2>&1
done
