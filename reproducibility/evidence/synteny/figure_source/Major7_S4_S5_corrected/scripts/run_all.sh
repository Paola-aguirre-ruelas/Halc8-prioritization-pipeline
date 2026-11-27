#!/usr/bin/env bash
# Regenerates the whole package from the read-only archive in ../synteny_figures_final.
# Requires python3 with matplotlib (figures) and biopython (verification only). Clinker is not run.
set -euo pipefail
cd "$(dirname "$0")"
python3 01_audit_sources.py > ../audit/source_audit_summary.txt
python3 02_build_annotations.py
python3 03_plot_figures.py
python3 04_verify.py
