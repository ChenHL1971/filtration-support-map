#!/usr/bin/env bash
# Workflow A — reproduce the published results from the files in this repository.
# Needs TWO downloads (third-party files that are not redistributed here; see raw/README.md):
#   raw/science.adp4753_tables_s1_to_s31.xlsx            Liu et al. 2025 supplementary tables (doi:10.1126/science.adp4753)
#   raw/benchmarks/panelapp_kidney_genes_2026-10-10.csv  PanelApp kidney panel genes (panel versions in raw/README.md)
# Everything else (scorecard copy, ClinGen snapshot, colocalization summaries, Titan maps, saved FinnGen results)
# is included. Workflow A verifies downstream analyses from archived intermediate results; it does not reconstruct
# the colocalization or the FinnGen analyses. Runtime about 10-15 minutes. Afterwards run: python3 scripts/verify_reproduction.py
set -euo pipefail
for f in raw/science.adp4753_tables_s1_to_s31.xlsx raw/benchmarks/panelapp_kidney_genes_2026-10-10.csv; do
  test -f "$f" || { echo "Missing $f (see raw/README.md)"; exit 1; }
done
mkdir -p manuscript figs/final

# Unpack the per-locus colocalization summaries (2,970 small files kept in one archive)
test -d out/aim1/sealed || tar -xzf out/coloc_per_locus.tar.gz

python3 scripts/build_gw_cache.py                 # out/gw_cache.pkl from data_derived/

# Aim 1: f_s from the saved per-locus colocalization summaries (out/aim1/sealed)
Rscript scripts/compute_fs.R
Rscript scripts/fs_fallback_p12.R                 # post hoc

# Titan consistency check
python3 scripts/titan_check.py
python3 scripts/titan_tests.py

# Aim 2
python3 scripts/aim2_models.py out/aim1           # primary analysis as implemented
python3 scripts/aim2_models_bun.py out/aim1
python3 scripts/aim2_secondary.py                 # learned models
python3 scripts/aim2_mgi.py                       # post hoc
python3 scripts/aim2_covered.py                   # post hoc
python3 scripts/aim2_p12.py                       # post hoc
python3 scripts/aim2_review3.py                   # post hoc
python3 scripts/aim2_innercv.py                   # post hoc
python3 scripts/aim2_within_locus.py              # post hoc, descriptive
Rscript scripts/review5_coloc_checks.R            # post hoc: direction-rule PP.H4, p12 range
Rscript scripts/review6_fallback.R                 # post hoc: fallback calibration checks
# review5_power.py (simulation) and review6_firth.py (needs FinnGen inputs) are run separately; see scripts/MANIFEST.md

# Data S1 and figures (Figure 5 uses the saved Aim 3 results in out/aim3_results.json)
python3 scripts/make_resource_table.py
python3 scripts/make_fig1.py
python3 scripts/make_fig_map.py
python3 scripts/make_figs.py
python3 scripts/make_graphical_abstract.py
