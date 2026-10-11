#!/usr/bin/env bash
# Workflow B — full reconstruction from the original public inputs (see raw/README.md for downloads).
# Steps marked [long] take hours. For a quick check of the published results use run_reproduce_results.sh (Workflow A).
set -euo pipefail
CORES=${CORES:-4}

# 0. Benchmarks and scorecard cache
python3 scripts/build_benchmarks.py out
python3 scripts/audit_gw.py raw/gw_scorecard.xlsx raw/science.adp4753_tables_s1_to_s31.xlsx out   # writes out/gw_cache.pkl from the Excel file
# (alternative without the Excel file: python3 scripts/build_gw_cache.py)

# 1. Colocalization (Aim 1)
# (re-creates out/aim1/sealed, out/aim1/tech and out/aim1_bun; the archived copies are in out/coloc_per_locus.tar.gz)
Rscript scripts/aim1_coloc_full.R raw/gwas raw/ldfull/full_ld_ref raw/science.adp4753_tables_s1_to_s31.xlsx out/aim1 "$CORES"   # [long]
Rscript scripts/collect_aim1_tech.R out/aim1
Rscript scripts/compute_fs.R
Rscript scripts/fs_fallback_p12.R                                   # post hoc
Rscript scripts/aim1_coloc_bun.R raw/gwas raw/ldfull/full_ld_ref raw/science.adp4753_tables_s1_to_s31.xlsx out/aim1_bun "$CORES"   # [long]
Rscript scripts/collect_aim1_tech.R out/aim1_bun
Rscript scripts/compute_fs_bun.R

# 2. Titan consistency check
python3 scripts/titan_prepare.py out
Rscript scripts/titan_map.R
python3 scripts/titan_check.py out/aim1/titan_check.csv
python3 scripts/titan_tests.py

# 3. Gene prioritization (Aim 2)
python3 scripts/aim2_models.py out/aim1                             # primary analysis as implemented
python3 scripts/aim2_models_bun.py out/aim1
python3 scripts/aim2_secondary.py
python3 scripts/aim2_mgi.py                                         # post hoc
python3 scripts/aim2_covered.py                                     # post hoc
python3 scripts/aim2_p12.py                                         # post hoc
python3 scripts/aim2_review3.py                                     # post hoc
python3 scripts/aim2_innercv.py                                     # post hoc
python3 scripts/aim2_within_locus.py                                # post hoc, descriptive

# 4. FinnGen (Aim 3)
Rscript scripts/aim3_ldscore.R
python3 scripts/aim3_finngen.py
python3 scripts/aim3_finngen_v2.py
DIRECT_ONLY=1 python3 scripts/aim3_finngen_v2.py
python3 scripts/aim3_balance.py                                     # post hoc

# 5. Resource table and figures (scripts write to manuscript/ and figs/)
mkdir -p manuscript figs/final
python3 scripts/make_resource_table.py
python3 scripts/make_fig1.py; python3 scripts/make_fig_map.py; python3 scripts/make_figs.py; python3 scripts/make_graphical_abstract.py
