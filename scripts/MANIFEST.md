# Script manifest

All scripts are run from the repository root. Scripts are released **unchanged** from the versions that produced the
reported results (SHA-256 checksums in `CHECKSUMS.sha256`), so that later edits cannot obscure the order in which
analyses were done. Some early docstrings use the phrase "nested cross-validation"; the procedure actually implemented
(λ tuned on the training chromosomes only) is described in the manuscript Methods and in Supplementary Table S1, item 7.

Status: **P** = prespecified (locked plan, 8 Oct 2026, or execution amendment made before unsealing);
**PH** = post hoc (logged in `docs/deviations.md`); **R** = added 11 Oct 2026 for this release to reproduce a step that was
originally run interactively (each verified to reproduce the original file exactly); **F** = feasibility audit (no results examined).

## 0. Inputs and feasibility audit
| Script | Status | Purpose | Main inputs → outputs |
|---|---|---|---|
| audit_supp.py | F | Structure of Liu et al. supplementary tables | raw/science…xlsx → out/feature_inventory.csv |
| audit_gw.py | F | Genome-wide scorecard structure; builds cache | raw/gw_scorecard.xlsx → out/gw_cache.pkl, gw_*.csv |
| audit_gwas.py | F | GWAS harmonization and credible-set coverage | raw/gwas → out/audit_gwas_summary.json |
| audit_ld.py | F | Classify LD crosswalk | raw/ld → out/crosswalk_*.csv |
| audit_finngen.py | F | FinnGen file structure only (no association values) | raw/finngen → out/audit_finngen_summary.json |
| local_helpers/*.R | F | Run on a local computer: LD reference extraction, UK Biobank LD, FinnGen and MGI filtering, LDlink crosswalk | see headers |
| build_benchmarks.py | R | ClinGen + PanelApp benchmark lists | raw/benchmarks → out/benchmark_*.csv |
| build_gw_cache.py | R | Rebuilds the scorecard cache from the redistributed text copy (identical to the cache used) | data_derived/ → out/gw_cache.pkl |

## 1. Colocalization (Aim 1)
| Script | Status | Purpose |
|---|---|---|
| pilot_select.py, coloc_pilot*.R, map_susie_s8.R | F | 10-locus technical pilot (PP.H4 sealed, not examined) |
| aim1_coloc_full.R | P | coloc-SuSiE / coloc.bf_bf / fallback coloc.abf at all 742 loci → out/aim1/sealed, out/aim1/tech (released in out/coloc_per_locus.tar.gz) |
| collect_aim1_tech.R | R | Collects per-locus files → out/aim1/aim1_signal_mapping.csv, aim1_technical_log.csv |
| compute_fs.R | P | Applies direction and coverage rules → out/aim1/f_s_signals.csv |
| aim1_coloc_bun.R, compute_fs_bun.R | P | BUN-based score (sensitivity) → out/aim1_bun |
| fs_fallback_p12.R | PH | Fallback PP.H4 with p12 = 5e-6 → out/aim1/f_s_signals_p12.csv |

## 2. Titan consistency check
| Script | Status | Purpose |
|---|---|---|
| titan_prepare.py | R | Titan ST1 export; region-based Titan→signal map |
| titan_map.R | P | Lead-variant Titan→signal map → out/titan_to_s8_map.csv |
| titan_check.py | R | Joins map with f_s → out/aim1/titan_check.csv |
| titan_tests.py | P/PH | Lead-variant (P), region-based and locus-level (PH) tests → out/titan_tests.json |

## 3. Gene prioritization (Aim 2)
| Script | Status | Purpose |
|---|---|---|
| aim2_models.py | P | **Primary analysis as implemented**, Models A–D, prespecified sensitivity and alternative benchmarks → out/aim1/aim2_results.json, gene_scores.csv |
| aim2_models_bun.py | P | BUN-based sensitivity → aim2_results_bun.json |
| aim2_secondary.py | P | Learned models (L1 logistic, XGBoost) → aim2_secondary.json |
| aim2_mgi.py | PH | MGI renal-function benchmark → aim2_mgi.json |
| aim2_covered.py | PH | Genes with a filtration score only → aim2_covered.json |
| aim2_p12.py | PH | Primary analysis with fallback p12 = 5e-6 → aim2_p12.json |
| aim2_review3.py | PH | Fold-specific standardization; no direction rule → aim2_review3.json |
| aim2_innercv.py | PH | **Inner chromosome-grouped CV for λ** → aim2_innercv.json |
| aim2_within_locus.py | PH | **Within-locus variation of S_g (descriptive)** → aim2_within_locus.json |
| review5_coloc_checks.R | PH | PP.H4 of direction-rule signals; class counts under p12 1e-6 to 1e-4 (all signals) → out/review5_coloc_checks.json |
| review6_fallback.R | PH | Classes by estimation method and credible-set size; p12 varied for fallback signals only → out/review6_fallback.json |

## 4. FinnGen (Aim 3)
| Script | Status | Purpose |
|---|---|---|
| aim3_ldscore.R | P | Signal leads and LD scores → out/aim3_signal_leads.csv |
| aim3_finngen.py | P | Propensity-stratified permutation tests → out/aim3_results.json |
| aim3_finngen_v2.py | P | Adds covariate-adjusted logistic OR; `DIRECT_ONLY=1` for directly matched signals → aim3_results_v2.json, aim3_results_direct.json |
| aim3_balance.py | PH | Covariate balance within strata → aim3_balance_dialysis.csv, aim3_strata_dialysis.csv |
| review5_power.py | PH | Design-matched power simulation for the dialysis comparison → out/review5_power.json |
| review6_firth.py | PH | Firth-penalized dialysis model and Newcombe interval (needs Workflow B inputs) → out/review6_firth.json |

## 5. Figures, tables and manuscript
| Script | Purpose |
|---|---|
| verify_reproduction.py | Compares regenerated outputs with the files behind the manuscript (Workflow A) |
| make_resource_table.py | Data S1 (CSV; the XLSX adds the README sheet) |
| make_fig1.py, make_fig_map.py, make_figs.py, make_graphical_abstract.py, make_ga_pptx.py | Figures and graphical abstract |
| build_docx.py, build_iscience.py | Manuscript assembly from the authors' working document (not needed to reproduce results) |
