# Filtration-support map of kidney function genetic signals

Code and derived data for the manuscript

> **A Filtration-Support Map of Kidney Function Genetic Signals: Consistency and Incremental Value for Gene Prioritization.**
> Hui-Lan Chen, Mei-Jy Jeng, Yi-Hsuan Tang, Chien-Hung Lin, Chao-Ting Teng, Hsin-Hui Wang. Manuscript prepared for submission to *iScience*.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23294824.svg)](https://doi.org/10.5281/zenodo.23294824)

Archived versions on Zenodo: version 1.1.0, https://doi.org/10.5281/zenodo.23294824 (version 1.1.1 adds the post hoc analyses of review rounds 5–6)

The study re-analyzes the Kidney Disease Genetic Scorecard of Liu et al. (Science 2025). For each of 2,051 fine-mapped
eGFRcrea signals it estimates **filtration support (f_s)**, defined as directionally concordant colocalization between
creatinine-based and cystatin C-based eGFR (coloc-SuSiE; fallback coloc.abf within the credible set). It then tests
(i) consistency with published creatinine-specific and GFR-associated loci (Titan et al. 2026), (ii) whether adding
filtration support to the scorecard improves recovery of established kidney disease genes (prespecified primary
analysis), and (iii) enrichment of FinnGen R13 kidney outcome associations.

f_s is an operational colocalization measure. It is not a direct measure of glomerular filtration and not the
probability that a candidate gene is causal.

## The resource

`data_S1/Supplementary_Data_1_filtration_support_map.xlsx` (also as CSV) lists all 2,051 signals, including the 370
inconclusive ones, with f_s, classification, estimation method and linked candidate genes. Field definitions and
recommended use are in the README sheet and in `docs/data_dictionary.md`.

## Two ways to reproduce the results

| | Workflow A — reproduce published results | Workflow B — full reconstruction |
|---|---|---|
| Script | `run_reproduce_results.sh`, then `python3 scripts/verify_reproduction.py` | `run_full_reconstruction.sh` |
| Starts from | intermediate files saved in this repository | original public inputs |
| Downloads needed | **two files**: Liu et al. supplementary tables and the PanelApp kidney-panel gene list (see `raw/README.md`) | all inputs in `raw/README.md` (GWAS, LD reference, FinnGen, Titan tables, MGI, scorecard) |
| Runtime | about 10–15 minutes | several hours (colocalization at 742 loci) |
| Reruns | f_s classes; Titan tests; all Aim 2 analyses (primary, sensitivity, learned models, post hoc); Data S1; all figures | everything, including colocalization, signal mapping, LD crosswalk and FinnGen tests |
| Does not rerun | colocalization itself (uses the saved per-locus summaries in `out/coloc_per_locus.tar.gz`, unpacked automatically); FinnGen tests (Figure 5 uses saved `out/aim3_results*.json`) | — |

**Scope of Workflow A.** Workflow A verifies the reproducibility of downstream analyses from archived intermediate
results; it does not independently reconstruct the original colocalization or FinnGen analyses. Those require
Workflow B, which has not been re-run end to end for this release.

Workflow A was tested on clean copies of this repository: all 26 regenerated result files and figures were
byte-identical to the files behind the manuscript, except that in some runs one file showed the known difference
listed below.
`scripts/verify_reproduction.py` performs this comparison and exits with status 1 if any other file differs or is
missing, so it can be used in automated checks.

### Which analysis needs what

| Analysis | Script(s) | Workflow A | Extra inputs for Workflow B |
|---|---|---|---|
| Colocalization and signal mapping | aim1_coloc_full.R, collect_aim1_tech.R | saved results used | eGFRcrea and eGFRcys GWAS, 1000 Genomes EUR LD |
| f_s and classes (287 / 132 / 1,262 / 370) | compute_fs.R | ✓ | — |
| Titan check | titan_prepare.py, titan_map.R, titan_check.py, titan_tests.py | ✓ (tests) | Titan tables, LD reference |
| Primary Aim 2 (ΔAUPRC −0.0003) and prespecified sensitivity analyses | aim2_models.py, aim2_models_bun.py | ✓ | BUN GWAS for the BUN-based score |
| Learned models | aim2_secondary.py | ✓ | — |
| Post hoc: inner-CV (+0.0033), covered genes, p12, fold standardization, MGI, within-locus | aim2_innercv.py, aim2_covered.py, aim2_p12.py, aim2_review3.py, aim2_mgi.py, aim2_within_locus.py | ✓ | — |
| FinnGen (dialysis P = 0.46; adjusted OR 0.88) and covariate balance | aim3_ldscore.R, aim3_finngen.py, aim3_finngen_v2.py, aim3_balance.py | saved results only | FinnGen R13 summary statistics, eGFRcrea GWAS, LD reference, gene annotation |
| Data S1 and figures | make_resource_table.py, make_fig*.py, make_graphical_abstract.py | ✓ | — |

## Repository layout

```
README.md                    this file
LICENSE                      MIT (code); data licences below
CITATION.cff
requirements.txt             Python packages and versions used
environment.yml              conda environment (Python 3.13, R 4.3.3, coloc 6.0.3, susieR 0.12.35)
config/parameters.yaml       reference copy of parameters (NOT read by the code; see below)
run_reproduce_results.sh     Workflow A
run_full_reconstruction.sh   Workflow B
scripts/                     analysis scripts, unchanged from the versions that produced the results
  MANIFEST.md                what each script does, prespecified vs post hoc, inputs and outputs
  local_helpers/             helpers run on a local computer (LD extraction, FinnGen/MGI filtering, LDlink)
data_derived/                text copy of the Liu et al. genome-wide scorecard (rebuilds out/gw_cache.pkl)
raw/README.md                where to download each input
raw/benchmarks/              ClinGen snapshot used for the benchmark (included; PanelApp must be downloaded)
out/                         intermediate files and results; per-locus colocalization summaries are in
                             out/coloc_per_locus.tar.gz (unpacks to out/aim1/sealed, out/aim1/tech, out/aim1_bun/)
data_S1/                     the filtration-support map (Data S1)
docs/
  deviations.md              all amendments and post hoc analyses (Supplementary Table S1)
  data_dictionary.md         Data S1 field definitions
  expected_outputs.sha256    checksums of the files behind the manuscript (used by verify_reproduction.py)
  known_differences.sha256   exact checksum accepted for the one documented known difference
CHECKSUMS.sha256             SHA-256 of every released file
```

`raw/` and `out/` keep the directory names used inside the scripts, so every script runs unchanged from the
repository root.

## Parameters

Analysis parameters are set inside each script (for example `MAXV` and `FLANK` in `aim1_coloc_full.R`, the λ grid
and bootstrap count in `aim2_models.py`). The scripts were deliberately left unchanged so that they are identical to
the versions that produced the reported results. `config/parameters.yaml` collects the same values in one place for
reference only; the code does not read it, and editing it has no effect.

## Prespecified versus post hoc

The analysis plan was locked on 8 October 2026. Execution amendments made before any colocalization result,
benchmark performance or FinnGen association was examined, and all post hoc analyses, are listed in
`docs/deviations.md`. The primary analysis as implemented tuned λ on the training chromosomes without an inner
cross-validation loop (the plan said "nested"; some early docstrings still use that word); the inner-CV analysis is a
separate post hoc script. Post hoc analyses were added as new scripts and did not modify the primary scripts.
`scripts/MANIFEST.md` marks each script as prespecified, post hoc, feasibility, or added for this release.
The analysis plan was locked internally by the authors; it was not publicly registered before the analysis.
All departures from it are listed in `docs/deviations.md`.

## Known differences

- `out/aim1/aim2_mgi.json`: across reruns, the upper 95% bootstrap bound of the L1-model difference on the MGI
  benchmark was either the saved −0.0003 or −0.0002 (small run-to-run variation in the L1 fit). This value is not
  reported in the manuscript; all reported MGI values are reproduced exactly. The file is accepted if it is identical
  to the saved file or matches the exact SHA-256 of the −0.0002 version recorded in
  `docs/known_differences.sha256`; any other change is reported as unexpected.
- Figures are byte-identical only with the matplotlib version in `requirements.txt`.

## Files, sources and licences

| Files | Source | Licence / terms |
|---|---|---|
| `scripts/`, run scripts | this study | MIT |
| `data_S1/`, `out/` (derived results, except the benchmark lists below) | this study | CC BY 4.0 |
| `data_derived/liu2025_genomewide_scorecard.csv.gz` | Liu et al. 2025, Figshare 26299093 | CC BY 4.0 (redistributed unchanged in content) |
| `raw/benchmarks/clingen_gdv_2026-10-09.csv` | ClinGen gene–disease validity | CC0 |
| `out/benchmark_*.csv` (gene-symbol lists) | derived from ClinGen, PanelApp, Groopman et al. 2019 and MGI | gene symbols only; the source terms apply to the underlying curations |
| inputs listed in `raw/README.md`, including the PanelApp gene list and the LDlink r² crosswalk | original sources | not redistributed (download instructions, access dates and SHA-256 of the files used are given) |

Only public summary-level data were used; no individual-level data were accessed and no attempt was made to identify
participants. The repository is about 25 MB and about 140 files; the 2,970 per-locus colocalization summaries are kept in a
single archive (`out/coloc_per_locus.tar.gz`). No single file exceeds 10 MB.

## Contact

Lead contact: Hsin-Hui Wang, M.D., Ph.D. (hhwang@vghtpe.gov.tw), Taipei Veterans General Hospital and National Yang Ming Chiao Tung University, Taipei, Taiwan.
