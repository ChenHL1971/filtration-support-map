# Input data

The analysis uses only public summary-level data. Workflow A needs the Liu et al. supplementary tables and the
PanelApp gene list. Apart from the ClinGen snapshot in `raw/benchmarks/` (CC0, included), third-party input files are
not redistributed; download them from the original sources and place them at the paths shown.
The GWAS files used for the reported analysis had been filtered on a local computer to the 742 merged-locus windows
(`out/kidney_loci_windows_hg19.txt`); the full summary-statistics files can be used directly because the scripts
subset each locus themselves.

| Path | Content | Source |
|---|---|---|
| raw/science.adp4753_tables_s1_to_s31.xlsx | Liu et al. 2025 supplementary tables (Table S8: 2,051 fine-mapped signals) | Liu et al., Science 2025, supplementary materials |
| raw/gw_scorecard.xlsx | Genome-wide Kidney Disease Genetic Scorecard (all variant–gene pairs, 32 features); optional, a text copy is in data_derived/ | Figshare 26299093 (CC BY 4.0) |
| raw/gwas/liu_crea.gz | eGFRcrea GWAS, European ancestry (Liu et al.) | Figshare 26299093 |
| raw/gwas/cys.gz, raw/gwas/bun.gz | eGFRcys and BUN GWAS (Stanzick et al. 2021) | https://ckdgen.imbi.uni-freiburg.de/datasets |
| raw/ldfull/full_ld_ref.{bed,bim,fam} | 1000 Genomes Phase 3 EUR, restricted to the 742 loci (`scripts/local_helpers/extract_ld_ref_full.R`) | MAGMA g1000_eur, https://ctg.cncr.nl/software/magma |
| raw/ukbld/ | UK Biobank LD for the 10 pilot loci (`scripts/local_helpers/extract_ukb_ld.R`) | https://alkesgroup.broadinstitute.org/UKBB_LD/ |
| raw/finngen/sumstats/fg_{endpoint}_kidneyloci_tsv.gz | FinnGen R13 summary statistics filtered to kidney-locus rsIDs (`scripts/local_helpers/filter_finngen.R`) | https://www.finngen.fi/en/access_results |
| raw/titan/mmc3.xlsx | Titan et al. 2026 Supplementary Tables (ST1, ST3) | Titan et al., Kidney Int 2026 |
| raw/benchmarks/clingen_gdv_2026-10-09.csv (included) | ClinGen gene–disease validity curations (downloaded 9 Oct 2026) | https://search.clinicalgenome.org/kb/gene-validity |
| raw/benchmarks/panelapp_kidney_genes_2026-10-10.csv | PanelApp kidney panel genes (downloaded 10 Oct 2026; panels and versions below) | https://panelapp.genomicsengland.co.uk |
| raw/annot/gene_annot_jun10.txt | GRCh37 gene annotation | PoPS, https://github.com/FinucaneLab/pops |
| raw/mgi/ | MGI phenotype reports (`scripts/local_helpers/filter_mgi_genepheno.R`) | https://www.informatics.jax.org |

| raw/ld/ld_crosswalk_results.csv | r² (1000 Genomes EUR, GRCh37) between 429 scorecard lead variants and nearby credible-set variants; used only in the feasibility audit (`audit_ld.py`) | LDlink, regenerated with `scripts/local_helpers/ld_crosswalk.R` and a personal LDlink token (no token is stored here) |

## PanelApp file

One row per gene per panel, with columns `panel_id, panel, version, gene, confidence` (confidence 3 = green,
2 = amber). Panels 145, 1317, 154 and 521 (adrenal and renal cancer) are downloaded but excluded by the scripts.

| panel_id | Panel | Version |
|---|---|---|
| 234 | CAKUT | 1.182 |
| 1317 | Congenital adrenal hyperplasia | 1.600 |
| 145 | Congenital adrenal hypoplasia | 5.600 |
| 283 | Cystic kidney disease | 9.300 |
| 487 | Cystic renal disease | 13.700 |
| 99 | Haematuria | 2.200 |
| 521 | Inherited renal cancer | 1.310 |
| 83 | Membranoproliferative glomerulonephritis including C3 glomerulopathy | 3.150 |
| 18 | Monogenic nephrogenic diabetes insipidus | 1.120 |
| 149 | Nephrocalcinosis or nephrolithiasis | 6.300 |
| 106 | Proteinuric renal disease | 6.220 |
| 154 | Renal cancer pertinent cancer susceptibility | 1.300 |
| 725 | Renal ciliopathies | 5.200 |
| 292 | Renal tubulopathies | 6.100 |
| 548 | Tubulointerstitial kidney disease | 3.350 |
| 156 | Unexplained kidney failure in young people | 1.125 |
| 678 | Unexplained young onset end-stage renal disease | 13.490 |
| 1537 | Unexplained young onset end-stage renal disease - additional genes | 2.700 |

Gene lists from these panel versions can be retrieved from the PanelApp website or API
(`/api/v1/panels/<panel_id>/?version=<version>`). The resulting green-gene list used in the analysis is in
`out/benchmark_panelapp_green.csv`.

## SHA-256 of third-party files used

| File | SHA-256 |
|---|---|
| raw/benchmarks/clingen_gdv_2026-10-09.csv (included) | 4df3aac3871b750c243a27f09a417d0704aa26d924d34b330cab6680686d2f2b |
| raw/benchmarks/panelapp_kidney_genes_2026-10-10.csv | e0154f3c551e1f4c9512a42ee37aee631415600a93187f7200d4a87595da5b0b |
| raw/ld/ld_crosswalk_results.csv | 50b3acaa8daebb73522e71f242395293210e1cdbf555f0cc4d55732ac89160d3 |

Benchmark gene lists derived from these sources are included in `out/` (`benchmark_*.csv`). The Groopman et al. 2019
625-gene list (`out/benchmark_groopman625.csv`) was transcribed from the supplementary appendix of Groopman et al.,
N Engl J Med 2019, and the MGI renal-function list (`out/benchmark_mgi_renal_function.csv`) from the MGI filtering helper.
