# Protocol amendments and deviations

Exported from Supplementary Table S1 of the manuscript. The analysis plan (Plan-5) was locked internally on 8 October 2026.

**Supplementary Table S1. Deviations from the prespecified analysis plan**

| # | Deviation | Phase | Effect |
|---|---|---|---|
| 1 | Pilot allele-alignment error (strand flips tested before allele swaps) corrected; UK Biobank LD briefly considered as primary, 1000 Genomes EUR kept as primary | Amendment: pilot, after lock and before unsealing | None on full analysis; UK Biobank LD kept as pilot sensitivity check |
| 2 | coloc.bf_bf rule added for loci where eGFRcys had no SuSiE credible set | Amendment: pilot, after lock and before unsealing | Defines f_s for these loci |
| 3 | Variant cap per locus lowered from 8,000 to 4,000 because of memory limits | Amendment: full run, before unsealing | 9 loci exceeded the cap and were filtered by MAF; 8 of these also needed the top-\|z\| cap |
| 4 | One line of colocalization output was printed to a run log during the first run; log moved to the sealed folder and printing disabled | Amendment: full run, before unsealing | No rule was changed after this exposure |
| 5 | Titan creatinine-specific loci: 54 regions in Titan Supplementary Table S1 versus 52 loci in the text and Supplementary Table S3; lead-variant mapping used all 54 regions, and the region-based sensitivity analysis is reported for the 52-locus list and for all 54 | Amendment: before unsealing | The two extra regions map to no signal by lead variant (n = 11); region-based: 71 vs. 946 (52 loci) and 73 vs. 946 signals (54 regions) |
| 6 | Gene-level Score_A defined as the maximum over variant–gene pairs (not specified in the plan) | Amendment: before unsealing | Defines the Model A comparator |
| 7 | Aim 2 cross-validation: the plan described nested cross-validation; λ, the only tuned parameter, was tuned on training-chromosome AUPRC without an inner cross-validation loop | Amendment: implementation, before unsealing | Held-out chromosome never used for tuning; Methods describe the procedure as implemented |
| 8 | Operational rules for Titan-to-signal mapping and for the FinnGen matching covariates written out in detail | Amendment: before Titan classes and FinnGen associations were opened | Rules fixed before the corresponding data were examined |
| 9 | MGI renal-function benchmark added | Addition after unsealing, after pre-submission review | Additional benchmark; reported as secondary |
| 10 | Titan comparison with one value per locus | Post hoc, after pre-submission review | P = 4.2 × 10⁻⁶; conclusion unchanged |
| 11 | Aim 2 restricted to genes with a filtration score | Post hoc, after pre-submission review | +0.0053 (−0.007 to 0.018); conclusion unchanged |
| 12 | Covariate balance diagnostics for the FinnGen permutation test | Post hoc, after pre-submission review | Reported in Supplementary Table S2 |
| 13 | Fallback coloc.abf had used default p₁₂ = 10⁻⁵; fallback PP.H4 recomputed with p₁₂ = 5 × 10⁻⁶ | Post hoc, after pre-submission technical audit | 287 → 278 filtration-consistent signals; Titan P = 3.5 × 10⁻⁶; Aim 2 −0.0003 (−0.007 to 0.005); conclusions unchanged |
| 14 | Fold-specific standardization of scores; f_s without the direction rule | Post hoc, after third pre-submission review | −0.0004 (−0.008 to 0.005) and −0.0004 (−0.006 to 0.005); Titan P = 1.7 × 10⁻⁶ without direction rule; conclusions unchanged |
| 15 | λ selected by chromosome-grouped inner cross-validation; fallback diagnostics (PP.H4 by method, class changes under alternative prior, fallback signals in multi-signal loci) | Post hoc, after fourth pre-submission review | +0.0033 (0.0004 to 0.009), below 0.03; 9 consistent→intermediate and 19 intermediate→low; 842/871 fallback signals in multi-signal loci; conclusions unchanged |
| 16 | Genome-wide map figure and illustrative loci (complete Köttgen 2010 list of seven suspected creatinine loci plus UMOD and SHROOM3) | Post hoc, descriptive | 3 low support, 1 intermediate, 3 filtration-consistent among Köttgen loci; no statistical test |
| 17 | Within-locus variation in gene-level filtration support (S_g) among candidate genes | Post hoc, descriptive, after fifth pre-submission review | 183/456 loci with identical S_g; 27% of S_g variance within loci (Score_A 71%); expected benchmark top-1 at 81 loci 12.4 (Score_A 24.5; random 12.8); no model refitted |
| 18 | Design-matched power simulation for the FinnGen dialysis comparison; PP.H4 of direction-rule signals; class counts with p₁₂ from 10⁻⁶ to 10⁻⁴ for all signals | Post hoc, after fifth pre-submission review | Power 14%, 37%, 56%, 77%, 87% at OR 1.5–3.5; 187/192 direction-rule signals had PP.H4 < 0.2; 193–410 filtration-consistent; conclusions unchanged |
| 19 | Firth-penalized dialysis model and Newcombe interval for the difference in hit proportions; class distribution by estimation method and credible-set size; p₁₂ varied for fallback signals only | Post hoc, after sixth pre-submission review | Firth OR 0.89 (0.30–2.39); difference −0.04 percentage points (−3.3 to 4.6); 13.1% of fallback vs. 21.4% of directly matched signals filtration-consistent; 259–331 filtration-consistent; conclusions unchanged |
