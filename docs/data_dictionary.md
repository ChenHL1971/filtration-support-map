# Data dictionary: Data S1 (filtration-support map)

Same content as the README sheet of Supplementary_Data_1_filtration_support_map.xlsx.

| Item | Description |
|---|---|
| Interpretation and recommended use | This dataset provides signal-level filtration-support annotations for 2,051 fine-mapped eGFRcrea genetic signals, of which 1,681 have classifiable scores. f_s represents direction-filtered colocalization evidence for a shared causal variant between eGFRcrea and eGFRcys. It is not a direct measurement of glomerular filtration and not a probability that a candidate gene is causal. |
| Low support | Low filtration support should not be interpreted as proof of creatinine-specific biology, because it may also reflect limited statistical power for eGFRcys, incomplete variant coverage, fallback estimation or uncertainty in colocalization. |
| Inconclusive records | Inconclusive signals (n = 370) are retained so that users can see that a signal exists but could not be scored; inconclusive_reason gives the cause. They should be treated as missing information, not as low support. |
| Linked candidate genes | Linked genes and weights are taken from the original Kidney Disease Genetic Scorecard (Liu et al., Science 2025) and are not independently validated causal genes. This resource complements, rather than replaces, existing variant-to-gene evidence. |
| How to query | Search by lead_variant_rsid, by chromosome and position_GRCh37, by locus_id, or by gene symbol in linked_candidate_genes_link_weight. A locus can contain several signals with different classifications. |
| locus_id | Merged eGFRcrea locus from Liu et al. (2025) Table S8 |
| signal_id | Fine-mapped CARMA credible set (locus\|credible set) |
| lead_variant_rsid | Variant with the highest CARMA posterior inclusion probability in the credible set |
| chromosome, position_GRCh37 | Lead-variant coordinates (GRCh37/hg19) |
| credible_set_size | Number of variants in the Liu et al. credible set |
| credible_set_coverage | Fraction of credible-set variants available in the eGFRcys GWAS and LD reference (>= 0.80 required) |
| estimation_method | Procedure attempted for each signal. Direct match (coloc-SuSiE; 831 signals attempted, 810 scored); fallback coloc.abf restricted to the credible set (879 attempted, 871 scored); not estimated (341, coverage < 0.80). 'Attempted' does not mean 'scored': see classification and inconclusive_reason. Direct-match and fallback estimates use different procedures and priors and should be distinguished when interpreting individual signals. |
| PP_H4 | Posterior probability of a shared causal variant for eGFRcrea and eGFRcys |
| direction_concordant | Yes / No (f_s set to 0) for scored signals; 'Yes/No (signal not scored)' for signals failing the coverage rule after colocalization; 'Not determined' when the lead variant was not among analyzed variants; 'Not estimated' when colocalization was not run (coverage < 0.80) |
| filtration_support_f_s | Primary filtration-support score. Equals PP_H4 when lead-variant effects on eGFRcrea and eGFRcys are concordant, and is set to 0 by study rule when they are discordant, even if PP_H4 is not zero. Blank for inconclusive signals. f_s is a direction-filtered operational annotation, not a measured GFR. |
| f_s_fallback_prior_sensitivity | f_s with fallback PP.H4 recomputed at p12 = 5e-6 (post hoc sensitivity) |
| classification | Filtration-consistent (f_s >= 0.8), Intermediate (0.2 to < 0.8), Low support (< 0.2), Inconclusive |
| inconclusive_reason | Why f_s could not be estimated |
| n_linked_genes, linked_candidate_genes_link_weight | Candidate genes linked to the signal in the Liu et al. scorecard, with normalized link weights (sum = 1 per signal) |
| Missing values | Blank cells mean not applicable or not estimated: PP_H4 is blank when colocalization was not run; f_s columns are blank for inconclusive signals; inconclusive_reason is blank for scored signals; linked-gene columns are blank when the scorecard links no gene to the signal |
