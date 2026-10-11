"""Post hoc descriptive analysis (second-round iScience pre-submission review, 2026-10-11):
does gene-level filtration support S_g vary between candidate genes within the same locus?
Uses the primary S_g (link-weighted mean f_s) and Score_A exactly as in aim2_models.py; no new modelling."""
import json
import numpy as np, pandas as pd
gw = pd.read_pickle("out/gw_cache.pkl")
gw = gw[gw["Variant RSID"].notna() & gw["Target gene symbol"].notna()].copy()
gw["Tartet Score"] = pd.to_numeric(gw["Tartet Score"], errors="coerce").fillna(0)
s8 = pd.read_excel("raw/science.adp4753_tables_s1_to_s31.xlsx", sheet_name="S8 Finemap.Multi", header=1, dtype=str)
s8 = s8[s8["RSID"].notna() & s8["Merged GWAS loci"].astype(str).str.startswith("chr")]
s8["signal_id"] = s8["Merged GWAS loci"] + "|" + s8["CARMA Credible Set"].astype(str)
L = s8[["signal_id", "Merged GWAS loci", "RSID"]].merge(gw[["Variant RSID", "Target gene symbol", "Tartet Score"]],
                                                        left_on="RSID", right_on="Variant RSID")
W = L.groupby(["signal_id", "Merged GWAS loci", "Target gene symbol"])["Tartet Score"].sum().reset_index(name="w")
W = W[W.w > 0]
F = pd.read_csv("out/aim1/f_s_signals.csv")
G = pd.read_csv("out/aim1/gene_scores.csv").set_index("Target gene symbol")
# candidate genes of a locus = genes linked to any scored signal in that locus
Wsc = W.merge(F[["signal_id", "f_s"]], on="signal_id")
Wsc = Wsc[Wsc.f_s.notna()]
LG = Wsc[["Merged GWAS loci", "Target gene symbol"]].drop_duplicates()
LG = LG.join(G[["S_g", "scoreA", "benchmark"]], on="Target gene symbol")
LG = LG[LG.S_g.notna()]
TOL = 1e-6
rows = []
for loc, d in LG.groupby("Merged GWAS loci"):
    rows.append(dict(locus=loc, n_genes=len(d), sd_S=d.S_g.std(ddof=0), range_S=d.S_g.max() - d.S_g.min(),
                     range_A=d.scoreA.max() - d.scoreA.min(), n_bench=int(d.benchmark.sum()),
                     nsig=Wsc[Wsc["Merged GWAS loci"] == loc].signal_id.nunique()))
R = pd.DataFrame(rows)
M = R[R.n_genes >= 2]
ident = M.range_S < TOL
near = M.range_S < 0.05
# variance decomposition of S_g across locus-gene rows (loci with >=2 genes)
LGm = LG[LG["Merged GWAS loci"].isin(M.locus)]
grand = LGm.S_g.mean()
within = LGm.groupby("Merged GWAS loci").S_g.transform(lambda x: x - x.mean())
frac_within = float((within ** 2).sum() / ((LGm.S_g - grand) ** 2).sum())
LGm_A = LGm.scoreA; withinA = LGm.groupby("Merged GWAS loci").scoreA.transform(lambda x: x - x.mean())
frac_within_A = float((withinA ** 2).sum() / ((LGm_A - LGm_A.mean()) ** 2).sum())
# loci containing >=1 benchmark gene and >=2 genes: can S_g single out the benchmark gene?
B = []
for loc, d in LGm.groupby("Merged GWAS loci"):
    if d.benchmark.sum() == 0: continue
    top = d.S_g.max(); bmax = d[d.benchmark].S_g.max()
    topA = d.scoreA.max()
    B.append(dict(locus=loc, bench_strict_top=bool(bmax >= top - TOL and (d[~d.benchmark].S_g < top - TOL).all()),
                  bench_strict_top_A=bool(d[d.benchmark].scoreA.max() >= topA - TOL and (d[~d.benchmark].scoreA < topA - TOL).all()),
                  expected_random=d.benchmark.sum() / len(d),
                  # tie-aware: probability that the top-ranked gene is a benchmark gene, ties broken at random
                  top1_tie_S=float(d[d.S_g >= top - TOL].benchmark.mean()),
                  top1_tie_A=float(d[d.scoreA >= topA - TOL].benchmark.mean()),
                  tied_with_nonbench=bool(((d[~d.benchmark].S_g - bmax).abs() < TOL).any())))
B = pd.DataFrame(B)
single = R[R.n_genes >= 2]
out = {
 "loci_with_scored_candidate_genes": int(len(R)),
 "loci_with_>=2_candidate_genes": int(len(M)),
 "median_candidate_genes_per_locus(>=2)": float(M.n_genes.median()),
 "IQR_candidate_genes": M.n_genes.quantile([.25, .75]).tolist(),
 "loci_identical_S_g(range<1e-6)": int(ident.sum()), "pct_identical": round(100 * ident.mean(), 1),
 "loci_range_S_g<0.05": int(near.sum()), "pct_range<0.05": round(100 * near.mean(), 1),
 "median_within_locus_SD_S_g": round(float(M.sd_S.median()), 4),
 "median_within_locus_range_S_g": round(float(M.range_S.median()), 4),
 "IQR_range_S_g": [round(v, 4) for v in M.range_S.quantile([.25, .75]).tolist()],
 "single_signal_loci(>=2 genes)": int((M.nsig == 1).sum()),
 "single_signal_loci_identical": int(((M.nsig == 1) & ident).sum()),
 "multi_signal_loci(>=2 genes)": int((M.nsig > 1).sum()),
 "multi_signal_loci_identical": int(((M.nsig > 1) & ident).sum()),
 "fraction_variance_within_locus_S_g": round(frac_within, 3),
 "fraction_variance_within_locus_ScoreA": round(frac_within_A, 3),
 "loci_identical_S_g_but_ScoreA_varies": int((ident & (M.range_A > 0)).sum()),
 "benchmark_loci(>=2 genes)": int(len(B)),
 "benchmark_loci_bench_gene_uniquely_top_S_g": int(B.bench_strict_top.sum()) if len(B) else 0,
 "benchmark_loci_bench_gene_uniquely_top_ScoreA": int(B.bench_strict_top_A.sum()),
 "benchmark_loci_expected_top_if_random_untied": round(float(B.expected_random.sum()), 1),
 "benchmark_loci_expected_top1_tie_broken_S_g": round(float(B.top1_tie_S.sum()), 1),
 "benchmark_loci_expected_top1_tie_broken_ScoreA": round(float(B.top1_tie_A.sum()), 1),
 "variance_decomposition_rows(locus-gene, unweighted, same rows for S_g and Score_A)": int(len(LGm)),
 "benchmark_loci_bench_tied_with_nonbench": int(B.tied_with_nonbench.sum()) if len(B) else 0,
}
json.dump(out, open("out/aim1/aim2_within_locus.json", "w"), indent=1)
R.to_csv("out/aim1/aim2_within_locus_by_locus.csv", index=False)
print(json.dumps(out, indent=1))
