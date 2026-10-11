"""Aim 2: does signal-level filtration specificity add to Liu's scorecard? (Plan-5, rules fixed before unsealing)
Gene-level Score_A = max over the gene's variant-gene pairs of Liu's Prioritization Score minus the benchmark-overlapping
feature (Gene Monogenic). Model D also removes the four existing eGFRcys/BUN features.
"""
import json, sys
import numpy as np, pandas as pd
from sklearn.metrics import average_precision_score as ap

OUT = sys.argv[1]
rng = np.random.default_rng(20261008)
gw = pd.read_pickle("out/gw_cache.pkl")
gw = gw[gw["Variant RSID"].notna() & gw["Target gene symbol"].notna()].copy()
num = ["Prioritization Score", "Gene Monogenic", "eGFRcys", "BUN", "Gene CDS eGFRcys", "Gene CDS BUN", "Tartet Score"]
for c in num: gw[c] = pd.to_numeric(gw[c], errors="coerce").fillna(0)
gw["sA"] = gw["Prioritization Score"] - gw["Gene Monogenic"]
gw["sD"] = gw["sA"] - gw[["eGFRcys", "BUN", "Gene CDS eGFRcys", "Gene CDS BUN"]].sum(axis=1)
G = gw.groupby("Target gene symbol").agg(scoreA=("sA", "max"), scoreD=("sD", "max"),
                                         chrom=("Chromosome", lambda x: x.mode().iloc[0]))
prim = set(pd.read_csv("out/benchmark_primary_union.csv").gene)
groop = set(pd.read_csv("out/benchmark_groopman625.csv").gene)
cl = pd.read_csv("out/benchmark_clingen_kidney_all_curations.csv")
cl_mod = set(cl[cl.CLASSIFICATION.isin(["Definitive", "Strong", "Moderate"])]["GENE SYMBOL"])
pa = pd.read_csv("raw/benchmarks/panelapp_kidney_genes_2026-10-10.csv")
pa = pa[~pa.panel_id.isin([145, 1317, 154, 521])]
pa_amb = set(pa[pa.confidence.astype(str).isin(["2", "3"])].gene)
ext = cl_mod | pa_amb

# ---- signal -> gene link weights (Liu target score), S8 credible-set variants ----
s8 = pd.read_excel("raw/science.adp4753_tables_s1_to_s31.xlsx", sheet_name="S8 Finemap.Multi", header=1, dtype=str)
s8 = s8[s8["RSID"].notna() & s8["Merged GWAS loci"].astype(str).str.startswith("chr")]
s8["signal_id"] = s8["Merged GWAS loci"] + "|" + s8["CARMA Credible Set"].astype(str)
L = s8[["signal_id", "Merged GWAS loci", "RSID"]].merge(gw[["Variant RSID", "Target gene symbol", "Tartet Score"]],
                                                        left_on="RSID", right_on="Variant RSID")
W = L.groupby(["signal_id", "Merged GWAS loci", "Target gene symbol"])["Tartet Score"].sum().reset_index(name="w")
W = W[W.w > 0]
W["w"] = W.w / W.groupby("signal_id").w.transform("sum")
F = pd.read_csv("out/aim1/f_s_signals.csv")

def gene_spec(fs, agg="mean", equal=False):
    d = W.merge(fs[["signal_id", "f_s"]], on="signal_id")
    d = d[d.f_s.notna()].copy()
    if equal: d["w"] = 1.0
    if agg == "max": s = d.groupby("Target gene symbol").f_s.max()
    else: s = d.assign(wf=d.w * d.f_s).groupby("Target gene symbol").agg(a=("wf", "sum"), b=("w", "sum")).eval("a / b")
    return s.reindex(G.index)

def z(x):
    x = pd.Series(x, index=G.index, dtype=float); m, sd = x.mean(skipna=True), x.std(skipna=True)
    return ((x - m) / sd).fillna(0.0)

chroms = sorted(G.chrom.unique())
GRID = np.round(np.arange(0, 2.01, 0.1), 1)

def heldout(y, base, spec):
    """Nested leave-one-chromosome-out: lambda chosen on other chromosomes, scored on the held-out one."""
    pred = pd.Series(np.nan, index=G.index); lam = {}
    for c in chroms:
        tr = G.chrom != c
        best = max(GRID, key=lambda l: ap(y[tr], (base + l * spec)[tr]))
        lam[c] = best; pred[~tr] = (base + best * spec)[~tr]
    return pred, lam

def boot_diff(y, p1, p0, B=2000):
    d = []; idx = {c: np.where(G.chrom.values == c)[0] for c in chroms}
    for _ in range(B):
        cs = rng.choice(chroms, len(chroms), replace=True); ii = np.concatenate([idx[c] for c in cs])
        yy = y.values[ii]
        if yy.sum() == 0: continue
        d.append(ap(yy, p1.values[ii]) - ap(yy, p0.values[ii]))
    return np.percentile(d, [2.5, 97.5]).round(4).tolist()

def across(label_set, spec, tag):
    y = pd.Series(G.index.isin(label_set).astype(int), index=G.index)
    A = z(G.scoreA); D = z(G.scoreD); S = z(spec)
    pC, lamC = heldout(y, A, S); pB, lamB = heldout(y, D, S)
    r = {"benchmark": tag, "positives": int(y.sum()), "genes": int(len(y)),
         "AUPRC_A": round(ap(y, A), 4), "AUPRC_C_heldout": round(ap(y, pC), 4),
         "AUPRC_D": round(ap(y, D), 4), "AUPRC_B_heldout": round(ap(y, pB), 4),
         "AUPRC_C_lambda1": round(ap(y, A + S), 4), "baseline_prevalence": round(y.mean(), 4),
         "lambda_C_by_chrom": sorted(set(lamC.values())), "lambda_C_median": float(np.median(list(lamC.values())))}
    r["diff_C_minus_A"] = round(r["AUPRC_C_heldout"] - r["AUPRC_A"], 4)
    r["diff_C_minus_A_CI95"] = boot_diff(y, pC, A)
    r["diff_D_minus_A"] = round(r["AUPRC_D"] - r["AUPRC_A"], 4)
    r["diff_B_minus_A"] = round(r["AUPRC_B_heldout"] - r["AUPRC_A"], 4)
    r["diff_C_lambda1_minus_A"] = round(r["AUPRC_C_lambda1"] - r["AUPRC_A"], 4)
    return r, pC, A

res = {}
spec_main = gene_spec(F)
res["coverage"] = {"genes_total": int(len(G)), "genes_with_S_g": int(spec_main.notna().sum()),
                   "benchmark_genes_with_S_g": int(spec_main[G.index.isin(prim)].notna().sum()),
                   "signals_classifiable": int(F.f_s.notna().sum()), "signals_linked": int(W.signal_id.nunique()),
                   "signals_classifiable_and_linked": int(F[F.f_s.notna()].signal_id.isin(W.signal_id).sum())}
r, pC, A = across(prim, spec_main, "primary (ClinGen Def/Strong + PanelApp green)")
res["primary"] = r
# sensitivity analyses
Fd = F.copy(); Fd.loc[Fd.method == "fallback_abf", "f_s"] = np.nan
res["sens_direct_matched_only"] = across(prim, gene_spec(Fd), "primary; directly matched signals only")[0]
res["sens_max_aggregation"] = across(prim, gene_spec(F, agg="max"), "primary; max aggregation")[0]
res["sens_equal_weights"] = across(prim, gene_spec(F, equal=True), "primary; equal link weights")[0]
res["sens_extended_benchmark"] = across(ext, spec_main, "ClinGen +Moderate, PanelApp +amber")[0]
res["secondary_groopman"] = across(groop, spec_main, "Groopman 2019 (625)")[0]

# ---- within-locus top-1 ----
lamC = res["primary"]["lambda_C_by_chrom"]
_, lam_by_chr = heldout(pd.Series(G.index.isin(prim).astype(int), index=G.index), z(G.scoreA), z(spec_main))
cand = W.groupby("Merged GWAS loci").agg(n_sig=("signal_id", "nunique"), genes=("Target gene symbol", lambda x: sorted(set(x))))
nsig_total = s8.groupby("Merged GWAS loci").signal_id.nunique()
cand["n_sig_total"] = nsig_total.reindex(cand.index)
elig = cand[(cand.n_sig >= 2) & (cand.n_sig == cand.n_sig_total) & (cand.genes.map(len) >= 2)]
elig = elig[elig.genes.map(lambda g: any(x in prim for x in g))]
Az, Sz = z(G.scoreA), z(spec_main)

def top1(genes, score):
    s = score[genes]; m = s.max(); tops = s[s == m].index
    return np.mean([g in prim for g in tops])  # expected top-1 under random tie-breaking

rows = []
for loc, rr in elig.iterrows():
    g = rr.genes; c = G.chrom[g[0]]; l = lam_by_chr[c]
    rows.append({"locus": loc, "n_genes": len(g), "top1_A": top1(g, Az), "top1_C": top1(g, Az + l * Sz)})
WL = pd.DataFrame(rows)
res["within_locus"] = {"eligible_loci_with_benchmark_gene": int(len(WL)),
                       "top1_A": round(WL.top1_A.mean(), 3) if len(WL) else None,
                       "top1_C": round(WL.top1_C.mean(), 3) if len(WL) else None,
                       "loci_C_better": int((WL.top1_C > WL.top1_A).sum()), "loci_C_worse": int((WL.top1_C < WL.top1_A).sum())}
WL.to_csv(f"{OUT}/within_locus_top1.csv", index=False)
pd.DataFrame({"scoreA": G.scoreA, "S_g": spec_main, "pred_C_heldout": pC, "benchmark": G.index.isin(prim)}).to_csv(f"{OUT}/gene_scores.csv")
json.dump(res, open(f"{OUT}/aim2_results.json", "w"), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
