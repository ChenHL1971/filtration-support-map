"""Aim 2 secondary (prespecified): learned models in chromosome-blocked nested CV; coding-variant benchmark (within-locus)."""
import json
import numpy as np, pandas as pd
from sklearn.metrics import average_precision_score as ap
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
import xgboost as xgb

rng = np.random.default_rng(20261008)
gw = pd.read_pickle("out/gw_cache.pkl")
gw = gw[gw["Variant RSID"].notna() & gw["Target gene symbol"].notna()].copy()
feats = ["eGFRcys", "BUN", "TSS1Kb", "CredibleSet", "snATACPeak", "bASA", "snASA", "HMMReg",
         "Gene CDS", "Gene CDS pLoF", "Gene CDS Missense", "Gene CDS eGFRcys", "Gene CDS BUN", "Gene CDS CredibleSet",
         "Gene ExWAS", "Gene Monogenic", "Gene MousePhynotype", "Gene PheWAS", "Gene MAGMA", "Gene DEPICT",
         "ASE of Tubule", "ASE of Glomerlus", "eQTL", "Coloc with tubule ASE", "Coloc with glomerlus ASE",
         "Coloc with kidney eQTL", "Open4Gene link", "Same Top Cell Type", "meQTL", "Coaccess", "ABC", "cS2G"]
for c in feats + ["Tartet Score"]: gw[c] = pd.to_numeric(gw[c], errors="coerce").fillna(0)
X32 = gw.groupby("Target gene symbol")[feats].max()
chrom = gw.groupby("Target gene symbol")["Chromosome"].agg(lambda x: x.mode().iloc[0]).reindex(X32.index)
gs = pd.read_csv("out/aim1/gene_scores.csv", index_col=0).reindex(X32.index)
S = gs.S_g; Sz = ((S - S.mean()) / S.std()).fillna(0); Smiss = S.isna().astype(int)
prim = set(pd.read_csv("out/benchmark_primary_union.csv").gene)
y = pd.Series(X32.index.isin(prim).astype(int), index=X32.index)
chroms = sorted(chrom.unique())
res = {}

def nested(model_fn, X):
    pred = pd.Series(np.nan, index=X.index)
    for c in chroms:
        tr, te = chrom != c, chrom == c
        m = model_fn(X[tr], y[tr], chrom[tr]); pred[te] = m(X[te])
    return pred

def l1_fn(Xtr, ytr, gtr):
    best, bc = -1, None
    for C in [0.01, 0.03, 0.1, 0.3, 1.0]:
        sc = []
        for a, b in GroupKFold(5).split(Xtr, ytr, gtr):
            m = LogisticRegression(penalty="l1", solver="liblinear", C=C, class_weight="balanced").fit(Xtr.iloc[a], ytr.iloc[a])
            sc.append(ap(ytr.iloc[b], m.predict_proba(Xtr.iloc[b])[:, 1]))
        if np.mean(sc) > best: best, bc = np.mean(sc), C
    m = LogisticRegression(penalty="l1", solver="liblinear", C=bc, class_weight="balanced").fit(Xtr, ytr)
    return lambda Xte: m.predict_proba(Xte)[:, 1]

def xgb_fn(Xtr, ytr, gtr):
    m = xgb.XGBClassifier(n_estimators=200, max_depth=3, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8,
                          scale_pos_weight=(1 - ytr.mean()) / ytr.mean(), eval_metric="aucpr", random_state=1, n_jobs=1)
    m.fit(Xtr, ytr); return lambda Xte: m.predict_proba(Xte)[:, 1]

def boot(p1, p0, B=1000):
    idx = {c: np.where(chrom.values == c)[0] for c in chroms}; d = []
    for _ in range(B):
        ii = np.concatenate([idx[c] for c in rng.choice(chroms, len(chroms))])
        if y.values[ii].sum() == 0: continue
        d.append(ap(y.values[ii], p1.values[ii]) - ap(y.values[ii], p0.values[ii]))
    return np.percentile(d, [2.5, 97.5]).round(4).tolist()

base = X32.drop(columns=["Gene Monogenic"])              # benchmark-overlapping feature removed
withS = base.assign(S_g=Sz, S_missing=Smiss)
for name, fn in [("L1_logistic", l1_fn), ("XGBoost", xgb_fn)]:
    p0 = nested(fn, base); p1 = nested(fn, withS)
    res[name] = {"AUPRC_without_S": round(ap(y, p0), 4), "AUPRC_with_S": round(ap(y, p1), 4),
                 "diff": round(ap(y, p1) - ap(y, p0), 4), "diff_CI95": boot(p1, p0)}
m = xgb.XGBClassifier(n_estimators=200, max_depth=3, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8,
                      scale_pos_weight=(1 - y.mean()) / y.mean(), random_state=1, n_jobs=1).fit(withS, y)
shap = m.get_booster().predict(xgb.DMatrix(withS), pred_contribs=True)[:, :-1]
imp = pd.Series(np.abs(shap).mean(0), index=withS.columns).sort_values(ascending=False)
res["XGBoost_mean_abs_SHAP_top10"] = imp.head(10).round(4).to_dict()
res["XGBoost_SHAP_rank_S_g"] = int(list(imp.index).index("S_g") + 1)

# ---- coding-variant benchmark (within-locus); coding features removed from every model ----
coding = ["Gene CDS", "Gene CDS pLoF", "Gene CDS Missense", "Gene CDS eGFRcys", "Gene CDS BUN", "Gene CDS CredibleSet"]
gw["sNC"] = gw["Prioritization Score"].astype(float) - gw[coding].sum(axis=1)
scoreNC = gw.groupby("Target gene symbol").sNC.max().reindex(X32.index)
Az = (scoreNC - scoreNC.mean()) / scoreNC.std()
pos_cod = set(X32.index[X32["Gene CDS CredibleSet"] == 1])
s8 = pd.read_excel("raw/science.adp4753_tables_s1_to_s31.xlsx", sheet_name="S8 Finemap.Multi", header=1, dtype=str)
s8 = s8[s8["RSID"].notna() & s8["Merged GWAS loci"].astype(str).str.startswith("chr")]
s8["signal_id"] = s8["Merged GWAS loci"] + "|" + s8["CARMA Credible Set"].astype(str)
L = s8[["signal_id", "Merged GWAS loci", "RSID"]].merge(gw[gw["Tartet Score"] > 0][["Variant RSID", "Target gene symbol"]], left_on="RSID", right_on="Variant RSID")
cand = L.groupby("Merged GWAS loci").agg(n_sig=("signal_id", "nunique"), genes=("Target gene symbol", lambda x: sorted(set(x))))
cand["n_sig_total"] = s8.groupby("Merged GWAS loci").signal_id.nunique().reindex(cand.index)
el = cand[(cand.n_sig >= 2) & (cand.n_sig == cand.n_sig_total) & (cand.genes.map(len) >= 2)]
el = el[el.genes.map(lambda g: any(x in pos_cod for x in g))]
def top1(g, sc):
    s = sc[g]; t = s[s == s.max()].index; return np.mean([x in pos_cod for x in t])
rows = [{"locus": l, "A": top1(r.genes, Az), "C_l1": top1(r.genes, Az + Sz), "C_l04": top1(r.genes, Az + 0.4 * Sz)} for l, r in el.iterrows()]
W = pd.DataFrame(rows)
res["coding_benchmark_within_locus"] = {"eligible_loci": int(len(W)), "top1_A": round(W.A.mean(), 3),
    "top1_C_lambda0.4": round(W.C_l04.mean(), 3), "top1_C_lambda1": round(W.C_l1.mean(), 3),
    "loci_better_l1": int((W.C_l1 > W.A).sum()), "loci_worse_l1": int((W.C_l1 < W.A).sum())}
json.dump(res, open("out/aim1/aim2_secondary.json", "w"), indent=1)
print(json.dumps(res, indent=1))
