"""Post hoc diagnostics after review: covariate balance, strata overlap, hit counts."""
"""Aim 3: signal-level FinnGen enrichment (rules fixed 2026-10-10 before FinnGen associations were opened)."""
import json, sys
import numpy as np, pandas as pd
from scipy.stats import norm
from sklearn.linear_model import LogisticRegression

rng = np.random.default_rng(20261008)
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
lead = pd.read_csv("out/aim3_signal_leads.csv", dtype={"chr": str})
F = pd.read_csv("out/aim1/f_s_signals.csv")
d = lead.merge(F[["signal_id", "f_s"]], on="signal_id")
d["cls"] = np.where(d.f_s >= 0.8, "filtration", np.where(d.f_s < 0.2, "low", None))
d = d[d.cls.notna()].copy()
import os
if os.environ.get("DIRECT_ONLY")=="1":
    meth = F.set_index("signal_id").method
    d = d[d.signal_id.map(meth) != "fallback_abf"].copy()

# eGFRcrea effect at lead (Liu European), aligned to eGFRcrea-lowering allele
liu = pd.read_csv("raw/gwas/liu_crea.gz", sep="\t", usecols=["MarkerName", "REF", "ALT", "MAF", "Zscore", "BETA"])
liu = liu.drop_duplicates("MarkerName").rename(columns={"MarkerName": "rsid"})
d = d.merge(liu, on="rsid", how="left")
d["risk_allele"] = np.where(d.Zscore < 0, d.ALT.str.upper(), d.REF.str.upper())
d["other_allele"] = np.where(d.Zscore < 0, d.REF.str.upper(), d.ALT.str.upper())
d["absz"] = d.Zscore.abs(); d["abs_beta"] = d.BETA.abs()

# gene density and TSS distance (PoPS hg19 annotation); candidate genes from scorecard links
ga = pd.read_csv("raw/annot/gene_annot_jun10.txt", sep="\t", dtype={"CHR": str})
def tss_feats(r):
    g = ga[ga.CHR == r.chr]; dist = (g.TSS - r.pos).abs()
    return pd.Series({"gene_density": int((dist <= 5e5).sum()), "tss_dist": float(dist.min()) if len(dist) else np.nan})
d[["gene_density", "tss_dist"]] = d.apply(tss_feats, axis=1)
gw = pd.read_pickle("out/gw_cache.pkl")
s8 = pd.read_excel("raw/science.adp4753_tables_s1_to_s31.xlsx", sheet_name="S8 Finemap.Multi", header=1, dtype=str)
s8 = s8[s8["RSID"].notna() & s8["Merged GWAS loci"].astype(str).str.startswith("chr")]
s8["signal_id"] = s8["Merged GWAS loci"] + "|" + s8["CARMA Credible Set"].astype(str)
gw["ts"] = pd.to_numeric(gw["Tartet Score"], errors="coerce").fillna(0)
L = s8[["signal_id", "RSID"]].merge(gw[gw.ts > 0][["Variant RSID", "Target gene symbol"]], left_on="RSID", right_on="Variant RSID")
d["n_genes"] = d.signal_id.map(L.groupby("signal_id")["Target gene symbol"].nunique()).fillna(0)

res = {"class_counts_before_pruning": d.cls.value_counts().to_dict(), "missing_eGFRcrea_effect": int(d.Zscore.isna().sum())}
d = d[d.Zscore.notna()]
# LD pruning: descending |z|, keep if >= 1 Mb from all kept on the same chromosome
kept = []; pos_by_chr = {}
for i, r in d.sort_values("absz", ascending=False).iterrows():
    ps = pos_by_chr.setdefault(r.chr, [])
    if all(abs(r.pos - p) >= 1e6 for p in ps): kept.append(i); ps.append(r.pos)
d = d.loc[kept].copy()
res["class_counts_after_pruning"] = d.cls.value_counts().to_dict()

def load_fg(ep):
    f = pd.read_csv(f"raw/finngen/sumstats/fg_{ep}_kidneyloci_tsv.gz", sep="\t", usecols=["ref", "alt", "rsids", "pval", "beta", "sebeta"])
    f = f.assign(rsid=f.rsids.str.split(",")).explode("rsid")
    return f[f.rsid.isin(d.rsid)]

def align(f):
    m = d[["signal_id", "rsid", "risk_allele", "other_allele"]].merge(f, on="rsid")
    pal = (m.risk_allele + m.other_allele).isin(["AT", "TA", "CG", "GC"])
    same = (m.alt == m.risk_allele) & (m.ref == m.other_allele); swap = (m.alt == m.other_allele) & (m.ref == m.risk_allele)
    fs = ~pal & (m.alt.map(COMP) == m.risk_allele) & (m.ref.map(COMP) == m.other_allele)
    fw = ~pal & (m.alt.map(COMP) == m.other_allele) & (m.ref.map(COMP) == m.risk_allele)
    sgn = np.select([same | fs, swap | fw], [1, -1], 0)
    m = m[sgn != 0].assign(b=lambda x: x.beta * sgn[sgn != 0])
    return m.drop_duplicates("signal_id")[["signal_id", "b", "sebeta", "pval"]]

covs = ["absz", "MAF", "ldscore", "gene_density", "tss_dist", "cs_size", "n_genes"]

def balance(ep):
    m = d.merge(align(load_fg(ep)), on="signal_id")
    X = m[covs].copy(); X = X.fillna(X.median()); Xs = (X - X.mean()) / X.std()
    y = (m.cls == "filtration").astype(int).values
    ps = LogisticRegression(max_iter=1000).fit(Xs, y).predict_proba(Xs)[:, 1]
    strata = pd.qcut(ps, 5, labels=False, duplicates="drop")
    hit = ((m.pval < 0.05) & (m.b > 0)).astype(int).values
    rows=[]
    for c in covs:
        x=X[c].values; sd=np.sqrt((x[y==1].var(ddof=1)+x[y==0].var(ddof=1))/2)
        raw=(x[y==1].mean()-x[y==0].mean())/sd
        # stratified SMD: weight strata by number of filtration signals
        num=0; w=0
        for s in np.unique(strata):
            ii=strata==s; a=x[ii&(y==1)]; b=x[ii&(y==0)]
            if len(a) and len(b): num+=len(a)*(a.mean()-b.mean()); w+=len(a)
        rows.append({"covariate":c,"mean_filtration":round(x[y==1].mean(),3),"mean_low":round(x[y==0].mean(),3),"SMD_raw":round(raw,3),"SMD_within_strata":round(num/w/sd,3)})
    q=pd.crosstab(strata,y); q.columns=["low","filtration"]
    q["ps_range"]=[f"{ps[strata==s].min():.3f}-{ps[strata==s].max():.3f}" for s in q.index]
    return pd.DataFrame(rows), q, (int(hit[y==1].sum()),int((y==1).sum()),int(hit[y==0].sum()),int((y==0).sum()))
res={}
for ep in ["N14_DIALYSIS","Z21_KIDNEY_TRANSPLANT_STATUS","N14_CHRONKIDNEYDIS","N14_ACUTERENFAIL","H7_VITROPACIT"]:
    B,Q,H=balance(ep); res[ep]=H
    if ep=="N14_DIALYSIS": print(B.to_string()); print(Q.to_string()); B.to_csv("out/aim3_balance_dialysis.csv",index=False); Q.to_csv("out/aim3_strata_dialysis.csv")
print(res)
json.dump(res,open("out/aim3_hit_counts.json","w"))
