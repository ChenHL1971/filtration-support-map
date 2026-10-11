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
def analyse(ep, B=10000):
    m = d.merge(align(load_fg(ep)), on="signal_id")
    X = m[covs].copy(); X = X.fillna(X.median()); X = (X - X.mean()) / X.std()
    y = (m.cls == "filtration").astype(int).values
    ps = LogisticRegression(max_iter=1000).fit(X, y).predict_proba(X)[:, 1]
    strata = pd.qcut(ps, 5, labels=False, duplicates="drop")
    hit = ((m.pval < 0.05) & (m.b > 0)).astype(int).values
    stat = hit[y == 1].mean() - hit[y == 0].mean()
    null = np.empty(B)
    for k in range(B):
        yy = y.copy()
        for s in np.unique(strata):
            ii = np.where(strata == s)[0]; yy[ii] = rng.permutation(y[ii])
        null[k] = hit[yy == 1].mean() - hit[yy == 0].mean()
    p = (1 + (null >= stat).sum()) / (B + 1)
    out = {"n_filtration": int(y.sum()), "n_low": int((1 - y).sum()),
           "prop_hit_filtration": round(hit[y == 1].mean(), 4), "prop_hit_low": round(hit[y == 0].mean(), 4),
           "difference": round(stat, 4), "perm_p_one_sided": round(p, 5)}
    for c, lab in [(1, "filtration"), (0, "low")]:
        b = m.b.values[y == c]; out[f"mean_logOR_{lab}"] = round(b.mean(), 4)
        se = b.std(ddof=1) / np.sqrt(len(b)); out[f"mean_logOR_{lab}_CI95"] = [round(b.mean() - 1.96 * se, 4), round(b.mean() + 1.96 * se, 4)]
        x = m.abs_beta.values[y == c]; w = 1 / m.sebeta.values[y == c] ** 2
        slope = np.sum(w * x * b) / np.sum(w * x * x); sse = np.sqrt(1 / np.sum(w * x * x))
        out[f"slope_logOR_per_eGFRcrea_{lab}"] = [round(slope, 3), round(sse, 3)]
    return out

for ep in ["N14_DIALYSIS", "Z21_KIDNEY_TRANSPLANT_STATUS", "N14_CHRONKIDNEYDIS", "N14_ACUTERENFAIL", "H7_VITROPACIT"]:
    res[ep] = analyse(ep)
sec = ["Z21_KIDNEY_TRANSPLANT_STATUS", "N14_CHRONKIDNEYDIS", "N14_ACUTERENFAIL"]
ps = np.array([res[e]["perm_p_one_sided"] for e in sec]); o = np.argsort(ps); q = np.empty(3)
q[o] = np.minimum.accumulate((ps[o] * 3 / np.arange(1, 4))[::-1])[::-1]
for e, v in zip(sec, q): res[e]["BH_q"] = round(min(v, 1), 5)
json.dump(res, open("out/aim3_results.json", "w"), indent=1)
print(json.dumps(res, indent=1))
