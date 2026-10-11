import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from sklearn.metrics import average_precision_score as ap
src = open("scripts/aim2_secondary.py").read().split("res = {}")[0]
exec(src)  # X32, chrom, Sz, Smiss, chroms, gw, feats, y(primary)
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
import xgboost as xgb
rng = np.random.default_rng(20261008)
mgi = set(pd.read_csv("out/benchmark_mgi_renal_function.csv").gene)
prim = set(pd.read_csv("out/benchmark_primary_union.csv").gene)
GRID = np.round(np.arange(0, 2.01, 0.1), 1)
def z(s): return ((s - s.mean()) / s.std()).fillna(0)
def boot(yv, p1, p0, B=1000):
    idx = {c: np.where(chrom.values == c)[0] for c in chroms}; d = []
    for _ in range(B):
        ii = np.concatenate([idx[c] for c in rng.choice(chroms, len(chroms))])
        if yv.values[ii].sum() == 0: continue
        d.append(ap(yv.values[ii], p1.values[ii]) - ap(yv.values[ii], p0.values[ii]))
    return np.percentile(d, [2.5, 97.5]).round(4).tolist()
def heldout(yv, base, spec):
    pred = pd.Series(np.nan, index=base.index)
    for c in chroms:
        tr = chrom != c; best = max(GRID, key=lambda l: ap(yv[tr], (base + l * spec)[tr])); pred[~tr] = (base + best * spec)[~tr]
    return pred
def score_minus(cols):
    s = gw[feats].sum(axis=1) - gw[cols].sum(axis=1)
    return s.groupby(gw["Target gene symbol"]).max().reindex(X32.index)
def nested(fn, X, yv):
    pred = pd.Series(np.nan, index=X.index)
    for c in chroms:
        tr, te = chrom != c, chrom == c; pred[te] = fn(X[tr], yv[tr], chrom[tr])(X[te])
    return pred
def l1_fn(Xtr, ytr, gtr):
    best, bc = -1, None
    for C in [0.01, 0.03, 0.1, 0.3, 1.0]:
        sc = [ap(ytr.iloc[b], LogisticRegression(penalty="l1", solver="liblinear", C=C, class_weight="balanced").fit(Xtr.iloc[a], ytr.iloc[a]).predict_proba(Xtr.iloc[b])[:, 1]) for a, b in GroupKFold(5).split(Xtr, ytr, gtr)]
        if np.mean(sc) > best: best, bc = np.mean(sc), C
    m = LogisticRegression(penalty="l1", solver="liblinear", C=bc, class_weight="balanced").fit(Xtr, ytr); return lambda X: m.predict_proba(X)[:, 1]
res = {}
for tag, labels, drop in [("MGI renal function (mouse feature removed)", mgi, ["Gene MousePhynotype"]),
                          ("primary, mouse feature also removed (exploratory)", prim, ["Gene Monogenic", "Gene MousePhynotype"])]:
    yv = pd.Series(X32.index.isin(labels).astype(int), index=X32.index)
    A = z(score_minus(drop)); pC = heldout(yv, A, Sz)
    base = X32.drop(columns=drop)
    l0 = nested(l1_fn, base, yv); l1 = nested(l1_fn, base.assign(S_g=Sz, S_missing=Smiss), yv)
    res[tag] = {"positives": int(yv.sum()), "AUPRC_A_equal_weights": round(ap(yv, A), 4), "AUPRC_C_heldout": round(ap(yv, pC), 4),
                "diff_C_minus_A": round(ap(yv, pC) - ap(yv, A), 4), "diff_CI95": boot(yv, pC, A),
                "AUPRC_C_lambda1": round(ap(yv, A + Sz), 4),
                "L1_without_S": round(ap(yv, l0), 4), "L1_with_S": round(ap(yv, l1), 4), "L1_diff_CI95": boot(yv, l1, l0)}
json.dump(res, open("out/aim1/aim2_mgi.json", "w"), indent=1); print(json.dumps(res, indent=1))
