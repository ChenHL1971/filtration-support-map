"""Post hoc (review round 6): sparse-outcome sensitivity for the FinnGen dialysis comparison.
Firth-penalized logistic regression with the same covariates as the primary adjusted model (profile penalized-likelihood
95% CI), and a Newcombe hybrid-score 95% CI for the difference in hit proportions. Reuses the data preparation of
aim3_finngen_v2.py unchanged (executed up to the model definitions)."""
import json, numpy as np, pandas as pd
from scipy.stats import norm, chi2
src = open("scripts/aim3_finngen_v2.py").read(); exec(src[:src.index("covs = [")])
covs = ["absz", "MAF", "ldscore", "gene_density", "tss_dist", "cs_size", "n_genes"]
m = d.merge(align(load_fg("N14_DIALYSIS")), on="signal_id")
X = m[covs].copy(); X = X.fillna(X.median()); X = (X - X.mean()) / X.std()
y = (m.cls == "filtration").astype(int).values
hit = ((m.pval < 0.05) & (m.b > 0)).astype(int).values
Xa = np.column_stack([np.ones(len(y)), y, X.values])

def firth(Xm, yv, fixed=None, iters=200):
    """Firth logistic regression; fixed=(index, value) constrains one coefficient (for profile likelihood)."""
    b = np.zeros(Xm.shape[1]); free = np.ones(len(b), bool)
    if fixed: b[fixed[0]] = fixed[1]; free[fixed[0]] = False
    for _ in range(iters):
        p = 1 / (1 + np.exp(-Xm @ b)); W = p * (1 - p)
        I = Xm.T @ (Xm * W[:, None]); Ii = np.linalg.inv(I)
        h = W * np.einsum("ij,jk,ik->i", Xm, Ii, Xm)
        U = Xm.T @ (yv - p + h * (0.5 - p))
        Iff = I[np.ix_(free, free)]; step = np.linalg.solve(Iff, U[free])
        step = np.clip(step, -5, 5); b[free] += step
        if np.max(np.abs(step)) < 1e-10: break
    p = 1 / (1 + np.exp(-Xm @ b)); W = p * (1 - p); I = Xm.T @ (Xm * W[:, None])
    ll = np.sum(yv * np.log(p) + (1 - yv) * np.log(1 - p)) + 0.5 * np.linalg.slogdet(I)[1]
    return b, ll

b, ll = firth(Xa, hit)
crit = chi2.ppf(0.95, 1) / 2
def bound(direction):
    lo, hi = b[1], b[1] + direction * 6
    for _ in range(60):
        mid = (lo + hi) / 2; _, l = firth(Xa, hit, fixed=(1, mid))
        if ll - l > crit: hi = mid
        else: lo = mid
    return (lo + hi) / 2
ci = sorted([bound(-1), bound(1)])
_, l0 = firth(Xa, hit, fixed=(1, 0.0)); p_firth = float(chi2.sf(2 * (ll - l0), 1))

def wilson(x, n, z=1.96):
    p = x / n; c = (p + z * z / (2 * n)) / (1 + z * z / n); h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h
x1, n1, x0, n0 = int(hit[y == 1].sum()), int((y == 1).sum()), int(hit[y == 0].sum()), int((y == 0).sum())
p1, p0 = x1 / n1, x0 / n0; l1, u1 = wilson(x1, n1); l0_, u0 = wilson(x0, n0)
diff = p1 - p0
newc = [diff - np.sqrt((p1 - l1) ** 2 + (u0 - p0) ** 2), diff + np.sqrt((u1 - p1) ** 2 + (p0 - l0_) ** 2)]
res = {"hits": {"filtration": [x1, n1], "low": [x0, n0]}, "events_per_parameter": round(int(hit.sum()) / 8, 2),
       "firth_OR": round(float(np.exp(b[1])), 3), "firth_CI95_profile": [round(float(np.exp(v)), 3) for v in ci],
       "firth_p_LR": round(p_firth, 4),
       "difference_in_hit_proportion": round(diff, 4), "difference_CI95_newcombe": [round(float(v), 4) for v in newc]}
json.dump(res, open("out/review6_firth.json", "w"), indent=1); print(json.dumps(res, indent=1))
