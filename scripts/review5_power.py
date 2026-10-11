"""Post hoc (review round 5): detectable enrichment for the FinnGen dialysis comparison, simulated with the same design as
the primary test: 151 filtration-consistent vs. 399 low-support signals distributed over the observed propensity-score
quintiles, a one-sided permutation test of the difference in hit proportions with labels permuted within quintiles.
Baseline hit probability = observed overall low-support rate (16/399), applied in every stratum; filtration-consistent
signals have odds multiplied by OR. Signals are treated as independent after the 1-Mb pruning used in the analysis."""
import json, numpy as np
rng = np.random.default_rng(20261008)
strata = [(105, 5), (97, 13), (85, 25), (71, 39), (41, 69)]          # (low, filtration) per quintile
p0 = 16 / 399
lab = np.concatenate([np.r_[np.zeros(l), np.ones(f)] for l, f in strata]).astype(bool)
sid = np.concatenate([np.full(l + f, k) for k, (l, f) in enumerate(strata)])
idx = [np.where(sid == k)[0] for k in range(len(strata))]
def stat(hit, lab): return hit[lab].mean() - hit[~lab].mean()
def perm_p(hit, B=999):
    obs = stat(hit, lab); cnt = 0
    for _ in range(B):
        pl = lab.copy()
        for ii in idx: pl[ii] = rng.permutation(pl[ii])
        cnt += stat(hit, pl) >= obs
    return (1 + cnt) / (B + 1)
def power(OR, sims=1000):
    o = p0 / (1 - p0) * OR; p1 = o / (1 + o)
    prob = np.where(lab, p1, p0); rej = 0
    for _ in range(sims):
        hit = rng.random(len(lab)) < prob
        rej += perm_p(hit) < 0.05
    return rej / sims
res = {"design": {"n_filtration": int(lab.sum()), "n_low": int((~lab).sum()), "baseline_hit_rate": round(p0, 4),
                  "test": "one-sided, labels permuted within propensity quintiles (999 permutations per simulated dataset)",
                  "simulations_per_OR": 1000}}
res["power"] = {str(OR): power(OR) for OR in (1.5, 2.0, 2.5, 3.0, 3.5)}
json.dump(res, open("out/review5_power.json", "w"), indent=1); print(json.dumps(res, indent=1))
