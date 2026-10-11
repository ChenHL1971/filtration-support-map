"""Titan consistency tests (reproducible record). Lead-variant map uses all 54 ST1 creatinine-specific regions
(the 2 regions absent from ST3 do not map); region-based test reported for the 52-locus list (ST3) and all 54."""
import json, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
rng = np.random.default_rng(20261010)
F = pd.read_csv("out/aim1/f_s_signals.csv")
def hl(c, o): return float(np.median(np.subtract.outer(c, o)))
def test(c, o, boot=False):
    r = {"n_cr": int(len(c)), "n_gfr": int(len(o)), "med_cr": round(float(np.median(c)), 3), "med_gfr": round(float(np.median(o)), 3),
         "p_one_sided": float(mannwhitneyu(c, o, alternative="less").pvalue), "HL_shift": round(hl(c, o), 3)}
    if boot:
        b = [hl(rng.choice(c, len(c)), rng.choice(o, len(o))) for _ in range(5000)]
        r["HL_CI95_percentile_bootstrap_5000"] = [round(v, 2) for v in np.percentile(b, [2.5, 97.5])]
    return r
m = pd.read_csv("out/titan_to_s8_map.csv").merge(F[["signal_id", "f_s", "method"]], on="signal_id").dropna(subset=["f_s"])
cr = m["class"] == "Creatinine-specific"
res = {"lead_variant": test(m[cr].f_s.values, m[~cr].f_s.values, boot=True)}
d = m[m.method.isin(["overlap", "LD"])]; dc = d["class"] == "Creatinine-specific"
res["direct_only"] = test(d[dc].f_s.values, d[~dc].f_s.values)
cl = pd.read_csv("out/titan_to_s8_cluster_map.csv")
cl = cl[cl.cls.isin(["Creatinine-specific", "GFR-associated"])].merge(F[["signal_id", "f_s"]], left_on="sid", right_on="signal_id").dropna(subset=["f_s"])
drop = {"1:9429264_G_GT", "2:208463210_TAA_T"}   # in ST1 but not in the 52-locus ST3 list
for tag, x in [("region_52_loci", cl[~cl.marker.isin(drop)]), ("region_54_loci", cl)]:
    x = x.drop_duplicates(["cls", "signal_id"]); c = x.cls == "Creatinine-specific"
    res[tag] = test(x[c].f_s.values, x[~c].f_s.values)
json.dump(res, open("out/titan_tests.json", "w"), indent=1); print(json.dumps(res, indent=1))
