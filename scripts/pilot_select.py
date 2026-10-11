"""Aim 1 coloc pilot: locus selection by a fixed structural rule (no association values used).
Pool A: merged loci with >=2 S8 signals, all signals meeting the eGFRcys 80% coverage rule -> sample 8.
Pool B: single-signal loci meeting the rule -> sample 2. Seed fixed in advance."""
import sys, json
import numpy as np, pandas as pd
SUPP, GDIR, OUT = sys.argv[1:4]
SEED = 20261008
s8 = pd.read_excel(SUPP, sheet_name="S8 Finemap.Multi", header=1, dtype=str)
s8 = s8[s8["RSID"].notna() & s8["Merged GWAS loci"].astype(str).str.startswith("chr")].copy()
s8["signal_id"] = s8["Merged GWAS loci"] + "|" + s8["CARMA Credible Set"].astype(str)
cys = pd.read_csv(f"{GDIR}/cys.gz", sep="\t", usecols=["RSID"])
ok = s8.assign(ok=s8.RSID.isin(set(cys.RSID))).groupby("signal_id")["ok"].mean() >= 0.8
sig = pd.DataFrame({"ok": ok}); sig["locus"] = sig.index.str.split("|").str[0]
loc = sig.groupby("locus").agg(n=("ok", "size"), all_ok=("ok", "all"))
rng = np.random.default_rng(SEED)
A = sorted(loc[(loc.n >= 2) & loc.all_ok].index); B = sorted(loc[(loc.n == 1) & loc.all_ok].index)
pick = list(rng.choice(A, 8, replace=False)) + list(rng.choice(B, 2, replace=False))
m = pd.Series(pick).str.extract(r"chr(\w+)_(\d+)_(\d+)"); m.columns = ["chr", "start", "end"]
m.insert(0, "locus", pick); m["n_signals"] = [int(loc.n[p]) for p in pick]
m.to_csv(f"{OUT}/coloc_pilot_loci.csv", index=False)
liu = pd.read_csv(f"{GDIR}/liu_crea.gz", sep="\t", usecols=["MarkerName", "CHR", "POS"], dtype={"CHR": str})
rs = []
for _, r in m.iterrows():
    sub = liu[(liu.CHR == r.chr) & liu.POS.between(int(r.start), int(r.end))]
    rs.append(sub.MarkerName)
rs = pd.concat(rs).drop_duplicates()
rs.to_csv(f"{OUT}/coloc_pilot_rsids.txt", index=False, header=False)
print(json.dumps({"poolA": len(A), "poolB": len(B), "seed": SEED, "rsids": int(len(rs))}, indent=1)); print(m)
