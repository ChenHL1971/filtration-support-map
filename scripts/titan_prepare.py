"""Added 2026-10-11 for the code release (reconstructs steps originally run interactively; verified identical to the files used).
(1) Export Titan et al. Supplementary Table ST1 to CSV. (2) Region-based Titan-to-S8 map: every S8 signal with at least one
credible-set variant inside a Titan associated region (cluster_start-cluster_end, hg19).
Usage: python3 scripts/titan_prepare.py [outdir=out]"""
import sys, pandas as pd
OUT = sys.argv[1] if len(sys.argv) > 1 else "out"
t = pd.read_excel("raw/titan/mmc3.xlsx", sheet_name="ST1 - Associated regions")
t.to_csv(f"{OUT}/titan_ST1.csv", index=False)
s8 = pd.read_excel("raw/science.adp4753_tables_s1_to_s31.xlsx", sheet_name="S8 Finemap.Multi", header=1, dtype=str)
s8 = s8[s8.RSID.notna() & s8["Merged GWAS loci"].astype(str).str.startswith("chr")].copy()
s8["sid"] = s8["Merged GWAS loci"] + "|" + s8["CARMA Credible Set"].astype(str)
s8["pos"] = s8.Position.astype(int)
rows = []
for _, r in t.iterrows():
    m = s8[(s8.Chromosome.astype(str) == str(r.CHROMOSSOME)) & s8.pos.between(r.cluster_start, r.cluster_end)]
    rows += [(r.Marker, r.Classification, s) for s in m.sid.unique()]
pd.DataFrame(rows, columns=["marker", "cls", "sid"]).drop_duplicates().to_csv(f"{OUT}/titan_to_s8_cluster_map.csv", index=False)
