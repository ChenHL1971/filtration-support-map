"""Added 2026-10-11 for the code release (reconstructs a step originally run interactively; verified identical).
Joins the lead-variant Titan map (titan_map.R) with f_s for creatinine-specific and GFR-associated loci.
Usage: python3 scripts/titan_check.py [outfile=out/aim1/titan_check.csv]"""
import sys, pandas as pd
OUT = sys.argv[1] if len(sys.argv) > 1 else "out/aim1/titan_check.csv"
m = pd.read_csv("out/titan_to_s8_map.csv"); F = pd.read_csv("out/aim1/f_s_signals.csv")
x = m[m.signal_id.notna() & m["class"].isin(["Creatinine-specific", "GFR-associated"])].copy()
x["grp"] = x["class"].map({"Creatinine-specific": "creat", "GFR-associated": "gfr"})
x.merge(F[["signal_id", "f_s", "method"]], on="signal_id", how="left")[["marker", "class", "signal_id", "how", "grp", "f_s", "method"]].to_csv(OUT, index=False)
