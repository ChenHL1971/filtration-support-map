"""Added for the code release. Rebuilds out/gw_cache.pkl (the cache read by the Aim 2/3 scripts) from the redistributed
copy of the Liu et al. 2025 genome-wide Kidney Disease Genetic Scorecard (Figshare 26299093, CC BY 4.0), without
downloading the 66-MB Excel file. The CSV holds the same 367,308 variant-gene rows and 46 columns, all as text, exactly
as read by scripts/audit_gw.py (pd.read_excel(..., dtype=str)); the rebuilt cache is identical to the one used.
Usage: python3 scripts/build_gw_cache.py"""
import pandas as pd
gw = pd.read_csv("data_derived/liu2025_genomewide_scorecard.csv.gz", dtype=str, keep_default_na=False, na_filter=False)
gw.to_pickle("out/gw_cache.pkl")
print(gw.shape)
