"""Added 2026-10-11 for the code release (reconstructs a step originally run interactively; verified identical).
Primary benchmark = ClinGen Definitive/Strong genes from the five kidney GCEPs (download 2026-10-09)
UNION PanelApp green (confidence 3) genes from kidney panels (download 2026-10-10), excluding panels 145, 1317, 154, 521
(adrenal and renal cancer). Groopman 2019 and MGI lists are provided as processed files (see raw/README.md).
Usage: python3 scripts/build_benchmarks.py [outdir=out]"""
import sys, pandas as pd
OUT = sys.argv[1] if len(sys.argv) > 1 else "out"
cl = pd.read_csv("raw/benchmarks/clingen_gdv_2026-10-09.csv", skiprows=4)
cl = cl[~cl["GENE SYMBOL"].astype(str).str.startswith("+")]
GCEP = ["Kidney Cystic and Ciliopathy Disorders", "Tubulopathy", "Glomerulopathy", "Complement-Mediated Kidney Diseases",
        "Congenital Anomalies of the Kidney and Urinary Tract"]
k = cl[cl.GCEP.isin([g + " Gene Curation Expert Panel" for g in GCEP])]
k.to_csv(f"{OUT}/benchmark_clingen_kidney_all_curations.csv", index=False)
cg = sorted(set(k[k.CLASSIFICATION.isin(["Definitive", "Strong"])]["GENE SYMBOL"]))
pd.DataFrame({"gene": cg}).to_csv(f"{OUT}/benchmark_clingen_kidney_primary.csv", index=False)
pa = pd.read_csv("raw/benchmarks/panelapp_kidney_genes_2026-10-10.csv")
pg = sorted(set(pa[~pa.panel_id.isin([145, 1317, 154, 521]) & (pa.confidence.astype(str) == "3")].gene))
pd.DataFrame({"gene": pg}).to_csv(f"{OUT}/benchmark_panelapp_green.csv", index=False)
u = sorted(set(cg) | set(pg))
pd.DataFrame({"gene": u, "clingen": [g in cg for g in u], "panelapp_green": [g in pg for g in u]}).to_csv(
    f"{OUT}/benchmark_primary_union.csv", index=False)
