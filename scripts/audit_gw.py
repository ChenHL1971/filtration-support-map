"""Feasibility audit, part 2: genome-wide scorecard (Figshare, all variant-gene pairs).

Structure, coverage and implementability only. No benchmark performance,
no FinnGen or outcome data.
"""
import sys
import json
import pandas as pd

GW, SUPP, OUT = sys.argv[1], sys.argv[2], sys.argv[3]

import os
CACHE = f"{OUT}/gw_cache.pkl"
if os.path.exists(CACHE):
    gw = pd.read_pickle(CACHE)
else:
    gw = pd.read_excel(GW, sheet_name=0, dtype=str)
    gw = gw.loc[:, ~gw.columns.astype(str).str.startswith("Unnamed")]
    gw.to_pickle(CACHE)
gw = gw[gw["Variant RSID"].notna()].copy()

variant_feats = ["eGFRcys", "BUN", "TSS1Kb", "CredibleSet", "snATACPeak", "bASA", "snASA", "HMMReg"]
gene_feats = ["Gene CDS", "Gene CDS pLoF", "Gene CDS Missense", "Gene CDS eGFRcys", "Gene CDS BUN",
              "Gene CDS CredibleSet", "Gene ExWAS", "Gene Monogenic", "Gene MousePhynotype",
              "Gene PheWAS", "Gene MAGMA", "Gene DEPICT"]
target_feats = ["ASE of Tubule", "ASE of Glomerlus", "eQTL", "Coloc with tubule ASE",
                "Coloc with glomerlus ASE", "Coloc with kidney eQTL", "Open4Gene link",
                "Same Top Cell Type", "meQTL", "Coaccess", "ABC", "cS2G"]
feats = variant_feats + gene_feats + target_feats
res = {"missing_feature_columns": [c for c in feats if c not in gw.columns]}

num = gw[feats + ["Variant Score", "Gene Score", "Tartet Score", "Prioritization Score"]].apply(pd.to_numeric, errors="coerce")

inv = []
for c in feats:
    v = num[c]
    inv.append({"feature": c,
                "group": "variant" if c in variant_feats else ("gene" if c in gene_feats else "target"),
                "values": ",".join(sorted(map(str, v.dropna().unique().astype(int)))),
                "missing_pct": round(100 * v.isna().mean(), 2),
                "pct_one": round(100 * (v == 1).mean(), 2)})
pd.DataFrame(inv).to_csv(f"{OUT}/gw_feature_inventory.csv", index=False)
res["feature_inventory"] = inv

res["score_reconstruction_match_fraction"] = {
    "variant": round(float((num[variant_feats].sum(axis=1) == num["Variant Score"]).mean()), 4),
    "gene": round(float((num[gene_feats].sum(axis=1) == num["Gene Score"]).mean()), 4),
    "target": round(float((num[target_feats].sum(axis=1) == num["Tartet Score"]).mean()), 4),
    "total": round(float((num[feats].sum(axis=1) == num["Prioritization Score"]).mean()), 4),
}
res["rows_variant_gene_pairs"] = int(len(gw))
res["unique_variants"] = int(gw["Variant RSID"].nunique())
res["unique_genes"] = int(gw["Target gene symbol"].nunique())
res["unique_loci"] = int(gw["GWAS loci (Lead variant)"].nunique())
res["prioritization_score_range"] = [float(num["Prioritization Score"].min()), float(num["Prioritization Score"].max())]
res["pct_pairs_target_score_0"] = round(100 * float((num["Tartet Score"] == 0).mean()), 1)
gpv = gw.groupby("Variant RSID")["Target gene symbol"].nunique()
res["genes_per_variant"] = {k: round(float(v), 2) for k, v in gpv.describe().items()}
res["pct_variants_linked_to_>=2_genes"] = round(100 * float((gpv >= 2).mean()), 1)

# eGFRcys/BUN: are the 4 biomarker columns separable and non-constant?
res["biomarker_columns_pct_one"] = {c: round(100 * float((num[c] == 1).mean()), 1)
                                    for c in ["eGFRcys", "BUN", "Gene CDS eGFRcys", "Gene CDS BUN"]}

# ---- link to multi-ancestry credible-set signals (S8) via rsID (same paper, same dbSNP) ----
s8 = pd.read_excel(SUPP, sheet_name="S8 Finemap.Multi", header=1, dtype=str)
s8 = s8[s8["RSID"].notna() & s8["Merged GWAS loci"].astype(str).str.startswith("chr")].copy()
s8["signal_id"] = s8["Merged GWAS loci"] + "|" + s8["CARMA Credible Set"].astype(str)
v2sig = s8.groupby("RSID")["signal_id"].agg(lambda x: ";".join(sorted(set(x))))
gw["signal_ids"] = gw["Variant RSID"].map(v2sig)
total_signals = s8["signal_id"].nunique()
res["multi_cs_signals_total"] = int(total_signals)
res["pct_pairs_with_cs_signal"] = round(100 * float(gw["signal_ids"].notna().mean()), 2)

lk = gw[gw["signal_ids"].notna()].copy()
lk = lk.assign(signal_id=lk["signal_ids"].str.split(";")).explode("signal_id")
lk["locus"] = lk["signal_id"].str.split("|").str[0]
lk_t = lk.assign(tscore=num.loc[lk.index, "Tartet Score"].values)
sig_any = lk.groupby("signal_id")["Target gene symbol"].nunique()
sig_tgt = lk_t[lk_t["tscore"] > 0].groupby("signal_id")["Target gene symbol"].nunique()
res["signals_with_>=1_gene_pair"] = int(len(sig_any))
res["pct_signals_with_>=1_gene_pair"] = round(100 * len(sig_any) / total_signals, 1)
res["signals_with_>=1_gene_link_target_score>0"] = int(len(sig_tgt))
res["pct_signals_with_link_target_score>0"] = round(100 * len(sig_tgt) / total_signals, 1)

# within-locus eligibility (Plan-5): >=2 signals, >=2 linked genes, every signal has >=1 link (target score > 0)
lt = lk_t[lk_t["tscore"] > 0]
per = lt.groupby("locus").agg(n_sig_linked=("signal_id", "nunique"), n_gene=("Target gene symbol", "nunique"))
tot = s8.groupby("Merged GWAS loci")["signal_id"].nunique()
per["n_sig_total"] = tot.reindex(per.index)
strict = per[(per["n_sig_linked"] >= 2) & (per["n_gene"] >= 2) & (per["n_sig_linked"] == per["n_sig_total"])]
relaxed = per[(per["n_sig_linked"] >= 2) & (per["n_gene"] >= 2)]
res["within_locus_eligible_loci_strict"] = int(len(strict))
res["within_locus_eligible_loci_relaxed"] = int(len(relaxed))
res["loci_with_>=2_signals_total"] = int((tot >= 2).sum())
per.to_csv(f"{OUT}/gw_locus_signal_gene_counts.csv")

# signals whose linked genes differ (needed for specificity to change within-locus order)
sig_gene_sets = lt.groupby(["locus", "signal_id"])["Target gene symbol"].agg(frozenset)
differ = 0
for locus, grp in sig_gene_sets.groupby(level=0):
    sets = list(grp.values)
    if len(sets) >= 2 and len(set(sets)) >= 2:
        differ += 1
res["loci_where_signals_link_to_different_gene_sets"] = int(differ)

with open(f"{OUT}/audit_gw_summary.json", "w") as fh:
    json.dump(res, fh, indent=2, default=str)
print(json.dumps({k: v for k, v in res.items() if k != "feature_inventory"}, indent=2, default=str))
