"""Feasibility audit, part 1: Liu et al. 2025 supplementary tables only.

Checks data structure, coverage and implementability. It never loads any
benchmark gene list beyond Liu's own S22 (needed to document leakage) and
never touches FinnGen or any outcome data.
"""
import sys
import json
import pandas as pd
import numpy as np

XLSX = sys.argv[1]
OUT = sys.argv[2]

def load(sheet, header=1):
    df = pd.read_excel(XLSX, sheet_name=sheet, header=header, dtype=str)
    df = df.loc[:, ~df.columns.astype(str).str.startswith("Unnamed")]
    return df

res = {}

# ---------- S23: regulatory scorecard (variant-gene pairs, score >= 10) ----------
s23 = load("S23 Reg.Variants")
s23 = s23[s23["Variant RSID"].notna()].copy()
variant_feats = ["eGFRcys", "BUN", "TSS1Kb", "CredibleSet", "snATACPeak", "bASA", "snASA", "HMMReg"]
gene_feats = ["Gene CDS", "Gene CDS pLoF", "Gene CDS Missense", "Gene CDS eGFRcys", "Gene CDS BUN",
              "Gene CDS CredibleSet", "Gene ExWAS", "Gene Monogenic", "Gene MousePhynotype",
              "Gene PheWAS", "Gene MAGMA", "Gene DEPICT"]
target_feats = ["Tubular ASE", "Glomerular ASE", "eQTL", "Coloc with tubular ASE",
                "Coloc with glomerular ASE", "Coloc with kidney eQTL", "Open4Gene link",
                "Same Top Cell Type", "meQTL", "Coaccess", "ABC", "cS2G"]
all_feats = variant_feats + gene_feats + target_feats
missing_cols = [c for c in all_feats if c not in s23.columns]
res["s23_missing_feature_columns"] = missing_cols

num = s23[all_feats + ["Variant Score", "Gene Score", "Tartet Score", "Prioritization Score"]].apply(pd.to_numeric, errors="coerce")
feat_inventory = []
for c in all_feats:
    v = num[c]
    feat_inventory.append({
        "feature": c,
        "group": "variant" if c in variant_feats else ("gene" if c in gene_feats else "target"),
        "values": ",".join(sorted(map(str, v.dropna().unique().astype(int)))),
        "missing_pct": round(100 * v.isna().mean(), 2),
        "pct_one": round(100 * (v == 1).mean(), 1),
    })
pd.DataFrame(feat_inventory).to_csv(f"{OUT}/feature_inventory.csv", index=False)

# score reconstruction
recon = {
    "variant": (num[variant_feats].sum(axis=1) == num["Variant Score"]).mean(),
    "gene": (num[gene_feats].sum(axis=1) == num["Gene Score"]).mean(),
    "target": (num[target_feats].sum(axis=1) == num["Tartet Score"]).mean(),
    "total": (num[all_feats].sum(axis=1) == num["Prioritization Score"]).mean(),
}
res["score_reconstruction_match_fraction"] = {k: round(float(v), 4) for k, v in recon.items()}
res["s23_rows"] = int(len(s23))
res["s23_unique_variants"] = int(s23["Variant RSID"].nunique())
res["s23_unique_genes"] = int(s23["Target gene symbol"].nunique())
res["s23_unique_loci"] = int(s23["GWAS loci (Lead variant)"].nunique())
res["s23_min_prioritization_score"] = float(num["Prioritization Score"].min())
genes_per_variant = s23.groupby("Variant RSID")["Target gene symbol"].nunique()
res["s23_genes_per_variant"] = genes_per_variant.describe().round(2).to_dict()
res["s23_pct_variants_linked_to_>=2_genes"] = round(100 * float((genes_per_variant >= 2).mean()), 1)
res["s23_target_score_distribution"] = num["Tartet Score"].value_counts().sort_index().astype(int).to_dict()

# ---------- S4: loci ----------
s4 = load("S4 Ind.Loci")
s4 = s4[s4["Lead variant"].notna() & s4["Independent loci"].astype(str).str.startswith("chr")].copy()
res["s4_loci"] = int(len(s4))
loc = s4["Independent loci"].str.extract(r"chr(\w+):(\d+)-(\d+)")
s4["chr"] = loc[0]; s4["start"] = loc[1].astype(int); s4["end"] = loc[2].astype(int)
s4["pos"] = pd.to_numeric(s4["Position of lead variant"], errors="coerce")

# ---------- S7/S8: credible sets ----------
def load_cs(sheet):
    d = load(sheet)
    d = d[d["RSID"].notna() & d["Merged GWAS loci"].astype(str).str.startswith("chr")].copy()
    m = d["Merged GWAS loci"].str.extract(r"chr(\w+)_(\d+)_(\d+)")
    d["chr"] = m[0]; d["mstart"] = m[1].astype(int); d["mend"] = m[2].astype(int)
    d["pos"] = pd.to_numeric(d["Position"], errors="coerce").astype("Int64")
    d["pip"] = pd.to_numeric(d["CARMA posterior inclusion probability (PIP)"], errors="coerce")
    d["cs"] = d["CARMA Credible Set"].astype(str)
    d["signal_id"] = d["Merged GWAS loci"] + "|" + d["cs"]
    return d

cs_out = {}
for sheet, tag in [("S8 Finemap.Multi", "multi"), ("S7 Finemap.EUR", "eur")]:
    d = load_cs(sheet)
    sig_per_locus = d.groupby("Merged GWAS loci")["signal_id"].nunique()
    vars_per_sig = d.groupby("signal_id")["RSID"].nunique()
    cs_out[tag] = {
        "variants": int(len(d)),
        "merged_loci": int(d["Merged GWAS loci"].nunique()),
        "signals": int(d["signal_id"].nunique()),
        "signals_per_locus": sig_per_locus.value_counts().sort_index().astype(int).to_dict(),
        "loci_with_>=2_signals": int((sig_per_locus >= 2).sum()),
        "variants_per_signal_median": float(vars_per_sig.median()),
        "single_variant_signals": int((vars_per_sig == 1).sum()),
        "cs_labels_sample": sorted(d["cs"].unique())[:8],
    }
    if tag == "multi":
        cs_multi = d
    else:
        cs_eur = d
res["credible_sets"] = cs_out

# ---------- Crosswalk: S4 lead variant -> multi-ancestry credible set ----------
def crosswalk(cs, label):
    rows = []
    cs_key = cs.assign(key=cs["chr"] + ":" + cs["pos"].astype(str) + ":" + cs["Allele 1 (Effect allele)"].str.upper() + ":" + cs["Allele 2"].str.upper())
    keyset = {}
    for k, sid in zip(cs_key["key"], cs_key["signal_id"]):
        keyset.setdefault(k, set()).add(sid)
    for _, r in s4.iterrows():
        a1 = str(r["Alternative allele (Effect allele)"]).upper(); a2 = str(r["Reference allele"]).upper()
        if pd.isna(r["pos"]):
            rows.append((r["Lead variant"], "build_allele_mismatch", "")); continue
        k1 = f'{r["chr"]}:{int(r["pos"])}:{a1}:{a2}'; k2 = f'{r["chr"]}:{int(r["pos"])}:{a2}:{a1}'
        hits = keyset.get(k1, set()) | keyset.get(k2, set())
        if hits:
            rows.append((r["Lead variant"], "matched_in_cs" if len(hits) == 1 else "ambiguous", ";".join(sorted(hits))))
            continue
        # position present in CS table but alleles differ?
        same_pos = cs[(cs["chr"] == r["chr"]) & (cs["pos"] == int(r["pos"]))]
        if len(same_pos):
            rows.append((r["Lead variant"], "build_allele_mismatch", "")); continue
        near = cs[(cs["chr"] == r["chr"]) & ((cs["pos"] - int(r["pos"])).abs() <= 500_000)]
        if len(near) == 0:
            rows.append((r["Lead variant"], "no_credible_set_within_500kb", ""))
        else:
            # LD not computable in this environment (no LD reference yet)
            rows.append((r["Lead variant"], "near_cs_needs_LD_check", ";".join(sorted(near["signal_id"].unique()))))
    cw = pd.DataFrame(rows, columns=["lead_variant", "status", "signal_ids"])
    cw.to_csv(f"{OUT}/crosswalk_s4_to_{label}.csv", index=False)
    return cw["status"].value_counts().astype(int).to_dict()

res["crosswalk_s4_to_multi_cs"] = crosswalk(cs_multi, "multi")
res["crosswalk_s4_to_eur_cs"] = crosswalk(cs_eur, "eur")

# ---------- Signal-to-gene coverage (S23 is truncated at score >= 10) ----------
v2sig = cs_multi.groupby("RSID")["signal_id"].agg(lambda x: ";".join(sorted(set(x))))
s23["signal_ids"] = s23["Variant RSID"].map(v2sig)
res["s23_pct_rows_with_multi_cs_signal"] = round(100 * float(s23["signal_ids"].notna().mean()), 1)
res["s23_pct_rows_flagged_CredibleSet"] = round(100 * float((num["CredibleSet"] == 1).mean()), 1)
linked = s23[s23["signal_ids"].notna()].copy()
linked = linked.assign(signal_id=linked["signal_ids"].str.split(";")).explode("signal_id")
linked["locus"] = linked["signal_id"].str.split("|").str[0]
sig_genes = linked.groupby("signal_id")["Target gene symbol"].nunique()
res["signals_with_>=1_gene_link_in_S23"] = int(len(sig_genes))
res["pct_multi_signals_with_link_in_S23"] = round(100 * len(sig_genes) / cs_out["multi"]["signals"], 1)
per_locus = linked.groupby("locus").agg(n_sig=("signal_id", "nunique"), n_gene=("Target gene symbol", "nunique"))
all_sig_per_locus = cs_multi.groupby("Merged GWAS loci")["signal_id"].nunique()
per_locus["n_sig_total"] = all_sig_per_locus.reindex(per_locus.index)
eligible = per_locus[(per_locus["n_sig"] >= 2) & (per_locus["n_gene"] >= 2) & (per_locus["n_sig"] == per_locus["n_sig_total"])]
res["within_locus_eligible_loci_lower_bound_S23"] = int(len(eligible))
res["within_locus_eligible_loci_relaxed_(>=2 linked signals, >=2 genes)"] = int(((per_locus["n_sig"] >= 2) & (per_locus["n_gene"] >= 2)).sum())

# ---------- Leakage: S22 nephropathy genes used to set the score cutoff ----------
s22 = load("S22 Neph.Genes")
s22 = s22[s22["Gene"].notna() & (s22["Gene"] != "Gene")]
res["s22_nephropathy_genes_used_for_cutoff"] = sorted(s22["Gene"].astype(str).str.strip().unique().tolist())

with open(f"{OUT}/audit_supp_summary.json", "w") as fh:
    json.dump(res, fh, indent=2, default=str)
print(json.dumps(res, indent=2, default=str))
