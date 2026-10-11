"""Supplementary Data 1: signal-level filtration-support map (all 2,051 Liu et al. S8 signals)."""
import numpy as np, pandas as pd
F = pd.read_csv("out/aim1/f_s_signals.csv")
L = pd.read_csv("out/aim3_signal_leads.csv", dtype={"chr": str})[["signal_id", "chr", "pos", "rsid"]]
P = pd.read_csv("out/aim1/f_s_signals_p12.csv")[["signal_id", "f_s_p12"]]
gw = pd.read_pickle("out/gw_cache.pkl"); gw = gw[gw["Variant RSID"].notna() & gw["Target gene symbol"].notna()]
gw["ts"] = pd.to_numeric(gw["Tartet Score"], errors="coerce").fillna(0)
s8 = pd.read_excel("raw/science.adp4753_tables_s1_to_s31.xlsx", sheet_name="S8 Finemap.Multi", header=1, dtype=str)
s8 = s8[s8["RSID"].notna() & s8["Merged GWAS loci"].astype(str).str.startswith("chr")]
s8["signal_id"] = s8["Merged GWAS loci"] + "|" + s8["CARMA Credible Set"].astype(str)
W = s8[["signal_id", "RSID"]].merge(gw[["Variant RSID", "Target gene symbol", "ts"]], left_on="RSID", right_on="Variant RSID")
W = W.groupby(["signal_id", "Target gene symbol"]).ts.sum().reset_index(); W = W[W.ts > 0]
W["w"] = W.ts / W.groupby("signal_id").ts.transform("sum")
genes = W.sort_values(["signal_id", "w"], ascending=[True, False]).groupby("signal_id").apply(
    lambda d: "; ".join(f"{g} ({w:.2f})" for g, w in zip(d["Target gene symbol"], d["w"])))
D = F.merge(L, on="signal_id", how="left").merge(P, on="signal_id", how="left")
D["linked_genes_weight"] = D.signal_id.map(genes).fillna("")
D["n_linked_genes"] = D.signal_id.map(W.groupby("signal_id").size()).fillna(0).astype(int)
def cls(r):
    if pd.isna(r.f_s): return "Inconclusive"
    return "Filtration-consistent" if r.f_s >= 0.8 else ("Low support" if r.f_s < 0.2 else "Intermediate")
D["classification"] = D.apply(cls, axis=1)
meth = {"overlap": "coloc-SuSiE (direct match: posterior overlap)", "LD": "coloc-SuSiE (direct match: lead-variant LD)",
        "fallback_abf": "coloc.abf within credible set (fallback)", "missing_coverage": "Not estimated"}
D["estimation_method"] = D.method.map(meth)
def reason(r):
    if r.classification != "Inconclusive": return ""
    if r.method == "missing_coverage": return "<80% of credible-set variants available"
    if r.note == "coverage < 0.8": return "<80% of credible-set variants available after mapping"
    return "Effect direction could not be determined"
D["inconclusive_reason"] = D.apply(reason, axis=1)
D["direction_concordant"] = D.same_dir.map({True: "Yes", False: "No (f_s set to 0)"}).fillna("Not determined")
D.loc[(D.method == "missing_coverage"), "direction_concordant"] = "Not estimated"
D.loc[D.f_s.isna() & (D.same_dir == False), "direction_concordant"] = "No (signal not scored)"
D.loc[D.f_s.isna() & (D.same_dir == True), "direction_concordant"] = "Yes (signal not scored)"
out = pd.DataFrame({"locus_id": D.locus, "signal_id": D.signal_id, "lead_variant_rsid": D.rsid, "chromosome": D.chr,
    "position_GRCh37": D.pos.astype("Int64"), "credible_set_size": D.cs_size, "credible_set_coverage": D.cs_coverage.round(3),
    "estimation_method": D.estimation_method, "PP_H4": D.pp4.round(4), "direction_concordant": D.direction_concordant,
    "filtration_support_f_s": D.f_s.round(4), "f_s_fallback_prior_sensitivity": D.f_s_p12.round(4),
    "classification": D.classification, "inconclusive_reason": D.inconclusive_reason,
    "n_linked_genes": D.n_linked_genes, "linked_candidate_genes_link_weight": D.linked_genes_weight})
out = out.sort_values(["chromosome", "position_GRCh37"], key=lambda s: pd.to_numeric(s, errors="coerce"))
out.to_csv("manuscript/Supplementary_Data_1_filtration_support_map.csv", index=False)
print(out.classification.value_counts()); print(out.isna().sum()[out.isna().sum() > 0]); print(len(out))
