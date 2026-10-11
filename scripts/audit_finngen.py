"""Feasibility audit, part 5: FinnGen R13 filtered files - structure only.
Reads only variant identifiers, alleles and allele frequency; counts missing beta/SE/P
without printing or summarizing their values. No association results are examined."""
import sys, json
import numpy as np, pandas as pd
SUPP, GDIR, FDIR, OUT = sys.argv[1:5]
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
flip = lambda s: s.fillna("").map(lambda a: "".join(COMP.get(x, x) for x in a))
EPS = ["N14_DIALYSIS", "Z21_KIDNEY_TRANSPLANT_STATUS", "N14_CHRONKIDNEYDIS", "N14_ACUTERENFAIL", "H7_VITROPACIT"]
liu = pd.read_csv(f"{GDIR}/liu_crea.gz", sep="\t", usecols=["MarkerName", "REF", "ALT", "MAF"]).rename(columns={"MarkerName": "rsid"})
liu["REF"] = liu.REF.str.upper(); liu["ALT"] = liu.ALT.str.upper()
liu = liu.drop_duplicates("rsid")
s8 = pd.read_excel(SUPP, sheet_name="S8 Finemap.Multi", header=1, dtype=str)
s8 = s8[s8["RSID"].notna() & s8["Merged GWAS loci"].astype(str).str.startswith("chr")].copy()
s8["signal_id"] = s8["Merged GWAS loci"] + "|" + s8["CARMA Credible Set"].astype(str)
res = {}; sets = {}
for ep in EPS:
    f = pd.read_csv(f"{FDIR}/fg_{ep}_kidneyloci_tsv.gz", sep="\t", dtype={"#chrom": str},
                    usecols=["#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"])
    r = {"rows": int(len(f)),
         "pct_missing_beta_se_p": round(100 * float(f[["beta", "sebeta", "pval"]].isna().any(axis=1).mean()), 3),
         "multi_rsid_rows": int(f.rsids.str.contains(",", na=False).sum())}
    f = f.drop(columns=["pval", "beta", "sebeta"])
    f = f.assign(rsid=f.rsids.str.split(",")).explode("rsid")
    f = f[f.rsid.isin(liu.rsid)]
    r["duplicated_rsid_rows"] = int(f.rsid.duplicated().sum())
    f = f.drop_duplicates("rsid")
    m = f.merge(liu, on="rsid")
    pal = ((m.REF.isin(["A", "T"]) & m.ALT.isin(["A", "T"])) | (m.REF.isin(["C", "G"]) & m.ALT.isin(["C", "G"]))) & m.MAF.between(0.42, 0.58)
    same = (m.alt == m.ALT) & (m.ref == m.REF); swap = (m.alt == m.REF) & (m.ref == m.ALT)
    strand = ((flip(m.alt) == m.ALT) & (flip(m.ref) == m.REF)) | ((flip(m.alt) == m.REF) & (flip(m.ref) == m.ALT))
    m["status"] = np.select([(same | swap) & ~pal, (same | swap) & pal, strand & ~pal], ["aligned", "palindromic_dropped", "strand_flip_aligned"], "allele_mismatch")
    r["rsid_matched_to_liu"] = int(len(m))
    r["pct_of_liu_window_rsids"] = round(100 * len(m) / len(liu), 1)
    r["harmonization_pct"] = m.status.value_counts(normalize=True).mul(100).round(2).to_dict()
    # allele-frequency sanity on aligned SNPs (EUR vs Finnish; large gaps flag misalignment)
    a = m[m.status == "aligned"].copy()
    liu_af_alt = np.where(a.alt == a.ALT, np.nan, np.nan)  # Liu gives MAF only; compare minor-allele freq
    fg_maf = np.minimum(a.af_alt, 1 - a.af_alt)
    r["pct_aligned_maf_diff_gt_0.2"] = round(100 * float((np.abs(fg_maf - a.MAF) > 0.2).mean()), 2)
    usable = set(m.loc[m.status.isin(["aligned", "strand_flip_aligned"]), "rsid"])
    per = s8.assign(ok=s8.RSID.isin(usable)).groupby("signal_id")["ok"].mean()
    r["S8_signals"] = int(per.size)
    r["S8_signals_meeting_80pct_rule"] = int((per >= 0.8).sum())
    r["pct_S8_signals_meeting_80pct_rule"] = round(100 * float((per >= 0.8).mean()), 1)
    r["S8_signals_with_any_usable_variant"] = int((per > 0).sum())
    res[ep] = r; sets[ep] = usable
base = sets[EPS[0]]
res["usable_variant_set_identical_across_endpoints"] = {ep: int(len(sets[ep] ^ base)) for ep in EPS[1:]}
json.dump(res, open(f"{OUT}/audit_finngen_summary.json", "w"), indent=2)
print(json.dumps(res, indent=2))
