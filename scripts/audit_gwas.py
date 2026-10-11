"""Feasibility audit, part 3: filtered GWAS files (Liu eGFRcrea EUR, CKDGen eGFRcys, BUN).

Data readiness only: build agreement, allele harmonization, credible-set coverage
(the prespecified >= 80% quality rule), and numeric precision. No colocalization,
no f_s, no benchmark or outcome data.
"""
import sys, json
import numpy as np
import pandas as pd

SUPP, GDIR, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
res = {}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

def flip(a):
    return "".join(COMP.get(x, x) for x in a)

liu = pd.read_csv(f"{GDIR}/liu_crea.gz", sep="\t", usecols=["MarkerName", "CHR", "POS", "REF", "ALT", "MAF", "Zscore", "BETA", "SE", "N"], dtype={"CHR": str})
liu = liu.rename(columns={"MarkerName": "RSID"})
liu["REF"] = liu["REF"].str.upper(); liu["ALT"] = liu["ALT"].str.upper()

def load_ckd(fn):
    d = pd.read_csv(f"{GDIR}/{fn}", sep="\t", usecols=["Allele1", "Allele2", "n", "Freq1", "Effect", "StdErr", "P.value", "chr", "pos", "StdErr.GC", "P.value.GC", "RSID", "Direction"], dtype={"chr": str})
    d["Allele1"] = d["Allele1"].str.upper(); d["Allele2"] = d["Allele2"].str.upper()
    return d

cys = load_ckd("cys.gz"); bun = load_ckd("bun.gz")
res["rows_after_filter"] = {"liu_eGFRcrea": len(liu), "ckdgen_eGFRcys": len(cys), "ckdgen_BUN": len(bun)}

# ---- precision ----
for name, d in [("eGFRcys", cys), ("BUN", bun)]:
    eff_dec = d["Effect"].astype(str).str.replace("-", "")
    res[f"{name}_pct_effect_exactly_0"] = round(100 * float((d["Effect"] == 0).mean()), 2)
    res[f"{name}_pct_SE_one_significant_digit"] = round(100 * float(d["StdErr"].astype(str).str.fullmatch(r"\d(\.0+)?e-\d+|0\.0*\d").mean()), 1)
    res[f"{name}_pct_P_zero_or_missing"] = round(100 * float((d["P.value"].isna() | (d["P.value"] == 0)).mean()), 3)

# ---- credible sets (S8 multi-ancestry; S7 EUR) ----
def load_cs(sheet):
    d = pd.read_excel(SUPP, sheet_name=sheet, header=1, dtype=str)
    d = d[d["RSID"].notna() & d["Merged GWAS loci"].astype(str).str.startswith("chr")].copy()
    d["signal_id"] = d["Merged GWAS loci"] + "|" + d["CARMA Credible Set"].astype(str)
    d["pos"] = pd.to_numeric(d["Position"], errors="coerce")
    d["chr"] = d["Chromosome"].astype(str)
    d["A1"] = d["Allele 1 (Effect allele)"].str.upper(); d["A2"] = d["Allele 2"].str.upper()
    return d

cs8 = load_cs("S8 Finemap.Multi"); cs7 = load_cs("S7 Finemap.EUR")

# ---- build agreement: positions by rsID ----
def build_check(cs, d, poscol, chrcol, label):
    dd = d[["RSID", chrcol, poscol]].drop_duplicates("RSID").rename(columns={chrcol: "d_chr", poscol: "d_pos"})
    m = cs.merge(dd, on="RSID", how="inner")
    same = (m["pos"] == m["d_pos"]) & (m["chr"] == m["d_chr"].astype(str))
    return {"rsids_compared": int(len(m)), "pct_same_chr_pos": round(100 * float(same.mean()), 2)}

res["build_check_vs_S8_positions"] = {
    "liu_eGFRcrea": build_check(cs8, liu, "POS", "CHR", "liu"),
    "ckdgen_eGFRcys": build_check(cs8, cys, "pos", "chr", "cys"),
    "ckdgen_BUN": build_check(cs8, bun, "pos", "chr", "bun"),
}

# ---- harmonization: match CKDGen to Liu by chr:pos + alleles ----
def harmonize(ckd, label):
    a = liu[["CHR", "POS", "REF", "ALT", "MAF"]].rename(columns={"CHR": "chr", "POS": "pos"})
    b = ckd[["chr", "pos", "Allele1", "Allele2"]]
    m = a.merge(b, on=["chr", "pos"], how="left")
    pal = (m["REF"].isin(["A", "T"]) & m["ALT"].isin(["A", "T"]) | m["REF"].isin(["C", "G"]) & m["ALT"].isin(["C", "G"])) & m["MAF"].between(0.42, 0.58)  # prespecified: drop only ambiguous palindromic SNPs
    same = (m["Allele1"] == m["ALT"]) & (m["Allele2"] == m["REF"])
    swap = (m["Allele1"] == m["REF"]) & (m["Allele2"] == m["ALT"])
    strand = (m["Allele1"].fillna("").map(flip) == m["ALT"]) & (m["Allele2"].fillna("").map(flip) == m["REF"]) | \
             (m["Allele1"].fillna("").map(flip) == m["REF"]) & (m["Allele2"].fillna("").map(flip) == m["ALT"])
    present = m["Allele1"].notna()
    status = np.select([~present, (same | swap) & ~pal, (same | swap) & pal, strand & ~pal],
                       ["absent", "aligned", "palindromic_dropped", "strand_flip_aligned"], "allele_mismatch")
    m["status"] = status
    m = m.drop_duplicates(["chr", "pos", "REF", "ALT"])
    return m, m["status"].value_counts(normalize=True).mul(100).round(2).to_dict()

mc, res["harmonization_liu_vs_eGFRcys_pct"] = harmonize(cys, "cys")
mb, res["harmonization_liu_vs_BUN_pct"] = harmonize(bun, "bun")

# ---- credible-set coverage (prespecified quality rule: >= 80% of CS variants usable in eGFRcys) ----
def cs_coverage(cs, m, label):
    usable = m[m["status"].isin(["aligned", "strand_flip_aligned"])][["chr", "pos"]].drop_duplicates()
    usable["ok"] = True
    c = cs.merge(usable, on=["chr", "pos"], how="left")
    c["ok"] = c["ok"].fillna(False).astype(bool)
    in_liu = cs.merge(liu[["CHR", "POS"]].rename(columns={"CHR": "chr", "POS": "pos"}).drop_duplicates(), on=["chr", "pos"], how="inner")
    per = c.groupby("signal_id")["ok"].mean()
    return {
        "signals": int(per.size),
        "signals_meeting_80pct_rule": int((per >= 0.8).sum()),
        "pct_signals_meeting_80pct_rule": round(100 * float((per >= 0.8).mean()), 1),
        "cs_variants_present_in_liu_EUR_file_pct": round(100 * len(in_liu) / len(cs), 1),
    }

res["cs_coverage_S8_multi_vs_eGFRcys"] = cs_coverage(cs8, mc, "cys")
res["cs_coverage_S8_multi_vs_BUN"] = cs_coverage(cs8, mb, "bun")
res["cs_coverage_S7_eur_vs_eGFRcys"] = cs_coverage(cs7, mc, "cys")

# ---- z-score recoverability from P + sign (workaround for rounded Effect/SE) ----
for name, d in [("eGFRcys", cys), ("BUN", bun)]:
    sign_known = (d["Effect"] != 0) | d["Direction"].astype(str).str.fullmatch(r"\++|-+")
    res[f"{name}_pct_z_recoverable_from_P_and_sign"] = round(100 * float(sign_known.mean()), 2)

with open(f"{OUT}/audit_gwas_summary.json", "w") as fh:
    json.dump(res, fh, indent=2, default=str)
print(json.dumps(res, indent=2, default=str))
