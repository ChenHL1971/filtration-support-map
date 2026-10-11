"""Feasibility audit, part 4: classify LD crosswalk (prespecified rules).
LD-rescued: exactly one signal with r2 >= 0.8 to the lead; ambiguous: >1; else unmatched with reason."""
import sys, json
import pandas as pd
LD, CW, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
d = pd.read_csv(LD, dtype=str); d["note"] = d["note"].fillna("")
assert not d["note"].str.startswith("error").any()
d["r2"] = pd.to_numeric(d["r2"], errors="coerce")
rows = []
for lead, g in d.groupby("lead_variant"):
    if (g["note"] == "lead_not_in_1000G").all():
        rows.append((lead, "unmatched", "LD_not_computable_lead_not_in_1000G", "", float("nan"))); continue
    per = g.groupby("signal_id")["r2"].max()
    q = per[per >= 0.8]
    uncomputable_sigs = per[per.isna()].index.tolist()
    if len(q) == 1:
        rows.append((lead, "LD_rescued", "", q.index[0], float(q.iloc[0])))
    elif len(q) > 1:
        rows.append((lead, "ambiguous", "", ";".join(q.sort_values(ascending=False).index), float(q.max())))
    elif per.notna().sum() == 0:
        rows.append((lead, "unmatched", "LD_not_computable_cs_not_in_1000G", "", float("nan")))
    else:
        reason = "LD_below_0.8" + ("_some_signals_not_computable" if uncomputable_sigs else "")
        rows.append((lead, "unmatched", reason, per.idxmax(), float(per.max())))
ld = pd.DataFrame(rows, columns=["lead_variant", "status", "reason", "signal_ids", "max_r2"])
cw = pd.read_csv(CW, dtype=str)
base = cw[cw["status"] != "near_cs_needs_LD_check"].copy()
base["reason"] = ""; base["max_r2"] = float("nan")
base.loc[base["status"] == "matched_in_cs", "status"] = "exact_match"
base.loc[base["status"] == "no_credible_set_within_500kb", ["status", "reason"]] = ["unmatched", "no_credible_set_within_500kb"]
final = pd.concat([base, ld], ignore_index=True)
final.to_csv(f"{OUT}/crosswalk_s4_to_multi_final.csv", index=False)
res = {"leads_queried": int(len(ld)),
       "ld_status": ld["status"].value_counts().to_dict(),
       "ld_unmatched_reasons": ld[ld.status == "unmatched"]["reason"].value_counts().to_dict(),
       "rescued_r2_median": round(float(ld[ld.status == "LD_rescued"]["max_r2"].median()), 3),
       "below_0.8_max_r2_quartiles": ld[ld.reason.str.startswith("LD_below")]["max_r2"].quantile([.25, .5, .75]).round(3).to_dict(),
       "final_crosswalk": final["status"].value_counts().to_dict(),
       "final_unmatched_reasons": final[final.status == "unmatched"]["reason"].value_counts().to_dict(),
       "total_loci": int(len(final))}
json.dump(res, open(f"{OUT}/audit_ld_summary.json", "w"), indent=2)
print(json.dumps(res, indent=2))
