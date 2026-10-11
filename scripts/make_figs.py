import json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#5f5e5a", "axes.labelcolor": "#2c2c2a", "xtick.color": "#5f5e5a", "ytick.color": "#5f5e5a"})
BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#2c2c2a", "#8a8980", "#e4e3dd"
def save(fig, name):
    fig.savefig(f"figs/{name}.png", dpi=300, bbox_inches="tight"); fig.savefig(f"figs/{name}.tif", dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"}); plt.close(fig)

# ---- Figure 1: f_s distribution and Titan check ----
F = pd.read_csv("out/aim1/f_s_signals.csv"); T = pd.read_csv("out/aim1/titan_check.csv")
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1.4, 1]})
v = F.f_s.dropna(); bins = np.linspace(0, 1, 21)
ax[0].hist(v, bins=bins, color=BLUE, edgecolor="white", linewidth=1.0)
ax[0].axvline(0.8, color=MUTED, ls="--", lw=1); ax[0].axvline(0.2, color=MUTED, ls="--", lw=1)
ax[0].set_xlabel("Filtration-support score (f$_s$)"); ax[0].set_ylabel("Signals")
ax[0].set_title(f"a  Signal-level filtration support (n = {len(v):,})", loc="left", fontsize=9, color=INK)
ax[0].text(0.81, ax[0].get_ylim()[1]*0.85, f"≥ 0.8: {(v>=0.8).sum()}", color=INK, fontsize=8)
ax[0].text(0.21, ax[0].get_ylim()[1]*0.85, f"< 0.2: {(v<0.2).sum():,}", color=INK, fontsize=8)
ax[0].yaxis.grid(True, color=GRID, lw=0.6); ax[0].set_axisbelow(True)
rng = np.random.default_rng(1)
for i, (g, col, lab) in enumerate([("gfr", BLUE, "GFR-associated"), ("creat", ORANGE, "Creatinine-specific")]):
    y = T[T.grp == g].f_s.dropna().values
    ax[1].boxplot(y, positions=[i], widths=0.5, showfliers=False, medianprops={"color": INK, "lw": 1.5},
                  boxprops={"color": MUTED}, whiskerprops={"color": MUTED}, capprops={"color": MUTED})
    ax[1].scatter(i + rng.uniform(-0.15, 0.15, len(y)), y, s=12, color=col, alpha=0.75, edgecolor="white", linewidth=0.4, zorder=3)
ax[1].set_xticks([0, 1]); ax[1].set_xticklabels([f"GFR-\nassociated\n(n = {T[T.grp=='gfr'].f_s.notna().sum()})", f"Creatinine-\nspecific\n(n = {T[T.grp=='creat'].f_s.notna().sum()})"])
ax[1].set_ylabel("f$_s$"); ax[1].set_ylim(-0.05, 1.12)
ax[1].set_title("b  Titan et al. loci", loc="left", fontsize=9, color=INK)
ax[1].text(0.5, 1.06, "Wilcoxon P = 3.8 × 10$^{-6}$", ha="center", fontsize=8, color=INK)
save(fig, "Fig2_fs_titan")

# ---- Figure 2: incremental value forest plot ----
r = json.load(open("out/aim1/aim2_results.json")); s = json.load(open("out/aim1/aim2_secondary.json"))
b = json.load(open("out/aim1/aim2_results_bun.json")); m = json.load(open("out/aim1/aim2_mgi.json")); ic = json.load(open("out/aim1/aim2_innercv.json"))["inner_cv_lambda"]
mg = m["MGI renal function (mouse feature removed)"]
rows = [("Primary (ClinGen + PanelApp)", r["primary"]["diff_C_minus_A"], r["primary"]["diff_C_minus_A_CI95"], "Primary"),
        ("Inner-CV λ selection (post hoc)", ic["diff"], ic["CI95"], "Sensitivity"),
        ("Directly matched signals only", r["sens_direct_matched_only"]["diff_C_minus_A"], r["sens_direct_matched_only"]["diff_C_minus_A_CI95"], "Sensitivity"),
        ("Maximum aggregation", r["sens_max_aggregation"]["diff_C_minus_A"], r["sens_max_aggregation"]["diff_C_minus_A_CI95"], "Sensitivity"),
        ("Equal link weights", r["sens_equal_weights"]["diff_C_minus_A"], r["sens_equal_weights"]["diff_C_minus_A_CI95"], "Sensitivity"),
        ("BUN-based f$_s$", b["primary"]["diff_C_minus_A"], b["primary"]["diff_C_minus_A_CI95"], "Sensitivity"),
        ("Extended benchmark (+Moderate, +amber)", r["sens_extended_benchmark"]["diff_C_minus_A"], r["sens_extended_benchmark"]["diff_C_minus_A_CI95"], "Benchmark"),
        ("Groopman 2019 (625 genes)", r["secondary_groopman"]["diff_C_minus_A"], r["secondary_groopman"]["diff_C_minus_A_CI95"], "Benchmark"),
        ("MGI renal-function phenotypes", mg["diff_C_minus_A"], mg["diff_CI95"], "Benchmark"),
        ("L1 logistic, with vs without S$_g$", s["L1_logistic"]["diff"], s["L1_logistic"]["diff_CI95"], "Learned"),
        ("XGBoost, with vs without S$_g$", s["XGBoost"]["diff"], s["XGBoost"]["diff_CI95"], "Learned")]
fig, ax = plt.subplots(figsize=(6.8, 4.2))
for i, (lab, d, ci, grp) in enumerate(rows[::-1]):
    col = BLUE if grp == "Primary" else INK
    ax.plot(ci, [i, i], color=col, lw=1.6); ax.scatter([d], [i], s=36 if grp == "Primary" else 22, color=col, zorder=3)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([x[0] for x in rows[::-1]])
ax.axvline(0, color=MUTED, lw=1); ax.axvline(0.03, color=ORANGE, ls="--", lw=1.2)
ax.text(0.0305, len(rows) - 0.4, "Prespecified\nmeaningful gain (0.03)", color=INK, fontsize=7.5, va="top")
ax.set_xlabel("Held-out AUPRC difference, with minus without filtration support (95% CI)")
ax.xaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True); ax.set_xlim(-0.02, 0.05)
save(fig, "Fig3_incremental_value")

# ---- Figure 3: FinnGen outcome enrichment ----
a = json.load(open("out/aim3_results.json"))
eps = [("N14_DIALYSIS", "Dialysis\n(primary)"), ("Z21_KIDNEY_TRANSPLANT_STATUS", "Kidney\ntransplant"),
       ("N14_CHRONKIDNEYDIS", "CKD"), ("N14_ACUTERENFAIL", "AKI"), ("H7_VITROPACIT", "Negative control\n(vitreous opacities)")]
fig, ax = plt.subplots(figsize=(7.2, 3.4))
x = np.arange(len(eps)); w = 0.36
f = [100 * a[e]["prop_hit_filtration"] for e, _ in eps]; l = [100 * a[e]["prop_hit_low"] for e, _ in eps]
ax.bar(x - w/2 - 0.01, f, w, color=BLUE, label=f"Filtration-consistent (f$_s$ ≥ 0.8; n = {a['N14_DIALYSIS']['n_filtration']})")
ax.bar(x + w/2 + 0.01, l, w, color=ORANGE, label=f"Low filtration support (f$_s$ < 0.2; n = {a['N14_DIALYSIS']['n_low']})")
for i, (e, _) in enumerate(eps):
    p = a[e]["perm_p_one_sided"]; ax.text(i, max(f[i], l[i]) + 0.8, (f"P = {p:.3f}" if p < 0.1 else f"P = {p:.2f}"), ha="center", fontsize=7.5, color=INK)
ax.axhline(2.5, color=MUTED, ls=":", lw=1)
ax.set_xticks(x); ax.set_xticklabels([n for _, n in eps]); ax.set_ylabel("Signals with P < 0.05\nin risk direction (%)")
ax.legend(frameon=False, fontsize=7.5, loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2); ax.yaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
ax.set_ylim(0, max(f + l) + 4)
save(fig, "Fig4_finngen")
print("ok")
