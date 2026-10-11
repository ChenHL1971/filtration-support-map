"""Figure: genome-wide filtration-support map (a) and predefined illustrative loci (b)."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#2c2c2a", "#8a8980", "#e4e3dd"
GREY = "#a9a8a0"
d = pd.read_csv("manuscript/Supplementary_Data_1_filtration_support_map.csv")
d["chr"] = pd.to_numeric(d.chromosome, errors="coerce"); d = d[d.chr.notna()].copy(); d["chr"] = d.chr.astype(int)
L = d.groupby("chr").position_GRCh37.max(); off = L.cumsum().shift(fill_value=0); gap = 2e7
d["x"] = d.position_GRCh37 + d.chr.map(off) + (d.chr - 1) * gap
col = {"Filtration-consistent": BLUE, "Intermediate": GREY, "Low support": ORANGE}
fig = plt.figure(figsize=(7.2, 6.2))
ax = fig.add_axes([0.08, 0.58, 0.9, 0.33])
for c in range(1, 23):
    if c % 2 == 0:
        lo = off[c] + (c - 1) * gap - gap / 2; ax.axvspan(lo, lo + L[c] + gap, color="#f3f2ee", lw=0, zorder=0)
for k in ["Low support", "Intermediate", "Filtration-consistent"]:
    s = d[d.classification == k]
    ax.scatter(s.x, s.filtration_support_f_s, s=7, color=col[k], lw=0, alpha=0.85, zorder=2,
               label=f"{k} ({len(s):,})")
inc = d[d.classification == "Inconclusive"]
ax.scatter(inc.x, np.full(len(inc), -0.09), marker="|", s=14, color=MUTED, lw=0.6, zorder=2, label=f"Inconclusive ({len(inc)})")
ax.axhline(0.8, color=MUTED, ls="--", lw=0.8); ax.axhline(0.2, color=MUTED, ls="--", lw=0.8)
ax.set_xticks([off[c] + (c - 1) * gap + L[c] / 2 for c in range(1, 23)]); ax.set_xticklabels([str(c) if c < 20 or c % 2 == 0 else "" for c in range(1, 23)])
ax.set_xlim(d.x.min() - gap, d.x.max() + gap); ax.set_ylim(-0.15, 1.05); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_xlabel("Chromosome (lead-variant position, GRCh37)"); ax.set_ylabel("Filtration support (f$_s$)")
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
ax.set_title("a  Filtration-support map of 2,051 eGFRcrea signals", loc="left", fontsize=9, color=INK)
ax.legend(frameon=False, fontsize=7, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.16), handletextpad=0.3, columnspacing=1.2, markerscale=1.8)
# panel b
genes = [("UMOD", "Classic kidney function locus"), ("SHROOM3", "Classic kidney function locus"),
         ("CPS1", "Köttgen 2010: creatinine"), ("SLC22A2", "Köttgen 2010: creatinine"), ("TMEM60", "Köttgen 2010: creatinine"),
         ("WDR37", "Köttgen 2010: creatinine"), ("SLC6A13", "Köttgen 2010: creatinine"), ("WDR72", "Köttgen 2010: creatinine"),
         ("BCAS3", "Köttgen 2010: creatinine")]
bx = fig.add_axes([0.27, 0.07, 0.53, 0.38]); s = d.linked_candidate_genes_link_weight.fillna("")
rng = np.random.default_rng(1)
for i, (g, grp) in enumerate(genes):
    y = len(genes) - 1 - i
    m = d[s.str.contains("(?:^|; )" + g + r" \(")]
    sc = m[m.classification != "Inconclusive"]
    bx.scatter(sc.filtration_support_f_s, y + rng.uniform(-0.15, 0.15, len(sc)), s=28, color=[col[k] for k in sc.classification],
               edgecolor="white", lw=0.6, zorder=3)
    ni = (m.classification == "Inconclusive").sum()
    bx.text(1.06, y, f"{len(sc)} scored" + (f", {ni} inconclusive" if ni else ""), va="center", fontsize=6.8, color=MUTED)
bx.set_yticks(range(len(genes))); bx.set_yticklabels([g for g, _ in genes][::-1], fontstyle="italic")
bx.axvline(0.8, color=MUTED, ls="--", lw=0.8); bx.axvline(0.2, color=MUTED, ls="--", lw=0.8)
bx.axhspan(len(genes) - 2.5, len(genes) - 0.5, color="#f3f2ee", lw=0, zorder=0)
fig.text(0.075, 0.07 + 0.38 * (len(genes) - 1.5 + 0.5) / len(genes), "Classic\nfiltration loci", va="center", ha="center", fontsize=7, color=INK)
bx.text(-0.42, 3, "Köttgen 2010\nsuspected\ncreatinine\nloci", va="center", ha="center", fontsize=7, color=INK)
bx.set_xlim(-0.05, 1.05); bx.set_ylim(-0.5, len(genes) - 0.5); bx.set_xlabel("Filtration support (f$_s$) of signals linked to the gene")
for sp in ["top", "right"]: bx.spines[sp].set_visible(False)
bx.xaxis.grid(True, color=GRID, lw=0.6); bx.set_axisbelow(True)
fig.text(0.02, 0.475, "b  Illustrative loci (complete K\u00f6ttgen 2010 list; post hoc)", fontsize=9, color=INK)
for ext in ["png", "tif"]:
    fig.savefig(f"figs/FigNew_map.{ext}", dpi=300, facecolor="white", **({"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tif" else {}))
