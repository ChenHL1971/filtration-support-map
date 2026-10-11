import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams.update({"font.family": "DejaVu Sans"})
BLUE, ORANGE, INK, MUTED, GREY = "#2a78d6", "#eb6834", "#2c2c2a", "#8a8980", "#a9a8a0"
d = pd.read_csv("manuscript/Supplementary_Data_1_filtration_support_map.csv")
d["chr"] = pd.to_numeric(d.chromosome, errors="coerce"); d = d[d.chr.notna() & d.filtration_support_f_s.notna()].copy()
L = d.groupby("chr").position_GRCh37.max(); off = L.cumsum().shift(fill_value=0); d["x"] = d.position_GRCh37 + d.chr.map(off) + d.chr * 2e7
fig = plt.figure(figsize=(5, 5)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
def box(x, y, w, h, fc, ec):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=2", fc=fc, ec=ec, lw=1.4))
def arrow(x0, y0, x1, y1, c=MUTED):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=16, lw=1.6, color=c))
ax.text(50, 95, "eGFRcrea genetic signal:\nfiltration or creatinine handling?", ha="center", va="center", fontsize=13, fontweight="bold", color=INK)
box(8, 76, 84, 9, "#f3f4f6", GREY)
ax.text(50, 80.5, "2,051 fine-mapped signals  ×  colocalization with eGFRcys", ha="center", va="center", fontsize=9.8, color=INK)
arrow(50, 76, 50, 69.5)
m = fig.add_axes([0.17, 0.43, 0.76, 0.22])
col = {"Filtration-consistent": BLUE, "Intermediate": GREY, "Low support": ORANGE}
for k, c in col.items():
    s = d[d.classification == k]; m.scatter(s.x, s.filtration_support_f_s, s=3, color=c, lw=0)
m.set_xticks([]); m.set_yticks([0, 1]); m.tick_params(labelsize=9)
for sp in ["top", "right"]: m.spines[sp].set_visible(False)
m.set_ylabel("Filtration support", fontsize=9.5)
m.set_title("Filtration-support map (1,681 scored)", fontsize=11, color=INK, pad=4)
arrow(30, 40, 24, 33, BLUE); arrow(70, 40, 76, 33, ORANGE)
box(3, 4, 44, 28, "#e8f1fc", BLUE); box(53, 4, 44, 28, "#fdeee7", ORANGE)
ax.text(25, 27.5, "Interpretation", ha="center", fontsize=12.5, fontweight="bold", color=BLUE)
ax.text(25, 15, "Separates creatinine-\nspecific from GFR-\nassociated loci\n(0.002 vs. 0.392)", ha="center", va="center", fontsize=10.5, color=INK, linespacing=1.3)
ax.text(75, 27.5, "Gene prioritization", ha="center", fontsize=12.5, fontweight="bold", color=ORANGE)
ax.text(75, 15, "No gain over a\nmultiomic scorecard\n(ΔAUPRC −0.0003);\nno dialysis enrichment", ha="center", va="center", fontsize=10.5, color=INK, linespacing=1.3)
fig.savefig("figs/final/Graphical_abstract.tif", dpi=240, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
fig.savefig("figs/final/Graphical_abstract.png", dpi=240, facecolor="white")
