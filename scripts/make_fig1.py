"""Figure 1: study design."""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5})
INK, MUTED = "#1f2328", "#57606a"
BLUE, BLUE_BG = "#2a78d6", "#e8f1fc"
ORANGE, ORANGE_BG = "#eb6834", "#fdeee7"
GRAY_BG, GRAY_LN = "#f3f4f6", "#9aa1a9"

fig = plt.figure(figsize=(7.2, 5.0))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

def box(x, y, w, h, title, body, fc, ec, tcol=INK, fs_body=6.8, fs_title=8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2",
                                fc=fc, ec=ec, lw=1.0))
    ax.text(x + 1.6, y + h - 2.2, title, ha="left", va="top", fontsize=fs_title, fontweight="bold", color=tcol)
    ax.text(x + 1.6, y + h - 6.8, body, ha="left", va="top", fontsize=fs_body, color=INK, linespacing=1.35)

def arrow(x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=10,
                                 lw=1.1, color=MUTED, shrinkA=0, shrinkB=0))

# ---- inputs ----
ax.text(2, 98.5, "Inputs (public summary data)", fontsize=8.5, fontweight="bold", color=MUTED, va="top")
iy, ih = 75, 19
box(2, iy, 30.5, ih, "Kidney Disease Genetic Scorecard",
    "Liu et al., Science 2025\n742 merged loci · 2,051 fine-mapped\nsignals (Table S8) · 8,209 genes\nScore_A, variant–gene links", GRAY_BG, GRAY_LN, fs_title=7.2)
box(34.75, iy, 30.5, ih, "Kidney-function GWAS",
    "eGFRcrea · eGFRcys · BUN (European)\nSummary statistics, hg19\nLD: 1000G EUR (primary);\nUK Biobank (sensitivity)", GRAY_BG, GRAY_LN)
box(67.5, iy, 30.5, ih, "External resources",
    "Benchmarks: ClinGen + PanelApp\n(primary); Groopman 2019; MGI\nConsistency check: Titan loci\nOutcomes: FinnGen R13", GRAY_BG, GRAY_LN)

# ---- aims ----
ay, ah, aw = 21, 47, 30.5
xs = [2, 34.75, 67.5]
box(xs[0], ay, aw, ah, "Aim 1  Filtration support, f_s",
    "Per signal: coloc-SuSiE,\neGFRcrea vs eGFRcys\n\nMap credible set → S8 signal\n(overlap ≥50% or lead r² ≥0.8;\nfallback flagged)\n\nf_s = max PP.H4 if effects\nconcordant; 0 if opposite\n\nExternal consistency check:\nTitan creatinine-specific vs\nGFR-associated loci",
    BLUE_BG, BLUE, tcol=BLUE)
box(xs[1], ay, aw, ah, "Aim 2  Incremental value",
    "Gene score S_g = link-weighted\nmean f_s over linked signals\n\nScore_C = Z(Score_A) + λ·Z(S_g)\nλ tuned on training chromosomes\n(leave-one-chromosome-out CV)\nModels A–D\n\nOutcome: AUPRC vs benchmark\n(154 positive genes)\nMinimum meaningful gain\nΔAUPRC = 0.03",
    BLUE_BG, BLUE, tcol=BLUE)
box(xs[2], ay, aw, ah, "Aim 3  Clinical endpoints",
    "Signals: f_s ≥ 0.8 vs f_s < 0.2\n1-Mb pruning by |z|\n\nFinnGen R13 lead-variant\nassociation (risk-aligned)\nPrimary: dialysis\nSecondary: CKD, AKI,\ntransplant (BH-adjusted)\nNegative control endpoint\n\nNull: propensity-quintile\npermutation (10,000)",
    ORANGE_BG, ORANGE, tcol=ORANGE)

# arrows inputs -> aims, aim -> aim
for x in xs:
    arrow(x + aw / 2, iy, x + aw / 2, ay + ah + 0.3)
arrow(xs[0] + aw, ay + ah / 2, xs[1] - 0.3, ay + ah / 2)
arrow(xs[1] + aw, ay + ah / 2, xs[2] - 0.3, ay + ah / 2)
ax.text((xs[0] + aw + xs[1]) / 2, ay + ah / 2 + 1.5, "f_s", ha="center", fontsize=6.5, color=MUTED)
ax.text((xs[1] + aw + xs[2]) / 2, ay + ah / 2 + 1.5, "f_s", ha="center", fontsize=6.5, color=MUTED)

# ---- outputs ----
oy, oh = 2, 12
outs = ["1,681 / 2,051 signals classifiable\n287 with f_s ≥ 0.8",
        "Score_C vs Score_A\nprimary, sensitivity and\nalternative-benchmark analyses",
        "Enrichment of FinnGen kidney\nendpoint associations"]
for x, t in zip(xs, outs):
    arrow(x + aw / 2, ay, x + aw / 2, oy + oh + 0.3)
    ax.add_patch(FancyBboxPatch((x, oy), aw, oh, boxstyle="round,pad=0,rounding_size=1.2", fc="white", ec=GRAY_LN, lw=0.8, ls="--"))
    ax.text(x + aw / 2, oy + oh / 2, t, ha="center", va="center", fontsize=6.8, color=INK, linespacing=1.35)

for ext in ["png", "tif"]:
    fig.savefig(f"figs/Fig1_study_design.{ext}", dpi=300, facecolor="white",
                **({"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tif" else {}))
