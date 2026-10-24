# -*- coding: utf-8 -*-
"""Fig.1 study flowchart — publication style (English, PDF vector + 600 dpi PNG)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams["font.family"] = "Arial"

OUT = r"D:\GEO-limma,肌少症1\18-Fig1\output"

# palette
C_SP_FC, C_SP_EC = "#FDECEA", "#E64B35"   # sarcopenia line
C_HP_FC, C_HP_EC = "#E8F4F8", "#4DBBD5"   # H. pylori line
C_SH_FC, C_SH_EC = "#F2F2F2", "#7F7F7F"   # shared steps
C_HU_FC, C_HU_EC = "#FFF6DC", "#C9A227"   # hub genes (highlight)
C_VA_FC, C_VA_EC = "#EAF3EA", "#3C8C50"   # validation band
C_TXT = "#1A1A1A"
C_ARR = "#4D4D4D"

fig, ax = plt.subplots(figsize=(7.4, 5.9), dpi=600)
ax.set_xlim(0, 100)
ax.set_ylim(-3.2, 78)
ax.axis("off")


def box(cx, cy, w, h, text, fc, ec, fs=6.8, lw=1.0, bold_first=False, rounding=0.9):
    ax.add_patch(FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle=f"round,pad=0,rounding_size={rounding}",
        linewidth=lw, edgecolor=ec, facecolor=fc, zorder=2))
    lines = text.split("\n")
    if bold_first and len(lines) > 1:
        ax.text(cx, cy + 0.30 * h, lines[0], ha="center", va="center",
                fontsize=fs + 0.6, fontweight="bold", color=C_TXT, zorder=3)
        ax.text(cx, cy - 0.18 * h, "\n".join(lines[1:]), ha="center", va="center",
                fontsize=fs, color=C_TXT, zorder=3, linespacing=1.25)
    else:
        ax.text(cx, cy, text, ha="center", va="center",
                fontsize=fs, color=C_TXT, zorder=3, linespacing=1.3)


def arrow(x1, y1, x2, y2, lw=1.1, color=C_ARR, rad=0.0):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=9,
        linewidth=lw, color=color, zorder=1,
        connectionstyle=f"arc3,rad={rad}", shrinkA=1.5, shrinkB=1.5))


def section(label, y):
    ax.text(1.5, y, label, ha="left", va="center", fontsize=6.2, style="italic",
            color="#8A8A8A", rotation=90)
    ax.plot([4.6, 4.6], [y - 8.5, y + 8.5], color="#BFBFBF", lw=0.9, zorder=1)


# ---------------- I. Data collection ----------------
section("Data collection", 62)
box(30, 68, 34, 9.5,
    "Sarcopenia (SP) line\nGSE1428 (GPL96)\nVastus lateralis, Young 10 vs Old 12",
    C_SP_FC, C_SP_EC, bold_first=True)
box(70, 68, 34, 9.5,
    "H. pylori (HP) line\nGSE27411 (GPL6255)\nGastric mucosa, Control 6 vs HP 6",
    C_HP_FC, C_HP_EC, bold_first=True)

box(30, 56.5, 30, 8,
    "limma differential analysis\n1,153 DEGs\n(605 Up / 548 Down)", C_SP_FC, C_SP_EC)
box(70, 56.5, 30, 8,
    "limma differential analysis\n2,430 DEGs\n(1,257 Up / 1,173 Down)", C_HP_FC, C_HP_EC)

arrow(30, 63.2, 30, 60.6)
arrow(70, 63.2, 70, 60.6)

# ---------------- II. Hub-gene screening ----------------
section("Hub-gene screening", 36)
box(50, 44.5, 36, 8,
    "Common DEGs\nVenn intersection: SP \u2229 HP = 118 genes",
    C_SH_FC, C_SH_EC, bold_first=True)
arrow(30, 52.4, 43, 48.6, rad=-0.12)
arrow(70, 52.4, 57, 48.6, rad=0.12)

box(50, 33, 30, 7,
    "WGCNA (SP & HP lines)\nkey module genes", C_SH_FC, C_SH_EC, bold_first=True)
arrow(50, 40.4, 50, 36.6)

# three ML methods
box(30, 23, 17.5, 6.5, "LASSO\nregression", C_SH_FC, C_SH_EC)
box(50, 23, 17.5, 6.5, "SVM-RFE", C_SH_FC, C_SH_EC)
box(70, 23, 17.5, 6.5, "Random\nforest", C_SH_FC, C_SH_EC)
arrow(41, 29.4, 32, 26.4, rad=-0.1)
arrow(50, 29.4, 50, 26.4)
arrow(59, 29.4, 68, 26.4, rad=0.1)

box(50, 13, 44, 8,
    "Hub genes\nC1QA  ·  NFKBIA  ·  DTL  ·  CTNNBIP1",
    C_HU_FC, C_HU_EC, bold_first=True, lw=1.6)
arrow(30, 19.6, 42, 17.1, rad=-0.12)
arrow(50, 19.6, 50, 17.1)
arrow(70, 19.6, 58, 17.1, rad=0.12)

# ---------------- III. Validation & functional analysis ----------------
section("Validation & functional analysis", 5.0)
vb = [
    (15.5, "External ROC validation\nGSE111016 · GSE5081\nGSE25941 · GPL20301"),
    (38.5, "Functional analysis\nGO & KEGG enrichment\nGSVA (SP & HP lines)"),
    (61.5, "Immune infiltration\nssGSEA + gene\u2013immune-cell\ncorrelation (SP & HP lines)"),
    (84.5, "Multi-omics validation\nscRNA: SKM atlas · GSE249874\nProtein: HPA (IHC)"),
]
for cx, txt in vb:
    box(cx, 1.2, 22.5, 7.4, txt, C_VA_FC, C_VA_EC, fs=6.2)
# fan-out from hub box bottom edge
for tx, cx in zip((36, 44, 56, 64), (15.5, 38.5, 61.5, 84.5)):
    arrow(tx, 9.0, cx, 5.0, rad=(cx - tx) / 55.0)

fig.savefig(OUT + r"\Fig1_flowchart.pdf", bbox_inches="tight")
fig.savefig(OUT + r"\Fig1_flowchart.png", bbox_inches="tight", dpi=600)
print("saved Fig1_flowchart.pdf/.png")
