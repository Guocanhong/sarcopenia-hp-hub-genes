# -*- coding: utf-8 -*-
"""Fig.1 flowchart - reference style (blue-outline white rounded boxes, top-down
spine with side columns, yellow ellipse for hub genes, bottom fan-out).
Ref: Apoptosis (2025) 30:1589 - GJB2 ovarian cancer paper Fig.1.
v5: two-disease flow made explicit - each line runs its own candidate pool and
its own ML screening (LASSO / SVM-RFE / RF), then both lanes converge on the
hub-gene ellipse (intersection of the two lines' ML candidates)."""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse, FancyArrowPatch

BORDER = '#2E75B6'
TEXT   = '#17375E'
ARROW  = '#2E75B6'
HUB_FC = '#FFE699'
HUB_EC = '#C55A11'
HUB_TX = '#C00000'
UP_TX  = '#C00000'
DN_TX  = '#2E75B6'

fig, ax = plt.subplots(figsize=(7.4, 7.5))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis('off')

def box(cx, cy, w, h, text, fs=6.3, bold_first=False, fc='white', ec=BORDER):
    ax.add_patch(FancyBboxPatch((cx-w/2, cy-h/2), w, h,
        boxstyle='round,pad=0.25,rounding_size=0.9',
        fc=fc, ec=ec, lw=1.3, zorder=2))
    lines = text.split('\n')
    if bold_first:
        ax.text(cx, cy + h*0.28, lines[0], ha='center', va='center',
                fontsize=fs+0.4, fontweight='bold', color=TEXT, zorder=3)
        rest = '\n'.join(lines[1:])
        if rest:
            ax.text(cx, cy - h*0.16, rest, ha='center', va='center',
                    fontsize=fs, color=TEXT, zorder=3, linespacing=1.35)
    else:
        ax.text(cx, cy, text, ha='center', va='center', fontsize=fs,
                color=TEXT, zorder=3, linespacing=1.35)

def ellipse(cx, cy, w, h, l1, l2, fs1=6.6):
    ax.add_patch(Ellipse((cx, cy), w, h, fc=HUB_FC, ec=HUB_EC, lw=1.6, zorder=2))
    ax.text(cx, cy+h*0.20, l1, ha='center', va='center', fontsize=fs1,
            fontweight='bold', color=TEXT, zorder=3)
    ax.text(cx, cy-h*0.20, l2, ha='center', va='center', fontsize=7.0,
            fontstyle='italic', fontweight='bold', color=HUB_TX, zorder=3)

def arrow(x1, y1, x2, y2, rad=0.0, lw=1.2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>',
        mutation_scale=10, lw=lw, color=ARROW, zorder=1,
        connectionstyle=f'arc3,rad={rad}', shrinkA=1.5, shrinkB=1.5))

def dir_icons(cx, cy, up, down, fs=5.6):
    ax.text(cx-1.7, cy, up, ha='right', va='center', fontsize=fs,
            color=UP_TX, fontweight='bold', zorder=3)
    ax.text(cx+1.7, cy, down, ha='left', va='center', fontsize=fs,
            color=DN_TX, fontweight='bold', zorder=3)

# ---------------- Row 1-2: datasets -> limma ----------------
box(28, 94, 27, 8.2,
    'GSE1428\nSkeletal muscle (vastus lateralis)\nYoung n=10 / Old n=12', bold_first=True)
box(72, 94, 27, 8.2,
    'GSE27411\nGastric mucosa (atrophic gastritis)\nControl n=6 / H. pylori n=6', bold_first=True)
box(28, 83.5, 27, 8.0,
    'limma differential analysis\nDEGs n = 1,153 (SP line)\nnominal P < 0.05')
box(72, 83.5, 27, 8.0,
    'limma differential analysis\nDEGs n = 2,430 (HP line)\nnominal P < 0.05')
arrow(28, 89.4, 28, 88.0, lw=1.3)
dir_icons(28, 88.7, '\u25b2 605', '\u25bc 548')
arrow(72, 89.4, 72, 88.0, lw=1.3)
dir_icons(72, 88.7, '\u25b2 1,257', '\u25bc 1,173')

# ---------------- Row 3: common DEGs (Venn) + GO/KEGG directly below ----------------
box(50, 73.5, 30, 6.4, 'Common DEGs (Venn intersection)\nn = 118')
box(50, 65.8, 24, 5.4, 'GO & KEGG enrichment\n(118 common DEGs)', fs=5.7)
arrow(50, 70.3, 50, 68.7)
arrow(33, 79.1, 43, 76.2, rad=-0.06)
arrow(67, 79.1, 57, 76.2, rad=0.06)
dir_icons(30.3, 77.7, '\u25b2', '\u25bc', fs=5.4)
dir_icons(69.7, 77.7, '\u25b2', '\u25bc', fs=5.4)

# ---------------- GSVA beside each dataset (split per line) ----------------
box(10.5, 74.5, 17, 7.0, 'GSVA (SP line)\nMSigDB C2\nGSE1428', fs=5.3)
box(89.5, 74.5, 17, 7.0, 'GSVA (HP line)\nMSigDB C2\nGSE27411', fs=5.3)
arrow(14.5, 91.5, 10.8, 78.3, rad=0.06)   # GSE1428 -> GSVA SP
arrow(85.5, 91.5, 89.2, 78.3, rad=-0.06)  # GSE27411 -> GSVA HP

# ---------------- Row 4: WGCNA per dataset (outer columns) ----------------
box(18, 64.5, 27, 7.0,
    'WGCNA (GSE1428)\nSP-line key module genes', bold_first=True)
box(82, 64.5, 27, 7.0,
    'WGCNA (GSE27411)\nHP-line key module genes', bold_first=True)
arrow(25, 79.4, 20, 68.6, rad=0.03)
arrow(75, 79.4, 80, 68.6, rad=-0.03)

# ---------------- Row 5: per-line candidate pools ----------------
box(30, 54, 28, 6.0,
    'Candidate pool (SP line)\ncommon DEGs \u2229 WGCNA modules', bold_first=True, fs=5.9)
box(70, 54, 28, 6.0,
    'Candidate pool (HP line)\ncommon DEGs \u2229 WGCNA modules', bold_first=True, fs=5.9)
arrow(18, 60.0, 26, 57.3, rad=0.05)      # WGCNA SP -> SP pool
arrow(82, 60.0, 74, 57.3, rad=-0.05)     # WGCNA HP -> HP pool
arrow(37, 70.3, 30, 57.3, rad=0.05)      # Venn -> SP pool
arrow(63, 70.3, 70, 57.3, rad=-0.05)     # Venn -> HP pool

# ---------------- Row 6: per-line ML screening ----------------
box(30, 44, 28, 8.0,
    'ML feature selection (SP line)\n\u2460 LASSO  \u2461 SVM-RFE\n\u2462 Random forest',
    bold_first=True, fs=5.9)
box(70, 44, 28, 8.0,
    'ML feature selection (HP line)\n\u2460 LASSO  \u2461 SVM-RFE\n\u2462 Random forest',
    bold_first=True, fs=5.9)
arrow(30, 50.5, 30, 48.5)
arrow(70, 50.5, 70, 48.5)

# ---------------- hub ellipse (convergence of both lines) ----------------
ellipse(50, 32.5, 46, 10.0,
        'Hub genes (intersection of two lines)',
        'C1QA \u00b7 NFKBIA \u00b7 DTL \u00b7 CTNNBIP1', fs1=6.2)
arrow(28, 39.6, 36, 36.8, rad=-0.08)
arrow(72, 39.6, 64, 36.8, rad=0.08)

# ---------------- side columns (flank the hub ellipse, split per line) ----------------
box(11.5, 35.5, 20, 6.6,
    'ROC analysis (SP line)\nGSE1428 \u00b7 GSE111016\nGSE25941 \u00b7 GPL20301', fs=5.8)
box(11.5, 28, 20, 6.0,
    'ROC analysis (HP line)\nGSE27411 \u00b7 GSE5081', fs=5.8)
box(11.5, 21, 20, 5.6, 'Nomogram analysis')
box(88.5, 35.5, 20, 7.0,
    'Immune infiltration\nssGSEA (Bindea 24 signatures)\nGSE1428 \u00b7 GSE27411', fs=5.4)
arrow(21.7, 35.5, 28.3, 34.0, rad=0.06)
arrow(21.7, 28.0, 28.6, 30.8, rad=0.02)
arrow(21.7, 21.0, 30.8, 28.0, rad=-0.10)
arrow(78.3, 35.5, 71.7, 34.0, rad=-0.06)

# ---------------- bottom fan-out (multi-omics, readouts merged) ----------------
box(20, 10.5, 26, 9.5,
    'scRNA-seq (muscle)\nSKM ageing cell atlas\nOld 9 / Young 8 \u00b7 183,161 cells\n\u2192 Cell-type localization',
    fs=5.8)
box(50, 10.5, 26, 9.5,
    'scRNA-seq (stomach)\nGSE249874 \u00b7 Hp\u2212 3 / Hp+ 3 \u00b7 26,250 cells\n\u2192 Microenvironment\n   & gene localization',
    fs=5.8)
box(80, 10.5, 26, 9.5,
    'Protein level\nHuman Protein Atlas (IHC)\nmuscle / stomach\n\u2192 Protein corroboration',
    fs=5.8)
arrow(38, 27.3, 24, 15.9, rad=-0.18)
arrow(50, 27.3, 50, 15.7)
arrow(62, 27.3, 76, 15.9, rad=0.18)

out = 'output/Fig1_flowchart_v2'
fig.savefig(out + '.pdf', bbox_inches='tight')
fig.savefig(out + '.png', dpi=600, bbox_inches='tight')
print('saved', out)
