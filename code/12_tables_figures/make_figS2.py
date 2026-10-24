# -*- coding: utf-8 -*-
"""Supplementary Figure S2: GSE111016 per-gene ROC 2x2 (A-D, cohort header)"""
import os
from PIL import Image, ImageDraw, ImageFont

P = r"D:/GEO-limma,肌少症1/26-Figures/panels"
OUT = r"D:/GEO-limma,肌少症1/27-GSE60662/output"
cell, padg, hdr = 900, 24, 100
names = ["ext1_GSE111016_C1QA", "ext1_GSE111016_NFKBIA", "ext1_GSE111016_CTNNBIP1", "ext1_GSE111016_DTL_ND"]
LABELS = ["A", "B", "C", "D"]

W = 2 * cell + 3 * padg
H = hdr + 2 * cell + 3 * padg
canvas = Image.new("RGB", (W, H), "white")
d = ImageDraw.Draw(canvas)
fh = ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", 54)
fl = ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", 56)
d.text((padg, padg), "GSE111016 (SP vs Control, RNA-seq)", fill="black", font=fh)

y = padg + hdr
for r in range(2):
    x = padg
    for c in range(2):
        im = Image.open(os.path.join(P, names[r * 2 + c] + ".png")).convert("RGB").resize((cell, cell))
        dl = ImageDraw.Draw(im)
        letter = LABELS[r * 2 + c]
        bb = dl.textbbox((0, 0), letter, font=fl)
        lw, lh = bb[2] - bb[0], bb[3] - bb[1]
        dl.rectangle([14, 14, 14 + lw + 20, 14 + lh + 16], fill="white", outline="black", width=3)
        dl.text((24, 20), letter, fill="black", font=fl)
        canvas.paste(im, (x, y))
        x += cell + padg
    y += cell + padg

out_png = os.path.join(OUT, "Supplementary_Figure_S2_GSE111016.png")
canvas.save(out_png, dpi=(300, 300))
print("[OK]", out_png, canvas.size)
