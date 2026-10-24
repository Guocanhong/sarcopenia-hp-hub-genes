# 2-montage.py  HPA IHC 三基因×三组织拼图（Fig12）
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "pdf.fonttype": 42})

GENES = [
    ("C1QA", {"muscle": ("2350_7234_B_1_6.jpg", "HPA002350", "Myocytes: Not detected"),
              "stomach": ("2350_7234_A_1_1.jpg", "HPA002350", "Glandular cells: Not detected"),
              "spleen": ("2350_7234_B_7_4.jpg", "HPA002350", "Red pulp: Low (immune-rich)")}),
    ("NFKBIA", {"muscle": ("3815_10280_B_1_6.jpg", "CAB003815", "Myocytes: Medium"),
                "stomach": ("3815_10280_A_1_1.jpg", "CAB003815", "Glandular cells: Medium"),
                "spleen": ("3815_10280_B_7_4.jpg", "CAB003815", "Red/white pulp: Low")}),
    ("CTNNBIP1", {"muscle": ("67089_166489_B_1_6.jpg", "HPA067089", "Myocytes: Not detected"),
                  "stomach": ("67089_166489_A_1_1.jpg", "HPA067089", "Glandular cells: Not detected"),
                  "spleen": ("67089_166489_B_7_4.jpg", "HPA067089", "Red/white pulp: Not detected")}),
    ("DTL", {"muscle": ("28016_59362_B_1_6.jpg", "HPA028016", "Myocytes: Not detected"),
             "stomach": ("28016_59362_A_1_1.jpg", "HPA028016", "Glandular cells: Not detected"),
             "spleen": ("28016_59362_B_7_4.jpg", "HPA028016", "Red/white pulp: Not detected")}),
]
TISSUES = [("muscle", "Skeletal muscle"), ("stomach", "Stomach"), ("spleen", "Spleen")]
PANEL = ["A", "B", "C", "D"]

fig, axes = plt.subplots(4, 3, figsize=(6.0, 8.4))
for j, (tk, tl) in enumerate(TISSUES):
    axes[0, j].set_title(tl, fontsize=9, fontweight="bold", pad=6)
for i, (g, cells) in enumerate(GENES):
    axes[i, 0].annotate(PANEL[i], xy=(-0.28, 0.5), xycoords="axes fraction",
                        fontsize=12, fontweight="bold", ha="center", va="center")
    for j, (tk, tl) in enumerate(TISSUES):
        ax = axes[i, j]
        fn, ab, note = cells[tk]
        img = mpimg.imread(f"raw/{fn}")
        ax.imshow(img)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_linewidth(0.5); s.set_color("#888888")
        ax.set_xlabel(f"{note}\n({ab})", fontsize=6.2, labelpad=2)
        if j == 0:
            ax.annotate(g, xy=(-0.14, 0.5), xycoords="axes fraction", rotation=90,
                        fontsize=9, fontweight="bold", ha="center", va="center")
fig.suptitle("Protein expression of hub genes in human tissues (Human Protein Atlas)",
             fontsize=9, y=0.985)
fig.text(0.5, 0.005,
         "Immunohistochemistry images from the Human Protein Atlas (proteinatlas.org), Uhl\u00e9n et al., Science 2015 / Thul et al., Science 2017.",
         ha="center", fontsize=6, style="italic")
fig.tight_layout(rect=(0.02, 0.02, 1, 0.965))
fig.savefig("output/Fig12_HPA_IHC.pdf")
fig.savefig("output/Fig12_HPA_IHC.png", dpi=300)
print("DONE")
