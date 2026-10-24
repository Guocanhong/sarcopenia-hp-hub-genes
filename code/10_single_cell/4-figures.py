# 4-figures.py  SKM full atlas: publication figures
# Fig_SC1 composition stacked | Fig_SC2 Immune/MuSC box | Fig_SC3 DotPlot
# Fig_SC4 UMAP compartment | FeaturePlot x4 (C1QA/NFKBIA/CTNNBIP1/DTL)
# Memory-safe: CSR chunked scan (nnz 247M, cannot load whole X)
import numpy as np, pandas as pd, h5py, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "axes.linewidth": 0.6, "pdf.fonttype": 42})

H5 = "raw/SKM_human_pp_full.h5ad"
HUB = ["C1QA", "NFKBIA", "CTNNBIP1", "DTL"]
MARK = ["C1QB", "CD68", "MYH7", "MYH2"]
GENES = HUB + MARK

# ---------- metadata (backed, no X) ----------
a = ad.read_h5ad(H5, backed="r")
sym = a.var["SYMBOL"].astype(str).values if "SYMBOL" in a.var else np.array(a.var_names, dtype=object)
gpos = {g: int(np.where(sym == g)[0][0]) for g in GENES}
print("gene positions:", gpos)

obs = a.obs[["DonorID", "Age_bin", "DonorID", "annotation_level0"]].copy()
obs = a.obs[["DonorID", "Age_bin", "annotation_level0"]].copy()
L0 = obs["annotation_level0"].astype(str).values
age = obs["Age_bin"].astype(str).values
don = obs["DonorID"].astype(str).values
umap = np.asarray(a.obsm["X_umap"], dtype=np.float32)
n = umap.shape[0]
print("cells:", n)

# ---------- compartment / display mapping (same as 3-full_analysis.py) ----------
comp_map = {}
for x in ["MF-I","MF-II","MF-Isc(fg)","MF-IIsn(fg)","MF-IIsc(fg)","MF-Isn(fg)","Specialised MF"]:
    comp_map[x] = "Myofiber"
for x in ["Macrophage","Monocyte","Mast","Neutrophil","Eosinophil","T-cell","B-cell","B-plasma","NK-cell","cDC1","cDC2","pDC","Hyb"]:
    comp_map[x] = "Immune"
for x in ["FB","EnFB","PnFB"]: comp_map[x] = "Fibroblast"
for x in ["ArtEC","CapEC","VenEC","LymphEC"]: comp_map[x] = "Endothelium"
for x in ["SMC","Pericyte"]: comp_map[x] = "SMC_Pericyte"
comp_map.update({"MuSC":"MuSC","Tenocyte":"Tenocyte","Mesothelium":"Other",
                 "mSchwann":"Schwann","nmSchwann":"Schwann","Adipocyte":"Adipocyte","RBC":"Other"})
COMP_ORDER = ["Myofiber","Immune","Fibroblast","Endothelium","SMC_Pericyte",
              "MuSC","Schwann","Tenocyte","Adipocyte","Other"]
comp = np.array([comp_map.get(x, "Other") for x in L0])

disp_map = dict(comp_map)
disp_map.update({x: x for x in ["Macrophage","Monocyte","Mast","Neutrophil","Eosinophil",
                                "T-cell","B-cell","B-plasma","NK-cell","cDC1","cDC2","pDC","Hyb"]})
DISP_ORDER = ["Myofiber","MuSC","Fibroblast","Endothelium","SMC/Pericyte" if False else "SMC_Pericyte",
              "Schwann","Tenocyte","Adipocyte",
              "Macrophage","Monocyte","cDC1","cDC2","pDC","Mast","Neutrophil","Eosinophil",
              "T-cell","NK-cell","B-cell","B-plasma"]
disp = np.array([disp_map.get(x, "Other") for x in L0])
extra = sorted(set(np.unique(disp)) - set(DISP_ORDER + ["Other"]))
if extra:
    print("extra display cats (counted, not in DotPlot):", extra)
disp_cat = pd.Categorical(disp, categories=DISP_ORDER + ["Other"] + extra)
dcodes = disp_cat.codes.astype(np.int64)

# ---------- chunked CSR scan ----------
f = h5py.File(H5, "r")
X = f["X"]
data, indices, indptr = X["data"], X["indices"], X["indptr"]
order = np.argsort([gpos[g] for g in GENES])   # sorted-position k -> GENES[order[k]]
t = np.array([gpos[g] for g in GENES])[order]
G = np.zeros((n, len(GENES)), dtype=np.float32)
n_disp = len(disp_cat.categories)
sum_d = np.zeros((n_disp, len(GENES))); pos_d = np.zeros((n_disp, len(GENES))); cnt_d = np.zeros(n_disp)

CH = 15000
for s in range(0, n, CH):
    e = min(s + CH, n)
    p0, p1 = int(indptr[s]), int(indptr[e])
    dseg = data[p0:p1]; iseg = indices[p0:p1]
    pos = np.searchsorted(t, iseg); np.clip(pos, 0, len(t) - 1, out=pos)
    hit = t[pos] == iseg
    cnt_d += np.bincount(dcodes[s:e], minlength=n_disp)
    if hit.any():
        rows_local = np.repeat(np.arange(s, e), np.diff(indptr[s:e + 1]))[hit]
        gc = pos[hit]
        G[rows_local, gc] = dseg[hit]
        w = dseg[hit]
        rr = rows_local - s
        for k in range(len(GENES)):
            m = gc == k
            if m.any():
                sum_d[:, k] += np.bincount(dcodes[s:e][rr[m]], weights=w[m], minlength=n_disp)
                pos_d[:, k] += np.bincount(dcodes[s:e][rr[m]], weights=(w[m] > 0).astype(float), minlength=n_disp)
    if (s // CH) % 3 == 0:
        print(f"  scan {e}/{n}", flush=True)
f.close()
inv = np.argsort(order)   # 把列从 sorted-position 顺序重排回 GENES 顺序
G = np.ascontiguousarray(G[:, inv])
sum_d = sum_d[:, inv]; pos_d = pos_d[:, inv]

mean_d = np.divide(sum_d, cnt_d[:, None], out=np.zeros_like(sum_d), where=cnt_d[:, None] > 0)
pct_d = np.divide(pos_d, cnt_d[:, None], out=np.zeros_like(pos_d), where=cnt_d[:, None] > 0) * 100
np.save("output/full_expr_per_cell.npy", G)

# ---------- Fig_SC1 composition stacked ----------
donors = obs.groupby("DonorID").agg(age=("Age_bin", "first")).reset_index()
donors["age"] = donors["age"].astype(str)
old_ids = sorted(donors.loc[donors.age == "old", "DonorID"])
young_ids = sorted(donors.loc[donors.age == "young", "DonorID"])
order = old_ids + young_ids
ct = pd.crosstab(don, comp)
frac = ct.div(ct.sum(1), axis=0).loc[order] * 100
COMP_COLORS = {"Myofiber":"#B0B7C0","Immune":"#D64545","Fibroblast":"#7FA8D9","Endothelium":"#F2A541",
               "SMC_Pericyte":"#8E6C9E","MuSC":"#2E86AB","Schwann":"#C9B458","Tenocyte":"#67B99A",
               "Adipocyte":"#E8C4B8","Other":"#DDDDDD"}
fig, ax = plt.subplots(figsize=(7.0, 3.4))
x = np.arange(len(order)); bottom = np.zeros(len(order))
for c in COMP_ORDER:
    v = frac[c].values if c in frac else np.zeros(len(order))
    ax.bar(x, v, bottom=bottom, color=COMP_COLORS[c], width=0.75, label=c, edgecolor="white", linewidth=0.2)
    bottom += v
xt = list(old_ids) + list(young_ids)
ax.set_xticks(x); ax.set_xticklabels(xt, rotation=90, fontsize=6.5)
for i, ag in enumerate(["old"] * len(old_ids) + ["young"] * len(young_ids)):
    ax.text(i, 103, ag[0].upper(), ha="center", fontsize=8, fontweight="bold",
            color="#D64545" if ag == "old" else "#2E86AB")
ax.set_ylabel("Cell/nuclei fraction (%)"); ax.set_ylim(0, 108)
ax.legend(ncol=2, fontsize=6.5, frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5))
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("output/Fig_SC1_Composition_stacked.pdf"); fig.savefig("output/Fig_SC1_Composition_stacked.png", dpi=300); plt.close(fig)
print("Fig_SC1 done")

# ---------- Fig_SC2 Immune & MuSC per donor ----------
stats = pd.read_csv("output/full_composition_stats.csv").set_index("compartment")
fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.6))
for ax, c in zip(axes, ["Immune", "MuSC"]):
    ov = frac.loc[old_ids, c].values; yg = frac.loc[young_ids, c].values
    bp = ax.boxplot([yg, ov], widths=0.5, patch_artist=True, medianprops=dict(color="black", linewidth=1))
    for patch, col in zip(bp["boxes"], ["#2E86AB", "#D64545"]):
        patch.set_facecolor(col); patch.set_alpha(0.35); patch.set_edgecolor("black"); patch.set_linewidth(0.6)
    rng = np.random.default_rng(1)
    ax.scatter(rng.normal(1, 0.05, len(yg)), yg, s=14, color="#2E86AB", edgecolor="black", linewidth=0.3, zorder=3)
    ax.scatter(rng.normal(2, 0.05, len(ov)), ov, s=14, color="#D64545", edgecolor="black", linewidth=0.3, zorder=3)
    p = stats.loc[c, "p"]
    ax.set_title(f"{c}\nMWU p = {p:.2f}" + (" (ns)" if p > 0.05 else ""), fontsize=8)
    ax.set_xticks([1, 2]); ax.set_xticklabels(["Young", "Old"], fontsize=8)
    ax.set_ylabel(f"% of all cells")
    ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("output/Fig_SC2_Immune_MuSC_boxplot.pdf"); fig.savefig("output/Fig_SC2_Immune_MuSC_boxplot.png", dpi=300); plt.close(fig)
print("Fig_SC2 done")

# ---------- Fig_SC3 DotPlot ----------
cats = [c for c in DISP_ORDER if cnt_d[disp_cat.categories.get_loc(c)] > 0]
ci = [disp_cat.categories.get_loc(c) for c in cats]
M = mean_d[ci, :len(HUB)]; P = pct_d[ci, :len(HUB)]
fig, ax = plt.subplots(figsize=(7.0, 2.6))
for k in range(len(HUB)):
    scaled = M[:, k] / max(M[:, k].max(), 1e-9)
    sc = ax.scatter(np.arange(len(cats)), np.full(len(cats), k, dtype=float),
                    s=P[:, k] * 1.9, c=scaled, cmap="YlOrRd",
                    vmin=0, vmax=1, edgecolor="#555555", linewidth=0.4)
ax.set_xticks(np.arange(len(cats))); ax.set_xticklabels(cats, rotation=90, fontsize=7)
ax.set_yticks(np.arange(len(HUB))); ax.set_yticklabels(HUB, fontsize=8)
ax.set_xlim(-0.7, len(cats) - 0.3); ax.invert_yaxis()
cb = fig.colorbar(sc, ax=ax, shrink=0.85, pad=0.01)
cb.set_label("Mean expression\n(scaled per gene)", fontsize=7); cb.ax.tick_params(labelsize=6)
for sz in [25, 50, 75, 100]:
    ax.scatter([], [], s=sz * 1.9, facecolor="white", edgecolor="#555555", label=f"{sz}%")
ax.legend(scatterpoints=1, frameon=False, fontsize=6.5, loc="center left",
          bbox_to_anchor=(1.12, 0.5), title="% expressing", title_fontsize=7)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("output/Fig_SC3_DotPlot_hubgenes.pdf"); fig.savefig("output/Fig_SC3_DotPlot_hubgenes.png", dpi=300); plt.close(fig)
print("Fig_SC3 done")

# ---------- Fig_SC4 UMAP compartment ----------
cmask = {c: COMP_COLORS[c] for c in COMP_ORDER}
fig, ax = plt.subplots(figsize=(4.6, 3.6))
for c in COMP_ORDER:
    m = comp == c
    ax.scatter(umap[m, 0], umap[m, 1], s=0.25, color=cmask[c], label=f"{c} (n={m.sum():,})",
               rasterized=True, linewidths=0)
ax.set_xticks([]); ax.set_yticks([]); ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
ax.legend(markerscale=12, fontsize=6.5, frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5))
fig.tight_layout(); fig.savefig("output/Fig_SC4_UMAP_compartment.pdf"); fig.savefig("output/Fig_SC4_UMAP_compartment.png", dpi=300); plt.close(fig)
print("Fig_SC4 done")

# ---------- FeaturePlots ----------
for k, g in enumerate(HUB):
    v = G[:, GENES.index(g)]
    pos = v > 0
    fig, ax = plt.subplots(figsize=(3.4, 3.0))
    ax.scatter(umap[~pos, 0], umap[~pos, 1], s=0.25, color="#E8E8E8", rasterized=True, linewidths=0)
    if pos.any():
        sc = ax.scatter(umap[pos, 0], umap[pos, 1], s=0.6, c=v[pos], cmap="magma",
                        rasterized=True, linewidths=0)
        cb = fig.colorbar(sc, ax=ax, shrink=0.8, pad=0.02)
        cb.set_label("Expression", fontsize=7); cb.ax.tick_params(labelsize=6)
    ax.set_title(f"{g}  (pos {pos.mean()*100:.1f}%, max {v.max():.1f})", fontsize=8)
    ax.set_xticks([]); ax.set_yticks([]); ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    fig.tight_layout(); fig.savefig(f"output/FeaturePlot_{g}.pdf"); fig.savefig(f"output/FeaturePlot_{g}.png", dpi=300); plt.close(fig)
    print(f"FeaturePlot {g} done")

# CSV for records
rec = pd.DataFrame(mean_d[ci, :], index=cats, columns=GENES)
rec_pct = pd.DataFrame(pct_d[ci, :], index=cats, columns=GENES)
rec.to_csv("output/fig_dotplot_mean.csv"); rec_pct.to_csv("output/fig_dotplot_pct.csv")
print("ALL DONE")
