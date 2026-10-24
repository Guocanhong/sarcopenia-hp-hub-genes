# -*- coding: utf-8 -*-
# 3-full_analysis.py  SKM full atlas: composition old vs young + hub gene localization
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import mannwhitneyu

a = ad.read_h5ad("raw/SKM_human_pp_full.h5ad")
print("loaded:", a.shape)
# var_names -> SYMBOL（去重）
if not a.var_names.equals(pd.Index(a.var["SYMBOL"].astype(str))):
    sym = a.var["SYMBOL"].astype(str).copy()
    sym = sym.where(~sym.duplicated(keep="first"), sym + "_dup" + pd.Series(range(len(sym))).astype(str))
    a.var_names = sym.values
a.var_names_make_unique()

# 重组到功能 compartment（level0 -> compartment）
L0 = a.obs["annotation_level0"].astype(str)
comp_map = {}
for x in ["MF-I","MF-II","MF-Isc(fg)","MF-IIsn(fg)","MF-IIsc(fg)","MF-Isn(fg)","MF-IIsn(fg)","Specialised MF"]:
    comp_map[x] = "Myofiber"
for x in ["Macrophage","Monocyte","Mast","Neutrophil","Eosinophil","T-cell","B-cell","B-plasma","NK-cell","cDC1","cDC2","pDC","Hyb"]:
    comp_map[x] = "Immune"
for x in ["FB","EnFB","PnFB"]: comp_map[x] = "Fibroblast"
for x in ["ArtEC","CapEC","VenEC","LymphEC"]: comp_map[x] = "Endothelium"
for x in ["SMC","Pericyte"]: comp_map[x] = "SMC_Pericyte"
comp_map["MuSC"] = "MuSC"
comp_map["Tenocyte"] = "Tenocyte"; comp_map["Mesothelium"] = "Other"
for x in ["mSchwann","nmSchwann"]: comp_map[x] = "Schwann"
comp_map["Adipocyte"] = "Adipocyte"; comp_map["RBC"] = "Other"
comp = L0.map(lambda x: comp_map.get(x, "Other"))
a.obs["compartment"] = pd.Categorical(comp)

HUB = ["C1QA","NFKBIA","CTNNBIP1","DTL"]
MARK = ["C1QB","CD68","LYZ","PTPRC","MYH7","MYH2"]
gi = {g: a.var_names.get_loc(g) for g in HUB + MARK}
print("gene index ok")

# ---------- 1. 细胞组成 old vs young（per donor） ----------
obs = a.obs[["DonorID","Age_bin","Age_group","compartment","Sex"]].copy()
donor_meta = obs.groupby("DonorID").agg(age=("Age_bin","first"), agegrp=("Age_group","first"), sex=("Sex","first"))
ct = pd.crosstab(obs["DonorID"], obs["compartment"])
frac = ct.div(ct.sum(1), axis=0) * 100
frac = frac.join(donor_meta)
frac.to_csv("output/full_composition_per_donor.csv")
print("\n=== per-donor compartment % ===")
print(frac.round(2).to_string())

res = []
for c in [x for x in frac.columns if x not in ("age","agegrp","sex")]:
    o = frac[frac.age=="old"][c]; y = frac[frac.age=="young"][c]
    u,p = mannwhitneyu(o, y, alternative="two-sided")
    res.append(dict(compartment=c, old_median=o.median(), young_median=y.median(),
                    old_n=len(o), young_n=len(y), p=p))
st = pd.DataFrame(res); st["p_adj"] = __import__("numpy").minimum(st["p"]*len(st), 1.0)
st.to_csv("output/full_composition_stats.csv", index=False)
print("\n=== composition stats (MWU old vs young) ===")
print(st.round(4).sort_values("p").to_string(index=False))

# ---------- 2. hub 基因按 compartment / 主要细胞类型表达 ----------
X = a.X
def col(g): return np.asarray(X[:, gi[g]].todense()).ravel()
G = {g: col(g) for g in HUB + MARK}
gd = pd.DataFrame(G); gd["compartment"] = comp.values
t = gd.groupby("compartment").agg(n=("compartment","size"),
        **{f"{g}_mean": (g,"mean") for g in HUB+MARK},
        **{f"{g}_pct": (g, lambda s,g=g: (s>0).mean()*100) for g in HUB+MARK})
t.to_csv("output/full_hubgene_by_compartment.csv")
print("\n=== hub genes by compartment ===")
print(t[["n"]+[f"{g}_mean" for g in HUB]+[f"{g}_pct" for g in HUB]].round(3).to_string())

# 免疫精细类型 level0（Macrophage/Mono/T/NK/DC...）C1QA
gd["L0"] = L0.values
t2 = gd.groupby("L0").agg(n=("L0","size"), **{f"{g}_mean": (g,"mean") for g in HUB+MARK})
t2.to_csv("output/full_hubgene_by_level0.csv")
print("\n=== C1QA by level0 (immune rows) ===")
imm = ["Macrophage","Monocyte","Mast","Neutrophil","Eosinophil","T-cell","B-cell","B-plasma","NK-cell","cDC1","cDC2","pDC","Hyb"]
print(t2.loc[imm, ["n","C1QA_mean","CD68_mean","LYZ_mean"]].round(3).to_string())

# ---------- 3. 相关性 & 巨噬细胞内 hub 基因 old vs young ----------
from scipy.stats import spearmanr
for m in ["C1QB","CD68","LYZ"]:
    r,p = spearmanr(gd["C1QA"], gd[m]); print(f"Spearman C1QA x {m}: rho={r:.3f}")
sel = (comp.values == "Immune")
mp = gd.loc[sel].copy(); mp["L0"] = L0.values[sel]
mpd = mp[mp.L0 == "Macrophage"].copy()
mpd["donor"] = obs["DonorID"].values[mpd.index.to_numpy()]
mpd["age"] = obs["Age_bin"].values[mpd.index.to_numpy()]
mg = mpd.groupby(["donor","age"])[HUB].mean().reset_index()
mg.to_csv("output/full_mphi_hubgene_per_donor.csv", index=False)
print("\n=== hub genes in Macrophages per donor ===")
print(mg.round(3).to_string(index=False))
for g in HUB:
    o = mg[mg.age=="old"][g]; y = mg[mg.age=="young"][g]
    u,p = mannwhitneyu(o,y,alternative="two-sided")
    print(f"{g}: old {o.median():.3f} vs young {y.median():.3f}, MWU p={p:.3f}")

# ---------- 4. 免疫 compartment 比例 old vs young（核心图数据已在 1 中） ----------
print("\nDONE")
