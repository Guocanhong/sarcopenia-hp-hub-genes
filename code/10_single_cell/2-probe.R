# 2-probe.R  GSE249874 胃炎层 Hp+ vs Hp- 试探（Seurat + harmony）
# 输入: filtered/ (10x: matrix.mtx.gz, barcodes.tsv, features.tsv) + sample_map.csv
# 输出: output/probe_*.csv
suppressMessages({library(Seurat); library(harmony); library(Matrix); library(dplyr); library(data.table)})

# filtered/ 目录需预先备好: matrix.mtx.gz / barcodes.tsv.gz / features.tsv.gz（真实 gzip 内容）

m <- Read10X("filtered")
obj <- CreateSeuratObject(counts = m, min.cells = 3)
smap <- fread("sample_map.csv")
bc <- colnames(obj)
sidx <- match(sub("^.*-", "", bc), sub("^.*-", "", smap$barcode))
obj$sample <- smap$sample[sidx]
obj$group <- smap$group[sidx]
print(table(obj$sample, obj$group))

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
print(dim(obj))

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- RunHarmony(obj, group.by.vars = "sample", verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
cat("clusters:", nlevels(Idents(obj)), "\n")

# 注释
DefaultAssay(obj) <- "RNA"
markers <- c("EPCAM","KRT8","KRT18","KRT19","MUC5AC","TFF1",          # 上皮
             "PTPRC","CD68","C1QA","C1QB","LYZ","FCGR3A","S100A8",    # 髓系
             "CD3D","CD3E","CD8A","TRAC",                             # T
             "MS4A1","CD79A","MZB1","JCHAIN",                         # B/浆
             "KIT","TPSAB1",                                          # 肥大细胞
             "COL1A1","DCN","ACTA2","RGS5",                           # 基质
             "PECAM1","VWF")                                          # 内皮
av <- NULL  # 注释用各 cluster 的 marker 平均表达（下方直接计算）
mexp <- GetAssayData(obj, layer = "data")
gpos <- function(g) which(rownames(obj) == g)
ann_mat <- sapply(markers, function(g) { i <- gpos(g); if (length(i)) as.numeric(mexp[i,]) else rep(0, ncol(obj)) })
cf <- data.frame(cluster = Idents(obj), ann_mat, check.names = FALSE)
prof <- cf %>% group_by(cluster) %>% summarise(across(everything(), mean), .groups = "drop")
fwrite(prof, "output/probe_cluster_markers.csv")
print(as.data.frame(prof[, c("cluster","EPCAM","CD68","C1QA","LYZ","CD3D","MS4A1","MZB1","KIT","COL1A1","PECAM1")]))

# 自动注释规则（按平均表达最高者定，人工核对 CSV 后可改）
lbl <- setNames(rep("Unknown", nlevels(Idents(obj))), levels(Idents(obj)))
is_max <- function(cl, g) prof[[g]][prof$cluster == cl] > 0.5
for (cl in levels(Idents(obj))) {
  r <- prof[prof$cluster == cl, ]
  if (r$EPCAM > 1 || r$KRT18 > 1) lbl[cl] <- "Epithelial"
  else if (r$CD68 > 0.6 || r$C1QA > 0.5) lbl[cl] <- "Macrophage"
  else if (r$S100A8 > 1) lbl[cl] <- "Monocyte/Neutrophil"
  else if (r$CD3D > 0.5) lbl[cl] <- "T-cell"
  else if (r$MZB1 > 0.5 || r$JCHAIN > 0.5) lbl[cl] <- "Plasma"
  else if (r$MS4A1 > 0.3) lbl[cl] <- "B-cell"
  else if (r$KIT > 0.3) lbl[cl] <- "Mast"
  else if (r$PECAM1 > 0.5) lbl[cl] <- "Endothelial"
  else if (r$COL1A1 > 0.5 || r$ACTA2 > 0.5) lbl[cl] <- "Fibroblast/SMC"
  else lbl[cl] <- "Other"
}
print(lbl)
obj$celltype <- unname(lbl[as.character(Idents(obj))])
saveRDS(obj, "probe_seurat.rds")

# ===== 核心检验 =====
ct_tab <- table(obj$sample, obj$celltype)
ct_frac <- prop.table(ct_tab, margin = 1) * 100
fwrite(as.data.frame(ct_frac), "output/probe_celltype_frac_by_sample.csv")
cat("\n=== cell type % by sample ===\n"); print(round(ct_frac, 1))

# C1QA 按细胞类型 mean / pct
c1qa <- mexp[gpos("C1QA"), ]
st <- data.frame(celltype = obj$celltype, group = obj$group,
                 C1QA = as.numeric(c1qa)) %>%
  group_by(celltype, group) %>%
  summarise(mean = mean(C1QA), pct = mean(C1QA > 0) * 100, n = n(), .groups = "drop")
st$celltype <- factor(st$celltype)
fwrite(st, "output/probe_C1QA_by_celltype.csv")
cat("\n=== C1QA by celltype x group ===\n"); print(as.data.frame(st %>% arrange(celltype, group)), digits = 3)

# 巨噬细胞比例 Hp+ vs Hp-（per donor MWU）
mp <- data.frame(sample = obj$sample, group = obj$group, is_macro = obj$celltype == "Macrophage") %>%
  group_by(sample, group) %>% summarise(frac = mean(is_macro) * 100, .groups = "drop")
print(mp)
if (length(unique(mp$group)) == 2) {
  wt <- wilcox.test(frac ~ group, data = mp)
  cat("MWU macrophage fraction Hp+ vs Hp-: p =", format.pval(wt$p.value), "\n")
}
fwrite(mp, "output/probe_macro_frac_per_sample.csv")
cat("\nPROBE DONE\n")
