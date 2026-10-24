# 1-probe.R  GSE167186 试探：YouHM5 (GSM5098741, 20y) vs OldHM23 (GSM5098753, 90y)
# 数据源：per-GSM csv.gz（2021-02-21 同批上传，全系列一致版本；h5 为 2020/2023 混传不可用）
# 注意：不从脚本内 setwd 中文路径（Rscript 编码问题），由启动时的 CWD 决定工作目录
suppressMessages({
  library(Seurat); library(harmony); library(Matrix); library(data.table)
  library(dplyr); library(ggplot2); library(patchwork)
})

# ---- csv.gz -> 稀疏矩阵（分块读取，控内存；行=基因，列=细胞）----
read_counts_csv <- function(path, chunk_rows = 2000) {
  con <- gzfile(path, "rt"); head.line <- readLines(con, n = 1L); close(con)
  nbarcode <- length(strsplit(head.line, ",", fixed = TRUE)[[1]]) - 1L
  i.list <- list(); j.list <- list(); x.list <- list(); genes.all <- character(0)
  gi <- 0L; sk <- 1L; width <- NA_integer_
  repeat {
    blk <- tryCatch(fread(path, skip = sk, nrows = chunk_rows, header = FALSE,
                          sep = ",", showProgress = FALSE), error = function(e) NULL)
    if (is.null(blk) || nrow(blk) == 0) break
    if (is.na(width)) width <- ncol(blk) - 1L
    if (ncol(blk) - 1L != width) stop("column width mismatch at row ", gi + 1)
    genes.all <- c(genes.all, as.character(blk[[1]]))
    m <- as.matrix(as.data.frame(lapply(blk[, -1], as.numeric)))
    nz <- which(m != 0, arr.ind = TRUE)
    if (nrow(nz)) {
      i.list[[length(i.list) + 1L]] <- nz[, 1] + gi
      j.list[[length(j.list) + 1L]] <- nz[, 2]
      x.list[[length(x.list) + 1L]] <- m[nz]
    }
    gi <- gi + nrow(blk); sk <- sk + chunk_rows
  }
  M <- sparseMatrix(i = unlist(i.list), j = unlist(j.list), x = unlist(x.list),
                    dims = c(gi, width))
  barcodes <- strsplit(head.line, ",", fixed = TRUE)[[1]][-1]
  rownames(M) <- genes.all; colnames(M) <- barcodes
  cat("matrix:", gi, "genes x", width, "cells\n")
  M
}

files <- c(
  Young = "raw/GSM5098741_YouHM5_counts.csv.gz",
  Old   = "raw/GSM5098753_OldHM23_counts.csv.gz"
)

obj.list <- lapply(names(files), function(g) {
  cat("reading", g, "...\n")
  M <- read_counts_csv(files[[g]])
  o <- CreateSeuratObject(M, project = g, min.cells = 3, min.features = 200)
  o$group <- g; o$sample <- g
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
  o
})
names(obj.list) <- names(files)
print(sapply(obj.list, ncol))

obj <- merge(obj.list[[1]], obj.list[[2]])
obj <- JoinLayers(obj)
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
cat("After QC:", ncol(obj), "nuclei\n")

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- RunHarmony(obj, group.by.vars = "sample", verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
cat("Clusters:", length(unique(Idents(obj))), "\n")

# ---- 细胞类型标注（canonical marker 打分）----
markers <- list(
  Myonuclei   = c("MYH2","MYH7","ACTA1","TNNI2"),
  Satellite   = c("PAX7","MYF5","MKI67"),
  Macrophage  = c("C1QA","C1QB","C1QC","CD68","LYZ","FCGR3A"),
  DC          = c("FCER1A","CD1C"),
  Fibroblast  = c("DCN","COL1A1","PDGFRA"),
  Adipocyte   = c("ADIPOQ","PLIN1","LPL"),
  Endothelial = c("PECAM1","VWF","CLDN5"),
  T_cell      = c("CD3D","CD3E","CD8A"),
  NK          = c("NKG7","GNLY"),
  B_cell      = c("MS4A1","CD79A"),
  Schwann     = c("MBP","PLP1","NRXN1"),
  Tenocyte    = c("TNMD","SCX")
)
all.genes <- unique(unlist(markers))
avg <- AverageExpression(obj, features = all.genes, layer = "data", group.by = "seurat_clusters")$RNA
score.mat <- sapply(markers, function(g) colMeans(avg[intersect(g, rownames(avg)), , drop = FALSE]))
assign <- apply(score.mat, 1, function(x) names(which.max(x)))
cluster.annot <- setNames(assign, rownames(score.mat))
print(cluster.annot)
keys <- paste0("g", as.character(Idents(obj)))
celltype.vec <- unname(cluster.annot[keys])
if (any(is.na(celltype.vec))) stop("unannotated clusters: ", paste(unique(keys[is.na(celltype.vec)]), collapse=","))
obj$celltype <- factor(celltype.vec)

# ---- 核心检验：C1QA 阳性核是否落在巨噬细胞 ----
c1qa <- FetchData(obj, vars = c("C1QA","CD68","LYZ","C1QB","C1QC","celltype","group"))
c1qa.pos <- subset(c1qa, C1QA > 0)
res <- data.frame(
  C1QA_positive_nuclei      = nrow(c1qa.pos),
  of_which_Macrophage       = sum(c1qa.pos$celltype == "Macrophage"),
  pct_in_Macrophage         = round(100*mean(c1qa.pos$celltype == "Macrophage"), 1),
  mean_C1QA_in_Macrophage   = round(mean(c1qa.pos$C1QA[c1qa.pos$celltype == "Macrophage"]), 2),
  corr_C1QA_CD68            = round(cor(c1qa$C1QA, c1qa$CD68, method = "spearman"), 3),
  corr_C1QA_LYZ             = round(cor(c1qa$C1QA, c1qa$LYZ, method = "spearman"), 3)
)
print(res)
write.csv(res, "output/probe_C1QA_check.csv", row.names = FALSE)

# ---- 细胞构成 Old vs Young ----
comp <- obj@meta.data %>% count(group, celltype) %>% group_by(group) %>%
  mutate(pct = round(100*n/sum(n), 1)) %>% arrange(group, desc(pct))
print(as.data.frame(comp))
write.csv(comp, "output/probe_composition.csv", row.names = FALSE)

# ---- hub 基因各细胞类型平均表达 ----
hub <- c("C1QA","NFKBIA","DTL","CTNNBIP1")
hub.avg <- AverageExpression(obj, features = hub, layer = "data", group.by = "celltype")$RNA
hub.avg$celltype <- rownames(hub.avg)
write.csv(hub.avg, "output/probe_hubgene_by_celltype.csv", row.names = FALSE)
print(hub.avg)

# ---- 出图 ----
p1 <- DimPlot(obj, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("GSE167186 probe: YouHM5(20y) + OldHM23(90y)")
p2 <- DotPlot(obj, features = hub, group.by = "celltype") + coord_flip() + ggtitle("Hub genes by cell type")
p3 <- FeaturePlot(obj, features = "C1QA", order = TRUE)
p4 <- DotPlot(obj, features = c("C1QA","C1QB","C1QC","CD68","LYZ"), group.by = "group") + coord_flip() + ggtitle("Macrophage markers Young vs Old")
ggsave("output/probe_UMAP_celltype.pdf", p1, width = 9, height = 6)
ggsave("output/probe_DotPlot_hubgenes.pdf", p2, width = 7, height = 7)
ggsave("output/probe_FeaturePlot_C1QA.pdf", p3, width = 10, height = 5)
ggsave("output/probe_Macrophage_YoungvsOld.pdf", p4, width = 7, height = 4)
saveRDS(obj, "output/probe_seurat.rds")
cat("PROBE DONE\n")
