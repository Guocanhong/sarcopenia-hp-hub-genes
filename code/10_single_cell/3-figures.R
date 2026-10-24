# 3-figures.R  Hp 侧单细胞论文级图（基于 probe_seurat.rds）
suppressMessages({library(Seurat); library(ggplot2); library(dplyr); library(patchwork); library(data.table)})
obj <- readRDS("probe_seurat.rds")
obj$group <- factor(obj$group, levels = c("HpN", "HpP"))
ct_order <- c("Epithelial", "Macrophage", "Monocyte/Neutrophil", "T-cell", "Plasma",
              "B-cell", "Mast", "Fibroblast/SMC", "Endothelial", "Other")
ct_levels <- intersect(ct_order, unique(as.character(obj$celltype)))
obj$celltype <- factor(obj$celltype, levels = ct_levels)
pal <- setNames(c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
                  "#8491B4", "#91D1C2", "#DC9440", "#7E6148", "#B09C85")[seq_along(ct_levels)], ct_levels)

dir.create("output", showWarnings = FALSE)
# ragg 无法写中文路径 -> PNG 先写到 ASCII 临时目录，再由外部拷回
TMPD <- file.path(Sys.getenv("TEMP"), "scfigs"); dir.create(TMPD, showWarnings = FALSE)
save2 <- function(name, plot, w, h) {
  ggsave(file.path(TMPD, name), plot, width = w, height = h, dpi = 300)
  cat("PNG->", file.path(TMPD, name), "\n", sep = "")
}
th <- theme_classic(base_size = 12) + theme(axis.text = element_text(colour = "black"))

## Fig A: UMAP by celltype + by group
p1 <- DimPlot(obj, group.by = "celltype", cols = pal) + th + ggtitle("Cell types") + theme(legend.title = element_blank())
p2 <- DimPlot(obj, group.by = "group", cols = c(HpN = "#3C5488", HpP = "#E64B35")) + th + ggtitle("H. pylori status")
ggsave("output/Fig_SC5_UMAP_celltype_group.pdf", p1 + p2, width = 11, height = 4.6)
save2("Fig_SC5_UMAP_celltype_group.png", p1 + p2, 11, 4.6)

## Fig B: hub genes FeaturePlot (2x2) —— 逐图加主题，避免 patchwork & 完整主题的兼容问题
DefaultAssay(obj) <- "RNA"
hub <- intersect(c("C1QA", "NFKBIA", "CTNNBIP1", "DTL"), rownames(obj))
pl <- lapply(hub, function(g) {
  q <- FeaturePlot(obj, features = g, order = TRUE, cols = c("#EEEEEE", "#B2182B"))
  q[[1]] + th + ggtitle(g) + theme(legend.position = "right",
                                   plot.title = element_text(face = "italic"))
})
fp <- patchwork::wrap_plots(pl, ncol = 2)
ggsave("output/Fig_SC6_FeaturePlot_hubgenes.pdf", fp, width = 8.5, height = 7)
save2("Fig_SC6_FeaturePlot_hubgenes.png", fp, 8.5, 7)

## Fig C: DotPlot 4 hub genes x celltype
dp <- DotPlot(obj, features = hub, group.by = "celltype", dot.scale = 9) +
  th + coord_flip() + xlab("") + ylab("Average expression") +
  scale_colour_gradient2(low = "#F7F7F7", mid = "#FDAE61", high = "#B2182B", midpoint = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
ggsave("output/Fig_SC7_DotPlot_hubgenes.pdf", dp, width = 6.8, height = 4.8)
save2("Fig_SC7_DotPlot_hubgenes.png", dp, 6.8, 4.8)
md <- dp$data; fwrite(md, "output/fig_dotplot_data.csv")

## Fig D: macrophage & plasma fraction per donor (Hp+ vs Hp-, 3v3)
frac <- fread("output/probe_celltype_frac_by_sample.csv")
fl <- melt(frac, id.vars = "sample", variable.name = "celltype", value.name = "pct")
fl$group <- ifelse(grepl("HP-N", fl$sample), "HpN", "HpP")
fl$celltype <- factor(fl$celltype, levels = ct_levels)
targets <- intersect(c("Macrophage", "Plasma", "T-cell", "Epithelial"), ct_levels)
sub <- fl[celltype %in% targets]
pd <- ggplot(sub, aes(group, pct, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = .6) +
  geom_jitter(width = .12, size = 2, aes(colour = sample)) +
  facet_wrap(~celltype, scales = "free_y") +
  scale_fill_manual(values = c(HpN = "#3C5488", HpP = "#E64B35"), guide = "none") +
  th + xlab("") + ylab("% of cells") + ggtitle("Hp- (n=3) vs Hp+ (n=3), per donor")
ggsave("output/Fig_SC8_celltype_frac_boxplot.pdf", pd, width = 9, height = 5.2)
save2("Fig_SC8_celltype_frac_boxplot.png", pd, 9, 5.2)

## 统计补充: C1QA+ macrophage per donor
mcells <- subset(obj, celltype == "Macrophage")
md2 <- data.frame(sample = mcells$sample, group = mcells$group,
                  C1QA = as.numeric(GetAssayData(mcells, layer = "data")["C1QA", ]))
st2 <- md2 %>% group_by(sample, group) %>%
  summarise(pct_C1QApos = mean(C1QA > 0) * 100, mean_C1QA = mean(C1QA), n = n(), .groups = "drop")
fwrite(st2, "output/probe_C1QApos_macro_per_sample.csv")
cat("\nC1QA+ macrophages per sample:\n"); print(as.data.frame(st2), digits = 3)
wt <- wilcox.test(pct_C1QApos ~ group, st2)
cat("MWU C1QA+ fraction in macrophages: p =", format.pval(wt$p.value), "\n")
cat("FIGS DONE\n")
