## =====================================================================
## 免疫浸润分析 - HP 集 (GSE27411, Control 6 vs HP 6)
## 方法: ssGSEA (Bindea 24 种免疫细胞签名, IOBR signature_collection)
## 输出: Fig11C 组间差异箱线图 | Fig11D 4基因x免疫细胞相关热图
## =====================================================================
options(width = 220, stringsAsFactors = FALSE)
rm(list = ls())

setwd('..')  # 到项目根: D:/GEO-limma,肌少症1

library(GSVA)      # 1.50.5
library(ggplot2)
library(ggpubr)
library(pheatmap)

## 1. 输入数据 --------------------------------------------------------------
expr <- readRDS('1-limma2/input/dat1_expr.Rds')      # 基因 x 样本
grp  <- readRDS('1-limma2/input/dat1_group.Rds')
kg   <- readRDS('4-Venn3/output/Key_gene.Rds')
genes_key <- as.character(unlist(kg))
cat('HP 表达矩阵:', nrow(expr), 'x', ncol(expr), '| 关键基因:', paste(genes_key, collapse = ', '), '\n')

## 2. Bindea 免疫细胞签名 -----------------------------------------------------
load('11-Immune2/input/signature_collection.rda')
sig_names <- grep('Bindea', names(signature_collection), value = TRUE)
sig_names <- setdiff(sig_names, c('SW480_cancer_cells_Bindea_et_al',
                                  'Normal_mucosa_Bindea_et_al',
                                  'Blood_vessels_Bindea_et_al',
                                  'Lymph_vessels_Bindea_et_al'))
sig_list <- signature_collection[sig_names]
cat('免疫细胞签名数:', length(sig_list), '\n')
cell_names <- sub('_Bindea_et_al$', '', sig_names)
names(sig_list) <- cell_names
for (i in seq_along(sig_list)) sig_list[[i]] <- unique(as.character(sig_list[[i]]))

# 签名基因与表达矩阵取交集（报告缺失）
sig_list <- lapply(sig_list, function(g) intersect(g, rownames(expr)))
sizes <- sapply(sig_list, length)
cat('各签名在表达矩阵中的基因数:\n'); print(sizes)
drop <- names(sizes)[sizes < 3]
if (length(drop) > 0) {
  cat('以下签名因平台基因覆盖不足(<3 基因)被剔除:', paste(drop, collapse = ', '), '\n')
  sig_list <- sig_list[!names(sig_list) %in% drop]
  cell_names <- setdiff(cell_names, drop)
}
stopifnot(length(sig_list) >= 15)

## 3. ssGSEA 计算浸润评分 ----------------------------------------------------
sp <- ssgseaParam(exprData = as.matrix(expr), geneSets = sig_list,
                  minSize = 3)
score <- gsva(sp, verbose = FALSE)
cat('ssGSEA 评分矩阵:', nrow(score), 'x', ncol(score), '\n')
write.csv(score, file = '11-Immune2/output/ImmuneScore_ssGSEA_HP.csv')

## 4. 组间差异 (Wilcoxon) ----------------------------------------------------
gvec <- as.character(grp$group)
if (is.null(names(gvec))) names(gvec) <- as.character(grp[[1]])
if (!is.null(rownames(grp))) names(gvec) <- rownames(grp)
gvec <- gvec[match(colnames(score), names(gvec))]
gvec <- factor(ifelse(grepl('HP|heli', gvec, ignore.case = TRUE), 'HP', 'Control'),
               levels = c('Control', 'HP'))
cat('分组:', paste(names(table(gvec)), table(gvec), collapse = ' | '), '\n')
stopifnot(!any(is.na(gvec)))

stat <- data.frame()
for (cc in cell_names) {
  p <- wilcox.test(score[cc, gvec == 'HP'], score[cc, gvec == 'Control'])$p.value
  m1 <- mean(score[cc, gvec == 'HP']); m0 <- mean(score[cc, gvec == 'Control'])
  stat <- rbind(stat, data.frame(Cell = cc, Control_mean = m0, HP_mean = m1,
                                 Diff_HP_minus_Control = m1 - m0,
                                 P = p))
}
stat$P.adj <- p.adjust(stat$P, 'BH')
write.csv(stat, file = '11-Immune2/output/Immune_groupdiff_HP.csv', row.names = FALSE)

star <- function(p) ifelse(p < 0.001, '***', ifelse(p < 0.01, '**', ifelse(p < 0.05, '*', 'ns')))
plt <- data.frame(Cell = rep(cell_names, each = ncol(score)),
                  Score = as.numeric(t(score)),
                  Group = rep(gvec, times = nrow(score)))
plt$Cell <- factor(plt$Cell, levels = cell_names[order(stat$P)])
p_fig11c <- ggboxplot(plt, x = 'Group', y = 'Score', fill = 'Group',
                      palette = c('#00A087', '#4DBBD5'), outlier.shape = 21,
                      facet.by = 'Cell', nrow = 4, ncol = 6, scales = 'free_y') +
  stat_compare_means(method = 'wilcox.test', label = 'p.signif', size = 3.4) +
  labs(title = sprintf('ssGSEA immune cell infiltration (Bindea signatures, n = %d) - GSE27411 (HP)', length(cell_names)),
       x = '', y = 'ssGSEA enrichment score') +
  theme(legend.position = 'top',
        strip.text = element_text(size = 8))
ggsave('11-Immune2/output/Fig11C_immune_boxplot_HP.pdf', p_fig11c,
       width = 12, height = 7)

## 5. 关键基因 x 免疫细胞 Spearman 相关 ----------------------------------------
rn <- intersect(genes_key, rownames(expr))
stopifnot(length(rn) == length(genes_key))
cor_mat <- matrix(NA_real_, length(genes_key), length(cell_names),
                  dimnames = list(genes_key, cell_names))
p_mat <- cor_mat
for (g in genes_key) for (cc in cell_names) {
  ct <- cor.test(as.numeric(expr[g, ]), score[cc, ], method = 'spearman')
  cor_mat[g, cc] <- as.numeric(ct$estimate)
  p_mat[g, cc]   <- ct$p.value
}
adj_p <- apply(p_mat, 2, p.adjust, 'BH')
stars <- matrix(ifelse(adj_p < 0.001, '***', ifelse(adj_p < 0.01, '**',
         ifelse(adj_p < 0.05, '*', 'ns'))), nrow = nrow(p_mat), dimnames = dimnames(p_mat))
cor_df <- data.frame(Gene = rep(rownames(cor_mat), times = ncol(cor_mat)),
                     Cell = rep(colnames(cor_mat), each = nrow(cor_mat)),
                     Spearman_r = as.numeric(cor_mat),
                     P = as.numeric(p_mat),
                     P.adj_BH_within_cell = as.numeric(adj_p),
                     Sig = as.vector(stars))
write.csv(cor_df, file = '11-Immune2/output/Keygene_immune_correlation_HP.csv', row.names = FALSE)

ph_mat <- t(cor_mat)  # 细胞 x 基因
num_lab <- matrix(paste0(sprintf('%.2f', ph_mat),
                         matrix(stars, nrow = nrow(ph_mat), ncol = ncol(ph_mat))),
                  nrow = nrow(ph_mat), dimnames = dimnames(ph_mat))
pdf('11-Immune2/output/Fig11D_keygene_immune_corr_HP.pdf', width = 10, height = 8)
pheatmap(ph_mat, display_numbers = num_lab, fontsize_number = 9,
         color = colorRampPalette(c('#3C5488', 'white', '#E64B35'))(100),
         cluster_rows = TRUE, cluster_cols = FALSE,
         main = 'Spearman correlation: key genes vs immune cells (HP)',
         border_color = 'grey80')
dev.off()

cat('HP 集免疫浸润完成: Fig11C / Fig11D / 3 个 CSV 已输出\n')
