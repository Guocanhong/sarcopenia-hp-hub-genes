# GSE25941 独立外部验证（Young vs Old, GPL570, vastus lateralis, n=36）
# 流程：解析 series matrix -> GPL570 注释映射 -> 探针合并 -> log2 -> ROC + 箱线图
rm(list = ls())
suppressPackageStartupMessages({library(pROC); library(ggplot2); library(ggpubr)})
# 运行方式：在 8-SP2 目录下 Rscript 1-External_ROC_GSE25941.R
setwd('..')  # 项目根

genes  <- c('C1QA', 'NFKBIA', 'DTL', 'CTNNBIP1')
outd   <- '8-SP2/output'; ind <- '8-SP2/input'
dir.create(outd, showWarnings = FALSE, recursive = TRUE)
dir.create(ind,  showWarnings = FALSE, recursive = TRUE)

## 1. 解析 series matrix -------------------------------------------------------
sm <- readLines(gzfile('external/GSE25941_series_matrix.txt.gz'))
get_line <- function(tag) {
  ln <- grep(paste0('^!', tag), sm, value = TRUE)[1]
  strsplit(ln, '\t')[[1]][-1] |> (\(x) gsub('"', '', trimws(x)))()
}
titles  <- get_line('Sample_title')          # Muscle_Young_Female_16.0 ...
gsm     <- get_line('Sample_geo_accession')  # GSM637513...
char_lines <- grep('^!Sample_characteristics_ch1', sm, value = TRUE)
age <- NULL
for (ln in char_lines) {
  v <- strsplit(ln, '\t')[[1]][-1] |> (\(x) gsub('"', '', trimws(x)))()
  if (any(grepl('^age:', v))) age <- sub('^age: ', '', v)
}
group   <- ifelse(age == 'Young', 'Control', 'SP')
cat('样本数:', length(gsm), '| Control(Young):', sum(group == 'Control'),
    '| SP(Old):', sum(group == 'SP'), '\n')

## 2. 表达矩阵（原始值 -> log2）------------------------------------------------
tab_start <- grep('^!series_matrix_table_begin', sm) + 1
tab_end   <- grep('^!series_matrix_table_end', sm) - 1
expr <- read.table(
  text = paste(sm[tab_start:tab_end], collapse = '\n'),
  sep = '\t', header = TRUE, quote = '', row.names = 1, check.names = FALSE
)
expr <- as.matrix(expr)
colnames(expr) <- gsub('"', '', colnames(expr))
rownames(expr) <- gsub('"', '', rownames(expr))
stopifnot(all(colnames(expr) == gsm))
cat('探针维度:', dim(expr), '\n')

## 3. GPL570 注释映射 ----------------------------------------------------------
ann_lines <- readLines(gzfile('external/GPL570.annot.gz'))
ann_lines <- ann_lines[grepl('^[0-9]', ann_lines)]           # 仅保留数据行
ann <- read.table(text = paste(ann_lines, collapse = '\n'),
                  sep = '\t', header = FALSE, quote = '', fill = TRUE,
                  comment.char = '', check.names = FALSE)
sym_map <- ann[, c(1, 3)]; colnames(sym_map) <- c('probe', 'symbol')
sym_map$symbol <- sub('///.*', '', sym_map$symbol)   # 多基因探针取第一个
sym_map <- sym_map[sym_map$symbol != '' & !is.na(sym_map$symbol), ]
sym_map <- sym_map[!duplicated(sym_map$probe), ]
hit <- match(rownames(expr), sym_map$probe)
sym <- sym_map$symbol[hit]
keep <- !is.na(sym)
expr <- expr[keep, ]; sym <- sym[keep]
cat('注释到符号:', length(sym), '/', length(hit), '\n')

# 同名基因多探针：取平均表达最高者
if (any(duplicated(sym))) {
  ord  <- order(sym, -rowMeans(expr))
  expr <- expr[ord, ]
  sym  <- sym[ord]
  first <- !duplicated(sym)
  expr <- expr[first, ]; sym <- sym[first]
  cat('重复合并后:', nrow(expr), '\n')
}
rownames(expr) <- sym

## 4. 目标基因 -----------------------------------------------------------------
miss <- setdiff(genes, rownames(expr))
if (length(miss)) stop('缺失基因: ', paste(miss, collapse = ', '))
for (g in genes) cat(g, '检出探针 OK | Old vs Young log2FC 均值差:',
                     round(mean(expr[g, group == 'SP']) - mean(expr[g, group == 'Control']), 3), '\n')
expr_log <- log2(expr + 1)

## 5. 保存 input Rds（沿用项目格式，方便复用其他脚本）---------------------------
dat_expr  <- expr_log[genes, gsm]
dat_group <- data.frame(sample = gsm, group = group, row.names = gsm)
saveRDS(dat_expr,  file.path(ind, 'dat_expr.Rds'))
saveRDS(dat_group, file.path(ind, 'dat_group.Rds'))
saveRDS(genes,     file.path(ind, 'Key_gene.Rds'))

## 6. ROC：4 基因逐个 + 合并图 --------------------------------------------------
auc_tab <- data.frame(); roc_list <- list()
dat <- data.frame(t(dat_expr), group = factor(group, levels = c('Control', 'SP')))
for (g in genes) {
  r <- roc(group ~ dat[[g]], data = dat, auc = TRUE, ci = TRUE, quiet = TRUE)
  roc_list[[g]] <- r
  auc_tab <- rbind(auc_tab, data.frame(
    Dataset = 'GSE25941 (external, Young vs Old)', Gene = g,
    AUC = round(as.numeric(r$auc), 3),
    CI_low = round(as.numeric(r$ci)[1], 3),
    CI_high = round(as.numeric(r$ci)[3], 3),
    Direction_in_SP = ifelse(mean(dat[[g]][dat$group == 'SP']) >
                             mean(dat[[g]][dat$group == 'Control']), 'up', 'down'),
    Wilcox_P = format.pval(wilcox.test(dat[[g]] ~ dat$group)$p.value, digits = 3)
  ))
}
print(auc_tab)
write.csv(auc_tab, file.path(outd, 'External_AUC_GSE25941.csv'), row.names = FALSE)

# 逐基因 ROC 图（与 ROC_C1QA.pdf 完全同款样式：同色 #E64B35、同布局同标注）
for (g in genes) {
  r <- roc_list[[g]]
  p <- ggroc(r, legacy.axes = TRUE, col = '#E64B35', linewidth = 0.9) +
    geom_abline(slope = 1, intercept = 0, color = 'grey', lty = 'dashed') +
    annotate('text', label = sprintf('%s: AUC = %.3f (%.3f-%.3f)',
             g, as.numeric(r$auc), r$ci[1], r$ci[3]),
             x = 0.4, y = 0.1, size = 3.5) +
    coord_fixed() + theme_bw() +
    labs(subtitle = paste0(g, ' | GSE25941 Young vs Old'),
         x = '1 - Specificity (FPR)', y = 'Sensitivity (TPR)')
  ggsave(file.path(outd, paste0('ROC_', g, '.pdf')), p, width = 5.5, height = 5.5)
}

# 合并 ROC 图
cols <- c('#E64B35', '#4DBBD5', '#00A087', '#3C5488')
roc_df <- do.call(rbind, lapply(genes, function(g) {
  r <- roc_list[[g]]
  cr <- coords(r, x = seq(0, 1, 0.005), input = 'fpr', ret = 'sensitivity',
               transpose = FALSE)
  data.frame(FPR = as.numeric(rownames(cr)), TPR = cr[, 1], Gene = g)
}))
p_all <- ggplot(roc_df, aes(FPR, TPR, color = Gene)) +
  geom_abline(slope = 1, intercept = 0, lty = 'dashed', color = 'grey') +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = cols, labels = paste0(
    genes, ' (AUC=', sapply(roc_list, \(r) round(as.numeric(r$auc), 3)), ')')) +
  coord_fixed() + theme_bw() +
  labs(title = 'External Validation: GSE25941 (vastus lateralis, Young vs Old)',
       subtitle = paste0('Control(Young) ', sum(group == 'Control'),
                         ' vs SP(Old) ', sum(group == 'SP')),
       x = '1 - Specificity (FPR)', y = 'Sensitivity (TPR)', color = '') +
  theme(legend.position = c(0.68, 0.22))
ggsave(file.path(outd, 'ROC_all_GSE25941.pdf'), p_all, width = 7, height = 6)

## 7. 小提琴图（log2 表达，Wilcoxon 检验，四基因合一图）------------------------
bx <- do.call(rbind, lapply(genes, function(g)
  data.frame(Gene = g, Expr = dat[[g]], Group = dat$group)))
bx$Group <- factor(bx$Group, levels = c('Control', 'SP'))
p_vio <- ggplot(bx, aes(Group, Expr, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.85, color = 'black', linewidth = 0.4,
              scale = 'width', position = position_dodge(0.9)) +
  geom_boxplot(width = 0.12, fill = 'white', alpha = 0.9,
               outlier.size = 0.8, linewidth = 0.3) +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', size = 3.2) +
  facet_wrap(~ Gene, scales = 'free_y', nrow = 1) +
  scale_fill_manual(values = c('Control' = '#00A087', 'SP' = '#E64B35')) +
  theme_bw() +
  labs(title = 'GSE25941 (Young vs Old): key gene expression',
       y = 'log2 expression', x = '', fill = '') +
  theme(legend.position = 'top')
ggsave(file.path(outd, 'Violin_4genes_GSE25941.pdf'), p_vio, width = 11, height = 3.8)

cat('完成：ROC 4 张 + 合并 ROC + 箱线图 + AUC 表 -> 8-SP2/output\n')
