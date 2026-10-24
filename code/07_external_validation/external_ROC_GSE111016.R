# =====================================================================
# 外部验证集 ROC — GSE111016 (Singapore Sarcopenia Study, RNA-seq)
# 训练集 GSE1428 的 4 关键基因: C1QA / NFKBIA / DTL / CTNNBIP1
# 输入: external/GSE111016_lcpm_sss.csv.gz (log2 CPM, ENSG 行名)
#       external/GSE111016_series_matrix.txt.gz (sarcopenia status)
# 输出: external/output/External_ROC.pdf / External_AUC.csv
# =====================================================================
rm(list = ls())
library(pROC)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)

genes <- c('C1QA', 'NFKBIA', 'DTL', 'CTNNBIP1')

# ---- 1. 读取表达矩阵 --------------------------------------------------------
lcpm <- read.csv(gzfile('GSE111016_lcpm_sss.csv.gz'), row.names = 1, check.names = FALSE)
cat('lcpm 维度:', dim(lcpm), '\n')

# ---- 2. 解析分组 ------------------------------------------------------------
sm <- readLines(gzfile('GSE111016_series_matrix.txt.gz'))
get_line <- function(tag) {
  ln <- grep(paste0('^!', tag), sm, value = TRUE)[1]
  parts <- strsplit(ln, '\t')[[1]][-1]
  gsub('"', '', parts)
}
titles <- get_line('Sample_title')          # "Sample 1 [sss]" ...
chars  <- sm[grepl('^!Sample_characteristics_ch1', sm)]
status <- NULL
for (ln in chars) {
  vals <- gsub('"', '', strsplit(ln, '\t')[[1]][-1])
  if (any(grepl('sarcopenia status', vals))) { status <- vals; break }
}
stopifnot(!is.null(status))
group  <- ifelse(grepl(': yes', status), 'SP', 'Control')

# 样本列名 "Sample 1" 与 titles 顺序一致
sample_id <- sub(' \\[.*\\]', '', titles)   # "Sample 1"
stopifnot(all(sample_id %in% colnames(lcpm)))
lcpm <- lcpm[, sample_id]
cat('分组:', paste(names(table(group)), table(group), collapse = ' '), '\n')

# ---- 3. ENSG -> SYMBOL ------------------------------------------------------
ens <- gsub('\\..*$', '', rownames(lcpm))   # 去掉版本号
map <- AnnotationDbi::select(org.Hs.eg.db,
                             keys = ens, keytype = 'ENSEMBL',
                             columns = 'SYMBOL')
map <- map[!duplicated(map$ENSEMBL) & !is.na(map$SYMBOL), ]
sym <- map$SYMBOL[match(ens, map$ENSEMBL)]
cat('注释到基因符号:', sum(!is.na(sym)), '/', length(sym), '\n')

# 去除未注释基因；重复符号保留平均表达最高的一行
keep <- !is.na(sym)
lcpm <- lcpm[keep, ]
sym  <- sym[keep]
if (anyDuplicated(sym)) {
  ord <- order(-rowMeans(lcpm, na.rm = TRUE))
  lcpm <- lcpm[ord, ]
  sym  <- sym[ord]
  first <- !duplicated(sym)
  lcpm <- lcpm[first, ]
  sym  <- sym[first]
}
rownames(lcpm) <- sym

# 目标基因是否存在
genes_found <- genes[genes %in% rownames(lcpm)]
cat('目标基因命中:', paste(genes, genes %in% rownames(lcpm)), '\n')
if (!all(genes %in% rownames(lcpm))) {
  cat('外部集中未检出:', paste(setdiff(genes, genes_found)), '（表达量过低/未检测，无法评估 ROC）\n')
}

# ---- 4. ROC 分析 ------------------------------------------------------------
# 单基因 ROC：各基因仅剔除自身 NA 样本
auc_tab <- data.frame()
roc_list <- list()
for (g in genes_found) {
  ok <- !is.na(lcpm[g, ])
  d1 <- data.frame(group = factor(group[ok], levels = c('Control', 'SP')),
                   x = as.numeric(lcpm[g, ok]))
  r <- roc(group ~ x, data = d1, auc = TRUE, ci = TRUE, quiet = TRUE)
  roc_list[[g]] <- r
  auc_tab <- rbind(auc_tab, data.frame(
    Dataset = 'GSE111016 (external)', Gene = g,
    AUC = round(as.numeric(r$auc), 3),
    CI_low = round(as.numeric(r$ci)[1], 3),
    CI_high = round(as.numeric(r$ci)[3], 3),
    Direction = r$direction, P = format.pval(r$auc, digits = 3)
  ))
}
# 联合模型（外部集内重新拟合 Logistic，仅作联合表现参考）
dat <- as.data.frame(t(lcpm[genes_found, ]))
dat$group <- factor(group, levels = c('Control', 'SP'))
dat <- na.omit(dat)
cat('联合模型样本数:', nrow(dat), '\n')
m <- glm(group ~ ., data = dat, family = 'binomial')
pr <- as.numeric(predict(m, type = 'response'))
r4 <- roc(dat$group, pr, auc = TRUE, ci = TRUE, quiet = TRUE)
auc_tab <- rbind(auc_tab, data.frame(
  Dataset = 'GSE111016 (external)',
  Gene = paste0(length(genes_found), '-gene logistic (refit)'),
  AUC = round(as.numeric(r4$auc), 3),
  CI_low = round(as.numeric(r4$ci)[1], 3),
  CI_high = round(as.numeric(r4$ci)[3], 3),
  Direction = '-', P = format.pval(r4$auc, digits = 3)
))
print(auc_tab)
write.csv(auc_tab, file = 'output/External_AUC.csv', row.names = FALSE)

# ---- 5. ROC 曲线图 ----------------------------------------------------------
dat_roc <- data.frame()
for (g in genes_found) {
  r <- roc_list[[g]]
  roc_tmp <- coords(r, x = seq(0, 1, 0.01), input = 'fpr',
                    ret = c('tpr'), transpose = FALSE)
  dat_roc <- rbind(dat_roc, data.frame(
    FPR = as.numeric(rownames(roc_tmp)), TPR = roc_tmp$tpr, Gene = g
  ))
}
cols <- c('#E64B35', '#4DBBD5', '#00A087')
p <- ggplot(dat_roc, aes(FPR, TPR, color = Gene)) +
  geom_abline(slope = 1, intercept = 0, lty = 'dashed', color = 'grey') +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = cols, labels = paste0(
    genes_found, ' (AUC=', sapply(roc_list, function(r) round(as.numeric(r$auc), 3)), ')'
  )) +
  coord_fixed() +
  theme_bw() +
  labs(
    title = 'External Validation: GSE111016 (Singapore Sarcopenia Study)',
    subtitle = paste0('Control ', sum(group == 'Control'), ' vs SP ', sum(group == 'SP'),
                      ' | DTL not detected in external set'),
    x = '1 - Specificity (FPR)', y = 'Sensitivity (TPR)', color = ''
  ) +
  theme(legend.position = c(0.65, 0.25))
ggsave(file = 'output/External_ROC.pdf', p, width = 6.5, height = 6)

cat('\n完成: output/External_ROC.pdf / output/External_AUC.csv\n')
