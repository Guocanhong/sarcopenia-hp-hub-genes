# =====================================================================
# 补充分析（最终版）：GSE111016 118 common DEGs 共识聚类分子亚型
# 距离度量经 3-compare_distance.R 对比（euclidean/pearson/spearman PAC），
# pearson 相关性距离 PAC 最低（k=2: 0.586），予以采用
# 输出: output/S3*.pdf + S3_cluster_data.Rds（Python 端转 PNG 拼图）
# =====================================================================
rm(list = ls())
set.seed(2026)
suppressMessages({
  library(ConsensusClusterPlus)
  library(ggplot2)
  library(AnnotationDbi); library(org.Hs.eg.db)
})

# ---- 1. 读取表达矩阵与分组 ----
lcpm <- read.csv(gzfile('GSE111016_lcpm_sss.csv.gz'), row.names = 1, check.names = FALSE)
sm <- readLines(gzfile('GSE111016_series_matrix.txt.gz'))
get_line <- function(tag) gsub('"', '', strsplit(grep(paste0('^!', tag), sm, value = TRUE)[1], '\t')[[1]][-1])
titles <- get_line('Sample_title')
chars <- sm[grepl('^!Sample_characteristics_ch1', sm)]
status <- NULL
for (ln in chars) {
  vals <- gsub('"', '', strsplit(ln, '\t')[[1]][-1])
  if (any(grepl('sarcopenia status', vals))) { status <- vals; break }
}
group <- ifelse(grepl(': yes', status), 'Sarcopenia', 'Control')
names(group) <- sub(' \\[.*\\]', '', titles)
lcpm <- lcpm[, names(group)]

# ---- 2. ENSG -> SYMBOL ----
ens <- gsub('\\..*$', '', rownames(lcpm))
map <- AnnotationDbi::select(org.Hs.eg.db, keys = ens, keytype = 'ENSEMBL', columns = 'SYMBOL')
map <- map[!duplicated(map$ENSEMBL) & !is.na(map$SYMBOL), ]
sym <- map$SYMBOL[match(ens, map$ENSEMBL)]
keep <- !is.na(sym); l2 <- lcpm[keep, ]; s2 <- sym[keep]
if (anyDuplicated(s2)) {
  o <- order(-rowMeans(l2)); l2 <- l2[o, ]; s2 <- s2[o]
  f <- !duplicated(s2); l2 <- l2[f, ]; s2 <- s2[f]
}
rownames(l2) <- s2

# ---- 3. 118 common DEGs 交集 + NA 清理 ----
common <- read.csv('../2-Venn/output/SP-HP.CSV', row.names = 1)$x
l2 <- l2[, colSums(is.na(l2)) == 0]
hit <- intersect(common, rownames(l2))
hit <- hit[rowSums(is.na(l2[hit, ])) == 0]
drop_smp <- setdiff(names(group), colnames(l2))
cat('剔除样本:', drop_smp, '| 有效:', ncol(l2), '例 | common DEGs:', length(hit), '\n')
g39 <- group[colnames(l2)]
mat <- t(scale(t(l2[hit, ])))

# ---- 4. ConsensusClusterPlus（pearson 距离, k=2..6, 1000 次） ----
ccp <- ConsensusClusterPlus(mat, maxK = 6, reps = 1000, pItem = 0.8, pFeature = 0.8,
                            clusterAlg = 'hc', distance = 'pearson',
                            innerLinkage = 'ward.D2', finalLinkage = 'ward.D2',
                            title = 'output/ccp_final', plot = 'pdf', seed = 2026)
kvals <- 2:6
pac <- sapply(kvals, function(k) {
  cm <- ccp[[k]]$consensusMatrix
  sum(cm > 0.1 & cm < 0.9) / (nrow(cm)^2 - nrow(cm))
})
cat('PAC:\n'); print(data.frame(k = kvals, PAC = round(pac, 4)))

# ---- 5. 面板 A：consensus matrix k=2（经典 CCP 输出样式，仿 Zhao et al. Fig.8A：
#           顶部树状图 + 亚型颜色条 + 右上图例 + 白→深蓝热图，按 consensusTree 排序） ----
cl2 <- ccp[[2]]$consensusClass
cl2 <- factor(cl2, labels = c('subtype I', 'subtype II'))
tree2 <- ccp[[2]]$consensusTree
# CCP 返回的 consensusTree 在 base plot.hclust 下报 C 层错误，
# 按 CCP 同样规则（1 - consensusMatrix, ward.D2）重建等价树
cmRaw <- ccp[[2]]$consensusMatrix
tree2 <- stats::hclust(stats::as.dist(1 - cmRaw), method = 'ward.D2')
tree2$labels <- rep('', nrow(cmRaw))   # 隐藏叶标签（as.dist 无维度名时 labels 为空，需显式置空）
ord <- tree2$order
cm2 <- cmRaw[ord, ord]
cl_ord <- cl2[ord]
nS <- length(ord)
col_cluster <- c('subtype I' = '#9ECAE1', 'subtype II' = '#2171B5')
ramp <- colorRampPalette(c('white', '#08306B'))(100)

pdf('output/S3A_consensus_matrix.pdf', width = 6.2, height = 5.8)
layout(rbind(c(1, 2), c(3, 4), c(5, 6)),
       heights = c(1.15, 0.22, 7.5), widths = c(6.3, 1.15))
# (1,1) 顶部树状图 + 标题
par(mar = c(0.2, 0.2, 3.0, 0))
plot(tree2, hang = -1, ann = FALSE, axes = FALSE, yaxt = 'n', xaxt = 'n')
title('consensus matrix  k = 2', cex.main = 1.2, font.main = 2, xpd = NA)
# (1,2) 图例：亚型颜色方块
par(mar = c(0, 0.4, 2.2, 0.4))
plot.new()
legend('topleft', legend = names(col_cluster), fill = unname(col_cluster),
       border = 'grey40', bty = 'o', cex = 0.95, xjust = 0, y.intersp = 1.1)
# (2,1) 亚型颜色条（与树状图/热图同 x 轴对齐）
par(mar = c(0, 0.2, 0, 0))
plot.new()
plot.window(xlim = c(0.5, nS + 0.5), ylim = c(0, 1))
seg <- rle(as.character(cl_ord))
xs <- cumsum(c(0.5, seg$lengths))
for (i in seq_along(seg$lengths)) {
  rect(xs[i], 0, xs[i + 1], 1, col = unname(col_cluster[seg$values[i]]), border = NA)
}
# (2,2) 右侧空格
par(mar = c(0, 0, 0, 0)); plot.new()
# (3,1) 共识矩阵热图（白→深蓝，无轴、无块内文字）
par(mar = c(0.2, 0.2, 0.2, 0.2))
image(x = 1:nS, y = 1:nS, t(cm2), col = ramp, zlim = c(0, 1),
      axes = FALSE, xlab = '', ylab = '')
box(col = 'grey70', lwd = 0.6)
# (3,2) 右下空格
par(mar = c(0, 0, 0, 0)); plot.new()
dev.off()

# ---- 6. 面板 B/C：CDF 曲线与 delta area ----
tri <- function(cm) cm[upper.tri(cm)]
cdf_dat <- do.call(rbind, lapply(kvals, function(k) {
  v <- sort(tri(ccp[[k]]$consensusMatrix))
  data.frame(k = k, x = v, y = seq_along(v) / length(v))
}))
auc <- sapply(kvals, function(k) mean(tri(ccp[[k]]$consensusMatrix)))
dK <- -diff(auc)   # 面积随 k 递减，取正向增量（同 CCP 惯例形态）
pB <- ggplot(cdf_dat, aes(x, y, colour = factor(k))) +
  geom_line(linewidth = 0.7) +
  scale_colour_brewer(palette = 'Spectral', name = 'k') +
  labs(x = 'consensus index', y = 'CDF', title = 'consensus CDF') +
  theme_bw(base_size = 12) + theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
ggsave('output/S3B_CDF.pdf', pB, width = 4.8, height = 4.2)
pC <- ggplot(data.frame(k = 3:6, dK = dK), aes(k, dK)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2.4) +
  scale_x_continuous(breaks = 2:6) +
  labs(x = 'k', y = expression(Delta*' area under CDF curve'),
       title = 'Delta area') +
  theme_bw(base_size = 12) + theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
ggsave('output/S3C_delta_area.pdf', pC, width = 4.8, height = 4.2)

# ---- 7. 面板 D：PCA ----
pca <- prcomp(t(mat))
df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], subtype = cl2, status = g39)
var <- pca$sdev^2 / sum(pca$sdev^2) * 100
pD <- ggplot(df, aes(PC1, PC2, colour = subtype, shape = status)) +
  geom_point(size = 3.4, alpha = 0.9) +
  stat_ellipse(aes(PC1, PC2, group = subtype), inherit.aes = FALSE,
               linetype = 2, linewidth = 0.4) +
  scale_colour_manual(values = c('subtype I' = '#E18727', 'subtype II' = '#20809E')) +
  labs(x = sprintf('PC1 (%.1f%%)', var[1]), y = sprintf('PC2 (%.1f%%)', var[2]),
       colour = 'Subtype', shape = 'Status') +
  theme_bw(base_size = 13)
ggsave('output/S3D_PCA.pdf', pD, width = 6.2, height = 4.8)

# ---- 8. 面板 E：亚型 x 疾病状态 ----
tb <- table(cl2, g39)
cat('--- crosstab ---\n'); print(tb)
ft <- fisher.test(tb); cat('Fisher P =', signif(ft$p.value, 3), '\n')
write.csv(data.frame(subtype = rownames(tb), tb, check.names = FALSE),
          'output/S3_crosstab_subtype_status.csv', row.names = FALSE)
tabdf <- as.data.frame(tb)
colnames(tabdf) <- c('Subtype', 'Status', 'n')
pE <- ggplot(tabdf, aes(Subtype, n, fill = Status)) +
  geom_col(position = 'dodge', width = 0.6, colour = 'black', linewidth = 0.3) +
  geom_text(aes(label = n), position = position_dodge(0.6), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c('Control' = '#8AB17D', 'Sarcopenia' = '#B6535B')) +
  labs(y = 'Number of samples', title = sprintf('Fisher exact P = %.3f', ft$p.value)) +
  theme_bw(base_size = 12) + theme(legend.position = 'top',
                                   plot.title = element_text(hjust = 0.5, face = 'bold'))
ggsave('output/S3E_crosstab.pdf', pE, width = 4.6, height = 4.2)

# ---- 9. 面板 F：hub 基因按亚型箱线图 ----
hubs <- intersect(c('C1QA', 'NFKBIA', 'DTL', 'CTNNBIP1'), rownames(l2))
long <- data.frame(
  gene = factor(rep(hubs, each = ncol(l2)), levels = hubs),
  expr = as.numeric(t(l2[hubs, ])),
  subtype = rep(cl2, times = length(hubs))
)
stat_rows <- do.call(rbind, lapply(hubs, function(g) {
  d <- long[long$gene == g, ]
  w <- suppressWarnings(wilcox.test(d$expr ~ d$subtype, exact = FALSE))
  data.frame(gene = g, P = signif(w$p.value, 3),
             med1 = median(d$expr[d$subtype == 'subtype I']),
             med2 = median(d$expr[d$subtype == 'subtype II']))
}))
cat('--- hub gene wilcoxon ---\n'); print(stat_rows)
write.csv(stat_rows, 'output/S3_hubgene_wilcox.csv', row.names = FALSE)
pF <- ggplot(long, aes(subtype, expr, fill = subtype)) +
  geom_boxplot(outlier.size = 1.2, linewidth = 0.35, width = 0.55) +
  geom_jitter(width = 0.12, size = 1.1, alpha = 0.6) +
  facet_wrap(~gene, scale = 'free_y', nrow = 1) +
  scale_fill_manual(values = c('subtype I' = '#E18727', 'subtype II' = '#20809E')) +
  labs(y = expression(log[2]*' CPM'), fill = 'Subtype') +
  theme_bw(base_size = 12) + theme(legend.position = 'top',
                                   axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank())
ggsave('output/S3F_hubgene_boxplot.pdf', pF, width = 9.5, height = 4.0)

# ---- 10. 存档 ----
saveRDS(list(cl2 = cl2, group = g39, samples = colnames(l2), pac = pac,
             crosstab = tb, fisher_p = ft$p.value, genes_used = hit,
             stat_rows = stat_rows, cm2 = cm2),
        'output/S3_cluster_data.Rds')
cat('ALL DONE\n')
