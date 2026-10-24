rm(list = ls())
library(tidyverse)
library(limma)
library(ggplot2)
dat_expr <- readRDS('input/dat_expr.Rds')
dat_group <- readRDS('input/dat_group.Rds')
# 过滤低表达基因
dat_expr <- dat_expr[
  apply(dat_expr, 1, function(x) sum(x > 0) > ncol(dat_expr) * 0.50),
]
ref_group <- 'Control'
test_group <- 'HP'


# 差异分析 --------------------------------------------------------------------

# 创建模型矩阵
mat_design <- model.matrix(~ 0 + factor(dat_group$group))
# 重新命名行名和列名
colnames(mat_design) <- levels(factor(dat_group$group))
rownames(mat_design) <- dat_group$sample

# 使用线性模型拟合表达数据
fit_limma <- lmFit(dat_expr, mat_design)

# 为组间差异创建对照
contrast <- paste0(test_group, '-', ref_group)
# 一定要实验组在前，对照组在后

# 创建线性对比矩阵
mat_contrast <- makeContrasts(contrasts = contrast, levels = mat_design)

# 对拟合的线性模型进行再次拟合
fit_limma <- contrasts.fit(fit_limma, mat_contrast)

# Bayes调整
fit_limma <- eBayes(fit_limma)

# 提取差异分析结果
dat_res_diff <- topTable(
  fit_limma,
  coef = contrast,
  n = Inf
) %>% 
  na.omit()

# saveRDS(dat_res_diff, file = 'dat_res_diff.Rds')

# 将为0的p值改成除0外最小的p值，以防止p值为0取不到对数
dat_res_diff$P.Value[dat_res_diff$P.Value == 0] <- min(dat_res_diff$P.Value[dat_res_diff$P.Value != 0])
dat_res_diff$adj.P.Val[dat_res_diff$adj.P.Val == 0] <- min(dat_res_diff$adj.P.Val[dat_res_diff$adj.P.Val != 0])

# logFC阈值
logFC_cutoff <- 0
# 筛选差异基因
dat_res_diff$change <- case_when(
  dat_res_diff$P.Value < 0.05 & dat_res_diff$logFC >= logFC_cutoff ~ 'Up',
  dat_res_diff$P.Value < 0.05 & dat_res_diff$logFC <= -logFC_cutoff ~ 'Down',
  .default = 'Not'
)
# 可根据实际情况选择筛选pvalue或padj
table(dat_res_diff$change)

# 差异基因
gene_DEGs <- dat_res_diff %>% 
  # 过滤不显著的基因
  dplyr::filter(change != 'Not') %>% 
  # 按照logFC排序
  dplyr::arrange(desc(logFC)) %>% 
  rownames()


# 火山图 ---------------------------------------------------------------------

# 火山图数据
dat_volcano <- dat_res_diff[, c('logFC', 'P.Value', 'change')]
# 加一列p值的负对数
dat_volcano$logp <- -log10(dat_volcano$P.Value)

# 因子化
dat_volcano$change <- factor(
  dat_volcano$change,
  levels = c('Up', 'Down', 'Not')
)

# 需要标注的基因
# gene_show <- 'TOP2A'
# dat_label <- dat_volcano[gene_show, ]

# 火山图
p_volcano <- ggplot(
  dat_volcano,
  aes(
    logFC,
    logp,
    color = change
  )
) +
  geom_point(
    alpha = 0.6
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    lty = 2
  ) +
  geom_vline(
    xintercept = c(-logFC_cutoff, logFC_cutoff),
    lty = 2
  ) +
# # 需要标签的基因
# geom_point(
#   data = dat_label,
#   color = 'black',
#   # 点型为空心圆点
#   shape = 1
# ) +
# # 基因标签
# ggrepel::geom_label_repel(
#   data = dat_label,
#   aes(label = rownames(dat_label)),
#   color = 'black'
# ) +
  theme_bw() +
  labs(
    x = 'LogFC',
    y = '-Log10(P Value)',
    # 或'-Log10(adj. P)'
    color = 'Change'
  ) +
  scale_x_continuous(
    limits = c(-3, 3)
  ) +
  scale_color_manual(
    values = c('#E64B35', '#4DBBD5', 'grey')
  )
ggsave(file = 'output/volcano1.pdf', p_volcano, width = 7, height = 6)

# 数据保存 --------------------------------------------------------------------

saveRDS(dat_res_diff, file = 'output/HP-diff_res.Rds')
saveRDS(gene_DEGs, file = 'output/HP-DEGs.Rds')

