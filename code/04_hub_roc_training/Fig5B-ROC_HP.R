rm(list = ls())
library(tidyverse)
library(pROC)
library(ggplot2)
dat_expr <- readRDS('input/dat1_expr.Rds')
dat_group <- readRDS('input/dat1_group.Rds')
gene_show <- readRDS('input/Key_gene.Rds')
ref_group <- 'Control'
test_group <- 'HP'
s_gene <- 'DTL'

dat_expr <- data.frame(dat_expr)

# ROC曲线 -------------------------------------------------------------------

# 提取基因
dat_expr_s_gene <- cbind(
  group = dat_group$group,
  t(dat_expr[s_gene, ])
) %>% 
  as.data.frame()
dat_expr_s_gene[, 2] <- as.numeric(dat_expr_s_gene[, 2])

# 因子化
dat_expr_s_gene$group <- factor(
  dat_expr_s_gene$group,
  levels = c(ref_group, test_group)
)

# ROC分析
res_roc_s_gene <- roc(
  as.formula(
    paste0('group ~ ', s_gene)
  ),
  data = dat_expr_s_gene,
  # 计算AUC和置信区间
  auc = T,
  ci = T
)

# ROC曲线
p_roc_s_gene <- ggroc(
  res_roc_s_gene,
  # 坐标轴格式更改（把Specificities改成1-Specificities）
  legacy.axes = T,
  col = '#4DBBD5'
) +
  # 添加对角线
  geom_abline(
    slope = 1,
    intercept = 0,
    color = 'grey',
    lty = 'dashed'
  ) +
  # 添加AUC标签
  annotate(
    'text',
    label = paste0(
      s_gene,
      ': AUC = ',
      round(res_roc_s_gene$auc, 3),
      ' (',
      round(res_roc_s_gene$ci[1], 3),
      '-',
      round(res_roc_s_gene$ci[3], 3),
      ')'
    ),
    x = 0.4,
    y = 0.1
  ) +
  # 固定xy比例
  coord_fixed() +
  theme_bw() +
  labs(
    subtitle = s_gene,
    x = '1 - Specificity (FPR)',
    y = 'Sensitivity (TPR)'
  )
ggsave(file = 'output/Fig6E.pdf', p_roc_s_gene, height = 6, width = 6)

# Fig5B其他基因的ROC曲线绘制方法相同，只需要将基因替换成对应基因即可
