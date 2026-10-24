rm(list = ls())
library(tidyverse)
library(pROC)
library(ggplot2)
dat_expr <- readRDS('output/Logistic_linear_pre_GSE18995.Rds')
dat_group <- readRDS('input/GSE18995_group.Rds')
s_gene <- 'linear_pre'
ref_group <- 'pre'
test_group <- 'post'

dat_expr_show <- dat_expr
# 因子化
dat_expr_show$group <- dat_group$group %>% 
  factor(
    levels = c(ref_group, test_group)
  )

# ROC分析
res_roc_s_gene <- roc(
  as.formula(
    paste0('group ~ ', s_gene)
  ),
  data = dat_expr_show,
  # 计算AUC和置信区间
  auc = T,
  ci = T
)

# 计算最佳cutoff值
dat_roc_cutoff_s_gene <- cutoff::roc(
  dat_expr_show[, s_gene],
  dat_expr_show$group
)

# ROC曲线
p_roc_s_gene <- ggroc(
  res_roc_s_gene,
  # 坐标轴格式更改（把Specificities改成1-Specificities）
  legacy.axes = T,
  # 线条颜色
  col = '#4DBBD5'
) +
  # 添加对角线
  geom_abline(
    # 斜率
    slope = 1,
    # 截距
    intercept = 0,
    # 颜色
    color = 'grey',
    # 线型
    lty = 'dashed'
  ) +
  # 曲线下面积
  # geom_area(
  #   fill = '#4DBBD5',
  #   alpha = 0.2
  # ) +
  # 展示最佳cutoff值点
  geom_point(
    data = dat_roc_cutoff_s_gene,
    aes(
      1 - specificity,
      sensitivity
    ),
    color = '#4DBBD5'
  ) +
  # 展示最佳cutoff值点标签
  ggrepel::geom_text_repel(
    data = dat_roc_cutoff_s_gene,
    aes(
      1 - specificity,
      sensitivity
    ),
    label = paste0(
      sprintf('%0.3f', dat_roc_cutoff_s_gene$cutoff),
      ' (',
      sprintf('%0.3f', 1 - dat_roc_cutoff_s_gene$specificity),
      ', ',
      sprintf('%0.3f', dat_roc_cutoff_s_gene$sensitivity),
      ')'
    ),
    color = '#4DBBD5'
  ) +
  # 添加AUC标签
  annotate(
    'text',
    label = paste0(
      s_gene,
      ': AUC = ',
      sprintf('%0.3f', res_roc_s_gene$auc),
      ' (',
      sprintf('%0.3f', res_roc_s_gene$ci[1]),
      '-',
      sprintf('%0.3f', res_roc_s_gene$ci[3]),
      ')'
    ),
    x = 0.4,
    y = 0.1,
    adj = 0
  ) +
  # 固定xy比例
  coord_fixed() +
  theme_bw() +
  labs(
    subtitle = s_gene,
    x = '1 - Specificity (FPR)',
    y = 'Sensitivity (TPR)'
  )
# 如果前面添加了曲线下面积就需要跑下面这句
# p_roc_s_gene$layers[[3]]$mapping <- p_roc_s_gene$layers[[1]]$mapping
ggsave(file = paste0('output/ROC_', s_gene, '_GSE18995.pdf'), p_roc_s_gene, height = 6, width = 6)

