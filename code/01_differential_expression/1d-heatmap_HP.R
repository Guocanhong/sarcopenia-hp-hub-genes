rm(list = ls())
library(tidyverse)
library(pheatmap)
dat_expr <-readRDS('input/dat_expr.Rds')
dat_group <- readRDS('input/dat_group.Rds')
ref_group <- 'Control'
test_group <- 'HP'
dat_res_diff <- readRDS('output/HP-diff_res.Rds')

# 表达热图 --------------------------------------------------------------------

# 提取表达谱
gene_show <- dat_res_diff %>% 
  dplyr::filter(change != 'Not') %>% 
  dplyr::arrange(desc(logFC)) %>% 
  .[c(1:50, (nrow(.)-49):nrow(.)), ] %>% 
  rownames() %>% 
  unique()
dat_expr_show <- dat_expr[gene_show, ]

# 样本注释信息
dat_anno <- data.frame(
  Group = dat_group$group,
  row.names = dat_group$sample
)
# 样本注释颜色
list_anno_color <- c('#4DBBD5', '#E64B35')
names(list_anno_color) <- c(ref_group, test_group)
list_anno_color <- list(Group = list_anno_color)

# 表达热图
pdf(file = 'output/Fig2D.pdf', width = 8, height = 9)
pheatmap(
  dat_expr_show,
  # 对行归一化
  scale = 'row',
  # 列注释
  annotation_col = dat_anno,
  # 列注释颜色
  annotation_colors = list_anno_color,
  # 热图颜色
  color = colorRampPalette(c('#4DBBD5', 'white', '#E64B35'))(50),
  # 颜色范围（这里的50要和上面括号中的50保持一致）
  breaks = c(seq(-3, 3, length = 50)),
  # 是否对行和列进行聚类
  cluster_cols = F,
  cluster_rows = T,
  labels_col = '',
  border_color = NA
)
dev.off()
