rm(list = ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)
dat_expr <- readRDS('input/dat_expr.Rds')
dat_group <- readRDS('input/dat_group.Rds')
gene_show <- readRDS('input/Key_gene.Rds')
ref_group <- 'Control'
test_group <- 'SP'

dat_expr <- data.frame(dat_expr)

# 分组比较图 -------------------------------------------------------------------

# 提取基因
dat_expr_show <- dat_expr[gene_show, ] %>% 
  t() %>% 
  as.data.frame()

# 加一列分组
dat_expr_show_long <- dat_expr_show
dat_expr_show_long$group <- dat_group$group
# 宽型数据转成长型数据
dat_expr_show_long <- reshape2::melt(
  dat_expr_show_long,
  variable.name = 'gene',
  value.name = 'expr'
)

# 因子化
dat_expr_show_long$gene <- factor(
  dat_expr_show_long$gene,
  levels = gene_show
)
dat_expr_show_long$group <- factor(
  dat_expr_show_long$group,
  levels = c(ref_group, test_group)
)

# 分组比较图
p_vlnplot <- ggplot(
  dat_expr_show_long,
  aes(
    gene,
    expr,
    fill = group
  )
) +
  # 小提琴图
  geom_violin() +
  # 添加显著性标志
  stat_compare_means(
    aes(group = group),
    # 显著性标志
    label = 'p.signif',
    # 统计方法
    method = 'wilcox.test',
    method.args = list(exact = FALSE),
    # 不是配对样本
    paired = F,
    # 星号
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c('***', '**', '*', 'ns')
    )
  ) +
  theme_bw() +
  labs(
    x = element_blank(),
    y = 'Expression',
    fill = 'Group'
  ) +
  scale_fill_manual(
    values = c('#4DBBD5', '#E64B35')
  )
ggsave(file = 'output/Fig5F.pdf', p_vlnplot, width = 10, height = 6)
