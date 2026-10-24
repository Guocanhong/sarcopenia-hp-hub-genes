rm(list = ls())
library(tidyverse)
library(ggvenn)
gene_show1 <- readRDS('output/gene_lasso.Rds')
gene_show2 <- readRDS('output/gene_svm.Rds')
gene_show3 <- readRDS('output/gene_rf.Rds')
gene1_name <- 'LASSO'
gene2_name <- 'SVM'
gene3_name <- 'RF'


# 交集韦恩图 ---------------------------------------------------------------

# 韦恩图
p_venn <- ggvenn(
  list(gene_show1, gene_show2, gene_show3) %>% 
    rlang::set_names(c(gene1_name, gene2_name, gene3_name)),
  # 下面的名称要与上面列表中的命名一致
  c(gene1_name, gene2_name, gene3_name),
  # 不展示比例
  show_percentage = F,
  # 透明度
  fill_alpha = 0.5,
  stroke_color = NA,
  fill_color = c('#4DBBD5', '#E64B35', '#3C5488')
)
ggsave(file = 'output/Fig4F.pdf', p_venn, width = 5, height = 5)


# 交集情况 ---------------------------------------------------------------

# 对基因取交集
gene_intersect <- Reduce(
  # 取交集
  intersect,
  list(
    gene_show1,
    gene_show2,
    gene_show3
  )
)
gene_intersect


# 数据保存 ---------------------------------------------------------------

saveRDS(gene_intersect, file = 'output/Key_gene.Rds')
