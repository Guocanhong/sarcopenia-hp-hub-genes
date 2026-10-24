rm(list = ls())
library(tidyverse)
library(randomForest)
dat_expr <- readRDS('input/dat_expr.Rds')
dat_group <- readRDS('input/dat_group.Rds')
gene_show <- readRDS('output/gene_svm.Rds')
ref_group <- 'Control'
test_group <- 'SP'


# RF ----------------------------------------------------------------------

# 提出基因并整理数据
dat_expr_show <- dat_expr[gene_show, ] %>% 
  t() %>% 
  as.data.frame() %>% 
  cbind(
    group = dat_group$group %>% 
      case_match(
        ref_group ~ 0,
        test_group ~ 1
      ),
    .
  )

# 种子
set.seed(2025)
# RF模型训练
res_rf <- randomForest(
  x = dat_expr_show[, -1],
  y = dat_expr_show[, 1],
  importance = T,
  # 决策树
  ntree = 1000,
  proximity = T
)

# 提取特征重要性
dat_rf_importance0 <- importance(
  res_rf,
  scale = F
) %>% 
  as.data.frame()
# 以平均节点纯度筛选变量
dat_rf_importance <- dat_rf_importance0[which(dat_rf_importance0$IncNodePurity > 0.2), ]
gene_show_rf <- rownames(dat_rf_importance)  

# 平均节点纯度
pdf(file = 'output/Fig4E.pdf', width = 6, height = 5)
varImpPlot(
  res_rf,
  type = 2,
  scale = F
)
dev.off()


# 数据保存 --------------------------------------------------------------------

saveRDS(gene_show_rf, file = 'output/gene_rf.Rds')

