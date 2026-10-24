rm(list = ls())
library(tidyverse)
library(randomForest)
dat_expr <- readRDS('input/dat1_expr.Rds')
dat_group <- readRDS('input/dat1_group.Rds')
# NOTE (archival correction 2026-10-11): in the primary analysis the random forest in the
# H. pylori line was applied to the FULL 14-gene candidate pool, not to the SVM-RFE-retained
# genes. Evidence: Supplementary Table S2 records rf_retained = 14 (all pool genes) in this
# line, including COL16A1 and CP, which the SVM-RFE step did not retain (svm_retained = 0).
# The earlier version of this archived script read output/gene_svm.Rds and would not
# reproduce the primary result; the input is therefore corrected to the full pool below
# (pool construction identical to 3e-svmrfe_HP.R). Threshold unchanged: IncNodePurity > 0.05.
gene_show2 <- intersect(readRDS('input/ModuleGenes1.Rds'), readRDS('input/ModuleGenes.Rds'))
gene_show1 <- readRDS('input/SP-HP.Rds')
gene_show <- intersect(gene_show1, gene_show2)

ref_group <- 'Control'
test_group <- 'HP'



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
dat_rf_importance <- dat_rf_importance0[which(dat_rf_importance0$IncNodePurity > 0.05), ]
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

saveRDS(gene_show_rf, file = 'output/gene1_rf.Rds')

