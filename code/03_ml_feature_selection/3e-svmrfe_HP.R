rm(list = ls())
library(tidyverse)
library(e1071)
dat_expr <- readRDS('input/dat_expr.Rds')
dat_group <- readRDS('input/dat_group.Rds')
gene_show2 <- intersect(readRDS('input/ModuleGenes1.Rds'), readRDS('input/ModuleGenes.Rds'))
gene_show1 <- readRDS('input/SP-HP.Rds')
gene_show <- intersect(gene_show1,gene_show2)
ref_group <- 'Control'
test_group <- 'SP'


# SVM ---------------------------------------------------------------------

# 导入一些关于SVM的函数
source('input/msvmRFE.R')

# 提出基因并整理数据
dat_expr_show <- dat_expr[gene_show, ] %>% 
  t() %>% 
  as.data.frame() %>% 
  cbind(
    group = dat_group$group %>% 
      # 因子化
      factor(
        levels = c(ref_group, test_group)
      ),
    .
  )

# 种子
set.seed(2025)
# 特征选择
svmRFE(
  dat_expr_show,
  # 递归迭代次数
  k = 10,
  halve.above = 100
)
# 用于交叉验证的折叠分配
res_folds_svm <- rep(
  1:10,
  len = nrow(dat_expr_show)
)[sample(nrow(dat_expr_show))]
res_folds_svm <- lapply(
  1:10,
  function(x) which(res_folds_svm == x)
)

# SVM分析
res_svm <- lapply(
  res_folds_svm,
  svmRFE.wrap,
  dat_expr_show,
  k = 10,
  # 递归迭代次数
  halve.above = 100
)
# 输出SVM结果
dat_top_feature_svm <- WriteFeatures(
  res_svm,
  dat_expr_show,
  save = F
)
rownames(dat_top_feature_svm) <- dat_top_feature_svm$FeatureName

# 评估SVM模型的泛化误差
res_svm_featsweep <- lapply(
  1:(ncol(dat_expr_show) - 1),
  FeatSweep.wrap,
  res_svm,
  dat_expr_show
)

# save(res_svm_featsweep, file = 'Rda/res_svm_featsweep.Rda')

info_svm <- min(prop.table(table(dat_expr_show[, 1])))
# 提取错误率
error_svm <- sapply(
  res_svm_featsweep,
  function(x) ifelse(is.null(x), NA, x$error)
)
# 提取SVM筛选到的基因
gene_show_svm <- dat_top_feature_svm[1:which.min(error_svm), ] %>% 
  rownames()

# 错误率
pdf(file = 'output/Fig4C.pdf', width = 6, height = 5)
PlotErrors(
  errors = error_svm,
  no.info = info_svm
)
dev.off()

# 正确率
pdf(file = 'output/Fig4D.pdf', width = 6, height = 5)
Plotaccuracy(
  1 - error_svm,
  no.info = info_svm
)
dev.off()


# 数据保存 --------------------------------------------------------------------

saveRDS(gene_show_svm, file = 'output/gene_svm.Rds')
