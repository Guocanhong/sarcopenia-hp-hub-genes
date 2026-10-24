rm(list = ls())
library(tidyverse)
library(glmnet)
dat_expr <- readRDS('input/dat1_expr.Rds')
dat_group <- readRDS('input/dat1_group.Rds')
gene_show <- readRDS('output/gene_svm.Rds')

ref_group <- 'Control'
test_group <- 'HP'
# 用于LASSO分析的基因不能超过样本数量的一半


# 诊断LASSO -----------------------------------------------------------------

# 提取基因
dat_expr_show <- dat_expr[gene_show, ] %>% 
  t() %>% 
  as.data.frame()

# 因子化
dat_group$group <- factor(
  dat_group$group,
  levels = c(ref_group, test_group)
)

# 种子
set.seed(2025)

# LASSO回归
res_fit_lasso <- glmnet(
  dat_expr_show,
  dat_group$group,
  # 结局为二分类
  family = 'binomial',
  # 弹性网络混合参数，1表示LASSO，0表示岭回归
  alpha = 1
)

# LASSO变量轨迹图
pdf(file = 'output/Fig4A.pdf', width = 6, height = 5)
plot(
  res_fit_lasso,
  # x轴为log(lambda)
  xvar = 'lambda',
  # 不展示基因标签
  label = F,
  # 坐标轴标签
  xlab = 'Log lambda',
  ylab = 'Coefficient'
)
dev.off()

# LASSO筛选
res_cv_lasso <- cv.glmnet(
  as.matrix(dat_expr_show),
  dat_group$group,
  # 十折交叉验证
  nfolds = 10,
  # 结局为二分类
  family = 'binomial'
)

# LASSO系数筛选图
pdf(file = 'output/Fig4B.pdf', width = 6, height = 5)
plot(
  res_cv_lasso,
  # 坐标轴标签
  xlab = 'Log lambda',
  ylab = 'Binomial Deviance'
)
dev.off()

# 提取LASSO筛选到的基因
gene_show_lasso <- coef(res_cv_lasso, s = 'lambda.min') %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  dplyr::filter(s1 != 0) %>% 
  rownames() %>% 
  .[-1]


# 数据保存 --------------------------------------------------------------------

saveRDS(gene_show_lasso, file = 'output/gene_lasso.Rds')

