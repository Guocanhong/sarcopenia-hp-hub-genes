rm(list = ls())
library(tidyverse)
library(rms)
library(forestplot)
dat_expr <- readRDS('input/dat_expr.Rds')
dat_group <- readRDS('input/dat_group.Rds')
gene_show <- readRDS('output/SPW.Rds')
ref_group <- 'Control'
test_group <- 'SP'

# 验证集里没有的基因可以去genecards上搜索别名进行替换
gene_show
gene_show <- c("TMEM80","PLK1","C1QA", "S100A11","FCGRT" ,"MAP3K8","PHACTR1","CTNNBIP1","CFDP1","CP","RALGPS1","COL16A1","DTL", "NFKBIA"   
              )

# 单因素Logistic回归 -----------------------------------------------------------

# 提取基因
dat_expr_show <- dat_expr[gene_show, ] %>% 
  t() %>% 
  as.data.frame()

# 因子化
dat_expr_show$group <- dat_group$group %>% 
  factor(
    levels = c(ref_group, test_group)
  )

# 单因素Logistic回归
# 遍历每一个基因
for(i in 1:length(gene_show)){
  # 初始化单因素Logistic回归结果
  if(i == 1) dat_log_uni <- data.frame()
  
  # Logistic回归
  res_fit_log_uni <- glm(
    as.formula(
      paste0('group ~ `', gene_show[i], '`')
    ),
    data = dat_expr_show,
    # 结局为二分类
    family = 'binomial'
  )
  
  # 计算系数
  res_sum_log_uni <- summary(res_fit_log_uni)
  dat_log_uni0 <- res_sum_log_uni$coefficients
  # 提取结果
  dat_log_uni0 <- cbind(
    dat_log_uni0,
    OR = exp(dat_log_uni0[, 'Estimate']),
    exp(confint(res_fit_log_uni))
  ) %>% 
    .[-1, ]
  
  # 整理结果
  dat_log_uni0 <- data.frame(
    gene = gene_show[i],
    OR = dat_log_uni0['OR'],
    OR1 = dat_log_uni0['2.5 %'],
    OR2 = dat_log_uni0['97.5 %'],
    pvalue = dat_log_uni0['Pr(>|z|)'],
    row.names = gene_show[i]
  )
  
  # 返回结果
  dat_log_uni <- rbind(dat_log_uni, dat_log_uni0)
  
  # 移除中间变量
  rm(i, res_fit_log_uni, res_sum_log_uni, dat_log_uni0)
}

# 输出单因素Logistic结果
dat_log_uni_output <- data.frame(
  Gene = dat_log_uni$gene,
  `Univariate OR (95% CI)` = paste0(
    sprintf('%0.2f', dat_log_uni$OR),
    ' (',
    sprintf('%0.2f', dat_log_uni$OR1),
    '-',
    sprintf('%0.2f', dat_log_uni$OR2),
    ')'
  ),
  `Univariate P Value` = case_when(
    dat_log_uni$pvalue < 0.001 ~ '*** < 0.001',
    dat_log_uni$pvalue < 0.01 ~ '** < 0.01',
    dat_log_uni$pvalue < 0.05 ~ '* < 0.05',
    dat_log_uni$pvalue >= 0.05 ~ '> 0.05'
  ),
  check.names = F
)
write.csv(dat_log_uni_output, file = 'output/Logistic_uni_res_GSE18995.csv', row.names = F)


# 多因素Logistic回归 -----------------------------------------------------------

# 单因素进入多因素的p值阈值
pvalue_cutoff <- 0.05
# 提取单因素显著的基因
gene_show_log <- dat_log_uni %>% 
  dplyr::filter(pvalue < pvalue_cutoff) %>% 
  .$gene

# VIF阈值
vif_cutoff <- 5
# VIF为方差膨胀因子Variance Inflation Factor
# VIF可以评估各基因之间的多重共线性
# 一般会以5或10为VIF的阈值，即若所有基因的VIF均小于5或10，则认为不存在多重共线性

# VIF_idc为是否要重新做多因素Logistic回归的指示变量
# 指示变量的初始值为T，强制进行第一次Logistic回归
vif_idc <- T
# 循环执行Logsitic回归，直到所有基因的VIF均小于阈值
while(vif_idc){
  # Logistic 回归
  res_fit_log_multi <- glm(
    as.formula(
      paste0('group ~ ', paste0(gene_show_log, collapse = ' + '))
    ),
    data = dat_expr_show,
    # 结局为二分类
    family = 'binomial'
  )
  # 计算VIF
  vif_log_show <- vif(res_fit_log_multi)
  # 检查是否所有基因的VIF均小于阈值
  vif_idc <- !all(vif_log_show < vif_cutoff)
  # 如果所有基因的VIF都小于阈值，终止循环
  if(!vif_idc) break
  # 如果存在至少一个基因的VIF值不小于阈值，则将VIF最大的基因去除，重新进入循环
  gene_show_log <- gene_show_log[- which.max(vif_log_show)]
}
# 检查所有基因的VIF值，此时所有剩余的基因即为单因素进入多因素的基因
vif_log_show

# 整理结果
res_sum_log_multi <- summary(res_fit_log_multi)
dat_log_multi <- data.frame(
  gene = gene_show_log,
  OR = exp(res_sum_log_multi$coefficients[-1, 'Estimate']),
  OR1 = exp(confint(res_fit_log_multi)[-1, '2.5 %']),
  OR2 = exp(confint(res_fit_log_multi)[-1, '97.5 %']),
  pvalue = res_sum_log_multi$coefficients[-1, 'Pr(>|z|)']
)

# 输出多因素Logistic结果
dat_log_multi_output <- data.frame(
  Gene = dat_log_multi$gene,
  `Multivariate OR (95% CI)` = paste0(
    sprintf('%0.2f', dat_log_multi$OR),
    ' (',
    sprintf('%0.2f', dat_log_multi$OR1),
    '-',
    sprintf('%0.2f', dat_log_multi$OR2),
    ')'
  ),
  `Multivariate P Value` = case_when(
    dat_log_multi$pvalue < 0.001 ~ '*** < 0.001',
    dat_log_multi$pvalue < 0.01 ~ '** < 0.01',
    dat_log_multi$pvalue < 0.05 ~ '* < 0.05',
    dat_log_multi$pvalue >= 0.05 ~ '> 0.05'
  ),
  check.names = F
)
write.csv(dat_log_multi_output, file = 'output/Logistic_multi_res_GSE18995.csv', row.names = F)

# Logistic回归系数
coef_log <- res_fit_log_multi$coef

# 提取显著基因的表达谱
dat_expr_log <- dat_expr_show[, gene_show_log]

# 预测值
linear_pre <- (as.matrix(dat_expr_log) %*% coef_log[-1]) %>% 
  as.data.frame()
colnames(linear_pre) <- 'linear_pre'
# 每个样本的预测值的计算方式为各样本中各基因的表达量乘以各基因的系数

# 合并单因素和多因素Logistic结果
dat_log_uni_multi_output <- full_join(
  dat_log_uni_output,
  dat_log_multi_output,
  by = 'Gene'
)
write.csv(dat_log_uni_multi_output, file = 'output/Logistic_res_GSE18995.csv', row.names = F, na = '')


# 森林图 ---------------------------------------------------------------------

# 森林图数据
dat_log_uni_forest <- rbind(
  colnames(dat_log_uni_output),
  dat_log_uni_output
)

# 森林图
pdf(file = paste0('output/forest_GSE18995.pdf'), height = 6, width = 8, onefile = F)
forestplot(
  # 森林图文字标签
  labeltext = as.matrix(dat_log_uni_forest),
  # OR值
  mean = c(NA, dat_log_uni$OR),
  # 置信区间
  lower = c(NA, dat_log_uni$OR1),
  upper = c(NA, dat_log_uni$OR2),
  # 森林图端点
  ci.vertices = T,
  # 图在表中列的位置
  graph.pos = 3,
  # 点型，可选参数有fpDrawCircleCI, fpDrawNormalCI, fpDrawDiamondCI, fpDrawPointCI
  fn.ci_norm = fpDrawCircleCI,
  # 点和线的颜色
  col = fpColors(
    box = ggsci::pal_npg()(2)[2],
    lines = 'black',
    zero = 'black'
  ),
  # 超出以下范围，端点用箭头表示
  clip = c(0, 5),
  # 点的大小
  boxsize = 0.15,
  # 参考线位置
  zero = 1,
  # 横坐标标签
  xlab = 'Odds Ratio',
  # 横坐标刻度
  # xticks = seq(0.5, 4, 3.5),
  # 三线图中线的行数
  hrzl_line = rlang::set_names(
    list(
      gpar(col = 'black'),
      gpar(col = 'black'),
      gpar(col = 'black')
    ),
    # 线分别在第1行、第2行、(特征数+2)行
    c(1, 2, nrow(dat_log_uni) + 2)
  )
)
dev.off()
# 多因素森林图同理


# 数据保存 --------------------------------------------------------------------

saveRDS(dat_log_uni, file = 'output/Logistic_uni_res_GSE18995.Rds')
saveRDS(dat_log_uni_output, file = 'output/Logistic_uni_output_GSE18995.Rds')
saveRDS(dat_log_multi, file = 'output/Logistic_multi_res_GSE18995.Rds')
saveRDS(dat_log_multi_output, file = 'output/Logistic_multi_output_GSE18995.Rds')
saveRDS(gene_show_log, file = 'output/gene_log_GSE18995.Rds')
saveRDS(coef_log, file = 'output/Logistic_coef_GSE18995.Rds')
saveRDS(linear_pre, file = 'output/Logistic_linear_pre_GSE18995.Rds')

# 仙桃
Logistic_input <- dat_expr_show
Logistic_input <- Logistic_input %>%
  mutate(outcome = if_else(group == "pre", 0, 1)) %>%
  select(outcome, everything(), -group)
write.csv(Logistic_input, file = 'output/二分类logistic_GSE18995.csv', row.names = FALSE, na = "")


a <- Logistic_uni_output_GSE18995[order(Logistic_uni_output_GSE18995$`Univariate P Value`),]
b <- a[1:9,]
saveRDS(b$Gene, file = 'output/Logistic_gene.Rds')
