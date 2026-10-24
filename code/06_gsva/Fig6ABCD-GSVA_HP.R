rm(list = ls())
library(tidyverse)
library(GSVA)
library(limma)
library(ggplot2)
dat_expr <- readRDS('input/dat1_expr.Rds')
dat_group <- readRDS('input/dat1_group.Rds')
gene_show <- readRDS('input/Key_gene.Rds')
ref_group <- 'Control'
test_group <- 'HP'



# GSVA分析 ------------------------------------------------------------------

# 整理单基因高低表达数据
#dat_group <- dat_group[dat_group$group == test_group, ]
#dat_expr <- dat_expr[, dat_group$sample]
#dat_group$group <- case_when(
#  as.numeric(dat_expr[s_gene, ]) >= median(as.numeric(dat_expr[s_gene, ])) ~ 'high',
#   as.numeric(dat_expr[s_gene, ]) < median(as.numeric(dat_expr[s_gene, ])) ~ 'low'
# )
#ref_group <- 'low'
#test_group <- 'high'

# 基因集
# GSVA一般用HALLMARK基因集
geneset_h <- GSEABase::getGmt('input/c2.all.v2025.1.Hs.symbols.gmt')

# GSVA分析
dat_res_gsva <- gsvaParam(
  as.matrix(dat_expr),
  geneset_h,
  # 富集得分计算方式为最大正负随机游走偏差之间的幅度差
  maxDiff = T,
  # 当输入的表达谱为Count(整数)时选择'Poisson'，其余情况均用'Gaussian'
  kcdf = 'Gaussian'
) %>% 
  gsva() %>% 
  as.data.frame()

# save(dat_res_gsva, file = 'Rda/dat_res_gsva.Rda')


# GSVA差异分析 ----------------------------------------------------------------

# 创建模型矩阵
mat_design <- model.matrix(~ 0 + factor(dat_group$group))
# 重新命名行名和列名
colnames(mat_design) <- levels(factor(dat_group$group))
rownames(mat_design) <- dat_group$sample

# 使用线性模型拟合表达数据
fit_limma <- lmFit(dat_res_gsva, mat_design)

# 为组间差异创建对照
contrast <- paste0(test_group, '-', ref_group)
# 一定要实验组在前，对照组在后

# 创建线性对比矩阵
mat_contrast <- makeContrasts(contrasts = contrast, levels = mat_design)

# 对拟合的线性模型进行再次拟合
fit_limma <- contrasts.fit(fit_limma, mat_contrast)

# Bayes调整
fit_limma <- eBayes(fit_limma)

# 提取差异分析结果
dat_res_gsva_diff <- topTable(
  fit_limma,
  coef = contrast,
  n = Inf
) %>% 
  na.omit()

# 将为0的p值改成除0外最小的p值，以防止p值为0取不到对数
dat_res_gsva_diff$P.Value[dat_res_gsva_diff$P.Value == 0] <- min(dat_res_gsva_diff$P.Value[dat_res_gsva_diff$P.Value != 0])
dat_res_gsva_diff$adj.P.Val[dat_res_gsva_diff$adj.P.Val == 0] <- min(dat_res_gsva_diff$adj.P.Val[dat_res_gsva_diff$adj.P.Val != 0])

# 筛选差异显著的通路
dat_res_gsva_diff_sig <- dat_res_gsva_diff %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::arrange(desc(t))
dat_res_gsva_diff_sig$t_group <- case_when(
  dat_res_gsva_diff_sig$t < 0 ~ 'Down',
  dat_res_gsva_diff_sig$t > 0 ~ 'Up'
)

# 如果通过过多，不方便可视化，可筛选上下调通路各top10
 dat_res_gsva_diff_sig <- lapply(
   c('Up', 'Down'),
   function(x){
     dat_res_gsva_diff_sig %>% 
      dplyr::filter(t_group == x) %>% 
       dplyr::arrange(desc(abs(t))) %>% 
       .[1:10, ] %>% 
       na.omit()
   }
 ) %>%
   do.call(rbind, .)

# 差异显著通路的富集矩阵
dat_res_gsva_sig <- dat_res_gsva[rownames(dat_res_gsva_diff_sig), ]


# GSVA双向柱状图 ---------------------------------------------------------------

# 因子化
dat_res_gsva_diff_sig$geneset <- factor(
  rownames(dat_res_gsva_diff_sig),
  levels = rownames(dat_res_gsva_diff_sig)
)

# 双向柱状图
p_gsva_barplot <- ggplot(
  dat_res_gsva_diff_sig,
  aes(
    t,
    geneset,
    fill = t_group
  )
) +
  geom_bar(
    stat = 'identity'
  ) +
  theme_bw() +
  guides(fill = 'none') +
  labs(
    y = element_blank()
  ) +
  scale_fill_manual(
    values = c('#4DBBD5', '#E64B35')
  )
ggsave(file = 'output/Fig7D11.pdf', p_gsva_barplot, width = 18, height = 8)

# Fig6BCD其他基因高低分组的GSVA双向柱状图绘制方法与Fig6A相同，只需要将基因替换成对应基因即可
