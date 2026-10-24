rm(list = ls())
library(tidyverse)
library(WGCNA)
dat_expr <- readRDS('input/dat_expr.Rds')
dat_group <- readRDS('input/dat_group.Rds')

ref_group <- 'Control'
test_group <- 'SP'


# 数据处理 --------------------------------------------------------------------

# 计算各基因的方差（或绝对中位差）
vars_all_gene <- apply(dat_expr, 1, mad)
# 如果以绝对中位差筛选基因，将var改成mad

# 选择方差（或绝对中位差）前n%的基因
n_vars <- 0.7
gene0 <- names(
  vars_all_gene[vars_all_gene >= quantile(vars_all_gene, probs = n_vars)]
)
# 具体提取前百分之多少，以提取后基因数量在5000个左右为宜
dat_expr <- dat_expr[gene0, ] %>% 
  t %>% 
  as.data.frame()

# 检查数据中是否存在缺失条目、权重低于阈值的条目和零方差的基因
res_gsg <- goodSamplesGenes(
  datExpr = dat_expr,
  verbose = 3
)

# 去除不好的基因和样本
if(res_gsg$allOK){
  print('all OK')
}else{
  if(sum(!res_gsg$goodGenes) == 0){
    print('genes OK')
  }else{
    print(
      paste0(
        'Removing gene(s): ',
        paste(colnames(dat_expr)[!res_gsg$goodGenes], collapse = ',')
      )
    )
  }
  if(sum(!res_gsg$goodSamples) == 0){
    print('Samples OK')
  }else{
    print(
      paste0(
        'Removing sample(s): ',
        paste(rownames(dat_expr)[!res_gsg$goodSamples], collapse = ',')
      )
    )
  }
  gene_good <- colnames(dat_expr)[res_gsg$goodGenes]
  sample_good <- rownames(dat_expr)[res_gsg$goodSamples]
  dat_expr <- dat_expr[sample_good, gene_good]
  print('Bad gene(s) and sample(s) have been removed')
}


# 样本聚类 --------------------------------------------------------------------

# 对样本聚类
clust_sample_tree <- hclust(
  dist(dat_expr),
  method = 'average'
)

# 样本聚类树
pdf(file = 'output/Fig3A.pdf', width = 7, height = 6)
plot(
  clust_sample_tree,
  xlab = NA,
  sub = NA,
  main = 'Sample Clustering to Detect Outliers'
)
dev.off()

# 因子化
dat_group$group <- factor(
  dat_group$group,
  levels = c(ref_group, test_group)
)

# 构建分组信息
dat_trait <- model.matrix(~ 0 + dat_group$group) %>% 
  as.data.frame()
colnames(dat_trait) <- levels(dat_group$group)
rownames(dat_trait) <- dat_group$sample

# 分组颜色
color_trait <- numbers2colors(
  dat_trait,
  color = c('grey90', '#E64B35')
)

# 样本聚类树
plotDendroAndColors(
  clust_sample_tree,
  color_trait,
  groupLabels = colnames(dat_trait),
  main = 'Sample Dendrogram and Trait Heatmap'
)
dev.off()


# 软阈值散点图 ------------------------------------------------------------------

# 允许多线程计算
enableWGCNAThreads()
# Mac电脑使用以下代码令计算机允许多线程
# allowWGCNAThreads()

# 无标度拓扑拟合指数的软阈值
power_vec <- 1:12

# 软阈值的无标度拓扑分析
r2_cutoff <- 0.95
res_sft <- pickSoftThreshold(
  dat_expr,
  # 软阈值
  powerVector = power_vec,
  # 无标度拟合指数
  RsquaredCut = r2_cutoff,
  verbose = 5
)

# save(res_sft, file = 'Rda/res_sft.Rda')

# 软阈值散点图
pdf(file = 'output/Fig3B.pdf', width = 8, height = 6)
# 散点图两个部分
par(mfrow = 1:2)
# 拟合指数与power值散点图
plot(
  # 软阈值
  res_sft$fitIndices$Power,
  # 无缩放拓扑模型拟合
  - sign(res_sft$fitIndices$slope) * res_sft$fitIndices$SFT.R.sq,
  xlab = 'Soft Threshold (power)',
  ylab = 'Scale Free Topology Model Fit, signed R^2',
  type = 'n',
  main = 'Scale Independence'
)
text(
  # 软阈值
  res_sft$fitIndices$Power,
  # 无缩放拓扑模型拟合
  - sign(res_sft$fitIndices$slope) * res_sft$fitIndices$SFT.R.sq,
  labels = power_vec,
  col = '#E64B35'
)
abline(
  # 无标度拟合指数
  h = r2_cutoff,
  col = '#E64B35'
)
# 平均连通性与power值散点图
plot(
  # 软阈值
  res_sft$fitIndices$Power,
  # 节点连接程度
  res_sft$fitIndices$mean.k.,
  xlab = 'Soft Threshold (power)',
  ylab = 'Mean Connectivity',
  type = 'n',
  main = 'Mean Connectivity'
)
text(
  # 软阈值
  res_sft$fitIndices$Power,
  # 节点连接程度
  res_sft$fitIndices$mean.k.,
  labels = power_vec,
  col = '#E64B35'
)
dev.off()


# 基因聚类 --------------------------------------------------------------------

# 最佳软阈值
soft_power <- res_sft$powerEstimate
soft_power
# 这里需要仔细看最佳软阈值是否合理，如果合理，下面这句是不用跑的
# soft_power <- sum(- sign(res_sft$fitIndices$slope) * res_sft$fitIndices$SFT.R.sq < r2_cutoff) + 1
# soft_power

# 根据最佳软阈值计算网络邻接性
mat_adjacency <- adjacency(
  dat_expr,
  power = soft_power
)

# 计算TOM矩阵（拓扑重叠矩阵，Topological Overlap Matrix）
mat_TOM <- TOMsimilarity(mat_adjacency)
# save(mat_TOM, file = 'Rda/mat_TOM.Rda')
# 计算TOM差异矩阵
mat_dissTOM <- 1 - mat_TOM
# TOM表示基因之间相似性，dissTOM表示基因之间的距离
# TOM越接近0 = 相似性越低 = dissTOM越接近1 = 距离越远

# 关闭多线程计算
disableWGCNAThreads()

# 对基因聚类
clust_gene_tree <- hclust(
  as.dist(mat_dissTOM),
  method = 'average'
)

# 基因聚类树
plot(
  clust_gene_tree,
  xlab = NA,
  sub = NA,
  main = 'Gene Clustering on TOM-Based Dissimilarity',
  labels = F
)
dev.off()

# 模块基因最小数量
min_cluster_size <- 25#25,50,100

# 动态剪切模块识别
dynamic_mod <- cutreeDynamic(
  dendro = clust_gene_tree,
  distM = mat_dissTOM,
  # 树的分裂程度
  deepSplit = 2,
  # 将PAM算法应用于树型结构
  pamRespectsDendro = F,
  minClusterSize = min_cluster_size
)
table(dynamic_mod)

# 模块颜色
color_dynamic <- labels2colors(dynamic_mod)
table(color_dynamic)

# 基因聚类树
plotDendroAndColors(
  clust_gene_tree,
  color_dynamic,
  groupLabels = 'Dynamic Tree Cut',
  dendroLabels = F,
  # 添加引导线
  addGuide = T,
  # 引导线高度比例
  guideHang = 0.05,
  main = 'Gene Dendrogram and Module Colors'
)
dev.off()


# 相似模块聚类和合并 ---------------------------------------------------------------

# 计算模块的特征基因
list_me <- moduleEigengenes(
  dat_expr,
  colors = color_dynamic
)

# 提取特征基因矩阵
dat_merged_me <- list_me$eigengenes

# 计算特征基因矩阵的相关性矩阵，即特征基因之间的差异
mat_dissme <- 1 - WGCNA::cor(dat_merged_me)

# 对特征基因聚类
clust_me_tree <- hclust(
  as.dist(mat_dissme),
  method = 'average'
)

# 合并相似模块的阈值
mediss_threshold <- 0.25

# 特征基因聚类树
plot(
  clust_me_tree,
  xlab = NA,
  sub = NA,
  main = 'Gene Dendrogram and Module Colors'
)
abline(
  h = mediss_threshold,
  col = '#E64B35'
)
dev.off()

# 合并相似模块
res_merge <- mergeCloseModules(
  dat_expr,
  color_dynamic,
  cutHeight = mediss_threshold,
  verbose = 3
)

# 新模块颜色
color_merged <- res_merge$colors

dat_merged_me_new <- res_merge$newMEs

# 动态模块颜色的基因聚类树
pdf(file = 'output/Fig3C.pdf', width = 6, height = 6)
plotDendroAndColors(
  # 新模块基因
  clust_gene_tree,
  data.frame(
    dynamic = color_dynamic,
    merged = color_merged
  ),
  groupLabels = c('Dynamic Tree Cut', 'Merged Dynamic'),
  dendroLabels = F,
  # 添加引导线
  addGuide = T,
  # 引导线高度比例
  guideHang = 0.05,
  main = 'Gene Dendrogram and Module Colors'
)
dev.off()

# 模块颜色
color_module <- color_merged
table(color_module)


# 模块基因热图 ------------------------------------------------------------------

# 计算相关性和p值
mat_module_trait_cor <- WGCNA::cor(
  dat_merged_me_new,
  dat_trait,
  # 统计方法
  method = 'spearman',
  use = 'p'
)
mat_module_trait_p <- corPvalueStudent(
  mat_module_trait_cor,
  nrow(dat_expr)
)

# 相关性和p值标注
mat_text <- paste0(
  round(mat_module_trait_cor, 2),
  '\n(',
  signif(mat_module_trait_p, 3),
  ')'
)
dim(mat_text) <- dim(mat_module_trait_cor)

# 热图
pdf(file = 'output/Fig3D.pdf', width = 8, height = 8)
labeledHeatmap(
  Matrix = mat_module_trait_cor,
  xLabels = colnames(dat_trait),
  yLabels = colnames(dat_merged_me_new),
  ySymbols = colnames(dat_merged_me_new),
  colors = colorRampPalette(colors = c('#4DBBD5', 'white', '#E64B35'))(50),
  textMatrix = mat_text,
  zlim = c(-1, 1),
  main = 'Module-Trait Relationships'
)
dev.off()


# 模块相关性热图 -----------------------------------------------------------------

# 计算相关性和p值
mat_module_module_cor <- WGCNA::cor(
  dat_merged_me_new,
  dat_merged_me_new,
  # 统计方法
  method = 'spearman',
  use = 'p'
)
mat_module_module_p <- corPvalueStudent(
  mat_module_module_cor,
  nrow(dat_expr)
)

# 相关性和p值标注
mat_text_module <- paste0(
  round(mat_module_module_cor, 2),
  '\n(',
  signif(mat_module_module_p, 3),
  ')'
)
dim(mat_text_module) <- dim(mat_module_module_cor)

# 热图
pdf(file = 'output/Fig3E.pdf', width = 8, height = 8)
labeledHeatmap(
  Matrix = mat_module_module_cor,
  xLabels = colnames(dat_merged_me_new),
  xSymbols = colnames(dat_merged_me_new),
  yLabels = colnames(dat_merged_me_new),
  ySymbols = colnames(dat_merged_me_new),
  colors = colorRampPalette(colors = c('#4DBBD5', 'white', '#E64B35'))(50),
  textMatrix = mat_text_module,
  zlim = c(-1, 1),
  main = 'Module-Module Relationships'
)
dev.off()


# 模块基因 --------------------------------------------------------------------

# 挑选模块，注意不要挑选到'grey'模块
# 注：'grey60','lightgrey'等看着好像是'grey'模块的，但是都不是'grey'模块，都可以选
# 可以挑选多个模块基因，最后合并


select_module <- c("black","pink","green","brown","red")
# 计算模块与基因的相关性和p值
mat_gene_module_cor <- WGCNA::cor(
  dat_expr,
  dat_merged_me_new,
  # 统计方法
  method = 'spearman',
  use = 'p'
)
mat_gene_module_p <- corPvalueStudent(
  mat_gene_module_cor,
  nrow(dat_expr)
)

# 计算分组与基因的相关性和p值
mat_gene_trait_cor <- WGCNA::cor(
  dat_expr,
  dat_trait,
  # 统计方法
  method = 'spearman',
  use = 'p'
)
mat_gene_trait_p <- corPvalueStudent(
  mat_gene_trait_cor,
  nrow(dat_expr)
)

# 选择模块基因
gene_module <- gene0[color_module %in% select_module]


# gene_module <- data.frame(gene_module)
# table(gene_module$gene_module)

# 散点图
# pdf(file = 'output/Fig3F.pdf', width = 8, height = 8)
# verboseScatterplot(
#   abs(mat_gene_module_cor[gene_module, paste0('ME', select_module)]),
#   abs(mat_gene_trait_cor[gene_module, 1]),
#   xlab = paste0('Module Membership in ',select_module,' Module'),
#   ylab = 'Gene Significance for Trait',
#   col = select_module,
#   main = 'Module membership vs Gene Significance\n'
# )
# dev.off()


# 保存数据 --------------------------------------------------------------------

saveRDS(gene_module, file = 'output/ModuleGenes.Rds')
