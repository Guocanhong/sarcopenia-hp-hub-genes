rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ggplot2)
library(enrichplot)
gene_show <- readRDS('output/Key_gene.Rds')

# GOKEGG分析 ----------------------------------------------------------------

# 基因名转化
dat_gene_entrez <- bitr(
  gene_show,
  fromType = 'SYMBOL',
  toType = 'ENTREZID',
  OrgDb = org.Hs.eg.db
)

# GO分析
res_go <- enrichGO(
  gene = dat_gene_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  # 富集类目，可选参数有'BP','CC','MF','ALL'
  ont = 'ALL',
  # p值校正方法，可选参数有'holm','hochberg','hommel','bonferroni','BH','BY','fdr','none'
  pAdjustMethod = 'BH',
  # p值阈值
  pvalueCutoff = 1,
  # q值阈值
  qvalueCutoff = 0.25,
  # 基因ID转成基因名
  readable = T
)

# KEGG分析
res_kegg <- enrichKEGG(
  gene = dat_gene_entrez$ENTREZID,
  organism = 'hsa',
  # p值校正方法，可选参数有'holm','hochberg','hommel','bonferroni','BH','BY','fdr','none'
  pAdjustMethod = 'BH',
  # p值阈值
  pvalueCutoff = 1,
  # q值阈值
  qvalueCutoff = 0.25
) %>%
  # 基因ID转成基因名
  setReadable(
    OrgDb = org.Hs.eg.db,
    keyType = 'ENTREZID'
  )

# save(res_go, res_kegg, file = 'Rda/res_gokegg.Rda')

# 提取GO分析结果
dat_res_go <- as.data.frame(res_go)
# BP、CC、MF各取top10条目（按p值）
dat_res_go_top10 <- rbind(
  dat_res_go %>% 
    dplyr::filter(ONTOLOGY == 'BP') %>% 
    .[1:10, ],
  dat_res_go %>% 
    dplyr::filter(ONTOLOGY == 'CC') %>% 
    .[1:10, ],
  dat_res_go %>% 
    dplyr::filter(ONTOLOGY == 'MF') %>% 
    .[1:10, ]
) %>% 
  na.omit()
# 提取KEGG分析结果
dat_res_kegg <- as.data.frame(res_kegg)
# KEGG取top10通路（按p值）
dat_res_kegg_top30 <- dat_res_kegg[1:30, ] %>% 
  na.omit()


# GO柱状图 -------------------------------------------------------------------

p_go_barplot <- barplot(
  res_go,
  # 横坐标，可选参数有'GeneRatio','Count'等
  x = 'GeneRatio',
  showCategory = 10,
  color = 'pvalue',
  split = 'ONTOLOGY',
  label_format = 40
) +
  facet_grid(
    ONTOLOGY ~ .,
    scales = 'free'
  ) +
  scale_fill_gradient(
    high = '#4DBBD5',
    low = '#E64B35',
    limits = c(0, 0.05)
  )
ggsave(file = 'output/Fig6a.pdf', p_go_barplot, width = 8, height = 12)


# KEGG气泡图 -----------------------------------------------------------------

p_kegg_dotplot <- dotplot(
  res_kegg,
  # 横坐标，可选参数有'GeneRatio','Count'等
  x = 'GeneRatio',
  showCategory = 30,
  color = 'pvalue',
  label_format = 40
) +
  scale_fill_gradient(
    high = '#4DBBD5',
    low = '#E64B35',
    limits = c(0, 0.05)
  )
ggsave(file = 'output/Fig26B.pdf', p_kegg_dotplot, width = 8, height = 12)
