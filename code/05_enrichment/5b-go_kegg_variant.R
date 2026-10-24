# 基于 SP∩HP 交集基因（118）重做 GO+KEGG 富集
# 口径：保留全量结果，作图前按未校正 p<0.05 筛选；p.adjust 在结果表中完整保留
rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(enrichplot)

# 读入交集基因 ----------------------------------------------------------------
venn <- read.csv('../2-Venn/output/SP-HP.CSV', row.names = 1)
gene_show <- as.character(venn$x)
cat('交集基因数:', length(gene_show), '\n')

# 基因名转化 ------------------------------------------------------------------
dat_gene_entrez <- bitr(
  gene_show,
  fromType = 'SYMBOL',
  toType = 'ENTREZID',
  OrgDb = org.Hs.eg.db
)
cat('成功映射 Entrez ID:', nrow(dat_gene_entrez), '/', length(gene_show), '\n')

# GO / KEGG 全量分析 ------------------------------------------------------------
res_go <- enrichGO(
  gene = dat_gene_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = 'ALL',
  pAdjustMethod = 'BH',
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)
res_kegg <- enrichKEGG(
  gene = dat_gene_entrez$ENTREZID,
  organism = 'hsa',
  pAdjustMethod = 'BH',
  pvalueCutoff = 1,
  qvalueCutoff = 1
) %>%
  setReadable(OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')

dat_res_go   <- as.data.frame(res_go)   %>% dplyr::filter(pvalue < 0.05) %>% arrange(pvalue)
dat_res_kegg <- as.data.frame(res_kegg) %>% dplyr::filter(pvalue < 0.05) %>% arrange(pvalue)
cat('GO  p<0.05 条目数:', nrow(dat_res_go), '\n')
cat('KEGG p<0.05 条目数:', nrow(dat_res_kegg), '\n')

# 结果表落盘（论文补充材料用，含 BH 校正 p 值）----------------------------------
write.csv(dat_res_go,   'output/GO_res_118genes.csv',   row.names = FALSE)
write.csv(dat_res_kegg, 'output/KEGG_res_118genes.csv', row.names = FALSE)

# GO 柱状图（BP/CC/MF 各 top10）-------------------------------------------------
if (nrow(dat_res_go) > 0) {
  go_top10 <- rbind(
    dat_res_go %>% dplyr::filter(ONTOLOGY == 'BP') %>% slice_head(n = 10),
    dat_res_go %>% dplyr::filter(ONTOLOGY == 'CC') %>% slice_head(n = 10),
    dat_res_go %>% dplyr::filter(ONTOLOGY == 'MF') %>% slice_head(n = 10)
  ) %>% na.omit()

  p_go_barplot <- ggplot(go_top10, aes(x = GeneRatio, y = reorder(Description, -pvalue))) +
    geom_col(aes(fill = pvalue), width = 0.7) +
    scale_fill_gradient(high = '#4DBBD5', low = '#E64B35') +
    facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
    labs(title = 'GO enrichment of SP-intersect-HP genes (n = 118)',
         x = 'GeneRatio', y = NULL, fill = 'p value') +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9))
  ggsave(file = 'output/Fig6A1.pdf', p_go_barplot, width = 8, height = 10)
} else {
  cat('GO 无 p<0.05 条目，未出图\n')
}

# KEGG 气泡图（p<0.05 中取 top30）-----------------------------------------------
if (nrow(dat_res_kegg) > 0) {
  kegg_top <- dat_res_kegg %>% slice_head(n = 30)
  p_kegg_dotplot <- ggplot(kegg_top,
                           aes(x = GeneRatio, y = reorder(Description, -pvalue))) +
    geom_point(aes(size = Count, color = pvalue)) +
    scale_color_gradient(high = '#4DBBD5', low = '#E64B35') +
    labs(title = 'KEGG enrichment of SP-intersect-HP genes (n = 118)',
         x = 'GeneRatio', y = NULL, color = 'p value', size = 'Count') +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9))
  ggsave(file = 'output/Fig6B1.pdf', p_kegg_dotplot, width = 8, height = 10)
} else {
  cat('KEGG 无 p<0.05 条目，未出图\n')
}

cat('完成\n')
