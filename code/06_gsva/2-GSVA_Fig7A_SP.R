# 方案A：SP vs Control 全组 GSVA（与 7-GSVA1 同款设计，替换原高低分组设计）
# 基因集：MSigDB C2 all（c2.all.v2025.1）；筛选：BH adj.P.Val<0.05，上下调各取 |t| top10
# 输出：output/Fig7A.pdf + output/GSVA_res_SP.csv（旧图 Fig7D.pdf 保留不删）
rm(list = ls())
library(tidyverse)
library(GSVA)
library(GSEABase)
library(limma)
library(ggplot2)

dat_expr  <- readRDS('input/dat_expr.Rds')
dat_group <- readRDS('input/dat_group.Rds')
ref_group  <- 'Control'
test_group <- 'SP'

# GSVA 分析 --------------------------------------------------------------------
geneset <- getGmt('input/c2.all.v2025.1.Hs.symbols.gmt')
dat_res_gsva <- gsvaParam(
  as.matrix(dat_expr),
  geneset,
  maxDiff = TRUE,
  kcdf = 'Gaussian'
) %>% gsva() %>% as.data.frame()

# GSVA 差异分析（SP - Control）---------------------------------------------------
mat_design <- model.matrix(~ 0 + factor(dat_group$group))
colnames(mat_design) <- levels(factor(dat_group$group))
rownames(mat_design) <- dat_group$sample

fit_limma <- lmFit(dat_res_gsva, mat_design)
contrast  <- paste0(test_group, '-', ref_group)
mat_contrast <- makeContrasts(contrasts = contrast, levels = mat_design)
fit_limma <- contrasts.fit(fit_limma, mat_contrast) %>% eBayes()

dat_res_gsva_diff <- topTable(fit_limma, coef = contrast, n = Inf) %>% na.omit()
dat_res_gsva_diff$P.Value[dat_res_gsva_diff$P.Value == 0] <-
  min(dat_res_gsva_diff$P.Value[dat_res_gsva_diff$P.Value != 0])
dat_res_gsva_diff$adj.P.Val[dat_res_gsva_diff$adj.P.Val == 0] <-
  min(dat_res_gsva_diff$adj.P.Val[dat_res_gsva_diff$adj.P.Val != 0])

write.csv(dat_res_gsva_diff, 'output/GSVA_res_SP.csv', row.names = TRUE)

dat_sig <- dat_res_gsva_diff %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::mutate(t_group = ifelse(t > 0, 'Up', 'Down'))
cat('SP: adj.P.Val<0.05 通路数 =', nrow(dat_sig),
    '(Up', sum(dat_sig$t_group == 'Up'), '/ Down', sum(dat_sig$t_group == 'Down'), ')\n')

dat_sig <- lapply(c('Up', 'Down'), function(x) {
  dat_sig %>% dplyr::filter(t_group == x) %>%
    dplyr::arrange(desc(abs(t))) %>% .[1:10, ]
}) %>% do.call(rbind, .) %>% na.omit()

# 双向柱状图 --------------------------------------------------------------------
p <- ggplot(dat_sig, aes(t, factor(rownames(dat_sig), levels = rev(rownames(dat_sig))), fill = t_group)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c(Down = '#4DBBD5', Up = '#E64B35')) +
  labs(title = paste0('GSVA differential pathways (', test_group, ' vs ', ref_group, ', n=',
                      sum(dat_group$group == test_group), '/', sum(dat_group$group == ref_group), ')'),
       subtitle = 'MSigDB C2, BH adj.P.Val < 0.05, top10 by |t|',
       x = 't statistic', y = NULL, fill = '') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8))
ggsave(file = 'output/Fig7A.pdf', p, width = 8, height = 8)
cat('完成: output/Fig7A.pdf / output/GSVA_res_SP.csv\n')
