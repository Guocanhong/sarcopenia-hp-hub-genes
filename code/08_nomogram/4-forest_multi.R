# 多因素（2基因）Logistic 森林图：NFKBIA + CTNNBIP1
# 输出：output/forest_SP_multi.pdf（与单因素 forest_SP.pdf 配对）
rm(list = ls())
library(forestplot)
library(ggsci)

dat_expr  <- readRDS('input/dat_expr.Rds')
dat_group <- readRDS('input/dat_group.Rds')
genes <- c('NFKBIA', 'CTNNBIP1')

dat <- as.data.frame(t(dat_expr[genes, ]))
dat$group <- factor(dat_group$group, levels = c('Control', 'SP'))

fit <- glm(group ~ NFKBIA + CTNNBIP1, data = dat, family = 'binomial')
s <- summary(fit)
ci <- exp(confint(fit))
dat_multi <- data.frame(
  gene   = genes,
  OR     = exp(s$coefficients[-1, 'Estimate']),
  OR1    = ci[-1, '2.5 %'],
  OR2    = ci[-1, '97.5 %'],
  pvalue = s$coefficients[-1, 'Pr(>|z|)'],
  row.names = genes
)
print(dat_multi)

dat_out <- data.frame(
  Gene = dat_multi$gene,
  `Multivariate OR (95% CI)` = paste0(
    sprintf('%0.2f', dat_multi$OR), ' (',
    sprintf('%0.2f', dat_multi$OR1), '-',
    sprintf('%0.2f', dat_multi$OR2), ')'
  ),
  `Multivariate P Value` = ifelse(
    dat_multi$pvalue < 0.001, '*** < 0.001',
    ifelse(dat_multi$pvalue < 0.01, '** < 0.01',
      ifelse(dat_multi$pvalue < 0.05, '* < 0.05', '> 0.05'))
  ),
  check.names = FALSE
)
write.csv(dat_out, 'output/Logistic_multi_res_SP_2gene.csv', row.names = FALSE)

# 森林图（与 forest_SP.pdf 同款样式）
lab <- rbind(colnames(dat_out), dat_out)
pdf(file = 'output/forest_SP_multi.pdf', height = 4.5, width = 8, onefile = FALSE)
forestplot(
  labeltext    = as.matrix(lab),
  mean         = c(NA, dat_multi$OR),
  lower        = c(NA, dat_multi$OR1),
  upper        = c(NA, dat_multi$OR2),
  ci.vertices  = TRUE,
  graph.pos    = 3,
  fn.ci_norm   = fpDrawCircleCI,
  col = fpColors(box = ggsci::pal_npg()(2)[2], lines = 'black', zero = 'black'),
  clip    = c(0.005, 100),
  xlog    = TRUE,
  boxsize = 0.15,
  zero    = 1,
  xlab    = 'Odds Ratio (log scale)',
  hrzl_lines = rlang::set_names(
    list(gpar(col = 'black'), gpar(col = 'black'), gpar(col = 'black')),
    c(1, 2, nrow(dat_multi) + 2)
  )
)
dev.off()
cat('完成: output/forest_SP_multi.pdf / output/Logistic_multi_res_SP_2gene.csv\n')
