# =====================================================================
# 第三疾病重叠: 置换零模型 (回应审稿 R2-M3 / R3 Tier2#7)
# 问题: 118 common DEGs 面板本身来自两个炎症条件, 与第三炎症疾病的
#       重叠是"构造性预期"; 超几何检验把基因当独立, 反保守。
# 两个零模型 (各 10000 次抽样):
#   Null1 表达匹配零模型: 从第三疾病 universe 中抽取与 pool118 在
#         该数据集 AveExpr 分位箱内计数完全相同的随机 118 基因集
#         -> 检验"任何可比可检出基因集"的重叠水平
#   Null2 炎症易感池零模型: 从 SP∪HP 训练线 DEG (P<0.05) 与第三
#         数据集 universe 的交集中抽取 118 基因 -> 直接检验
#         "面板只是任意炎症响应基因集合"这一构造性解释
# 输出: output/Third_disease_permutation_summary.csv
#       output/Third_disease_permutation_nulls.csv
# =====================================================================
rm(list = ls())
set.seed(2026)
suppressPackageStartupMessages({library(limma)})

outd <- 'output'
common118 <- as.character(read.csv('../2-Venn/output/SP-HP.CSV', row.names = 1, check.names = FALSE)$x)
stopifnot(length(common118) == 118)

## 训练线 DEG 并集 (炎症易感池): SP 与 HP 主线 limma, 名义 P<0.05
sp_tab <- readRDS('../1-limma/output/diff_res.Rds')
hp_tab <- readRDS('../1-limma2/output/HP-diff_res.Rds')
sp_deg <- rownames(sp_tab)[sp_tab$P.Value < 0.05]
hp_deg <- rownames(hp_tab)[hp_tab$P.Value < 0.05]
inflam_pool <- union(sp_deg, hp_deg)
cat('SP DEGs:', length(sp_deg), '| HP DEGs:', length(hp_deg),
    '| union (inflammation-permissive pool):', length(inflam_pool), '\n')

N_PERM <- 10000
contrasts <- list(
  list(file = 'GSE13355_Psoriasis_skin_PP-Control_all.csv', label = 'GSE13355_Psoriasis_skin', contrast = 'PP-Control'),
  list(file = 'GSE13355_Psoriasis_skin_PP-PN_all.csv',      label = 'GSE13355_Psoriasis_skin', contrast = 'PP-PN'),
  list(file = 'GSE50772_SLE_PBMC_SLE-Control_all.csv',      label = 'GSE50772_SLE_PBMC',       contrast = 'SLE-Control')
)

row_list <- list(); null_list <- list(); k <- 0
for (ct in contrasts) {
  tt <- read.csv(file.path(outd, ct$file), check.names = FALSE, stringsAsFactors = FALSE)
  rownames(tt) <- tt$gene
  universe <- tt$gene
  degs <- tt$gene[tt$change != 'Not']
  pool_in_u <- intersect(common118, universe)
  m <- length(pool_in_u); kd <- length(degs)
  q <- length(intersect(common118, degs))
  obs_fe <- q / (m * kd / length(universe))
  cat('\n=====', ct$label, ct$contrast, '| overlap', q, '/', m,
      '| observed FE =', round(obs_fe, 2), '\n')

  ## Null1: 表达分位箱匹配抽样
  av_pool <- tt[pool_in_u, 'AveExpr']
  brks <- quantile(av_pool, probs = seq(0, 1, length.out = 11))
  brks[1] <- -Inf; brks[11] <- Inf
  bin_of_pool <- cut(av_pool, breaks = brks, include.lowest = TRUE)
  n_per_bin <- table(bin_of_pool)
  cand <- setdiff(universe, common118)
  av_cand <- tt[cand, 'AveExpr']
  bin_cand <- cut(av_cand, breaks = brks, include.lowest = TRUE)
  cand_by_bin <- split(cand, bin_cand)
  o1 <- integer(N_PERM); fe1 <- numeric(N_PERM)
  for (i in seq_len(N_PERM)) {
    pick <- unlist(lapply(names(n_per_bin), function(b) {
      pool_b <- cand_by_bin[[b]]
      if (is.null(pool_b) || length(pool_b) == 0) return(character(0))
      sample(pool_b, min(n_per_bin[[b]], length(pool_b)))
    }), use.names = FALSE)
    o1[i] <- length(intersect(pick, degs))
    fe1[i] <- o1[i] / (m * kd / length(universe))
  }
  p1 <- (1 + sum(o1 >= q)) / (N_PERM + 1)

  ## Null2: 炎症易感池抽样 (union DEGs ∩ universe, 去掉 pool118 自身)
  pool2 <- setdiff(intersect(inflam_pool, universe), common118)
  n2 <- length(pool2)
  o2 <- if (n2 >= m) integer(N_PERM) else rep(NA_integer_, N_PERM)
  fe2 <- rep(NA_real_, N_PERM)
  if (n2 >= m) for (i in seq_len(N_PERM)) {
    pick <- sample(pool2, m)
    o2[i] <- length(intersect(pick, degs))
    fe2[i] <- o2[i] / (m * kd / length(universe))
  }
  p2 <- if (n2 >= m) (1 + sum(o2 >= q)) / (N_PERM + 1) else NA_real_

  cat(sprintf('Null1 (expr-matched): null overlap mean %.1f [95%% %d-%d], empirical P = %s\n',
              mean(o1), quantile(o1, .025), quantile(o1, .975), format.pval(p1, 2)))
  if (!is.na(p2)) cat(sprintf('Null2 (inflammation pool, n=%d): null overlap mean %.1f [95%% %d-%d], empirical P = %s\n',
              n2, mean(o2), quantile(o2, .025), quantile(o2, .975), format.pval(p2, 2)))

  k <- k + 1
  row_list[[k]] <- data.frame(
    dataset = ct$label, contrast = ct$contrast,
    pool_size = m, n_DEGs = kd, observed_overlap = q,
    hypergeom_P = signif(phyper(q - 1, m, length(universe) - m, kd, lower.tail = FALSE), 3),
    null1_mean_overlap = round(mean(o1), 1),
    null1_CI95 = paste0(quantile(o1, .025), '-', quantile(o1, .975)),
    null1_empirical_P = signif(p1, 3),
    null1_mean_FE = round(mean(fe1), 2),
    null2_pool_size = n2,
    null2_mean_overlap = if (!is.na(p2)) round(mean(o2), 1) else NA,
    null2_CI95 = if (!is.na(p2)) paste0(quantile(o2, .025), '-', quantile(o2, .975)) else NA,
    null2_empirical_P = if (!is.na(p2)) signif(p2, 3) else NA,
    null2_mean_FE = if (!is.na(p2)) round(mean(fe2), 2) else NA,
    n_perm = N_PERM, seed = 2026
  )
  null_list[[k]] <- data.frame(
    dataset = ct$label, contrast = ct$contrast,
    draw = c(1:N_PERM, 1:N_PERM),
    null = rep(c('expr_matched', 'inflammation_pool'), each = N_PERM),
    overlap = c(o1, o2)
  )
}

summary_tab <- do.call(rbind, row_list)
print(summary_tab, row.names = FALSE)
write.csv(summary_tab, file.path(outd, 'Third_disease_permutation_summary.csv'), row.names = FALSE)
write.csv(do.call(rbind, null_list), file.path(outd, 'Third_disease_permutation_nulls.csv'), row.names = FALSE)
saveRDS(list(summary = summary_tab, inflam_pool = inflam_pool), file.path(outd, 'results_permutation.Rds'))
cat('\nDONE\n')
