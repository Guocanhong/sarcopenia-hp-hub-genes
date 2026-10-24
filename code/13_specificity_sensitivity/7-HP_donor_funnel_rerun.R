# =====================================================================
# GSE27411 供者结构漏斗重跑 (回应审稿 R2-M2 / R3-M4)
# 问题: 主筛选把 12 活检当独立样本 (6 供者 x 2 部位), 供者水平检查
#       是事后而非筛选闸口。本脚本做三件事:
#  A. 供者结构尊重的重跑: 候选 = pool118 ∩ 供者校正 DEG (Model A,
#     duplicateCorrelation, 部位调整), 在供者聚合 3v3 数据 (6 样本)
#     上重跑 SVM-RFE 主筛选 (seed 2025, 镜像原漏斗), 报告四 hub 的名次
#  B. 精确枚举零模型: C(6,3)=20 种供者分组全枚举, 对每组做 3v3 limma,
#     统计四 hub 显著个数的零分布 -> 量化"3v3 完全分离"的偶然概率
#  C. WGCNA 供者级不可行性说明 (6 样本 < 最低要求, 仅记录)
# 输出: output/HP_donor_funnel_rerun_summary.csv
#       output/HP_donor_funnel_rfe_ranking.csv
#       output/HP_donor_funnel_enumeration.csv
# =====================================================================
rm(list = ls())
suppressMessages({library(limma); library(e1071); library(pROC)})
source('../input/msvmRFE.R')

expr <- readRDS('../1-limma2/input/dat1_expr.Rds')      # 18190 x 12, log2
grp  <- readRDS('../1-limma2/input/dat1_group.Rds')
pool <- as.character(read.csv('../2-Venn/output/SP-HP.CSV', row.names = 1, check.names = FALSE)$x)
orig <- readRDS('../1-limma2/output/HP-diff_res.Rds')
hubs <- c('C1QA','NFKBIA','DTL','CTNNBIP1')

ann <- data.frame(
  gsm   = c("GSM677730","GSM677731","GSM677732","GSM677733","GSM677734","GSM677735",
            "GSM677736","GSM677737","GSM677738","GSM677739","GSM677746","GSM677747"),
  donor = c("D23","D49","D23","D49","D8","D9","D30","D8","D9","D30","D35","D35"),
  site  = c("antrum","antrum","corpus","corpus","antrum","antrum","antrum",
            "corpus","corpus","corpus","antrum","corpus"),
  group = c("Control","Control","Control","Control","HP","HP","HP",
            "HP","HP","HP","Control","Control"),
  stringsAsFactors = FALSE)
stopifnot(setequal(ann$gsm, colnames(expr)))
expr <- expr[, ann$gsm]

## Model A: 12 biopsies, donor block + site adjust -----------------------------
designA <- model.matrix(~ ann$group + ann$site)
colnames(designA) <- c("Intercept","HPvsControl","sitecorpus")
corfit <- duplicateCorrelation(expr, designA, block = ann$donor)
fitA <- eBayes(lmFit(expr, designA, block = ann$donor, correlation = corfit$consensus))
ttA <- na.omit(topTable(fitA, coef = "HPvsControl", number = Inf))

## Model B: donor-aggregated 3v3 -----------------------------------------------
agg <- sapply(split(seq_len(ncol(expr)), ann$donor),
              function(ix) rowMeans(expr[, ix, drop = FALSE]))
dgrp <- ann$group[match(colnames(agg), ann$donor)]

cat('Model A consensus within-donor correlation:', round(corfit$consensus, 3), '\n')
cat('Pool118 ∩ Model-A DEGs (P<0.05):', sum(pool %in% rownames(ttA)[ttA$P.Value < 0.05]),
    '/', length(pool), '\n')
gate <- sapply(hubs, function(g) ttA[g, 'P.Value'] < 0.05)
cat('Hub survival at donor-corrected DEG gate:', paste(hubs[gate], collapse = ', '), '\n')

## A. donor-respecting SVM-RFE on 3v3 aggregated data ---------------------------
cand <- intersect(pool, rownames(ttA)[ttA$P.Value < 0.05])
mat <- agg[cand, , drop = FALSE]
dat <- data.frame(group = factor(dgrp, levels = c('Control','HP')), t(mat))
set.seed(2025)
rfe_full <- svmRFE(dat, k = 6, halve.above = 100)   # 6 samples -> 6-fold
best_first <- rev(as.integer(rfe_full))             # ranked.list is worst-first; reverse
rank_df <- data.frame(rank = seq_along(best_first), gene = rownames(mat)[best_first])
hub_rank <- rank_df[rank_df$gene %in% hubs, ]
cat('\nSVM-RFE ranking on donor-aggregated 3v3 data (top 20):\n')
print(head(rank_df, 20), row.names = FALSE)
cat('Hub ranks:\n'); print(hub_rank, row.names = FALSE)

## LASSO confirmatory on the same 3v3 data (cv.glmnet, seed 2025) ---------------
set.seed(2025)
suppressMessages(library(glmnet))
cvfit <- cv.glmnet(t(mat), dgrp, family = 'binomial', alpha = 1)
lco <- coef(cvfit, s = 'lambda.min')
lasso_genes <- rownames(lco)[as.numeric(lco) != 0]
lasso_genes <- setdiff(lasso_genes, '(Intercept)')
cat('\nLASSO (lambda.min) selected on 3v3 data:', length(lasso_genes),
    '| hubs in:', intersect(hubs, lasso_genes), '\n')

## B. exact enumeration of all 20 donor splits ----------------------------------
donors <- unique(ann$donor)
combos <- combn(donors, 3)               # choose 3 as "case" group, rest control
enum <- lapply(seq_len(ncol(combos)), function(j) {
  case_d <- combos[, j]
  g <- factor(ifelse(donors %in% case_d, 'HP', 'Control'), levels = c('Control','HP'))
  agg_p <- agg[, donors]
  designB <- model.matrix(~ g)
  fitB <- eBayes(lmFit(agg_p, designB))
  ttB <- na.omit(topTable(fitB, coef = 2, number = Inf))
  hub_p <- sapply(hubs, function(x) if (x %in% rownames(ttB)) ttB[x, 'P.Value'] else NA)
  data.frame(split = j, case_donors = paste(case_d, collapse = '+'),
             n_hub_sig_05 = sum(hub_p < 0.05, na.rm = TRUE),
             n_hub_sig_01 = sum(hub_p < 0.01, na.rm = TRUE),
             min_hub_P = signif(min(hub_p, na.rm = TRUE), 3),
             C1QA_P = signif(hub_p[1], 3), NFKBIA_P = signif(hub_p[2], 3),
             DTL_P = signif(hub_p[3], 3), CTNNBIP1_P = signif(hub_p[4], 3))
})
enum_tab <- do.call(rbind, enum)
true_split <- which(enum_tab$case_donors == 'D8+D9+D30')
cat('\nExact enumeration of 20 donor splits (3v3 limma):\n')
cat('  observed split =', true_split, '(case = D8+D9+D30)\n')
cat('  splits with all 4 hubs significant (P<0.05):',
    sum(enum_tab$n_hub_sig_05 == 4), '/ 20\n')
cat('  splits with >=3 hubs significant:',
    sum(enum_tab$n_hub_sig_05 >= 3), '/ 20\n')
cat('  median hubs significant across null splits:',
    median(enum_tab$n_hub_sig_05[-true_split]), '\n')
print(enum_tab, row.names = FALSE)

## summary -----------------------------------------------------------------------
summary_tab <- data.frame(
  item = c('pool_size', 'pool_survived_ModelA_gate', 'hub_survived_ModelA_gate',
           'hub_min_rank_in_RFE_3v3', 'hub_in_LASSO_3v3',
           'enum_splits_all4sig', 'enum_splits_ge3sig', 'enum_observed_split_rank_note'),
  value = c(length(pool),
            sum(pool %in% rownames(ttA)[ttA$P.Value < 0.05]),
            paste(hubs[gate], collapse = '|'),
            if (nrow(hub_rank)) paste(paste0(hub_rank$gene, '=', hub_rank$rank), collapse = '|') else 'none_ranked',
            paste(intersect(hubs, lasso_genes), collapse = '|'),
            sum(enum_tab$n_hub_sig_05 == 4), sum(enum_tab$n_hub_sig_05 >= 3),
            paste0('case=D8+D9+D30 is split #', true_split, ' of 20'))
)
write.csv(summary_tab, 'output/HP_donor_funnel_rerun_summary.csv', row.names = FALSE)
write.csv(rank_df, 'output/HP_donor_funnel_rfe_ranking.csv', row.names = FALSE)
write.csv(enum_tab, 'output/HP_donor_funnel_enumeration.csv', row.names = FALSE)
cat('\n[DONE]\n')
