# =====================================================================
# R3-m7 更正: LOOCV AUC 的置信区间方法
# 问题: 原用 pROC::ci.auc (DeLong) 于 out-of-fold 预测, 但 DeLong 假设观测独立,
#       交叉验证预测跨折共享训练子集, 违反独立性 (且区间 0.599-1.000 顶格可疑)。
# 更正: 200 次重复的分层五折 CV (seed=2026), 每次重复在原 22 样本上重新划分、
#       每样本恰被预测一次且必然在训练集之外; 取 200 个 CV AUC 的 2.5/50/97.5 百分位。
# 附:   复核 LOOCV 点估计 0.817 不变。
# =====================================================================
rm(list = ls())
suppressPackageStartupMessages(library(pROC))
set.seed(2026)

sp_expr <- as.matrix(readRDS('../1-limma/input/dat_expr.Rds'))
sp_grp  <- readRDS('../1-limma/input/dat_group.Rds')
stopifnot(identical(colnames(sp_expr), sp_grp$sample))
y <- as.integer(sp_grp$group == 'SP')
cat('n =', length(y), '| SP(older):', sum(y), '/ young:', sum(1 - y), '\n')
X <- t(sp_expr[c('NFKBIA', 'CTNNBIP1'), ])

## ---- 1. LOOCV 点估计 (复核) ------------------------------------------
lp <- numeric(length(y))
for (i in seq_along(y)) {
  fi <- suppressWarnings(glm(y ~ X, family = binomial, subset = -i))
  lp[i] <- predict(fi, newdata = data.frame(X)[i, , drop = FALSE], type = 'response')
}
auc_loo <- as.numeric(pROC::roc(y ~ lp, quiet = TRUE)$auc)
cat('LOOCV AUC (复核):', round(auc_loo, 4), '\n')

## ---- 2. 200 次重复分层五折 CV ----------------------------------------
B <- 200; k <- 5
aucs <- numeric(B); fails <- 0
idx_sp <- which(y == 1); idx_ct <- which(y == 0)
for (b in seq_len(B)) {
  fold <- integer(length(y))
  is <- sample(idx_sp); ic <- sample(idx_ct)
  fold[is] <- rep(1:k, length.out = length(is))
  fold[ic] <- rep(1:k, length.out = length(ic))
  pred <- rep(NA_real_, length(y)); ok <- TRUE
  for (kk in 1:k) {
    te <- which(fold == kk); tr <- which(fold != kk)
    fi <- try(suppressWarnings(glm(y ~ X, family = binomial, subset = tr)), silent = TRUE)
    if (inherits(fi, 'try-error')) { ok <- FALSE; break }
    pred[te] <- predict(fi, newdata = data.frame(X)[te, , drop = FALSE], type = 'response')
  }
  if (!ok || any(!is.finite(pred)) || length(unique(pred)) < 2) { fails <- fails + 1; next }
  aucs[b] <- as.numeric(pROC::roc(y ~ pred, quiet = TRUE)$auc)
}
valid <- aucs[is.finite(aucs)]
q <- quantile(valid, c(.025, .5, .975))
cat(sprintf('repeated 5-fold CV: B=%d 有效=%d 失败=%d\n', B, length(valid), fails))
cat(sprintf('median=%.3f  95%% percentile CI (%.3f-%.3f)\n', q[2], q[1], q[3]))

## ---- 3. 落盘 ---------------------------------------------------------
dir.create('output', showWarnings = FALSE)
write.csv(data.frame(
  method = 'repeated stratified 5-fold CV, 200 repeats, seed 2026',
  loocv_auc = round(auc_loo, 4),
  median_cv_auc = round(unname(q[2]), 4),
  ci_lower = round(unname(q[1]), 4),
  ci_upper = round(unname(q[3]), 4),
  n_valid = length(valid), n_failed = fails
), 'output/LOOCV_CI_corrected.csv', row.names = FALSE)
write.csv(data.frame(rep = seq_along(valid), cv_auc = round(valid, 4)),
          'output/LOOCV_CI_reps.csv', row.names = FALSE)
cat('DONE\n')
