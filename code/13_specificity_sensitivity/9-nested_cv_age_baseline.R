# =====================================================================
# 审稿响应分析: nomogram 全流程嵌套交叉验证 + age-only 基线对照
# 回应 (模拟评审 R1-M4 / R2-M4 / R2-M6 / R3-m9):
#   A. 嵌套 CV: 把 §2.7 预设定的单变量筛选 (4 hub 基因, univariable logistic
#      P < 0.05) 嵌套进每一折内重做, 得到全流程诚实性能估计;
#      LOOCV 与 200 次重复分层五折 CV (seed 2026) 两种口径。
#   B. age-only 基线: 在各外部年龄可及队列中, 以个体年龄为唯一预测变量
#      计算 AUC, 与 nomogram 锁定评分 AUC 并列 (R2-M6 的量化增量问题)。
#      GSE25941 的 age 字段为 Young/Old 类别 (与分组定义重合), GSE1428 的
#      对比即按年龄定义, 两者按设计同义反复, 不做数值比较, 仅记录。
# 输出: output/Nested_CV_age_baseline.csv + output/nested_cv_reps.csv
# =====================================================================
rm(list = ls())
suppressPackageStartupMessages({library(pROC); library(data.table); library(AnnotationDbi); library(org.Hs.eg.db)})
outd <- 'output'; dir.create(outd, showWarnings = FALSE, recursive = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a
G2 <- c('NFKBIA', 'CTNNBIP1'); G4 <- c('C1QA', 'NFKBIA', 'DTL', 'CTNNBIP1')
rows_out <- list(); add_row <- function(...) rows_out[[length(rows_out)+1]] <<- data.frame(...)

## ============ A. 嵌套 CV (GSE1428) =============================================
sp_expr <- as.matrix(readRDS('../1-limma/input/dat_expr.Rds'))
sp_grp  <- readRDS('../1-limma/input/dat_group.Rds')
stopifnot(identical(colnames(sp_expr), sp_grp$sample))
y <- as.integer(sp_grp$group == 'SP')
X4 <- t(sp_expr[G4, ])          # n x 4
n <- length(y); cat('GSE1428: n =', n, '| SP(older):', sum(y), '/ young:', sum(1-y), '\n')

# 全数据单变量筛选 (复核 §2.7 报告的 P 值)
uni_p <- sapply(G4, function(g) summary(glm(y ~ X4[, g], family = binomial))$coefficients[2, 4])
cat('全数据单变量 P:\n'); print(round(uni_p, 4))

screen_fit <- function(Xtr, ytr, xte) {
  ## 单折内: 单变量 logistic P<0.05 筛选 -> 多变量 logistic (失败回退 z 和分)
  ps <- sapply(colnames(Xtr), function(g)
    tryCatch(summary(glm(ytr ~ Xtr[, g], family = binomial))$coefficients[2, 4], error = function(e) 1))
  sel <- names(ps)[ps < 0.05]
  if (length(sel) == 0) return(list(pred = 0.5, sel = sel, mode = 'none-passed'))
  if (length(sel) == 1) {
    m <- mean(Xtr[, sel]); s <- sd(Xtr[, sel])
    return(list(pred = (xte[sel] - m) / max(s, 1e-9), sel = sel, mode = 'single-gene'))
  }
  dtr <- data.frame(y = ytr, Xtr[, sel, drop = FALSE])
  fi <- try(suppressWarnings(glm(y ~ ., family = binomial, data = dtr)), silent = TRUE)
  if (inherits(fi, 'try-error')) {
    m <- colMeans(Xtr[, sel, drop = FALSE]); s <- apply(Xtr[, sel, drop = FALSE], 2, sd)
    z <- (xte[sel] - m) / pmax(s, 1e-9)
    return(list(pred = sum(z), sel = sel, mode = 'glm-failed-zscore'))
  }
  nd <- as.data.frame(t(xte)); colnames(nd) <- names(xte)
  list(pred = as.numeric(predict(fi, newdata = nd[, sel, drop = FALSE],
                                 type = 'response')),
       sel = sel, mode = 'glm')
}

# ---- A1. 嵌套 LOOCV ----
pred_loo <- rep(NA_real_, n); sel_count <- list(); modes <- character(0)
for (i in seq_len(n)) {
  tr <- seq_len(n)[-i]
  r <- screen_fit(X4[tr, ], y[tr], X4[i, ])
  pred_loo[i] <- r$pred; modes[i] <- r$mode
  sel_count[[paste(r$sel, collapse = '+')]] <-
    (sel_count[[paste(r$sel, collapse = '+')]] %||% 0) + 1
}
auc_loo_nested <- as.numeric(pROC::roc(y ~ pred_loo, quiet = TRUE)$auc)
cat('\n[嵌套 LOOCV] AUC =', round(auc_loo_nested, 3),
    '| 折内模式:', paste(names(table(modes)), table(modes), collapse = '; '), '\n')
cat('折内选中基因组合分布:\n'); print(unlist(sel_count))

# ---- A2. 嵌套重复五折 CV (200 次, seed 2026) ----
set.seed(2026); B <- 200; k <- 5
aucs <- numeric(B); fails <- 0; mode_tab <- table(modes)
for (b in seq_len(B)) {
  fold <- integer(n)
  is <- sample(which(y == 1)); ic <- sample(which(y == 0))
  fold[is] <- rep(1:k, length.out = length(is))
  fold[ic] <- rep(1:k, length.out = length(ic))
  pred <- rep(NA_real_, n)
  for (kk in 1:k) {
    te <- which(fold == kk); tr <- which(fold != kk)
    pred[te] <- sapply(te, function(i) screen_fit(X4[tr, ], y[tr], X4[i, ])$pred)
  }
  if (any(!is.finite(pred)) || length(unique(pred)) < 2) { fails <- fails + 1; next }
  aucs[b] <- as.numeric(pROC::roc(y ~ pred, quiet = TRUE)$auc)
}
valid <- aucs[is.finite(aucs)]
q <- quantile(valid, c(.025, .5, .975))
cat(sprintf('[嵌套重复五折CV] B=%d 有效=%d 失败=%d | median=%.3f (95%%PI %.3f-%.3f)\n',
            B, length(valid), fails, q[2], q[1], q[3]))

add_row(setting = 'GSE1428 nested LOOCV (screen refit per fold)', n = n,
        estimate = round(auc_loo_nested, 3),
        detail = paste('modes:', paste(names(table(modes)), table(modes), collapse = '; ')))
add_row(setting = 'GSE1428 nested repeated 5-fold CV (200 reps, seed 2026)', n = n,
        estimate = round(unname(q[2]), 3),
        detail = sprintf('95%% percentile interval %.3f-%.3f, failed reps %d', q[1], q[3], fails))

# 基准 (固定 2 基因, 不嵌套) 复核
lp_fixed <- numeric(n)
X2 <- t(sp_expr[G2, ])
for (i in seq_len(n)) {
  fi <- suppressWarnings(glm(y ~ X2, family = binomial, subset = -i))
  lp_fixed[i] <- predict(fi, newdata = data.frame(X2)[i, , drop = FALSE], type = 'response')
}
auc_loo_fixed <- as.numeric(pROC::roc(y ~ lp_fixed, quiet = TRUE)$auc)
cat('[固定2基因LOOCV复核]', round(auc_loo_fixed, 3), '\n')
add_row(setting = 'GSE1428 LOOCV, fixed 2-gene structure (benchmark)', n = n,
        estimate = round(auc_loo_fixed, 3), detail = 'screen not nested (as published)')

## ============ B. age-only 基线 (外部队列) ======================================
coef_lock <- coef(glm(y ~ X2, family = binomial))[-1]; names(coef_lock) <- G2
cat('\n锁定系数:', paste(G2, round(coef_lock, 4), collapse = ', '), '\n')

score_roc <- function(e, grp, test, ref) {
  keep <- grp %in% c(test, ref)
  ee <- e[, keep, drop = FALSE]; gg <- grp[keep]
  zw <- sweep(sweep(ee, 1, rowMeans(ee, na.rm = TRUE), '-'), 1, apply(ee, 1, sd, na.rm = TRUE), '/')
  sc <- as.numeric(coef_lock[1] * zw[1, ] + coef_lock[2] * zw[2, ])
  r <- pROC::roc(gg, sc, levels = c(ref, test), direction = '<', quiet = TRUE)
  list(auc = as.numeric(r$auc), score = setNames(sc, colnames(ee)), gg = gg)
}
age_roc <- function(age, gg, test, ref) {
  ok <- is.finite(age)
  r <- pROC::roc(gg[ok], age[ok], levels = c(ref, test), direction = '<', quiet = TRUE)
  as.numeric(r$auc)
}

# ---- GSE167186 ----
cnt <- fread('../13-GPL20301/raw/bulk_counts.csv.gz', header = TRUE)
mat <- as.matrix(cnt[, -1]); rownames(mat) <- as.character(cnt[[1]])
con <- gzfile('../13-GPL20301/raw/GSE167186-GPL20301_series_matrix.txt.gz')
lines <- readLines(con); close(con)
pick <- function(pfx) {
  ln <- lines[grepl('^!Sample_characteristics_ch1', lines)]
  out <- rep(NA_character_, length(strsplit(ln[1], '\t')[[1]]) - 1)
  for (l in ln) {
    v <- gsub('^"|"$', '', strsplit(l, '\t', fixed = TRUE)[[1]][-1])
    hit <- grep(pfx, v); if (length(hit)) out[is.na(out)] <- sub(pfx, '', v[hit][is.na(out)])
  }
  out
}
titles <- gsub('^"|"$', '', strsplit(lines[grepl('^!Sample_title', lines)][1], '\t')[[1]][-1])
grp167 <- pick('^group: '); age167 <- suppressWarnings(as.numeric(pick('^age: ')))
lib <- colSums(mat); lcpm <- log2(sweep(mat, 2, lib, '/') * 1e6 + 1)
stopifnot(all(titles %in% colnames(lcpm)))
e167 <- lcpm[G2, titles]
for (ct in list(c('Sarcopenia', 'Old Healthy', 'Sarco_vs_OldMatched'),
                c('Sarcopenia+Old Healthy', 'Young Healthy', 'SarcoPlusOld_vs_Young'))) {
  g2 <- grp167
  for (g in strsplit(ct[1], '\\+')[[1]]) g2 <- ifelse(grp167 == g, 'TEST', g2)
  g2 <- ifelse(grp167 == ct[2], 'REF', g2)
  sr <- score_roc(e167, g2, 'TEST', 'REF')
  a <- age_roc(age167, g2, 'TEST', 'REF')
  ct_keep <- g2 %in% c('TEST', 'REF')
  ct_sp <- suppressWarnings(cor(sr$score, age167[ct_keep], method = 'spearman',
                                use = 'complete.obs'))
  cat(sprintf('\n[GSE167186 %s] nomogram AUC=%.3f | age-only AUC=%.3f | rho(score,age)=%.3f\n',
              ct[3], sr$auc, a, ct_sp))
  add_row(setting = paste0('GSE167186 ', ct[3], ' nomogram score'), n = sum(ct_keep),
          estimate = round(sr$auc, 3),
          detail = sprintf('Spearman rho(score, age) = %.3f', ct_sp))
  add_row(setting = paste0('GSE167186 ', ct[3], ' age-only baseline'), n = sum(ct_keep),
          estimate = round(a, 3), detail = 'individual age as sole predictor')
}

# ---- GSE111016 ----
lcpm <- read.csv(gzfile('../external/GSE111016_lcpm_sss.csv.gz'), row.names = 1, check.names = FALSE)
sm <- readLines(gzfile('../external/GSE111016_series_matrix.txt.gz'))
get_line <- function(tag) {
  ln <- grep(paste0('^!', tag), sm, value = TRUE)[1]
  strsplit(ln, '\t')[[1]][-1] |> (\(x) gsub('"', '', trimws(x)))()
}
titles16 <- get_line('Sample_title')
grp16 <- NULL; age16 <- NULL
for (ln in sm[grepl('^!Sample_characteristics_ch1', sm)]) {
  v <- strsplit(ln, '\t')[[1]][-1] |> (\(x) gsub('"', '', trimws(x)))()
  if (any(grepl('sarcopenia status', v)))
    grp16 <- ifelse(grepl(': yes', v), 'SP', 'Control')
  if (any(grepl('age \\(yr\\)', v)))
    age16 <- suppressWarnings(as.numeric(sub('.*age \\(yr\\): ', '', v)))
}
sid <- sub(' \\[.*\\]', '', titles16)
ens <- gsub('\\..*$', '', rownames(lcpm))
mp <- AnnotationDbi::select(org.Hs.eg.db, keys = ens, keytype = 'ENSEMBL', columns = 'SYMBOL')
mp <- mp[!duplicated(mp$ENSEMBL) & !is.na(mp$SYMBOL), ]
sym <- mp$SYMBOL[match(ens, mp$ENSEMBL)]
lcpm <- lcpm[!is.na(sym), ]; sym <- sym[!is.na(sym)]
o <- order(-rowMeans(lcpm, na.rm = TRUE), sym); lcpm <- lcpm[o, ]; sym <- sym[o]
f <- !duplicated(sym); lcpm <- lcpm[f, ]; sym <- sym[f]; rownames(lcpm) <- sym
e16 <- as.matrix(lcpm[G2, sid])
sr <- score_roc(e16, grp16, 'SP', 'Control')
a <- age_roc(age16, grp16, 'SP', 'Control')
ct_rho <- suppressWarnings(cor(sr$score, age16, method = 'spearman', use = 'complete.obs'))
cat(sprintf('\n[GSE111016 SP vs Control(older)] nomogram AUC=%.3f | age-only AUC=%.3f | rho=%.3f\n',
            sr$auc, a, ct_rho))
add_row(setting = 'GSE111016 SP_vs_Control_older nomogram score', n = ncol(e16),
        estimate = round(sr$auc, 3),
        detail = sprintf('Spearman rho(score, age) = %.3f', ct_rho))
add_row(setting = 'GSE111016 SP_vs_Control_older age-only baseline', n = length(grp16),
        estimate = round(a, 3), detail = 'individual age 65-79 as sole predictor')

# ---- 同义反复杂项记录 ----
add_row(setting = 'GSE25941 Old_vs_Young age-only baseline', n = NA, estimate = NA,
        detail = 'age field is categorical Young/Old and defines the contrast; tautological, not computed')
add_row(setting = 'GSE1428 training contrast age-only', n = n, estimate = NA,
        detail = 'contrast defined by donor age group by design; tautological, not computed')

res <- do.call(rbind, rows_out)
print(res, row.names = FALSE)
write.csv(res, file.path(outd, 'Nested_CV_age_baseline.csv'), row.names = FALSE)
write.csv(data.frame(rep = seq_along(valid), nested_cv_auc = round(valid, 4)),
          file.path(outd, 'nested_cv_reps.csv'), row.names = FALSE)
cat('\nDONE\n')
