# =====================================================================
# 方案③ 模型性能诚实化 — locked-model 外部转移 + 训练集 LOOCV
# 背景: 论文 nomogram (NFKBIA+CTNNBIP1, SP 线 GSE1428 训练) 的 C-index=0.950
#       为重代入估计 (特征筛选与 ROC 拟合同样本)。
# 交付:
#   A. GSE1428 训练集: apparent AUC + LOOCV AUC (逐折重拟合 2 基因 logistic);
#      4 基因模型 complete separation 记录在案
#   B. 锁定系数跨队列转移: 训练 beta 应用到 4 个外部队列
#      (各队列内 z 标准化后按锁定权重合成线性评分; 训练锚定 z 作敏感性)
#      队列: GSE167186 (Sarco vs Old 年龄匹配; Sarco+Old vs Young 论文主对比)
#            GSE111016 (SP vs Control, 老年) / GSE25941 (Old vs Young, GPL570)
#            GSE60662 (HP vs Control, 胃黏膜, 跨线检验)
# 输出: output/Locked_model_performance.csv + results_locked.Rds
# =====================================================================
rm(list = ls())
suppressPackageStartupMessages({library(pROC); library(data.table); library(AnnotationDbi); library(org.Hs.eg.db)})

outd <- 'output'; dir.create(outd, showWarnings = FALSE, recursive = TRUE)
G2 <- c('NFKBIA', 'CTNNBIP1')   # nomogram 基因
G4 <- c('C1QA', 'NFKBIA', 'DTL', 'CTNNBIP1')

## ============ A. 训练集 (GSE1428): apparent + LOOCV ============================
sp_expr <- as.matrix(readRDS('../1-limma/input/dat_expr.Rds'))
sp_grp  <- readRDS('../1-limma/input/dat_group.Rds')
stopifnot(identical(colnames(sp_expr), sp_grp$sample))
y <- as.integer(sp_grp$group == 'SP')
cat('GSE1428: n =', length(y), '| SP:', sum(y), '/ Control:', sum(1 - y), '\n')

X <- t(sp_expr[G2, ])   # n x 2
fit_full <- glm(y ~ X, family = binomial)
coef_lock <- coef(fit_full)[-1]; names(coef_lock) <- G2
auc_app <- as.numeric(pROC::roc(y ~ predict(fit_full, type = 'response'), quiet = TRUE)$auc)
cat('apparent AUC (2-gene):', round(auc_app, 3), '\n')
cat('锁定系数:', paste(G2, round(coef_lock, 4), collapse = ', '), '\n')

# 4-gene 模型: 预期 complete separation
X4 <- t(sp_expr[G4, ])
fit4 <- try(glm(y ~ X4, family = binomial), silent = TRUE)
sep_note <- if (inherits(fit4, 'try-error')) 'complete separation (non-convergent)' else
  if (any(!is.finite(coef(fit4))) || max(abs(coef(fit4)[-1])) > 10) 'complete separation (divergent coefficients)' else 'converged'
cat('4-gene coefficients:', paste(round(coef(fit4)[-1], 1), collapse = ', '), '\n')
cat('4-gene 模型状态:', sep_note, '\n')

# LOOCV (逐折重拟合)
set.seed(2026)
lp_loo <- numeric(length(y))
for (i in seq_along(y)) {
  fi <- glm(y ~ X, family = binomial, subset = -i)
  lp_loo[i] <- predict(fi, newdata = data.frame(X)[i, , drop = FALSE], type = 'response')
}
auc_loo <- as.numeric(pROC::roc(y ~ lp_loo, quiet = TRUE)$auc)
ci_loo  <- as.numeric(pROC::ci.auc(y ~ lp_loo, quiet = TRUE))
cat('LOOCV AUC (2-gene):', round(auc_loo, 3),
    sprintf('(%.3f-%.3f)', ci_loo[1], ci_loo[3]), '\n')

## ============ B. 外部队列加载函数 ==============================================
load_167186 <- function() {
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
  grp <- pick('^group: ')
  lib <- colSums(mat); lcpm <- log2(sweep(mat, 2, lib, '/') * 1e6 + 1)
  stopifnot(all(titles %in% colnames(lcpm)))
  e <- lcpm[G2, titles]
  list(e = e, grp = grp,
       contrasts = list(
         c('Sarcopenia', 'Old Healthy', 'Sarco_vs_OldMatched'),
         c('Sarcopenia+Old Healthy', 'Young Healthy', 'SarcoPlusOld_vs_Young')))
}
load_111016 <- function() {
  lcpm <- read.csv(gzfile('../external/GSE111016_lcpm_sss.csv.gz'), row.names = 1, check.names = FALSE)
  sm <- readLines(gzfile('../external/GSE111016_series_matrix.txt.gz'))
  get_line <- function(tag) {
    ln <- grep(paste0('^!', tag), sm, value = TRUE)[1]
    strsplit(ln, '\t')[[1]][-1] |> (\(x) gsub('"', '', trimws(x)))()
  }
  titles <- get_line('Sample_title')
  chars <- sm[grepl('^!Sample_characteristics_ch1', sm)]
  grp <- NULL
  for (ln in chars) {
    v <- strsplit(ln, '\t')[[1]][-1] |> (\(x) gsub('"', '', trimws(x)))()
    if (any(grepl('sarcopenia status', v))) grp <- ifelse(grepl(': yes', v), 'SP', 'Control')
  }
  sid <- sub(' \\[.*\\]', '', titles)
  ens <- gsub('\\..*$', '', rownames(lcpm))
  mp <- AnnotationDbi::select(org.Hs.eg.db, keys = ens, keytype = 'ENSEMBL',
                              columns = 'SYMBOL')
  mp <- mp[!duplicated(mp$ENSEMBL) & !is.na(mp$SYMBOL), ]
  sym <- mp$SYMBOL[match(ens, mp$ENSEMBL)]
  lcpm <- lcpm[!is.na(sym), ]; sym <- sym[!is.na(sym)]
  o <- order(-rowMeans(lcpm, na.rm = TRUE), sym); lcpm <- lcpm[o, ]; sym <- sym[o]
  f <- !duplicated(sym); lcpm <- lcpm[f, ]; sym <- sym[f]
  rownames(lcpm) <- sym
  e <- as.matrix(lcpm[G2, sid])
  list(e = e, grp = grp, contrasts = list(c('SP', 'Control', 'SP_vs_Control_older')))
}
load_25941 <- function() {
  sm <- readLines(gzfile('../external/GSE25941_series_matrix.txt.gz'))
  get_line <- function(tag) {
    ln <- grep(paste0('^!', tag), sm, value = TRUE)[1]
    strsplit(ln, '\t')[[1]][-1] |> (\(x) gsub('"', '', trimws(x)))()
  }
  gsm <- get_line('Sample_geo_accession')
  age <- NULL
  for (ln in sm[grepl('^!Sample_characteristics_ch1', sm)]) {
    v <- strsplit(ln, '\t')[[1]][-1] |> (\(x) gsub('"', '', trimws(x)))()
    if (any(grepl('^age:', v))) age <- sub('^age: ', '', v)
  }
  grp <- ifelse(age == 'Young', 'Control', 'SP')
  ts <- grep('^!series_matrix_table_begin', sm) + 1
  te <- grep('^!series_matrix_table_end', sm) - 1
  expr <- read.table(text = paste(sm[ts:te], collapse = '\n'), sep = '\t',
                     header = TRUE, quote = '', row.names = 1, check.names = FALSE)
  expr <- as.matrix(expr)
  colnames(expr) <- gsub('"', '', colnames(expr)); rownames(expr) <- gsub('"', '', rownames(expr))
  ann <- readLines(gzfile('../external/GPL570.annot.gz'))
  ann <- ann[grepl('^[0-9]', ann)]
  a <- read.table(text = paste(ann, collapse = '\n'), sep = '\t', header = FALSE,
                  quote = '', fill = TRUE, comment.char = '', check.names = FALSE)
  m <- a[, c(1, 3)]; colnames(m) <- c('probe', 'symbol')
  m$symbol <- sub('///.*', '', m$symbol)
  m <- m[m$symbol != '' & !is.na(m$symbol) & !duplicated(m$probe), ]
  hit <- match(rownames(expr), m$probe); sym <- m$symbol[hit]
  k <- !is.na(sym); expr <- expr[k, ]; sym <- sym[k]
  o <- order(sym, -rowMeans(expr)); expr <- expr[o, ]; sym <- sym[o]
  f <- !duplicated(sym); expr <- expr[f, ]; rownames(expr) <- sym[f]
  e <- log2(expr[G2, gsm] + 1)
  list(e = e, grp = grp, contrasts = list(c('SP', 'Control', 'Old_vs_Young')))
}
load_60662 <- function() {
  d <- readRDS('../27-GSE60662/output/GSE60662_analysis_input.Rds')
  e <- t(log2(as.matrix(d[, G2]) + 1))  # 2 x n, 行名 = G2
  grp <- as.character(d$group)
  list(e = e, grp = grp, contrasts = list(c('HP', 'Control', 'HP_vs_Control_gastric')))
}

## ============ C. 锁定系数转移 ==================================================
transfer <- function(e, grp, test, ref, label) {
  keep <- grp %in% c(test, ref)
  ee <- e[, keep, drop = FALSE]; gg <- grp[keep]
  zw <- sweep(ee, 1, rowMeans(ee, na.rm = TRUE), '-')  # 队列内 z
  zw <- sweep(zw, 1, apply(ee, 1, sd, na.rm = TRUE), '/')
  za <- sweep(ee, 1, colMeans(X), '-')           # 训练锚定 z
  za <- sweep(za, 1, apply(X, 2, sd), '/')
  sc_w <- as.numeric(coef_lock[1] * zw[1, ] + coef_lock[2] * zw[2, ])
  sc_a <- as.numeric(coef_lock[1] * za[1, ] + coef_lock[2] * za[2, ])
  gg_f <- factor(gg, levels = c(ref, test))
  make_roc <- function(s) {
    ok <- is.finite(s)
    r <- pROC::roc(response = gg_f[ok], predictor = s[ok], quiet = TRUE, ci = TRUE,
                   levels = c(ref, test), direction = '<')
    c(AUC = as.numeric(r$auc), low = as.numeric(r$ci)[1], high = as.numeric(r$ci)[3])
  }
  rw <- make_roc(sc_w); ra <- make_roc(sc_a)
  data.frame(cohort = label, n_test = sum(gg == test), n_ref = sum(gg == ref),
             AUC_withinZ = round(rw['AUC'], 3),
             CI_withinZ = sprintf('%.3f-%.3f', rw['low'], rw['high']),
             AUC_trainAnchorZ = round(ra['AUC'], 3),
             CI_trainAnchorZ = sprintf('%.3f-%.3f', ra['low'], ra['high']))
}

rows <- list()
d <- load_167186()
for (ct in d$contrasts) {
  g2 <- d$grp
  for (g in strsplit(ct[1], '\\+')[[1]]) g2 <- ifelse(d$grp == g, 'TEST', g2)
  g2 <- ifelse(d$grp == ct[2], 'REF', g2)
  rows[[ct[3]]] <- transfer(d$e, g2, 'TEST', 'REF', paste0('GSE167186 ', ct[3]))
}
d <- load_111016()
rows[[d$contrasts[[1]][3]]] <- transfer(d$e, d$grp, d$contrasts[[1]][1], d$contrasts[[1]][2],
                                        paste0('GSE111016 ', d$contrasts[[1]][3]))
d <- load_25941()
rows[[d$contrasts[[1]][3]]] <- transfer(d$e, d$grp, d$contrasts[[1]][1], d$contrasts[[1]][2],
                                        paste0('GSE25941 ', d$contrasts[[1]][3]))
d <- load_60662()
rows[[d$contrasts[[1]][3]]] <- transfer(d$e, d$grp, d$contrasts[[1]][1], d$contrasts[[1]][2],
                                        paste0('GSE60662 ', d$contrasts[[1]][3]))

ext_tab <- do.call(rbind, rows)
perf_tab <- rbind(
  data.frame(cohort = 'GSE1428 training (apparent, resubstitution)',
             n_test = sum(y), n_ref = sum(1 - y),
             AUC_withinZ = round(auc_app, 3), CI_withinZ = NA,
             AUC_trainAnchorZ = NA, CI_trainAnchorZ = NA),
  data.frame(cohort = 'GSE1428 training (LOOCV, refit per fold)',
             n_test = sum(y), n_ref = sum(1 - y),
             AUC_withinZ = round(auc_loo, 3),
             CI_withinZ = sprintf('%.3f-%.3f', ci_loo[1], ci_loo[3]),
             AUC_trainAnchorZ = NA, CI_trainAnchorZ = NA),
  ext_tab)
cat('\n===== Nomogram (NFKBIA + CTNNBIP1) 性能汇总 =====\n')
print(perf_tab, row.names = FALSE)
write.csv(perf_tab, file.path(outd, 'Locked_model_performance.csv'), row.names = FALSE)
saveRDS(list(perf = perf_tab, coef_lock = coef_lock,
             auc_app = auc_app, auc_loo = auc_loo, sep_note = sep_note),
        file.path(outd, 'results_locked.Rds'))
cat('\n锁定系数:', paste(G2, round(coef_lock, 4), collapse = ', '), '\n')
cat('4-gene:', sep_note, '| LOOCV AUC:', round(auc_loo, 3), '\nDONE\n')
