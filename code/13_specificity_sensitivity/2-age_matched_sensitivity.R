# =====================================================================
# 方案② 表型锚定敏感性分析 — GSE111016 (Singapore Sarcopenia Study)
# 问题: 四 hub 基因只是 "young vs old" 的 ageing 标志, 还是真的区分 sarcopenia 表型?
# 设计: 40 例老年男性 (65-79y) 股外侧肌, 20 sarcopenic vs 20 non-sarcopenic
#   (a) 全队列 limma: lcpm ~ status            (报告 logFC/P/方向/AUC)
#   (b) 年龄协变量调整: lcpm ~ age + status
#   (c) 1:1 最近邻年龄匹配子集 (caliper 3y)
#   (d) 118 common DEGs 池在本对比中的存活率 (超几何) — 量化 "签名 = 衰老?" 假设
# 输入: external/GSE111016_lcpm_sss.csv.gz (ENSG 行, log2 CPM)
#       external/GSE111016_series_matrix.txt.gz (status/age/sex)
# 输出: output/GSE111016_age_matched_*.csv + results_age_matched.Rds
# =====================================================================
rm(list = ls())
suppressPackageStartupMessages({library(limma); library(pROC); library(AnnotationDbi); library(org.Hs.eg.db)})

outd <- 'output'; dir.create(outd, showWarnings = FALSE, recursive = TRUE)
HUB <- c('C1QA', 'NFKBIA', 'DTL', 'CTNNBIP1')

## 0. 主线 SP 训练方向 (GSE1428 old vs young) -----------------------------------
sp_tab <- readRDS('../1-limma/output/diff_res.Rds')

## 1. 解析 series matrix --------------------------------------------------------
sm <- readLines(gzfile('../external/GSE111016_series_matrix.txt.gz'))
get_line <- function(tag) {
  ln <- grep(paste0('^!', tag), sm, value = TRUE)[1]
  strsplit(ln, '\t')[[1]][-1] |> (\(x) gsub('"', '', trimws(x)))()
}
titles <- get_line('Sample_title')
chars  <- sm[grepl('^!Sample_characteristics_ch1', sm)]
status <- age <- sex <- NULL
for (ln in chars) {
  v <- strsplit(ln, '\t')[[1]][-1] |> (\(x) gsub('"', '', trimws(x)))()
  if (any(grepl('sarcopenia status', v))) status <- ifelse(grepl(': yes', v), 'SP', 'Control')
  if (any(grepl('^age \\(yr\\):', v)))    age    <- as.numeric(sub('^age \\(yr\\): ', '', v))
  if (any(grepl('^Sex:', v)))             sex    <- sub('^Sex: ', '', v)
}
stopifnot(!is.null(status), !is.null(age))
sample_id <- sub(' \\[.*\\]', '', titles)

## 2. 表达矩阵 ------------------------------------------------------------------
lcpm <- read.csv(gzfile('../external/GSE111016_lcpm_sss.csv.gz'),
                 row.names = 1, check.names = FALSE)
lcpm <- lcpm[, sample_id]
cat('lcpm 维度:', dim(lcpm)[1], 'x', dim(lcpm)[2], '\n')

## 3. ENSG -> SYMBOL (同 1-External_ROC.R) ---------------------------------------
ens <- gsub('\\..*$', '', rownames(lcpm))
map <- AnnotationDbi::select(org.Hs.eg.db, keys = ens, keytype = 'ENSEMBL',
                             columns = 'SYMBOL')
map <- map[!duplicated(map$ENSEMBL) & !is.na(map$SYMBOL), ]
sym <- map$SYMBOL[match(ens, map$ENSEMBL)]
keep <- !is.na(sym)
lcpm <- lcpm[keep, ]; sym <- sym[keep]
if (anyDuplicated(sym)) {
  ord <- order(-rowMeans(lcpm, na.rm = TRUE)); lcpm <- lcpm[ord, ]; sym <- sym[ord]
  first <- !duplicated(sym); lcpm <- lcpm[first, ]; sym <- sym[first]
}
rownames(lcpm) <- sym
lcpm <- as.matrix(lcpm)
cat('注释后基因数:', nrow(lcpm), '\n')
cat('四 hub 基因检出:', paste(HUB, HUB %in% rownames(lcpm)), '\n')

## 4. 年龄平衡性 -----------------------------------------------------------------
cat('\n年龄分布: SP =', sprintf('%.1f±%.1f', mean(age[status=='SP']), sd(age[status=='SP'])),
    '| Control =', sprintf('%.1f±%.1f', mean(age[status=='Control']), sd(age[status=='Control'])),
    '| Wilcox P =', signif(wilcox.test(age ~ status)$p.value, 3), '\n')

## 5. 三种模型的逐基因结果 --------------------------------------------------------
st  <- factor(status, levels = c('Control', 'SP'))
age_c <- age - mean(age)
run_genes <- function(e, st_, label, form = c('status', 'status+age')) {
  do.call(rbind, lapply(HUB, function(g) {
    if (!g %in% rownames(e)) return(data.frame(model = label, gene = g,
                                               logFC = NA, P = NA, dir = 'not_detected',
                                               AUC = NA, AUC_low = NA, AUC_high = NA,
                                               dir_matches_training = NA))
    x <- as.numeric(e[g, ])
    if (length(form) == 1 && form == 'status') {
      d <- model.matrix(~ st_); cf <- 'st_SP'
    } else {
      d <- model.matrix(~ age_c + st_); cf <- 'st_SP'
    }
    fit <- lmFit(e[g, , drop = FALSE], d)
    fit <- eBayes(fit)
    tt <- topTable(fit, coef = cf, number = 1)
    r <- roc(st_ ~ x, quiet = TRUE, ci = TRUE)
    data.frame(model = label, gene = g,
               logFC = round(tt$logFC, 3), P = signif(tt$P.Value, 3),
               dir = ifelse(tt$logFC > 0, 'up', 'down'),
               AUC = round(as.numeric(r$auc), 3),
               AUC_low = round(as.numeric(r$ci)[1], 3),
               AUC_high = round(as.numeric(r$ci)[3], 3),
               dir_matches_training =
                 ifelse(tt$logFC > 0,
                        sp_tab[g, 'logFC'] > 0, sp_tab[g, 'logFC'] < 0))
  }))
}
res_full  <- run_genes(lcpm, st, 'full_cohort (n=40)', 'status')
res_adj   <- run_genes(lcpm, st, 'age_adjusted (n=40)', 'status+age')

## 1:1 年龄最近邻匹配 (greedy, caliper 3y, 无放回)
set.seed(42)
sp_i <- which(status == 'SP'); ct_i <- which(status == 'Control')
pair <- list(); used <- c()
for (i in sp_i) {
  cand <- ct_i[!(ct_i %in% used)]
  cand <- cand[abs(age[cand] - age[i]) <= 3]
  if (!length(cand)) next
  j <- cand[which.min(abs(age[cand] - age[i]))]
  pair <- c(pair, list(c(i, j))); used <- c(used, j)
}
m_idx <- unlist(pair)
cat('\n1:1 年龄匹配: ', length(pair), '对, |Δage| 中位数 =',
    median(sapply(pair, function(p) abs(age[p[1]] - age[p[2]]))), '岁\n')
st_m <- droplevels(st[m_idx])
res_matched <- run_genes(lcpm[, m_idx], st_m, 'age_matched_1to1', 'status')
if (length(pair)) cat('匹配对年龄差:', paste(sapply(pair, function(p) abs(age[p[1]]-age[p[2]])), collapse=', '), '\n')

hub_tab <- rbind(res_full, res_adj, res_matched)
cat('\n===== 四 hub 基因敏感性结果 =====\n'); print(hub_tab, row.names = FALSE)
write.csv(hub_tab, file.path(outd, 'GSE111016_age_matched_hub_genes.csv'), row.names = FALSE)

## 6. 118 common DEGs 在 "sarcopenic vs non-sarcopenic (老年组内)" 的存活率 -------
design <- model.matrix(~ age_c + st)
fit <- lmFit(lcpm, design); fit <- eBayes(fit)
coef_nm <- colnames(design)[ncol(design)]
tt_all <- topTable(fit, coef = coef_nm, n = Inf) |> na.omit()
tt_all$P.Value[tt_all$P.Value == 0] <- min(tt_all$P.Value[tt_all$P.Value != 0])
degs <- rownames(tt_all)[tt_all$P.Value < 0.05]
pool <- read.csv('../2-Venn/output/SP-HP.CSV', row.names = 1, check.names = FALSE)
pool118 <- as.character(pool$x)
pool_in <- intersect(pool118, rownames(tt_all))
ov <- intersect(pool118, degs)
U <- nrow(tt_all)
p_hyper <- phyper(length(ov) - 1, length(pool_in), U - length(pool_in), length(degs),
                  lower.tail = FALSE)
summary_row <- data.frame(
  contrast = 'GSE111016 sarcopenic vs non-sarcopenic (age-adjusted)',
  universe = U, n_DEGs = length(degs),
  pool118_tested = length(pool_in), overlap = length(ov),
  expected = round(length(pool_in) * length(degs) / U, 1),
  fold_enrichment = round(length(ov) / (length(pool_in) * length(degs) / U), 2),
  hypergeometric_P = signif(p_hyper, 3))
cat('\n===== 118 池在表型锚定对比中的存活 =====\n'); print(summary_row)
## 方向: 与训练 (young vs old) 方向的一致率
d_now <- ifelse(tt_all[pool_in, 'logFC'] > 0, 'up', 'down')
d_tr  <- ifelse(sp_tab[pool_in, 'logFC'] > 0, 'up', 'down')
cat('与训练方向一致率:', sum(d_now == d_tr), '/', length(pool_in), '=',
    round(mean(d_now == d_tr), 2), '| binomial P =',
    signif(binom.test(sum(d_now == d_tr), length(pool_in), 0.5)$p.value, 3), '\n')

## 7. 保存 -----------------------------------------------------------------------
write.csv(summary_row, file.path(outd, 'GSE111016_pool118_survival.csv'), row.names = FALSE)
write.csv(data.frame(gene = rownames(tt_all), tt_all, check.names = FALSE),
          file.path(outd, 'GSE111016_age_adjusted_all_genes.csv'), row.names = FALSE)
saveRDS(list(hub = hub_tab, pool_summary = summary_row,
             age_sp = list(mean = mean(age[status=='SP']), sd = sd(age[status=='SP'])),
             age_ct = list(mean = mean(age[status=='Control']), sd = sd(age[status=='Control']))),
        file.path(outd, 'results_age_matched.Rds'))
cat('\nDONE\n')
