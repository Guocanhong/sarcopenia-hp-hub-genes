# =====================================================================
# 第三疾病特异性对照分析 (Third-disease specificity control)
# 目的: 检验 118 common DEGs / 四 hub 基因的重叠是否为 "任意慢性炎症" 的通用产物
# 对照线 (均 GPL570, 注释与主流程一致):
#   A. GSE13355  银屑病受累皮肤 (PP, n=58) vs 正常皮肤 (n=64)   [组织层面, 主对照]
#      (敏感性: PP vs 未受累皮肤 PN, n=58, 排除人群背景差异)
#   B. GSE50772  SLE PBMC (n=61) vs Control PBMC (n=20)         [系统层面]
# 流程镜像主线: 注释映射 -> 同名基因合并(取表达最高探针) -> log2 ->
#               低表达过滤(>0 in >50% 样本) -> limma -> 名义 P<0.05 (无FC截断)
# 输出: output/ 下 DEG 全表 / 重叠明细 / 汇总 CSV + 控制台摘要
# =====================================================================
rm(list = ls())
suppressPackageStartupMessages({library(limma)})

outd <- 'output'; dir.create(outd, showWarnings = FALSE, recursive = TRUE)

## 0. 公共输入 -----------------------------------------------------------------
common <- read.csv('../2-Venn/output/SP-HP.CSV', row.names = 1, check.names = FALSE)
common118 <- as.character(common$x)   # 118 common DEGs (SP ∩ HP)
stopifnot(length(common118) == 118)
HUB <- c('C1QA', 'NFKBIA', 'DTL', 'CTNNBIP1')

## 主线两线的 full limma 表 (用于方向一致性判定)
sp_tab <- readRDS('../1-limma/output/diff_res.Rds')
hp_tab <- readRDS('../1-limma2/output/HP-diff_res.Rds')
dir_of <- function(tab, genes) {
  ifelse(tab[genes, 'logFC'] > 0, 'up', 'down')
}
sp_dir <- setNames(dir_of(sp_tab, common118), common118)
hp_dir <- setNames(dir_of(hp_tab, common118), common118)

## 1. GPL570 注释映射 (与 8-SP2/1-External_ROC_GSE25941.R 完全一致) ------------
ann_lines <- readLines(gzfile('../external/GPL570.annot.gz'))
ann_lines <- ann_lines[grepl('^[0-9]', ann_lines)]
ann <- read.table(text = paste(ann_lines, collapse = '\n'),
                  sep = '\t', header = FALSE, quote = '', fill = TRUE,
                  comment.char = '', check.names = FALSE)
sym_map <- ann[, c(1, 3)]; colnames(sym_map) <- c('probe', 'symbol')
sym_map$symbol <- sub('///.*', '', sym_map$symbol)
sym_map <- sym_map[sym_map$symbol != '' & !is.na(sym_map$symbol), ]
sym_map <- sym_map[!duplicated(sym_map$probe), ]

## 2. 单数据集处理函数 ----------------------------------------------------------
## log_transform=FALSE, low_filter=FALSE: 用于沉积值已是 z-score/对数中心化数据的数据集
##   (GSE13355 沉积值均值=0, 为逐样本标准化数据, 再取 log2 或按 ">0" 过滤均无意义)
run_line <- function(file, group_rule, label, test, ref,
                     log_transform = TRUE, low_filter = TRUE) {
  sm <- readLines(gzfile(file))
  get_line <- function(tag) {
    ln <- grep(paste0('^!', tag), sm, value = TRUE)[1]
    strsplit(ln, '\t')[[1]][-1] |> (\(x) gsub('"', '', trimws(x)))()
  }
  gsm  <- get_line('Sample_geo_accession')
  chr  <- get_line('Sample_characteristics_ch1')
  group <- group_rule(chr)
  names(group) <- gsm
  use_all <- group %in% c(test, ref)
  cat('\n=====', label, '|', test, 'vs', ref,
      '| n =', sum(group == test), 'vs', sum(group == ref), '\n')

  tab_start <- grep('^!series_matrix_table_begin', sm) + 1
  tab_end   <- grep('^!series_matrix_table_end', sm) - 1
  expr <- read.table(text = paste(sm[tab_start:tab_end], collapse = '\n'),
                     sep = '\t', header = TRUE, quote = '', row.names = 1,
                     check.names = FALSE)
  expr <- as.matrix(expr)
  colnames(expr) <- gsub('"', '', colnames(expr))
  rownames(expr) <- gsub('"', '', rownames(expr))
  stopifnot(all(colnames(expr) == gsm))
  cat('探针维度:', dim(expr)[1], 'x', dim(expr)[2], '\n')

  ## 注释映射 + 同名合并 (同主流程)
  hit  <- match(rownames(expr), sym_map$probe)
  sym  <- sym_map$symbol[hit]
  keep <- !is.na(sym)
  expr <- expr[keep, ]; sym <- sym[keep]
  ord  <- order(sym, -rowMeans(expr))
  expr <- expr[ord, ]; sym <- sym[ord]
  first <- !duplicated(sym)
  expr <- expr[first, ]; rownames(expr) <- sym[first]
  cat('注释到符号并合并后基因数:', nrow(expr), '\n')

  ## log2 + 低表达过滤 (镜像主线; 对 z-score 型沉积数据两步均跳过)
  el <- expr
  if (log_transform) el <- log2(el + 1)
  if (low_filter) {
    el <- el[apply(el, 1, function(x) sum(x > 0) > ncol(el) * 0.50), ]
  }
  universe <- rownames(el)
  cat('log_transform =', log_transform, '| low_filter =', low_filter,
      '| universe:', length(universe), '\n')

  e    <- el[, names(group)[use_all]]
  g    <- droplevels(factor(group[use_all], levels = c(ref, test)))
  design <- model.matrix(~ 0 + g); colnames(design) <- levels(g)
  contrast <- paste0(test, '-', ref)
  fit <- lmFit(e, design)
  fit <- contrasts.fit(fit, makeContrasts(contrasts = contrast, levels = design))
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = contrast, n = Inf) |> na.omit()
  tt$P.Value[tt$P.Value == 0] <- min(tt$P.Value[tt$P.Value != 0])
  tt$change <- ifelse(tt$P.Value < 0.05 & tt$logFC > 0, 'Up',
               ifelse(tt$P.Value < 0.05 & tt$logFC < 0, 'Down', 'Not'))
  degs <- rownames(tt)[tt$change != 'Not']
  cat('DEGs (名义 P<0.05):', length(degs),
      '| Up:', sum(tt$change == 'Up'), '| Down:', sum(tt$change == 'Down'), '\n')
  write.csv(data.frame(gene = rownames(tt), tt, check.names = FALSE),
            file.path(outd, paste0(label, '_', contrast, '_all.csv')),
            row.names = FALSE)

  ## 与 118 common DEGs 的重叠 + 超几何检验
  pool_in_u <- intersect(common118, universe)
  ov <- intersect(common118, degs)
  m  <- length(pool_in_u)
  k  <- length(degs)
  q  <- length(ov)
  U  <- length(universe)
  p_hyper <- phyper(q - 1, m, U - m, k, lower.tail = FALSE)
  expected <- m * k / U

  ## 方向一致性: 重叠基因在第三疾病中与 SP / HP 主线同方向的比例
  if (q > 0) {
    d3 <- ifelse(tt[ov, 'logFC'] > 0, 'up', 'down')
    same_sp <- mean(d3 == sp_dir[ov])
    same_hp <- mean(d3 == hp_dir[ov])
    p_sp <- binom.test(sum(d3 == sp_dir[ov]), q, 0.5)$p.value
    p_hp <- binom.test(sum(d3 == hp_dir[ov]), q, 0.5)$p.value
    dir_tab <- data.frame(
      gene = ov,
      dir_SP = unname(sp_dir[ov]), dir_HP = unname(hp_dir[ov]),
      logFC_third = round(tt[ov, 'logFC'], 3),
      P_third = signif(tt[ov, 'P.Value'], 3), dir_third = d3,
      same_as_SP = d3 == sp_dir[ov], same_as_HP = d3 == hp_dir[ov]
    )
    dir_tab <- dir_tab[order(-abs(dir_tab$logFC_third)), ]
    write.csv(dir_tab,
              file.path(outd, paste0(label, '_', contrast, '_overlap_genes.csv')),
              row.names = FALSE)
    cat('方向一致率 (vs SP):', sum(d3 == sp_dir[ov]), '/', q, '=',
        round(same_sp, 2), '| binomial P =', signif(p_sp, 3), '\n')
    cat('方向一致率 (vs HP):', sum(d3 == hp_dir[ov]), '/', q, '=',
        round(same_hp, 2), '| binomial P =', signif(p_hp, 3), '\n')
  } else {
    same_sp <- same_hp <- p_sp <- p_hp <- NA
  }

  summary_row <- data.frame(
    dataset = label, contrast = contrast,
    n_test = sum(group == test, na.rm = TRUE),
    n_ref  = sum(group == ref,  na.rm = TRUE),
    universe = U, n_DEGs = k,
    pool118_in_universe = m, overlap = q,
    expected_overlap = round(expected, 1),
    fold_enrichment = round(q / expected, 2),
    hypergeometric_P = signif(p_hyper, 3),
    concord_vs_SP = if (q > 0) paste0(round(same_sp, 2), ' (P=', signif(p_sp, 2), ')') else NA,
    concord_vs_HP = if (q > 0) paste0(round(same_hp, 2), ' (P=', signif(p_hp, 2), ')') else NA
  )
  print(summary_row)
  if (q > 0) cat('重叠基因:', paste(head(ov, 25), collapse = ', '),
                 if (q > 25) ' ...' else '', '\n')

  ## Hub 基因在第三疾病中的表现
  hub_tab <- do.call(rbind, lapply(intersect(HUB, rownames(tt)), function(gn) {
    data.frame(dataset = label, contrast = contrast, gene = gn,
               logFC = round(tt[gn, 'logFC'], 3),
               P = signif(tt[gn, 'P.Value'], 3),
               in_DEG_list = gn %in% degs,
               dir = ifelse(tt[gn, 'logFC'] > 0, 'up', 'down'))
  }))
  print(hub_tab)
  list(summary = summary_row, hub = hub_tab)
}

## 3. A 线: GSE13355 银屑病 ----------------------------------------------------
res <- list()
res$A1 <- run_line('input/GSE13355_series_matrix.txt.gz',
                   function(chr) {
                     g <- rep('Other', length(chr))
                     g[grepl('involved skin from cases', chr) &
                       !grepl('uninvolved', chr)] <- 'PP'
                     g[grepl('normal skin from controls', chr)] <- 'Control'
                     g
                   }, 'GSE13355_Psoriasis_skin', test = 'PP', ref = 'Control',
                   log_transform = FALSE, low_filter = FALSE)
## 敏感性: PP vs PN (同为银屑病患者, 排除人群背景)
res$A2 <- run_line('input/GSE13355_series_matrix.txt.gz',
                   function(chr) {
                     g <- rep('Other', length(chr))
                     g[grepl('involved skin from cases', chr) &
                       !grepl('uninvolved', chr)] <- 'PP'
                     g[grepl('uninvolved skin from cases', chr)] <- 'PN'
                     g
                   }, 'GSE13355_Psoriasis_skin', test = 'PP', ref = 'PN',
                   log_transform = FALSE, low_filter = FALSE)

## 4. B 线: GSE50772 SLE PBMC --------------------------------------------------
res$B <- run_line('input/GSE50772_series_matrix.txt.gz',
                  function(chr) ifelse(grepl('SLE$', chr), 'SLE',
                                ifelse(grepl('Control$', chr), 'Control', 'Other')),
                  'GSE50772_SLE_PBMC', test = 'SLE', ref = 'Control')

## 5. 汇总 ---------------------------------------------------------------------
ov_summary <- do.call(rbind, lapply(res, `[[`, 'summary'))
hub_all    <- do.call(rbind, lapply(res, `[[`, 'hub'))
write.csv(ov_summary, file.path(outd, 'Third_disease_overlap_summary.csv'),
          row.names = FALSE)
write.csv(hub_all, file.path(outd, 'Third_disease_hub_genes.csv'), row.names = FALSE)
cat('\n===== 总汇总 =====\n'); print(ov_summary)
cat('\n===== Hub 基因在第三疾病中的表现 =====\n'); print(hub_all)
saveRDS(list(summary = ov_summary, hub = hub_all), file.path(outd, 'results.Rds'))
cat('\nDONE\n')
