# =====================================================================
# Four-table data builder (Table 2 & Table 3 raw numbers)
# Table 2: 118 SP-and-HP intersection DEGs with logFC/P per line
# Table 3: 4 hub genes x 6 bulk cohorts (direction, AUC 95%CI, Wilcoxon P)
# Run: cd 17-Tables && Rscript 1-table_data.R   (relative paths only)
# =====================================================================
rm(list = ls())
suppressMessages({
  library(tidyverse)
  library(pROC)
})
dir.create('output', showWarnings = FALSE)
genes4 <- c('C1QA', 'NFKBIA', 'DTL', 'CTNNBIP1')

## ------------------------------------------------------------------ ##
## Table 2: 118 intersection genes                                    ##
## ------------------------------------------------------------------ ##
inter <- as.character(read.csv('../2-Venn/output/SP-HP.CSV', row.names = 1)$x)
cat('intersection genes:', length(inter), '\n')
sp <- readRDS('../1-limma/output/diff_res.Rds')
hp <- readRDS('../1-limma2/output/HP-diff_res.Rds')
stopifnot(all(inter %in% rownames(sp)), all(inter %in% rownames(hp)))

t2 <- data.frame(
  Gene      = inter,
  logFC_SP  = round(sp[inter, 'logFC'], 3),
  P_SP      = signif(sp[inter, 'P.Value'], 3),
  adjP_SP   = signif(sp[inter, 'adj.P.Val'], 3),
  Reg_SP    = sp[inter, 'change'],
  logFC_HP  = round(hp[inter, 'logFC'], 3),
  P_HP      = signif(hp[inter, 'P.Value'], 3),
  adjP_HP   = signif(hp[inter, 'adj.P.Val'], 3),
  Reg_HP    = hp[inter, 'change'],
  stringsAsFactors = FALSE
) %>% arrange(desc(logFC_SP))
write.csv(t2, 'output/Table2_118DEGs_data.csv', row.names = FALSE)
cat('Table 2 written:', nrow(t2), 'rows\n')

## ------------------------------------------------------------------ ##
## Table 3: uniform per-cohort analysis                               ##
## ------------------------------------------------------------------ ##
analyze_cohort <- function(expr, group, ref, test, dataset, comparison) {
  stopifnot(all(group %in% c(ref, test)))
  rows <- data.frame()
  for (g in genes4) {
    if (!g %in% rownames(expr)) {
      rows <- rbind(rows, data.frame(Dataset = dataset, Comparison = comparison,
        Gene = g, Direction = 'Not detected', AUC = NA, CI_low = NA, CI_high = NA,
        P_wilcox = NA, n_ref = sum(group == ref), n_test = sum(group == test)))
      next
    }
    d <- data.frame(group = factor(group, levels = c(ref, test)),
                    x = as.numeric(expr[g, ]))
    d <- d[!is.na(d$x), ]
    wt  <- wilcox.test(x ~ group, data = d, exact = FALSE)
    dir <- ifelse(median(d$x[d$group == test]) > median(d$x[d$group == ref]), 'Up', 'Down')
    r   <- roc(group ~ x, data = d, auc = TRUE, ci = TRUE, quiet = TRUE)
    rows <- rbind(rows, data.frame(
      Dataset = dataset, Comparison = comparison, Gene = g, Direction = dir,
      AUC = round(as.numeric(r$auc), 3),
      CI_low = round(as.numeric(r$ci)[1], 3),
      CI_high = round(as.numeric(r$ci)[3], 3),
      P_wilcox = signif(wt$p.value, 3),
      n_ref = sum(d$group == ref), n_test = sum(d$group == test)))
  }
  rows
}

load_pair <- function(dir, expr_name, group_name) {
  expr  <- readRDS(file.path(dir, expr_name))
  grp   <- readRDS(file.path(dir, group_name))
  expr  <- expr[, grp$sample]
  list(expr = expr, group = as.character(grp$group))
}

t3 <- data.frame()

# 1. GSE1428 SP training
p <- load_pair('../5-ROC/input', 'dat_expr.Rds', 'dat_group.Rds')
cat('GSE1428 groups:', paste(names(table(p$group)), table(p$group)), '\n')
t3 <- rbind(t3, analyze_cohort(p$expr, p$group, 'Control', 'SP',
                               'GSE1428 (training, SP)', 'SP vs Control'))

# 2. GSE27411 HP training
p <- load_pair('../5-ROC2/input', 'dat1_expr.Rds', 'dat1_group.Rds')
cat('GSE27411 groups:', paste(names(table(p$group)), table(p$group)), '\n')
t3 <- rbind(t3, analyze_cohort(p$expr, p$group, 'Control', 'HP',
                               'GSE27411 (training, HP)', 'HP vs Control'))

# 3. GSE5081 HP external
p <- load_pair('../8-HP/input', 'dat_expr.Rds', 'dat_group.Rds')
cat('GSE5081 groups:', paste(names(table(p$group)), table(p$group)), '\n')
t3 <- rbind(t3, analyze_cohort(p$expr, p$group, 'Control', 'HP',
                               'GSE5081 (external, HP)', 'HP vs Control'))

# 4. GSE25941 aging external (labels Control/SP = Young/Old, n=15/21)
p <- load_pair('../8-SP2/input', 'dat_expr.Rds', 'dat_group.Rds')
cat('GSE25941 groups:', paste(names(table(p$group)), table(p$group)), '\n')
stopifnot(sum(p$group == 'Control') == 15, sum(p$group == 'SP') == 21)
t3 <- rbind(t3, analyze_cohort(p$expr, p$group, 'Control', 'SP',
                               'GSE25941 (external, ageing)', 'Old vs Young'))

# 5. GSE111016 RNA-seq external (ENSG -> SYMBOL, DTL not detected)
lcpm <- read.csv(gzfile('../external/GSE111016_lcpm_sss.csv.gz'),
                 row.names = 1, check.names = FALSE)
sm <- readLines(gzfile('../external/GSE111016_series_matrix.txt.gz'))
get_line <- function(tag) {
  ln <- grep(paste0('^!', tag), sm, value = TRUE)[1]
  gsub('"', '', strsplit(ln, '\t')[[1]][-1])
}
titles <- get_line('Sample_title')
chars  <- sm[grepl('^!Sample_characteristics_ch1', sm)]
status <- NULL
for (ln in chars) {
  vals <- gsub('"', '', strsplit(ln, '\t')[[1]][-1])
  if (any(grepl('sarcopenia status', vals))) { status <- vals; break }
}
group111 <- ifelse(grepl(': yes', status), 'SP', 'Control')
lcpm <- lcpm[, sub(' \\[.*\\]', '', titles)]
ens  <- gsub('\\..*$', '', rownames(lcpm))
library(AnnotationDbi); library(org.Hs.eg.db)
map <- AnnotationDbi::select(org.Hs.eg.db, keys = ens, keytype = 'ENSEMBL',
                             columns = 'SYMBOL')
map <- map[!duplicated(map$ENSEMBL) & !is.na(map$SYMBOL), ]
sym <- map$SYMBOL[match(ens, map$ENSEMBL)]
keep <- !is.na(sym)
lcpm <- lcpm[keep, ]; sym <- sym[keep]
if (anyDuplicated(sym)) {
  ord   <- order(-rowMeans(lcpm, na.rm = TRUE))
  lcpm  <- lcpm[ord, ]; sym <- sym[ord]
  first <- !duplicated(sym)
  lcpm  <- lcpm[first, ]; sym <- sym[first]
}
rownames(lcpm) <- sym
cat('GSE111016 groups:', paste(names(table(group111)), table(group111)), '\n')
cat('GSE111016 hub genes found:', paste(genes4[genes4 %in% rownames(lcpm)]), '\n')
t3 <- rbind(t3, analyze_cohort(lcpm, group111, 'Control', 'SP',
                               'GSE111016 (external, SP RNA-seq)', 'SP vs Control'))

# 6. GPL20301 bulk (reuse module-13 outputs)
gc13 <- read.csv('../13-GPL20301/output/group_comparison_GPL20301.csv')
rocO <- read.csv('../13-GPL20301/output/ROC_GPL20301_Older_vs_Young.csv')
rocS <- read.csv('../13-GPL20301/output/ROC_GPL20301_Sarco_vs_Old.csv')
rocA <- read.csv('../13-GPL20301/output/ROC_GPL20301_Sarco_vs_AllControls.csv')
add_gpl <- function(rocdf, comparison) {
  data.frame(Dataset = 'GPL20301 (external, bulk RNA-seq)', Comparison = comparison,
             Gene = rocdf$gene, Direction = NA, AUC = round(rocdf$AUC, 3),
             CI_low = round(rocdf$AUC_low, 3), CI_high = round(rocdf$AUC_high, 3),
             P_wilcox = NA, n_ref = NA, n_test = NA)
}
t3 <- rbind(t3,
  add_gpl(rocO, 'Sarcopenia+Old vs Young'),
  add_gpl(rocS, 'Sarcopenia vs Old'),
  add_gpl(rocA, 'Sarcopenia vs All controls'))
# merge GPL20301 group-comparison P + means into the Older-vs-Young rows
idx <- t3$Dataset == 'GPL20301 (external, bulk RNA-seq)' &
       t3$Comparison == 'Sarcopenia+Old vs Young'
t3$Direction[idx] <- ifelse(
  gc13$mean_Old[match(t3$Gene[idx], gc13$gene)] >
  gc13$mean_Young[match(t3$Gene[idx], gc13$gene)], 'Up', 'Down')
t3$P_wilcox[idx] <- signif(gc13$p_Sarco_vs_Young[match(t3$Gene[idx], gc13$gene)], 3)
t3$n_ref[idx] <- 19; t3$n_test[idx] <- 53

t3 <- t3 %>% arrange(match(Dataset, unique(Dataset)),
                     match(Gene, genes4))
write.csv(t3, 'output/Table3_hubgene_stats_data.csv', row.names = FALSE)
cat('\nTable 3 written:', nrow(t3), 'rows\n')
print(t3, right = FALSE)
cat('\nAll done.\n')
