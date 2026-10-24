# 距离度量比较：euclidean vs pearson vs spearman（PAC 越低越好）
rm(list = ls())
set.seed(2026)
suppressMessages({
  library(ConsensusClusterPlus)
  library(AnnotationDbi); library(org.Hs.eg.db)
})
lcpm <- read.csv(gzfile('GSE111016_lcpm_sss.csv.gz'), row.names = 1, check.names = FALSE)
sm <- readLines(gzfile('GSE111016_series_matrix.txt.gz'))
get_line <- function(tag) gsub('"', '', strsplit(grep(paste0('^!', tag), sm, value = TRUE)[1], '\t')[[1]][-1])
titles <- get_line('Sample_title')
chars <- sm[grepl('^!Sample_characteristics_ch1', sm)]
status <- NULL
for (ln in chars) {
  vals <- gsub('"', '', strsplit(ln, '\t')[[1]][-1])
  if (any(grepl('sarcopenia status', vals))) { status <- vals; break }
}
group <- ifelse(grepl(': yes', status), 'Sarcopenia', 'Control')
names(group) <- sub(' \\[.*\\]', '', titles)
lcpm <- lcpm[, sub(' \\[.*\\]', '', titles)]

ens <- gsub('\\..*$', '', rownames(lcpm))
map <- AnnotationDbi::select(org.Hs.eg.db, keys = ens, keytype = 'ENSEMBL', columns = 'SYMBOL')
map <- map[!duplicated(map$ENSEMBL) & !is.na(map$SYMBOL), ]
sym <- map$SYMBOL[match(ens, map$ENSEMBL)]
keep <- !is.na(sym); l2 <- lcpm[keep, ]; s2 <- sym[keep]
if (anyDuplicated(s2)) {
  o <- order(-rowMeans(l2)); l2 <- l2[o, ]; s2 <- s2[o]
  f <- !duplicated(s2); l2 <- l2[f, ]; s2 <- s2[f]
}
rownames(l2) <- s2
common <- read.csv('../2-Venn/output/SP-HP.CSV', row.names = 1)$x
hit <- intersect(common, s2)
l2 <- l2[, colSums(is.na(l2)) == 0]
hit <- intersect(common, rownames(l2))
hit <- hit[rowSums(is.na(l2[hit, ])) == 0]
mat <- t(scale(t(l2[hit, ])))
cat('矩阵:', paste(dim(mat), collapse = 'x'), '\n')

for (dist in c('euclidean', 'pearson', 'spearman')) {
  d <- if (dist == 'euclidean') dist else dist
  set.seed(2026)
  ccp <- ConsensusClusterPlus(mat, maxK = 4, reps = 500, pItem = 0.8, pFeature = 0.8,
                              clusterAlg = 'hc', distance = dist,
                              innerLinkage = 'ward.D2', finalLinkage = 'ward.D2',
                              title = paste0('output/ccp_', dist), plot = 'pdf', seed = 2026)
  pac <- sapply(2:4, function(k) {
    cm <- ccp[[k]]$consensusMatrix
    sum(cm > 0.1 & cm < 0.9) / (nrow(cm)^2 - nrow(cm))
  })
  cl2 <- ccp[[2]]$consensusClass
  g39 <- group[colnames(mat)]
  tb <- table(Subtype = cl2, Status = g39)
  ft <- fisher.test(tb)$p.value
  cat('\n====', dist, '====\n')
  print(data.frame(k = 2:4, PAC = round(pac, 4)))
  print(tb)
  cat('Fisher P =', signif(ft, 3), '\n')
}
