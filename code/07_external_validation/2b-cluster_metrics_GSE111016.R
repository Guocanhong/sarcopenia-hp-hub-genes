# =====================================================================
# Consensus clustering selection metrics for GSE111016 (Methods 2.6):
# PAC (0.1-0.9 window), marginal decrease in area under the consensus CDF
# curve (delta area), and mean silhouette width on the consensus-based
# dissimilarity (1 - consensus matrix), for k = 2-6, at the primary
# parameters (reps=1000, pItem=pFeature=0.8, Pearson distance, ward.D2,
# seed 2026). Expected output (verified 2026-10-11, reproduced exactly):
#   n_samples 39 | genes_from_118 102
#   PAC k=2-6:   0.5924 0.6599 0.5655 0.4372 0.3630
#   dK  k=3-6:   0.1607 0.0908 0.0575 0.0387 (largest at k=3)
#   SIL k=2-6:   0.7600 0.4898 0.4136 0.4144 0.4575
# Complements 2-Consensus_subtype_GSE111016.R (which does not itself
# compute the silhouette / delta-area selection metrics).
# =====================================================================
suppressMessages({library(ConsensusClusterPlus); library(AnnotationDbi); library(org.Hs.eg.db); library(cluster)})
set.seed(2026)
lcpm <- read.csv(gzfile('external/GSE111016_lcpm_sss.csv.gz'), row.names = 1, check.names = FALSE)
sm <- readLines(gzfile('external/GSE111016_series_matrix.txt.gz'))
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
lcpm <- lcpm[, names(group)]
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
common <- read.csv('2-Venn/output/SP-HP.CSV', row.names = 1)$x
l2 <- l2[, colSums(is.na(l2)) == 0]
hit <- intersect(common, rownames(l2))
hit <- hit[rowSums(is.na(l2[hit, ])) == 0]
cat('n_samples:', ncol(l2), '| genes_from_118:', length(hit), '\n')
mat <- t(scale(t(l2[hit, ])))
ccp <- ConsensusClusterPlus(mat, maxK = 6, reps = 1000, pItem = 0.8, pFeature = 0.8,
                            clusterAlg = 'hc', distance = 'pearson',
                            innerLinkage = 'ward.D2', finalLinkage = 'ward.D2',
                            title = file.path(tempdir(), 'ccp_verify'), plot = 'pdf', seed = 2026)
kvals <- 2:6
pac <- sapply(kvals, function(k) {
  cm <- ccp[[k]]$consensusMatrix
  sum(cm > 0.1 & cm < 0.9) / (nrow(cm)^2 - nrow(cm))
})
auc <- sapply(kvals, function(k) mean(ccp[[k]]$consensusMatrix[upper.tri(ccp[[k]]$consensusMatrix)]))
dK <- -diff(auc)
sil <- sapply(kvals, function(k) {
  d <- as.dist(1 - ccp[[k]]$consensusMatrix)
  mean(silhouette(ccp[[k]]$consensusClass, d)[, 3])
})
cat('PAC_k2_6:', sprintf('%.4f', pac), '\n')
cat('dK_k3_6:', sprintf('%.4f', dK), '\n')
cat('SIL_k2_6:', sprintf('%.4f', sil), '\n')
