# 2-diag.R  诊断：基因行是否错位、MT 基因是否存在、免疫核真实分布
suppressMessages({library(data.table); library(Matrix)})
src <- '
read_counts_csv <- function(path, chunk_rows = 2000) {
  con <- gzfile(path, "rt"); head.line <- readLines(con, n = 1L); close(con)
  nbarcode <- length(strsplit(head.line, ",", fixed = TRUE)[[1]]) - 1L
  i.list <- list(); j.list <- list(); x.list <- list(); genes.all <- character(0)
  gi <- 0L; sk <- 1L; width <- NA_integer_
  repeat {
    blk <- tryCatch(fread(path, skip = sk, nrows = chunk_rows, header = FALSE,
                          sep = ",", showProgress = FALSE), error = function(e) NULL)
    if (is.null(blk) || nrow(blk) == 0) break
    if (is.na(width)) width <- ncol(blk) - 1L
    if (ncol(blk) - 1L != width) stop("column width mismatch at row ", gi + 1)
    genes.all <- c(genes.all, as.character(blk[[1]]))
    m <- as.matrix(as.data.frame(lapply(blk[, -1], as.numeric)))
    nz <- which(m != 0, arr.ind = TRUE)
    if (nrow(nz)) {
      i.list[[length(i.list) + 1L]] <- nz[, 1] + gi
      j.list[[length(j.list) + 1L]] <- nz[, 2]
      x.list[[length(x.list) + 1L]] <- m[nz]
    }
    gi <- gi + nrow(blk); sk <- sk + chunk_rows
  }
  M <- sparseMatrix(i = unlist(i.list), j = unlist(j.list), x = unlist(x.list), dims = c(gi, width))
  barcodes <- strsplit(head.line, ",", fixed = TRUE)[[1]][-1]
  rownames(M) <- genes.all; colnames(M) <- barcodes
  M
}
'
eval(parse(text = src))
M <- read_counts_csv("raw/GSM5098741_YouHM5_counts.csv.gz")

cat("== 基因总数:", nrow(M), " 细胞数:", ncol(M), "\n")
cat("== MT 基因存在:", sum(grepl("^MT-", rownames(M))), "个 ->", head(grep("^MT-", rownames(M), value=TRUE), 5), "\n")
chk <- c("C1QA","C1QB","C1QC","CD68","LYZ","FCGR3A","MYH7","MYH2","DCN","PECAM1","PTPRC","ACTB","GAPDH")
for (g in chk) {
  if (g %in% rownames(M)) {
    r <- M[g, ]
    cat(sprintf("%-8s nonzero_cells=%5d  sum=%9.0f  max=%5.0f\n", g, sum(r>0), sum(r), max(r)))
  } else cat(sprintf("%-8s NOT FOUND\n", g))
}
cat("\n== 成对 Spearman 相关（免疫基因应互高，与 MYH7 应≈0）==\n")
sub <- M[intersect(chk, rownames(M)), ]
cc <- cor(t(as.matrix(sub)), method = "spearman")
print(round(cc, 3))
cat("\n== 行名抽检（第1/100/8000/16555 行）==\n")
print(rownames(M)[c(1, 100, 8000, 16555)])
cat("\n== 值域抽检：整数？==\n")
print(table(M[sample(which(M@x>0), 20)]))
saveRDS(M, "output/diag_M_young.rds")
cat("DIAG DONE\n")
