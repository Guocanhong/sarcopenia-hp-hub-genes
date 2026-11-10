# 03: 下载 GWAS Catalog 最新 flat files，筛选含汇总统计的目标性状研究
suppressMessages({library(data.table)})
url <- "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/studies.tsv"
dest <- "data/gwascat_studies.tsv"
if (!file.exists(dest)) {
  cat("下载 studies.tsv ...\n")
  r <- tryCatch(download.file(url, dest, mode = "wb", method = "wininet"), error = function(e) e)
  if (inherits(r, "error")) {
    cat("wininet 失败:", conditionMessage(r), "改试 libcurl\n")
    r <- tryCatch(download.file(url, dest, mode = "wb", method = "libcurl"), error = function(e) e)
    if (inherits(r, "error")) stop(conditionMessage(r))
  }
  cat("完成\n")
}
cat("文件大小:", file.size(dest) / 1e6, "MB\n")
st <- fread(dest, sep = "\t", quote = "", showProgress = FALSE, encoding = "UTF-8")
cat("列名:", paste(names(st), collapse = " | "), "\n")
setnames(st, make.names(names(st)))
# 找汇总统计链接列
sscol <- grep("SUMMARY|SUMMARY_STATS|STATS", names(st), value = TRUE, ignore.case = TRUE)
cat("SS 列:", if (length(sscol)) paste(sscol, collapse = ", ") else "(无)", "\n")
n_ss <- if (length(sscol)) sum(!is.na(st[[sscol[1]]]) & st[[sscol[1]]] != "") else 0
cat("有 SS 链接的研究数:", n_ss, "\n")
fwrite(st, "data/gwascat_studies_clean.tsv", sep = "\t")
