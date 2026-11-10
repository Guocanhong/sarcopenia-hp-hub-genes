# 70-MR 环境检查:R 包 + 网络可达性
cat("=== R version ===\n"); print(R.version.string)

cat("\n=== 关键包检查 ===\n")
pkgs <- c("TwoSampleMR", "ieugwasr", "coloc", "MungeSumstats", "httr", "jsonlite", "data.table")
inst <- rownames(installed.packages())
for (p in pkgs) cat(sprintf("%-14s : %s\n", p, if (p %in% inst) as.character(packageVersion(p)) else "MISSING"))

cat("\n=== 网络可达性 ===\n")
if ("httr" %in% inst) {
  urls <- c(
    opengwas_api   = "https://api.opengwas.io/api/status",
    gwas_mrcieu    = "https://gwas.mrcieu.ac.uk/",
    gwas_catalog   = "https://www.ebi.ac.uk/gwas/api/search/status",
    download_ftp   = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
  )
  for (nm in names(urls)) {
    r <- tryCatch(httr::GET(urls[[nm]], httr::timeout(20)),
                  error = function(e) e)
    if (inherits(r, "error")) {
      cat(sprintf("%-14s : FAIL - %s\n", nm, conditionMessage(r)))
    } else {
      cat(sprintf("%-14s : HTTP %d\n", nm, httr::status_code(r)))
    }
  }
} else {
  cat("httr 未安装,网络检查跳过\n")
}

cat("\n=== OpenGWAS token ===\n")
tok <- Sys.getenv("OPENGWAS_JWT", unset = "")
cat(if (nzchar(tok)) "OPENGWAS_JWT 已设置\n" else "OPENGWAS_JWT 未设置\n")
