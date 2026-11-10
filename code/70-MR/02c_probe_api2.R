# 02c: 修正 URL 编码 + 找到正确的检索域
suppressMessages({library(httr); library(jsonlite)})
probe <- function(url) {
  r <- tryCatch(GET(url, timeout(45), add_headers(Accept = "application/json")), error = function(e) e)
  if (inherits(r, "error")) return(paste("ERR:", conditionMessage(r)))
  txt <- substr(content(r, "text", encoding = "UTF-8"), 1, 150)
  paste0("HTTP ", status_code(r), " | ", gsub("[\r\n]", " ", txt))
}
# 1) EBI search 域列表里找 gwas 相关
r <- tryCatch(GET("https://www.ebi.ac.uk/ebisearch/ws/rest", timeout(45)), error = function(e) e)
if (!inherits(r, "error") && status_code(r) == 200) {
  txt <- content(r, "text", encoding = "UTF-8")
  doms <- regmatches(txt, gregexpr('domainId="[^"]+"', txt))[[1]]
  cat("EBI 检索域(含 gwas):\n"); print(grep("gwas", doms, value = TRUE))
}
# 2) GWAS Catalog 各端点(全部 URLencode)
enc <- function(u) URLencode(u, reserved = TRUE)
cat("\nCatalog 端点探测:\n")
cat("ss_trait  :", probe(enc("https://www.ebi.ac.uk/gwas/rest/api/studies/search/summaryStatistics?trait=hand grip strength&size=3")), "\n")
cat("ss_dt     :", probe(enc("https://www.ebi.ac.uk/gwas/rest/api/studies/search/summaryStatistics?diseaseTrait=hand grip strength&size=3")), "\n")
cat("ssapi_q   :", probe(enc("https://www.ebi.ac.uk/gwas/summary-statistics/api/studies?q=hand grip strength&size=3")), "\n")
cat("studies_q :", probe(enc("https://www.ebi.ac.uk/gwas/rest/api/studies/search?q=hand grip strength&size=3")), "\n")
