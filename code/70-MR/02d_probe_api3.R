# 02d: 只对查询参数做 URL 编码
suppressMessages({library(httr); library(jsonlite)})
probe <- function(url) {
  r <- tryCatch(GET(url, timeout(45), add_headers(Accept = "application/json")), error = function(e) e)
  if (inherits(r, "error")) return(paste("ERR:", conditionMessage(r)))
  txt <- substr(content(r, "text", encoding = "UTF-8"), 1, 150)
  paste0("HTTP ", status_code(r), " | ", gsub("[\r\n]", " ", txt))
}
enc_param <- function(x) URLencode(x, reserved = TRUE)
cat("ss_trait  :", probe(paste0("https://www.ebi.ac.uk/gwas/rest/api/studies/search/summaryStatistics?trait=", enc_param("hand grip strength"), "&size=3")), "\n")
cat("ss_dt     :", probe(paste0("https://www.ebi.ac.uk/gwas/rest/api/studies/search/summaryStatistics?diseaseTrait=", enc_param("hand grip strength"), "&size=3")), "\n")
cat("ssapi_q   :", probe(paste0("https://www.ebi.ac.uk/gwas/summary-statistics/api/studies?q=", enc_param("hand grip strength"), "&size=3")), "\n")
cat("studies_q :", probe(paste0("https://www.ebi.ac.uk/gwas/rest/api/studies/search?q=", enc_param("hand grip strength"), "&size=3")), "\n")
