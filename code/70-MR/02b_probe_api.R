# 02b: 探测 GWAS Catalog summary statistics API 的可用端点/参数
suppressMessages({library(httr); library(jsonlite)})
probe <- function(url) {
  r <- tryCatch(GET(url, timeout(45), add_headers(Accept = "application/json")), error = function(e) e)
  if (inherits(r, "error")) return(paste("ERR:", conditionMessage(r)))
  txt <- substr(content(r, "text", encoding = "UTF-8"), 1, 120)
  paste0("HTTP ", status_code(r), " | ", txt)
}
base <- "https://www.ebi.ac.uk/gwas"
urls <- c(
  ss_trait   = paste0(base, "/rest/api/studies/search/summaryStatistics?trait=hand grip strength"),
  ss_dt      = paste0(base, "/rest/api/studies/search/summaryStatistics?diseaseTrait=hand grip strength"),
  ss_q       = paste0(base, "/rest/api/studies/search/summaryStatistics?q=hand grip strength&size=3"),
  studies_dt = paste0(base, "/rest/api/studies/search?diseaseTrait=hand grip strength&size=3"),
  studies_efo= paste0(base, "/rest/api/studies/search?efoTrait=hand grip strength&size=3"),
  ssapi      = "https://www.ebi.ac.uk/gwas/summary-statistics/api/studies?trait=hand grip strength&size=3",
  ssapi_q    = "https://www.ebi.ac.uk/gwas/summary-statistics/api/studies?q=hand grip strength&size=3",
  ebi_search = "https://www.ebi.ac.uk/ebisearch/ws/rest/gwas_studies?query=hand%20grip%20strength&size=3&format=json"
)
for (nm in names(urls)) {
  cat(sprintf("%-11s : %s\n", nm, probe(urls[[nm]])))
}
