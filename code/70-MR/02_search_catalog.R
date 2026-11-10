# 02: 通过 GWAS Catalog REST API 检索 Hp / 肌少症性状 GWAS 汇总数据
suppressMessages({library(httr); library(jsonlite)})

base <- "https://www.ebi.ac.uk/gwas/rest/api"
search_terms <- c("Helicobacter pylori", "hand grip strength", "appendicular lean mass",
                  "lean mass", "walking pace", "walking speed", "frailty", "sarcopenia",
                  "grip strength")

out <- data.frame()
for (q in search_terms) {
  cat("==== 检索:", q, "====\n")
  url <- sprintf("%s/studies/search/summaryStatistics?q=%s&size=50", base, URLencode(q, reserved = TRUE))
  r <- tryCatch(GET(url, timeout(60), add_headers(Accept = "application/json")), error = function(e) e)
  if (inherits(r, "error")) { cat("  FAIL:", conditionMessage(r), "\n"); next }
  if (status_code(r) != 200) { cat("  HTTP", status_code(r), "\n"); next }
  d <- fromJSON(content(r, "text", encoding = "UTF-8"), simplifyVector = FALSE)
  studies <- d$`_embedded`$studies
  if (is.null(studies)) { cat("  0 hits\n"); next }
  for (s in studies) {
    acc <- s$accessionId
    trait <- if (!is.null(s$diseaseTrait$trait)) s$diseaseTrait$trait else ""
    pub <- if (!is.null(s$publicationInfo$title)) substr(s$publicationInfo$title, 1, 60) else ""
    nss <- ""
    if (!is.null(s$initialSampleSize)) {
      nss <- if (!is.null(s$initialSampleSize$sampleDescription)) substr(s$initialSampleSize$sampleDescription, 1, 90) else ""
    }
    sslink <- s$`_links`$summaryStatisticsCollection$href
    if (is.null(sslink)) sslink <- s$`_links`$summaryStatistics$href
    if (is.null(sslink)) sslink <- ""
    out <- rbind(out, data.frame(query = q, accession = acc, trait = trait, sample = nss,
                                 pubtitle = pub, ss_href = ifelse(is.null(sslink), "", sslink),
                                 stringsAsFactors = FALSE))
    cat(sprintf("  %-14s %-45s %s\n", acc, substr(trait, 1, 45), substr(nss, 1, 60)))
  }
  Sys.sleep(1)
}
dir.create("output", showWarnings = FALSE)
write.csv(unique(out), "output/gwas_catalog_search.csv", row.names = FALSE)
cat("\n共", nrow(unique(out)), "条 -> output/gwas_catalog_search.csv\n")
