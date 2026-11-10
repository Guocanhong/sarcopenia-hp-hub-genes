# 03b: 探测 GWAS Catalog FTP 目录结构
suppressMessages({library(httr)})
get1 <- function(url) {
  r <- tryCatch(GET(url, timeout(60)), error = function(e) e)
  if (inherits(r, "error")) return(list(err = conditionMessage(r)))
  list(status = status_code(r), text = content(r, "text", encoding = "UTF-8"))
}
r1 <- get1("https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/")
if (!is.null(r1$err)) {
  cat("latest ERR:", r1$err, "\n")
} else {
  cat("latest HTTP", r1$status, "\n")
  links <- regmatches(r1$text, gregexpr('href="[^"]+"', r1$text))[[1]]
  links <- gsub('href="|"', "", links)
  cat("条目数:", length(links), "\n")
  print(head(links, 25))
}
# 也试 file-downloads 页提到的固定路径
r2 <- get1("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/")
cat("\nsummary_statistics/ HTTP:", if (is.null(r2$err)) r2$status else r2$err, "\n")
