Sys.setenv(OPENGWAS_JWT = trimws(readLines("data/token.txt", warn = FALSE)[1]))
suppressMessages(library(ieugwasr))
gi <- tryCatch(gwasinfo("ieu-b-4869"), error = function(e) e)
if (inherits(gi, "error")) {
  cat("FAIL:", conditionMessage(gi), "\n")
} else {
  cat("TOKEN_OK | trait:", gi$trait[1], "| n:", gi$n[1], "| year:", gi$year[1], "\n")
}
