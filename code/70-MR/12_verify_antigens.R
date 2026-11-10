source(file.path(getwd(), "code", "04_config.R"))
suppressMessages(library(ieugwasr))
ids <- paste0("ebi-a-GCST9000691", 0:6)
gi <- tryCatch(gwasinfo(id = ids), error = function(e) e)
if (inherits(gi, "error")) cat("fail:", conditionMessage(gi), "\n") else {
  d <- as.data.frame(gi)
  cat("columns:", paste(names(d), collapse=", "), "\n\n")
  keep <- intersect(c("id","trait","ncase","ncontrol","nsnp","unit","year","population"), names(d))
  print(d[, keep], right = FALSE)
}
