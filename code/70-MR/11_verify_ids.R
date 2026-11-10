source(file.path(getwd(), "code", "04_config.R"))
suppressMessages({library(ieugwasr); library(data.table)})
ids <- c("ukb-b-4903","ukb-b-10215","ukb-b-7478","ukb-b-4711",
         "ebi-a-GCST90000025","ebi-a-GCST90007526","ukb-d-23035_0_0",
         "ebi-a-GCST006674","ebi-a-GCST90018125")
gi <- tryCatch(gwasinfo(id = ids), error = function(e) e)
if (inherits(gi, "error")) { cat("batch fail:", conditionMessage(gi), "\n") } else {
  print(as.data.frame(gi)[, c("id","trait","ncase","ncontrol","nsnp","year")], right = FALSE)
}
