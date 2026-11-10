# 15_ukb_field_ids.R -- fetch gwasinfo metadata (notes contain UKB Field IDs)
source(file.path(getwd(), "code", "04_config.R"))
suppressMessages({library(ieugwasr); library(data.table)})
stopifnot("OPENGWAS_JWT missing" = nzchar(Sys.getenv("OPENGWAS_JWT")))

ids <- c("ukb-b-531", "ukb-b-10215", "ukb-b-7478", "ukb-b-4711",
         "ieu-b-4905", "ebi-a-GCST90000025", "ebi-a-GCST90007526",
         "ebi-a-GCST90006910")
info <- gwasinfo(ids)
for (i in ids) {
  row <- info[info$id == i, ]
  if (nrow(row) == 0) { cat("\n===", i, ": NOT FOUND ===\n"); next }
  cat("\n===", i, "===\n")
  for (f in c("trait", "ncase", "ncontrol", "sample_size", "population", "unit",
              "author", "year", "pmid", "note")) {
    v <- row[[f]]
    if (!is.null(v) && !is.na(v) && nzchar(as.character(v))) cat(" ", f, ":", as.character(v), "\n")
  }
}
