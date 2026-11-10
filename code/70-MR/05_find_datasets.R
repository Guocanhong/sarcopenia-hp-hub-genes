# 05_find_datasets.R -- search OpenGWAS datasets by keyword, pick ids
source(file.path(getwd(), "code", "04_config.R"))
suppressMessages({library(ieugwasr); library(data.table)})

stopifnot("OPENGWAS_JWT missing" = nzchar(Sys.getenv("OPENGWAS_JWT")))
cat("token loaded, nchar =", nchar(Sys.getenv("OPENGWAS_JWT")), "\n")

res <- tryCatch(gwasinfo(id = NULL), error = function(e) e)
if (inherits(res, "error")) stop("gwasinfo failed: ", conditionMessage(res))
if (!is.null(dim(res))) {
  res <- as.data.frame(res)
} else {
  # in case a named list is returned
  res <- rbindlist(lapply(res, function(x) as.data.frame(t(unlist(x)))), fill = TRUE)
}
if (!"id" %in% names(res)) res$id <- rownames(res)
res$trait_l <- tolower(as.character(res$trait))
cat("total datasets:", nrow(res), "\n")
fwrite(res, file.path(OUT, "opengwas_all_datasets.csv"))

kws <- c("helicobacter", "sarcopeni", "grip strength", "lean mass",
         "walking pace", "walking speed", "gait speed", "frailty",
         "grip", "appendicular", "muscle mass", "walk")
for (k in kws) {
  hit <- res[grepl(k, res$trait_l, fixed = FALSE),
             c("id", "trait", "ncase", "ncontrol", "nsnp", "population", "year")]
  hit <- hit[order(-as.numeric(hit$ncontrol) - ifelse(is.na(hit$ncase), 0, as.numeric(hit$ncase))), ]
  cat("\n---- ", k, " (", nrow(hit), ") ----\n", sep = "")
  if (nrow(hit)) print(head(as.data.frame(hit), 25), right = FALSE)
}
cat("\nfull list saved: output/opengwas_all_datasets.csv\n")
