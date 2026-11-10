Sys.setenv(OPENGWAS_JWT = trimws(readLines("data/token.txt", warn = FALSE)[1]))
suppressMessages(library(ieugwasr))
ok <- tryCatch({
  gi <- gwasinfo("ieu-b-4869")
  cat("OK: gwasinfo returned ->", substr(names(gi)[1], 1, 40), "\n")
  cat("trait:", gi$trait[1], "| n:", gi$n[1], "\n")
  TRUE
}, error = function(e) { cat("FAIL:", conditionMessage(e), "\n"); FALSE })
if (ok) {
  # 顺手搜一轮肌少症相关数据集，验证搜索权限
  res <- tryCatch(gwasinfo(), error = function(e) NULL)
  if (!is.null(res)) {
    cat("total datasets visible:", nrow(res), "\n")
    kw <- "grip|sarcopen|lean mass|appendicular|walking speed|gait|frailt|helicopter|helico"
    hit <- res[grepl(kw, res$trait, ignore.case = TRUE), c("id","trait","ncase","ncontrol","year")]
    hit <- hit[order(-ifelse(is.na(hit$ncase),0,hit$ncase)), ]
    cat("relevant hits:", nrow(hit), "\n")
    print(head(hit, 40))
  }
}
