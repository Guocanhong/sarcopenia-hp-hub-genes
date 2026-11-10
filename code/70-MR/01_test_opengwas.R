# 01: 测试 OpenGWAS 匿名访问 + 搜索 H. pylori / 肌少症相关 GWAS 数据集
suppressMessages({library(httr); library(jsonlite)})

cat("=== 1. API 状态 ===\n")
r <- tryCatch(GET("https://api.opengwas.io/api/status", timeout(30)), error = function(e) e)
if (inherits(r, "error")) { cat("status FAIL:", conditionMessage(r), "\n") } else {
  cat("HTTP", status_code(r), "\n")
  if (status_code(r) == 200) print(fromJSON(content(r, "text"))$time)
}

cat("\n=== 2. 匿名拉取全部数据集列表,过滤 helicobacter/sarcopenia/lean/grip ===\n")
r2 <- tryCatch(GET("https://api.opengwas.io/api/gwasinfo", timeout(60)), error = function(e) e)
if (inherits(r2, "error")) {
  cat("gwasinfo FAIL:", conditionMessage(r2), "\n")
} else {
  cat("gwasinfo HTTP", status_code(r2), "\n")
  if (status_code(r2) == 200) {
    info <- fromJSON(content(r2, "text", encoding = "UTF-8"))
    # info 是大列表,按 id 转 data.frame
    df <- do.call(rbind, lapply(names(info), function(id) {
      x <- info[[id]]
      data.frame(id = id,
                 trait = if (is.null(x$trait)) "" else as.character(x$trait),
                 n = if (is.null(x$ncase) && is.null(x$n)) "" else as.character(ifelse(is.null(x$n), x$nsnp, x$n)),
                 stringsAsFactors = FALSE)
    }))
    kw <- "helicobacter|sarcopen|lean mass|grip strength|walking pace"
    hits <- df[grepl(kw, df$trait, ignore.case = TRUE), ]
    cat("命中", nrow(hits), "个数据集:\n")
    print(hits, right = FALSE)
    write.csv(hits, "D:/GEO-limma,肌少症1/70-MR/output/gwas_hits_search.csv", row.names = FALSE)
    cat("\n已存 output/gwas_hits_search.csv\n")
  } else {
    cat("非 200,可能需要 token。响应头三行:\n")
    cat(substr(content(r2, "text"), 1, 200), "\n")
  }
}
