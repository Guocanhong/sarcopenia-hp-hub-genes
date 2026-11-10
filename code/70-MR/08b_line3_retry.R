# 08b_line3_retry.R -- Line 3 retry: grip/ALM exposures with CHUNKED outcome extraction
# Previous run failed because extract_outcome_data sends all SNPs in one request;
# large batches (176/687 SNPs) time out while 57 (walking) succeeded. Here we chunk.
source(file.path(getwd(), "code", "04_config.R"))
suppressMessages({library(TwoSampleMR); library(data.table)})
stopifnot("OPENGWAS_JWT missing" = nzchar(Sys.getenv("OPENGWAS_JWT")))

id_hp  <- c(selfreport = DATASETS$exposure_hp_selfreport,
            igg_pos    = DATASETS$exposure_hp_igg_pos)
id_exp <- c(grip_right = DATASETS$outcome_grip_right,
            alm        = DATASETS$outcome_alm)
CHUNK  <- 40   # walking worked with 57; use 40 to be safe
RETRY  <- 3

chunk_extract_outcome <- function(snps, outcome_id) {
  chunks <- split(snps, ceiling(seq_along(snps) / CHUNK))
  parts <- list()
  for (i in seq_along(chunks)) {
    ok <- FALSE
    for (a in seq_len(RETRY)) {
      d <- tryCatch(extract_outcome_data(snps = chunks[[i]], outcomes = outcome_id,
                                         proxies = TRUE),
                    error = function(e) NULL)
      if (!is.null(d) && nrow(d) > 0) { parts[[length(parts) + 1]] <- d; ok <- TRUE; break }
      Sys.sleep(5)
    }
    cat("    chunk", i, "/", length(chunks), if (ok) "ok\n" else "FAILED\n")
  }
  if (!length(parts)) return(NULL)
  rbindlist(parts, fill = TRUE)
}

res_list <- list()
for (en in names(id_exp)) {
  exp_dat <- tryCatch(extract_instruments(id_exp[[en]], p1 = P_EXP, clump = TRUE,
                                          r2 = CLUMP_R2, kb = CLUMP_KB),
                      error = function(e) { cat(en, "IV extract fail:", conditionMessage(e), "\n"); NULL })
  if (is.null(exp_dat) || nrow(exp_dat) < 3) { cat(en, ": IVs <3, skip\n"); next }
  cat(en, "IVs:", nrow(exp_dat), "\n")
  for (hn in names(id_hp)) {
    cat("  outcome:", hn, "\n")
    out_dat <- tryCatch(chunk_extract_outcome(exp_dat$SNP, id_hp[[hn]]), error = function(e) NULL)
    if (is.null(out_dat) || nrow(out_dat) == 0) next
    out_dat <- out_dat[!duplicated(out_dat$SNP), ]
    dat <- tryCatch(harmonise_data(exp_dat, out_dat, action = 2), error = function(e) NULL)
    if (is.null(dat) || nrow(dat) < 3) next
    key <- paste0(en, "_to_", hn)
    fwrite(as.data.frame(dat), file.path(OUT, paste0("line3_harmonised_", key, ".csv")))
    m <- tryCatch(mr(dat, method_list = c("mr_ivw_mre", "mr_ivw_fe",
                                          "mr_egger_regression", "mr_weighted_median")),
                  error = function(e) NULL)
    if (!is.null(m) && nrow(m) > 0) {
      m$exposure_trait <- en; m$outcome_trait <- hn
      res_list[[key]] <- m
      h <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
      if (!is.null(h) && nrow(h) > 0) fwrite(as.data.frame(h), file.path(OUT, paste0("line3_heterogeneity_", key, ".csv")))
      s <- tryCatch(directionality_test(dat), error = function(e) NULL)
      if (!is.null(s) && nrow(s) > 0) fwrite(as.data.frame(s), file.path(OUT, paste0("line3_steiger_", key, ".csv")))
      print(m[, c("method", "nsnp", "b", "se", "pval")])
    }
  }
}
if (length(res_list)) {
  fwrite(rbindlist(res_list, fill = TRUE), file.path(OUT, "line3b_reverse_mr_results.csv"))
  cat("line3b results saved\n")
} else cat("line3b: no results\n")
