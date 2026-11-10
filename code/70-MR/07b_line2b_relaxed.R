# 07b_line2b_relaxed.R -- Line 2 retry: C1QA/DTL/CTNNBIP1 with relaxed cis-eQTL threshold
# Rationale: eQTLGen cis-window datasets (~18k SNPs +/-1Mb); genome-wide 5e-8 leaves
# only 2 instruments for these genes. Standard cis-eQTL MR practice relaxes to 5e-6
# within the cis window (documented in Methods).
source(file.path(getwd(), "code", "04_config.R"))
suppressMessages({library(TwoSampleMR); library(data.table)})
stopifnot("OPENGWAS_JWT missing" = nzchar(Sys.getenv("OPENGWAS_JWT")))

genes <- c(C1QA     = "eqtl-a-ENSG00000173372",
           DTL      = "eqtl-a-ENSG00000143476",
           CTNNBIP1 = "eqtl-a-ENSG00000178585")

outcomes <- c(grip_right = DATASETS$outcome_grip_right,
              alm        = DATASETS$outcome_alm,
              walking    = DATASETS$outcome_walk)

res_list <- list(); iv_list <- list(); het_list <- list(); stg_list <- list()

for (g in names(genes)) {
  for (thr in c(5e-6, 1e-5)) {
    iv <- tryCatch(extract_instruments(genes[[g]], p1 = thr, clump = TRUE,
                                       r2 = CLUMP_R2, kb = CLUMP_KB),
                   error = function(e) { cat(g, thr, "extract fail:", conditionMessage(e), "\n"); NULL })
    iv <- if (is.null(iv)) NULL else iv[!duplicated(iv$SNP), ]
    n <- if (is.null(iv)) 0 else nrow(iv)
    cat(g, "threshold", thr, "-> IVs:", n, "\n")
    if (n >= 3) {
      iv$gene <- g; iv$threshold <- thr
      iv_list[[paste(g, thr)]] <- iv
      for (on in names(outcomes)) {
        od <- tryCatch(extract_outcome_data(snps = iv$SNP, outcomes = outcomes[[on]],
                                            proxies = TRUE),
                       error = function(e) { cat("  outcome fail:", on, conditionMessage(e), "\n"); NULL })
        if (is.null(od) || nrow(od) == 0) next
        dat <- tryCatch(harmonise_data(iv, od, action = 2), error = function(e) NULL)
        if (is.null(dat) || nrow(dat) < 3) next
        key <- paste0(g, "_to_", on)
        fwrite(as.data.frame(dat), file.path(OUT, paste0("line2b_harmonised_", key, "_p", thr, ".csv")))
        m <- tryCatch(mr(dat, method_list = c("mr_ivw_mre", "mr_ivw_fe",
                                              "mr_egger_regression", "mr_weighted_median")),
                      error = function(e) NULL)
        if (!is.null(m) && nrow(m) > 0) {
          m$gene <- g; m$outcome_trait <- on; m$threshold <- thr
          res_list[[paste(key, thr)]] <- m
          h <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
          if (!is.null(h) && nrow(h) > 0) { h$gene <- g; h$outcome_trait <- on; h$threshold <- thr
                                           het_list[[paste(key, thr)]] <- h }
          s <- tryCatch(directionality_test(dat), error = function(e) NULL)
          if (!is.null(s) && nrow(s) > 0) { s$gene <- g; s$outcome_trait <- on; s$threshold <- thr
                                           stg_list[[paste(key, thr)]] <- s }
          cat("  ", key, "thr", thr, "done\n")
          print(m[, c("method", "nsnp", "b", "se", "pval")])
        }
      }
      break  # got enough IVs at this threshold, stop relaxing
    }
  }
}
if (length(iv_list)) fwrite(rbindlist(iv_list, fill = TRUE), file.path(OUT, "line2b_instruments_relaxed.csv"))
if (length(res_list)) fwrite(rbindlist(res_list, fill = TRUE), file.path(OUT, "line2b_mr_results.csv"))
if (length(het_list)) fwrite(rbindlist(het_list, fill = TRUE), file.path(OUT, "line2b_heterogeneity.csv"))
if (length(stg_list)) fwrite(rbindlist(stg_list, fill = TRUE), file.path(OUT, "line2b_steiger.csv"))
cat("line2b saved:", length(res_list), "result blocks\n")
