# 08_line3_reverse_mr.R -- Line 3 (reverse): sarcopenia traits -> Hp infection
source(file.path(getwd(), "code", "04_config.R"))
suppressMessages({library(TwoSampleMR); library(data.table)})
stopifnot("OPENGWAS_JWT missing" = nzchar(Sys.getenv("OPENGWAS_JWT")))

id_hp <- c(selfreport = DATASETS$exposure_hp_selfreport,
           igg_pos    = DATASETS$exposure_hp_igg_pos)
id_exp <- c(grip_right = DATASETS$outcome_grip_right,
            alm        = DATASETS$outcome_alm,
            walking    = DATASETS$outcome_walk)

res_list <- list()
for (en in names(id_exp)) {
  exp_dat <- tryCatch(extract_instruments(id_exp[[en]], p1 = P_EXP, clump = TRUE,
                                          r2 = CLUMP_R2, kb = CLUMP_KB),
                      error = function(e) { cat(en, "IV extract fail:", conditionMessage(e), "\n"); NULL })
  if (is.null(exp_dat) || nrow(exp_dat) < 3) { cat(en, ": IVs <3, skip\n"); next }
  cat(en, "IVs:", nrow(exp_dat), "\n")
  for (hn in names(id_hp)) {
    out_dat <- tryCatch(extract_outcome_data(snps = exp_dat$SNP, outcomes = id_hp[[hn]],
                                             proxies = TRUE),
                        error = function(e) { cat("  outcome fail:", conditionMessage(e), "\n"); NULL })
    if (is.null(out_dat) || nrow(out_dat) == 0) next
    dat <- tryCatch(harmonise_data(exp_dat, out_dat, action = 2), error = function(e) NULL)
    if (is.null(dat) || nrow(dat) < 3) next
    key <- paste0(en, "_to_", hn)
    fwrite(as.data.frame(dat), file.path(OUT, paste0("line3_harmonised_", key, ".csv")))
    m <- tryCatch(mr(dat, method_list = c("mr_ivw_mre", "mr_ivw_fe",
                                          "mr_egger_regression", "mr_weighted_median")),
                  error = function(e) NULL)
    if (!is.null(m)) { m$exposure_trait <- en; m$outcome_trait <- hn
                       res_list[[key]] <- m; print(m[, c("method","nsnp","b","se","pval")]) }
  }
}
if (length(res_list)) {
  fwrite(rbindlist(res_list, fill = TRUE), file.path(OUT, "line3_reverse_mr_results.csv"))
  cat("line3 results saved\n")
} else cat("line3: no results\n")
