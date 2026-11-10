# 13_line1c_antigen_mr.R -- Exposure-spectrum MR: 6 Hp antigen antibody levels -> sarcopenia traits
# Same design as line1 main; p1=5e-6 consistent with main Hp exposures (small effective N).
source(file.path(getwd(), "code", "04_config.R"))
suppressMessages({library(TwoSampleMR); library(data.table)})
stopifnot("OPENGWAS_JWT missing" = nzchar(Sys.getenv("OPENGWAS_JWT")))

antigens <- c(CagA     = "ebi-a-GCST90006911",
              Catalase = "ebi-a-GCST90006912",
              GroEL    = "ebi-a-GCST90006913",
              OMP      = "ebi-a-GCST90006914",
              UREA     = "ebi-a-GCST90006915",
              VacA     = "ebi-a-GCST90006916")
outcomes <- c(grip_right = DATASETS$outcome_grip_right,
              alm        = DATASETS$outcome_alm,
              walking    = DATASETS$outcome_walk,
              lowgrip    = DATASETS$outcome_lowgrip)

res_list <- list()
for (an in names(antigens)) {
  exp_dat <- tryCatch(extract_instruments(antigens[[an]], p1 = 5e-6, clump = TRUE,
                                          r2 = CLUMP_R2, kb = CLUMP_KB),
                      error = function(e) { cat(an, "IV extract fail:", conditionMessage(e), "\n"); NULL })
  if (is.null(exp_dat) || nrow(exp_dat) < 3) { cat(an, ": IVs <3, skip\n"); next }
  cat(an, "IVs:", nrow(exp_dat), "\n")
  fwrite(as.data.frame(exp_dat), file.path(OUT, paste0("line1c_instruments_", an, ".csv")))
  for (on in names(outcomes)) {
    out_dat <- tryCatch(extract_outcome_data(snps = exp_dat$SNP, outcomes = outcomes[[on]],
                                             proxies = TRUE),
                        error = function(e) { cat("  outcome fail:", conditionMessage(e), "\n"); NULL })
    if (is.null(out_dat) || nrow(out_dat) == 0) next
    dat <- tryCatch(harmonise_data(exp_dat, out_dat, action = 2), error = function(e) NULL)
    if (is.null(dat) || nrow(dat) < 3) next
    key <- paste0(an, "_to_", on)
    fwrite(as.data.frame(dat), file.path(OUT, paste0("line1c_harmonised_", key, ".csv")))
    m <- tryCatch(mr(dat, method_list = c("mr_ivw_mre", "mr_ivw_fe",
                                          "mr_egger_regression", "mr_weighted_median")),
                  error = function(e) NULL)
    if (!is.null(m) && nrow(m) > 0) {
      m$exposure_trait <- an; m$outcome_trait <- on
      res_list[[key]] <- m
      h <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
      if (!is.null(h) && nrow(h) > 0) fwrite(as.data.frame(h), file.path(OUT, paste0("line1c_heterogeneity_", key, ".csv")))
      print(m[, c("method", "nsnp", "b", "se", "pval")])
    }
  }
}
if (length(res_list)) {
  fwrite(rbindlist(res_list, fill = TRUE), file.path(OUT, "line1c_antigen_mr_results.csv"))
  cat("line1c results saved\n")
} else cat("line1c: no results\n")
