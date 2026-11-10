# 06_line1_hp_mr.R -- Line 1: Hp infection -> sarcopenia traits
# Exposures: 3 Hp GWAS; Outcomes: 5 sarcopenia-related traits
source(file.path(getwd(), "code", "04_config.R"))
suppressMessages({library(TwoSampleMR); library(data.table)})
stopifnot("OPENGWAS_JWT missing" = nzchar(Sys.getenv("OPENGWAS_JWT")))

id_hp <- c(selfreport = DATASETS$exposure_hp_selfreport,
           igg_pos    = DATASETS$exposure_hp_igg_pos,
           igg_level  = DATASETS$exposure_hp_igg_level)
id_out <- c(grip_right = DATASETS$outcome_grip_right,
            grip_left  = DATASETS$outcome_grip_left,
            alm        = DATASETS$outcome_alm,
            walking    = DATASETS$outcome_walk,
            lowgrip60  = DATASETS$outcome_lowgrip)

res_all <- list(); har_all <- list(); sens_all <- list(); inst_all <- list()

for (hn in names(id_hp)) {
  cat("\n============ exposure:", hn, id_hp[[hn]], "============\n")
  exp_dat <- tryCatch(
    extract_instruments(id_hp[[hn]], p1 = P_EXP, clump = TRUE,
                        r2 = CLUMP_R2, kb = CLUMP_KB),
    error = function(e) { cat("IV extract fail:", conditionMessage(e), "\n"); NULL })
  if (is.null(exp_dat) || nrow(exp_dat) < 3) {
    cat("genome-wide IVs < 3, relax to 5e-6\n")
    exp_dat <- tryCatch(
      extract_instruments(id_hp[[hn]], p1 = 5e-6, clump = TRUE,
                          r2 = CLUMP_R2, kb = CLUMP_KB),
      error = function(e) { cat("relaxed extract fail:", conditionMessage(e), "\n"); NULL })
    if (!is.null(exp_dat) && nrow(exp_dat)) exp_dat$relaxed <- TRUE
  }
  if (is.null(exp_dat) || nrow(exp_dat) < 3) { cat("skip exposure:", hn, "\n"); next }
  exp_dat$exposure_id <- id_hp[[hn]]
  if (!is.null(exp_dat$relaxed) && any(exp_dat$relaxed %in% TRUE)) {
    exp_dat$p_threshold <- 5e-6
  } else {
    exp_dat$p_threshold <- P_EXP
  }
  inst_all[[hn]] <- exp_dat
  cat("IVs:", nrow(exp_dat), " F-stat min:",
      round(min(exp_dat$beta.exposure^2 / exp_dat$se.exposure^2), 1), "\n")

  for (on in names(id_out)) {
    cat("\n---- outcome:", on, id_out[[on]], "----\n")
    out_dat <- tryCatch(extract_outcome_data(snps = exp_dat$SNP,
                                             outcomes = id_out[[on]], proxies = TRUE),
                        error = function(e) { cat("outcome extract fail:", conditionMessage(e), "\n"); NULL })
    if (is.null(out_dat) || nrow(out_dat) == 0) next
    dat <- tryCatch(harmonise_data(exp_dat, out_dat, action = 2),
                    error = function(e) { cat("harmonise fail:", conditionMessage(e), "\n"); NULL })
    if (is.null(dat) || nrow(dat) < 3) { cat("harmonised <3 SNPs, skip\n"); next }
    key <- paste0(hn, "_to_", on)
    har_all[[key]] <- dat
    fwrite(as.data.frame(dat), file.path(OUT, paste0("line1_harmonised_", key, ".csv")))
    m <- tryCatch(mr(dat, method_list = c("mr_ivw_mre", "mr_ivw_fe",
                                          "mr_egger_regression", "mr_weighted_median",
                                          "mr_weighted_mode", "mr_wald_ratio")),
                  error = function(e) { cat("mr fail:", conditionMessage(e), "\n"); NULL })
    if (!is.null(m)) {
      m$exposure_trait <- hn; m$outcome_trait <- on
      res_all[[key]] <- m
      print(m[, c("method", "nsnp", "b", "se", "pval")])
    }
    sens_all[[key]] <- list(
      pleio   = tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL),
      het     = tryCatch(mr_heterogeneity(dat), error = function(e) NULL),
      leave1  = tryCatch(mr_leaveoneout(dat), error = function(e) NULL),
      steiger = tryCatch(directionality_test(dat), error = function(e) NULL)
    )
  }
}

if (length(inst_all)) {
  fwrite(rbindlist(lapply(names(inst_all), function(n)
    { d <- as.data.frame(inst_all[[n]]); d$exposure_trait <- n; d }), fill = TRUE),
    file.path(OUT, "line1_instruments_hp.csv"))
}
if (length(res_all)) {
  fwrite(rbindlist(res_all, fill = TRUE), file.path(OUT, "line1_mr_results.csv"))
  cat("\nline1 results written\n")
}
if (length(sens_all)) {
  pl <- rbindlist(lapply(names(sens_all), function(n)
    { x <- sens_all[[n]]$pleio; if (!is.null(x)) x$key <- n; x }), fill = TRUE)
  if (nrow(pl)) fwrite(pl, file.path(OUT, "line1_pleiotropy.csv"))
  ht <- rbindlist(lapply(names(sens_all), function(n)
    { x <- sens_all[[n]]$het; if (!is.null(x)) x$key <- n; x }), fill = TRUE)
  if (nrow(ht)) fwrite(ht, file.path(OUT, "line1_heterogeneity.csv"))
  st <- rbindlist(lapply(names(sens_all), function(n)
    { x <- sens_all[[n]]$steiger; if (!is.null(x)) x$key <- n; x }), fill = TRUE)
  if (nrow(st)) fwrite(st, file.path(OUT, "line1_steiger.csv"))
  lo <- rbindlist(lapply(names(sens_all), function(n)
    { x <- sens_all[[n]]$leave1; if (!is.null(x)) x$key <- n; x }), fill = TRUE)
  if (nrow(lo)) fwrite(lo, file.path(OUT, "line1_leaveoneout.csv"))
}

# plots for the primary pair (selfreport -> each outcome)
for (nm in names(har_all)) {
  if (!grepl("^selfreport_", nm)) next
  dat <- har_all[[nm]]
  tryCatch({
    p1 <- mr_scatter_plot(mr(dat), dat)
    ggsave(file.path(OUT, paste0("line1_scatter_", nm, ".png")), p1[[1]],
           width = 9, height = 7, dpi = 300)
    p4 <- mr_forest_plot(mr_singlesnp(dat))
    ggsave(file.path(OUT, paste0("line1_forest_", nm, ".png")), p4[[1]],
           width = 9, height = 7, dpi = 300)
    p3 <- mr_leaveoneout_plot(mr_leaveoneout(dat))
    ggsave(file.path(OUT, paste0("line1_loo_", nm, ".png")), p3[[1]],
           width = 9, height = 7, dpi = 300)
    cat(nm, ": plots saved\n")
  }, error = function(e) cat(nm, "plot fail:", conditionMessage(e), "\n"))
}
cat("\n=== line1 done ===\n")
