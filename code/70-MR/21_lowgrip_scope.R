# 21_lowgrip_scope.R -- complete the missing low-grip tests for C1QA / CTNNBIP1
# (relaxed 5e-6 cis sets; resolves the plan-B scope gap flagged in Results 3.4)
# Pipeline conventions identical to 07b_line2b_relaxed.R (clump params, action=2,
# proxies=TRUE, same method list) so results are directly comparable to Table 4.
source(file.path(getwd(), "code", "04_config.R"))
suppressMessages({library(TwoSampleMR); library(data.table)})
stopifnot("OPENGWAS_JWT missing" = nzchar(Sys.getenv("OPENGWAS_JWT")))

genes <- c(C1QA     = "eqtl-a-ENSG00000173372",
           CTNNBIP1 = "eqtl-a-ENSG00000178585")
out_lowgrip <- DATASETS$outcome_lowgrip   # ebi-a-GCST90007526 (low grip 60+, binary)

res_list <- list(); het_list <- list(); stg_list <- list()

for (g in names(genes)) {
  iv <- tryCatch(extract_instruments(genes[[g]], p1 = 5e-6, clump = TRUE,
                                     r2 = CLUMP_R2, kb = CLUMP_KB),
                 error = function(e) { cat(g, "extract fail:", conditionMessage(e), "\n"); NULL })
  if (is.null(iv)) next
  iv <- iv[!duplicated(iv$SNP), ]
  cat(g, "candidate IVs:", nrow(iv), "\n")
  od <- tryCatch(extract_outcome_data(snps = iv$SNP, outcomes = out_lowgrip, proxies = TRUE),
                 error = function(e) { cat("  outcome fail:", conditionMessage(e), "\n"); NULL })
  if (is.null(od) || nrow(od) == 0) next
  dat <- tryCatch(harmonise_data(iv, od, action = 2), error = function(e) NULL)
  if (is.null(dat)) next
  dat <- dat[dat$mr_keep, ]
  cat(g, "harmonised:", nrow(dat), "\n")
  if (nrow(dat) < 3) next
  key <- paste0(g, "_to_lowgrip")
  fwrite(as.data.frame(dat), file.path(OUT, paste0("line2b_harmonised_", key, "_p5e-06.csv")))
  m <- tryCatch(mr(dat, method_list = c("mr_ivw_mre", "mr_ivw_fe",
                                        "mr_egger_regression", "mr_weighted_median")),
                error = function(e) NULL)
  if (!is.null(m) && nrow(m) > 0) {
    m$gene <- g; m$outcome_trait <- "lowgrip"; m$threshold <- 5e-6
    res_list[[key]] <- m
    print(m[, c("method", "nsnp", "b", "se", "pval")])
  }
  h <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
  if (!is.null(h) && nrow(h) > 0) { h$gene <- g; h$outcome_trait <- "lowgrip"; h$threshold <- 5e-6
                                   het_list[[key]] <- h }
  s <- tryCatch(directionality_test(dat), error = function(e) NULL)
  if (!is.null(s) && nrow(s) > 0) { s$gene <- g; s$outcome_trait <- "lowgrip"; s$threshold <- 5e-6
                                   stg_list[[key]] <- s }
}
if (length(res_list)) fwrite(rbindlist(res_list, fill = TRUE), file.path(OUT, "line2b_lowgrip_mr_results.csv"))
if (length(het_list)) fwrite(rbindlist(het_list, fill = TRUE), file.path(OUT, "line2b_lowgrip_heterogeneity.csv"))
if (length(stg_list)) fwrite(rbindlist(stg_list, fill = TRUE), file.path(OUT, "line2b_lowgrip_steiger.csv"))
cat("lowgrip scope tests saved:", length(res_list), "blocks\n")
