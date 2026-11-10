# 07_line2_eqtl_mr.R -- Line 2: four hub genes' blood cis-eQTL -> sarcopenia traits
# Source: OpenGWAS eqtl-a-<ENSG> (eQTLGen blood cis-eQTL)
source(file.path(getwd(), "code", "04_config.R"))
suppressMessages({library(TwoSampleMR); library(data.table)})
stopifnot("OPENGWAS_JWT missing" = nzchar(Sys.getenv("OPENGWAS_JWT")))

# versioned ids verified via GTEx reference API; eqtl-a uses unversioned ENSG
GENE_ENSG <- list(
  C1QA     = "ENSG00000173372",
  NFKBIA   = "ENSG00000100906",
  DTL      = "ENSG00000143476",
  CTNNBIP1 = "ENSG00000178585"
)

id_out <- c(grip_right = DATASETS$outcome_grip_right,
            grip_left  = DATASETS$outcome_grip_left,
            alm        = DATASETS$outcome_alm,
            walking    = DATASETS$outcome_walk,
            lowgrip60  = DATASETS$outcome_lowgrip)

res_list <- list()
for (gene in names(GENE_ENSG)) {
  dsid <- paste0("eqtl-a-", GENE_ENSG[[gene]])
  cat("\n======== ", gene, ":", dsid, "========\n")
  inst <- tryCatch(extract_instruments(dsid, p1 = FDR_EQTL, clump = TRUE,
                                       r2 = CLUMP_R2, kb = CLUMP_KB),
                   error = function(e) { cat("  extract fail:", conditionMessage(e), "\n"); NULL })
  if (is.null(inst) || nrow(inst) < 2) {
    cat("  IVs <2 at 5e-8, relax to 1e-5 (cis window, noted in methods)\n")
    inst <- tryCatch(extract_instruments(dsid, p1 = 1e-5, clump = TRUE,
                                         r2 = CLUMP_R2, kb = CLUMP_KB),
                     error = function(e) { cat("  relaxed fail:", conditionMessage(e), "\n"); NULL })
    if (!is.null(inst) && nrow(inst)) inst$p_threshold <- 1e-5
  } else {
    inst$p_threshold <- FDR_EQTL
  }
  if (is.null(inst) || nrow(inst) < 2) { cat(gene, ": IVs <2, skip\n"); next }
  inst$exposure <- paste0(gene, "_blood_cis_eQTL")
  inst$id.exposure <- gene
  cat("  IVs:", nrow(inst), " F-stat min:",
      round(min(inst$beta.exposure^2 / inst$se.exposure^2), 1), "\n")
  fwrite(as.data.frame(inst), file.path(OUT, paste0("line2_instruments_", gene, ".csv")))

  for (nm in names(id_out)) {
    out_dat <- tryCatch(extract_outcome_data(snps = inst$SNP, outcomes = id_out[[nm]],
                                             proxies = TRUE),
                        error = function(e) { cat("  outcome fail:", nm, conditionMessage(e), "\n"); NULL })
    if (is.null(out_dat) || nrow(out_dat) == 0) next
    dat <- tryCatch(harmonise_data(inst, out_dat, action = 2), error = function(e) NULL)
    if (is.null(dat) || nrow(dat) < 2) next
    m <- tryCatch(mr(dat, method_list = if (nrow(dat) == 2)
                        c("mr_weighted_median", "mr_wald_ratio")
                      else c("mr_ivw_mre", "mr_ivw_fe", "mr_egger_regression",
                             "mr_weighted_median", "mr_wald_ratio")),
                  error = function(e) NULL)
    if (!is.null(m) && nrow(m) > 0) {
      m$gene <- gene; m$outcome_trait <- nm
      res_list[[paste(gene, nm)]] <- m
      print(m[, c("method", "nsnp", "b", "se", "pval")])
    }
    cat(gene, "->", nm, "done\n")
  }
}
if (length(res_list)) {
  fwrite(rbindlist(res_list, fill = TRUE), file.path(OUT, "line2_mr_results.csv"))
  cat("\nline2 results saved\n")
} else cat("\nline2: no results\n")
