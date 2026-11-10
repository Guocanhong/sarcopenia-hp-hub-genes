# 14_nfkbia_heterogeneity.R -- recompute Cochran's Q for line2 (gene -> traits)
# The line2 run saved mr results but no heterogeneity CSV. NFKBIA RE/FE SE
# differ ~3x (strong heterogeneity suspected), so Q must be reported.
# Also covers DTL (strict-threshold line2) whose heterogeneity is missing too.
source(file.path(getwd(), "code", "04_config.R"))
suppressMessages({library(TwoSampleMR); library(data.table)})
stopifnot("OPENGWAS_JWT missing" = nzchar(Sys.getenv("OPENGWAS_JWT")))

GENE_ENSG <- list(
  NFKBIA = "ENSG00000100906",
  DTL    = "ENSG00000143476"
)

id_out <- c(grip_right = DATASETS$outcome_grip_right,
            grip_left  = DATASETS$outcome_grip_left,
            alm        = DATASETS$outcome_alm,
            walking    = DATASETS$outcome_walk,
            lowgrip60  = DATASETS$outcome_lowgrip)

het_list <- list()
for (gene in names(GENE_ENSG)) {
  dsid <- paste0("eqtl-a-", GENE_ENSG[[gene]])
  cat("\n======== ", gene, ":", dsid, "========\n")
  inst <- tryCatch(extract_instruments(dsid, p1 = FDR_EQTL, clump = TRUE,
                                       r2 = CLUMP_R2, kb = CLUMP_KB),
                   error = function(e) { cat("  extract fail:", conditionMessage(e), "\n"); NULL })
  if (is.null(inst) || nrow(inst) < 2) { cat(gene, ": IVs <2 at 5e-8, skip\n"); next }
  inst$exposure <- paste0(gene, "_blood_cis_eQTL")
  inst$id.exposure <- gene
  cat("  IVs:", nrow(inst), "\n")

  for (nm in names(id_out)) {
    out_dat <- tryCatch(extract_outcome_data(snps = inst$SNP, outcomes = id_out[[nm]],
                                             proxies = TRUE),
                        error = function(e) { cat("  outcome fail:", nm, conditionMessage(e), "\n"); NULL })
    if (is.null(out_dat) || nrow(out_dat) == 0) next
    dat <- tryCatch(harmonise_data(inst, out_dat, action = 2), error = function(e) NULL)
    if (is.null(dat) || nrow(dat) < 3) { cat(gene, "->", nm, ": <3 SNPs, no Q\n"); next }
    h <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
    if (!is.null(h) && nrow(h) > 0) {
      h$gene <- gene; h$outcome_trait <- nm
      het_list[[paste(gene, nm)]] <- h
      print(h[, c("method", "Q", "Q_df", "Q_pval")])
    }
    cat(gene, "->", nm, "done\n")
  }
}
if (length(het_list)) {
  fwrite(rbindlist(het_list, fill = TRUE), file.path(OUT, "line2_heterogeneity_NFKBIA_DTL.csv"))
  cat("\nline2 heterogeneity (NFKBIA/DTL) saved\n")
} else cat("\nno heterogeneity results\n")
