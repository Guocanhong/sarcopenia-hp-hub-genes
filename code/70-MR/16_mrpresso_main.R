# 16_mrpresso_main.R -- MR-PRESSO robustness check for the primary pair:
# seropositivity (ebi-a-GCST90006910) -> ALM (ebi-a-GCST90000025), 11 IVs.
# Uses the FROZEN harmonised dataset from line1 (line1_harmonised_igg_pos_to_alm.csv)
# so the instrument set is identical to the main analysis. No OpenGWAS re-pull.
# English comments only (GBK pitfall). Output to D:/MRout/70-MR (ASCII path).
suppressMessages({library(data.table); library(MRPRESSO)})
OUT <- "D:/MRout/70-MR"
set.seed(20261105)  # reproducible permutation p-values

har <- fread(file.path(OUT, "line1_harmonised_igg_pos_to_alm.csv"))
cat("harmonised SNPs:", nrow(har), "\n")
stopifnot(nrow(har) == 11)

res <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                 SdOutcome = "se.outcome", SdExposure = "se.exposure",
                 OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                 data = as.data.frame(har), NbDistribution = 1000, SignifThresh = 0.05)
r <- res[["MR-PRESSO results"]]
cat("\n--- structure ---\n")
str(r, max.level = 2)

gt <- r[["Global Test"]]
ra <- r[["MR Analysis"]]
dt <- r[["Distortion Test"]]

# outlier info: either a dedicated table in results or appended columns in res[[1]]
ot <- NULL
if (!is.null(r[["Outlier Test"]])) {
  ot <- r[["Outlier Test"]]
} else if (is.data.frame(res[[1]]) &&
           all(c("Distortion Statistics", "Outlier P-value") %in% names(res[[1]]))) {
  ot <- res[[1]][, c("SNP", "Distortion Statistics", "Outlier P-value")]
  ot <- as.data.frame(ot)
}

cat("\n--- Global test ---\n"); print(gt)
cat("\n--- MR analysis (raw / outlier-corrected) ---\n"); str(ra, max.level = 2)
cat("\n--- Outlier test ---\n"); print(ot)
if (!is.null(dt)) { cat("\n--- Distortion test ---\n"); print(dt) }

# outlier SNPs at P < 0.05
out_snps <- character(0)
if (!is.null(ot) && is.data.frame(ot) && nrow(ot)) {
  out_rows <- ot[which(ot[["Outlier P-value"]] < 0.05), , drop = FALSE]
  out_snps <- out_rows[["SNP"]]
  cat("\noutlier SNPs (P<0.05):", paste(out_snps, collapse = ", "), "\n")
}

# extract raw / corrected estimates defensively; when no outliers are found the
# package omits the MR Analysis element and sets corrected = raw (package message),
# so fall back to the frozen line1 IVW estimates for the raw/corrected columns.
raw <- if (!is.null(ra[["Raw"]])) ra[["Raw"]] else NULL
corr <- if (!is.null(ra[["Outlier-corrected"]])) ra[["Outlier-corrected"]] else raw
gval <- function(x, nm) if (!is.null(x[[nm]])) as.numeric(x[[nm]]) else NA_real_

fr <- fread(file.path(OUT, "line1_mr_results.csv"))
fr <- fr[grepl("seropositivity", fr$exposure, ignore.case = TRUE) &
         grepl("alm", fr$outcome_trait, ignore.case = TRUE) &
         grepl("multiplicative random effects", fr$method), ]
ivw_b <- as.numeric(fr$b); ivw_se <- as.numeric(fr$se); ivw_p <- as.numeric(fr$pval)

if (is.null(raw) && nrow(fr) == 1) {
  raw <- list(Estimate = ivw_b, Sd = ivw_se, Pvalue = ivw_p)
  corr <- raw  # no outliers detected -> corrected equals raw
}

out <- data.frame(
  key = "igg_pos_to_alm",
  n_iv = nrow(har),
  nb_distribution = 1000L,
  seed = 20261105L,
  global_rssobs = gval(gt, "RSSobs"),
  global_pval   = gval(gt, "Pvalue"),
  raw_estimate  = gval(raw, "Estimate"),
  raw_se        = gval(raw, "Sd"),
  raw_pval      = gval(raw, "Pvalue"),
  n_outliers = length(out_snps),
  outlier_snps = paste(out_snps, collapse = ";"),
  corrected_estimate = gval(corr, "Estimate"),
  corrected_se = gval(corr, "Sd"),
  corrected_pval = gval(corr, "Pvalue"),
  distortion_pval = if (is.null(dt)) NA_real_ else gval(dt, "Pvalue")
)
fwrite(out, file.path(OUT, "line1_mrpresso_igg_pos_alm.csv"))
if (is.data.frame(ot)) fwrite(as.data.frame(ot), file.path(OUT, "line1_mrpresso_outlier_table.csv"))
cat("\n=== MR-PRESSO done, results written ===\n")
