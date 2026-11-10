# 37_mediation_test.R -- Supervisor-review item 1 (2026-11-07)
# Proximate mediation test: H. pylori seropositivity -> NFKBIA whole-blood expression
# Exposure  = frozen seropositivity IVs (line1_instruments_hp.csv, igg_pos, 11 SNPs, P<5e-6)
# Outcome   = eQTLGen blood cis-eQTL of NFKBIA (OpenGWAS eqtl-a-ENSG00000100906)
# NOTE: this is an exploratory test of the Hp -> NF-kB pathway link (proximate,
# expression in whole blood, not muscle); results feed Results 3.4 / Discussion 4.3.
source(file.path(getwd(), "code", "04_config.R"))
suppressMessages({library(TwoSampleMR); library(data.table)})
stopifnot("OPENGWAS_JWT missing" = nzchar(Sys.getenv("OPENGWAS_JWT")))

OUT_ID <- "eqtl-a-ENSG00000100906"   # NFKBIA blood cis-eQTL (eQTLGen, same source as line 2)

# 1. frozen exposure instruments
inst_all <- fread(file.path(OUT, "line1_instruments_hp.csv"))
inst <- inst_all[exposure_trait == "igg_pos"]
stopifnot(nrow(inst) == 11)
inst$exposure <- "Anti-helicobacter pylori IgG seropositivity (frozen IVs)"
inst$id.exposure <- "igg_pos"
cat("Frozen seropositivity IVs:", nrow(inst), "\n")

# 2. outcome data: NFKBIA eQTL summary stats at the exposure SNPs
out_dat <- extract_outcome_data(snps = inst$SNP, outcomes = OUT_ID, proxies = TRUE)
cat("Outcome rows retrieved:", nrow(out_dat), "\n")

# 3. harmonise + estimate
dat <- harmonise_data(inst, out_dat, action = 2)
cat("Harmonised SNPs:", sum(dat$mr_keep), "\n")
dat <- dat[dat$mr_keep, ]

res <- mr(dat, method_list = c("mr_ivw_mre", "mr_ivw_fe", "mr_egger_regression",
                               "mr_weighted_median", "mr_weighted_mode"))
res$exposure <- "H. pylori seropositivity"; res$outcome <- "NFKBIA whole-blood expression (eQTLGen)"
print(res[, c("method", "nsnp", "b", "se", "pval")])
fwrite(as.data.frame(res), file.path(OUT, "line4_mediation_iggpos_to_NFKBIA_mr.csv"))

# 4. heterogeneity
het <- mr_heterogeneity(dat)
print(het[, c("method", "Q", "Q_df", "Q_pval")])
fwrite(as.data.frame(het), file.path(OUT, "line4_mediation_iggpos_to_NFKBIA_heterogeneity.csv"))

# 5. Egger intercept (>=3 SNPs)
if (nrow(dat) >= 3) {
  int <- mr_pleiotropy_test(dat)
  print(int)
  fwrite(as.data.frame(int), file.path(OUT, "line4_mediation_iggpos_to_NFKBIA_egger_intercept.csv"))
}

# 6. Steiger directionality (binary exposure N=8,735; outcome per-SNP N varies)
dat$r.outcome <- get_r(dat$beta.outcome, dat$se.outcome, dat$samplesize.outcome)
dat$r.exposure <- get_r(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure)
st <- directionality_test(dat)
print(st)
fwrite(as.data.frame(st), file.path(OUT, "line4_mediation_iggpos_to_NFKBIA_steiger.csv"))

# 7. leave-one-out
loo <- mr_leaveoneout(dat)
fwrite(as.data.frame(loo), file.path(OUT, "line4_mediation_iggpos_to_NFKBIA_loo.csv"))

cat("\nDone. Results in", OUT, "\n")
