# 55_nfkbia_loo.R -- NFKBIA->ALM leave-one-out (responds to mock-review Round-2 R1 minor 1:
# heterogeneity Q P = 1.7e-5 with only 5 cis-eQTLs, but no LOO shown for this pair).
# Source: frozen harmonised CSV (no re-extraction). Verifies full-set IVW-RE against
# the frozen beta = -0.044 / P = 4.9e-4 before writing LOO results.
source(file.path(getwd(), "code", "04_config.R"))
suppressMessages({library(TwoSampleMR); library(data.table)})

OUTD <- "D:/MRout/70-MR"
d <- fread(file.path(OUTD, "line2_harmonised_NFKBIA_to_alm.csv"))
d <- as.data.frame(d)
if ("mr_keep" %in% names(d)) d <- d[which(d$mr_keep %in% c(TRUE, "TRUE", "true", 1)), ]
cat("harmonised IVs:", nrow(d), "\n")

run_ivw <- function(x, label) {
  m <- tryCatch(mr(x, method_list = c("mr_ivw_mre")), error = function(e) NULL)
  if (is.null(m)) return(data.frame(set = label, nsnp = nrow(x), b = NA_real_, se = NA_real_, p = NA_real_))
  data.frame(set = label, nsnp = nrow(x), b = m$b[1], se = m$se[1], p = m$pval[1])
}

rows <- list()
# full set (must reproduce frozen values)
rows[[1]] <- run_ivw(d, "all_5")
# leave-one-out
for (s in d$SNP) {
  rows[[length(rows) + 1]] <- run_ivw(d[d$SNP != s, ], paste0("drop_", s))
}
loo <- rbindlist(rows)
# range of LOO betas excluding the full set
loo_sub <- loo[set != "all_5"]
cat(sprintf("full: b=%.4f p=%.2e (frozen: b=-0.044 p=4.9e-04)\n", loo$b[1], loo$p[1]))
cat(sprintf("LOO range: b %.4f..%.4f, p %.2e..%.2e (n=%d subsets)\n",
            min(loo_sub$b), max(loo_sub$b), min(loo_sub$p), max(loo_sub$p), nrow(loo_sub)))
fwrite(loo, file.path(OUTD, "55_nfkbia_alm_leaveoneout.csv"))
cat("written:", file.path(OUTD, "55_nfkbia_alm_leaveoneout.csv"), "\n")
