# 53_egger_backfill.R -- Backfill MR-Egger intercepts + weighted-mode estimates
# Responds to mock-review R1-M4 (Egger intercepts missing where >=3 IVs) and
# R1-M2 (weak-instrument robustness: report weighted mode alongside IVW/WM).
# Sources: frozen harmonised CSVs (no re-extraction needed except NFKBIA->ALM).
source(file.path(getwd(), "code", "04_config.R"))
suppressMessages({library(TwoSampleMR); library(data.table)})

OUTD <- "D:/MRout/70-MR"
rows <- list()

# ---------- helper: run mr() wrapper on a filtered harmonised data.frame ----------
calc_pair <- function(line, pair, d) {
  d <- as.data.frame(d)
  if ("mr_keep" %in% names(d)) d <- d[which(d$mr_keep %in% c(TRUE, "TRUE", "true", 1)), ]
  n <- nrow(d)
  if (n < 3) return(NULL)
  meanF <- mean((d$beta.exposure / d$se.exposure)^2)
  m <- tryCatch(mr(d, method_list = c("mr_ivw_mre", "mr_egger_regression",
                                      "mr_weighted_median", "mr_weighted_mode")),
                error = function(e) { cat("  mr ERR:", pair, conditionMessage(e), "\n"); NULL })
  if (is.null(m)) return(NULL)
  pick <- function(nm, col) {
    i <- which(m$method == nm)
    if (!length(i)) return(NA_real_)
    m[[col]][i[1]]
  }
  # this TwoSampleMR version returns the Egger intercept as b_i/se_i/pval_i
  egr <- tryCatch(TwoSampleMR:::mr_egger_regression(d$beta.exposure, d$beta.outcome,
                                                    d$se.exposure, d$se.outcome),
                  error = function(e) NULL)
  get_i <- function(a) if (is.null(egr) || is.null(egr[[a]])) NA_real_ else as.numeric(egr[[a]])
  data.frame(line = line, pair = pair, nsnp = n, mean_F = round(meanF, 1),
             ivw_b = pick("Inverse variance weighted (multiplicative random effects)", "b"),
             ivw_p = pick("Inverse variance weighted (multiplicative random effects)", "pval"),
             egger_b = pick("MR Egger", "b"), egger_p = pick("MR Egger", "pval"),
             intercept_b = get_i("b_i"),
             intercept_se = get_i("se_i"),
             intercept_p = get_i("pval_i"),
             wm_p = pick("Weighted median", "pval"),
             mode_b = pick("Weighted mode", "b"), mode_se = pick("Weighted mode", "se"),
             mode_p = pick("Weighted mode", "pval"),
             stringsAsFactors = FALSE)
}

# ---------- 1) antigen pairs + igg-level pairs from frozen harmonised CSVs ----------
files <- c(
  sort(Sys.glob(file.path(OUTD, "line1c_harmonised_*.csv"))),
  sort(Sys.glob(file.path(OUTD, "line1_harmonised_igg_level_to_*.csv"))),
  sort(Sys.glob(file.path(OUTD, "line1_harmonised_igg_pos_to_*.csv"))),
  sort(Sys.glob(file.path(OUTD, "line2b_harmonised_*.csv")))
)
for (f in files) {
  fn <- basename(f)
  expo <- sub("^line1c_harmonised_(.+)_to_.+\\.csv$", "\\1", fn)
  if (grepl("^line1_harmonised_igg_level", fn)) expo <- "IgG level (titer)"
  if (grepl("^line1_harmonised_igg_pos", fn)) expo <- "Seropositivity"
  if (grepl("^line2b_harmonised", fn)) expo <- sub("^line2b_harmonised_(.+)_to_.+\\.csv$", "\\1", fn)
  outc <- sub("^line1c?_harmonised_.+_to_(.+)\\.csv$", "\\1", fn)
  if (grepl("^line2b_harmonised", fn)) outc <- sub("^line2b_harmonised_.+_to_(.+?)_p5e-06\\.csv$", "\\1", fn)
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d)) next
  r <- calc_pair("Antigen/Titer", paste0(expo, " -> ", outc), d)
  if (!is.null(r)) { rows[[length(rows) + 1]] <- r
    cat(sprintf("%-30s n=%2d  intP=%.3g  modeP=%.3g\n", paste0(expo, "->", outc),
                r$nsnp, r$intercept_p, r$mode_p)) }
}

# ---------- 2) NFKBIA -> all five outcomes (strict line; re-extract) ----------
cat("\n== NFKBIA re-extraction ==\n")
inst <- as.data.frame(fread(file.path(OUTD, "line2_instruments_NFKBIA.csv")))
inst$exposure <- "NFKBIA_blood_cis_eQTL"; inst$id.exposure <- "NFKBIA"
for (nm in c("alm", "grip_right", "grip_left", "walking", "lowgrip60")) {
  key <- switch(nm, walking = "walk", lowgrip60 = "lowgrip", nm)
  out <- tryCatch(extract_outcome_data(snps = inst$SNP,
                                       outcomes = DATASETS[[paste0("outcome_", key)]],
                                       proxies = TRUE),
                  error = function(e) { cat("extract fail:", nm, conditionMessage(e), "\n"); NULL })
  if (is.null(out) || nrow(out) == 0) next
  dat <- tryCatch(harmonise_data(inst, out, action = 2), error = function(e) NULL)
  if (is.null(dat) || nrow(dat) < 3) next
  val <- tryCatch(mr(dat, method_list = "mr_ivw_mre"), error = function(e) NULL)
  r <- calc_pair("Gene (strict)", paste0("NFKBIA -> ", nm), dat)
  if (!is.null(r)) {
    if (!is.null(val)) { r$validation_ivw_b <- val$b[1]; r$validation_ivw_se <- val$se[1]; r$validation_ivw_p <- val$pval[1] }
    rows[[length(rows) + 1]] <- r
  }
  fwrite(as.data.table(dat), file.path(OUTD, paste0("line2_harmonised_NFKBIA_to_", nm, "_backfill.csv")))
  cat(sprintf("NFKBIA -> %-10s n=%d  ivwP=%.3g  intP=%.3g\n", nm, r$nsnp, r$ivw_p, r$intercept_p))
}

res <- rbindlist(rows, fill = TRUE)
fwrite(res, file.path(OUTD, "53_egger_mode_backfill.csv"))
cat("\nsaved:", file.path(OUTD, "53_egger_mode_backfill.csv"), " rows:", nrow(res), "\n")
print(as.data.frame(res)[, c("pair", "nsnp", "intercept_b", "intercept_se", "intercept_p", "mode_p")])
