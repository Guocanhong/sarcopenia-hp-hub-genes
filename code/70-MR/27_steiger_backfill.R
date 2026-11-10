# 27_steiger_backfill.R -- Backfill Steiger directionality for pairs marked pending
# in Table S2: (a) 24 antigen-line pairs (line1c, Steiger never run), (b) the two
# reverse walking-pace pairs (line3, extraction failed at first pass and the retry
# script 08b omitted walking). Fully OFFLINE: reads the frozen harmonised CSVs,
# which already contain eaf/samplesize columns required by directionality_test().
# Mirrors 06_line1_hp_mr.R exactly (directionality_test on the harmonised dat).
source(file.path(getwd(), "code", "04_config.R"))
suppressMessages({library(TwoSampleMR); library(data.table)})

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

jobs <- list()   # key -> harmonised csv path
for (an in names(antigens)) for (on in names(outcomes)) {
  key <- paste0(an, "_to_", on)
  jobs[[key]] <- file.path(OUT, paste0("line1c_harmonised_", key, ".csv"))
}
# reverse walking pairs
for (hn in c("igg_pos", "selfreport")) {
  key <- paste0("walking_to_", hn)
  jobs[[key]] <- file.path(OUT, paste0("line3_harmonised_", key, ".csv"))
}

res_list <- list()
for (key in names(jobs)) {
  f <- jobs[[key]]
  if (!file.exists(f)) { cat(key, ": harmonised file missing, skip\n"); next }
  dat <- fread(f)
  if (!"mr_keep" %in% names(dat)) dat$mr_keep <- TRUE
  n_keep <- sum(dat$mr_keep)
  s <- tryCatch(directionality_test(dat), error = function(e) {
    cat(key, ": steiger error -", conditionMessage(e), "\n"); NULL })
  if (is.null(s) || nrow(s) == 0) next
  s$pair <- key; s$n_snps_used <- n_keep
  res_list[[key]] <- s
  # per-pair file follows the existing naming convention consumed by 18/19
  out_name <- if (grepl("^walking_to_", key)) paste0("line3_steiger_", key, ".csv")
              else paste0("line1c_steiger_", key, ".csv")
  fwrite(as.data.frame(s), file.path(OUT, out_name))
  cat(sprintf("%-28s correct=%s  steiger_pval=%.3g  (n=%d)\n",
              key, s$correct_causal_direction, s$steiger_pval, n_keep))
}
if (length(res_list)) {
  fwrite(rbindlist(res_list, fill = TRUE),
         file.path(OUT, "steiger_backfill_2026-11-07.csv"))
  cat("backfill summary saved: steiger_backfill_2026-11-07.csv\n")
} else cat("no steiger results computed\n")
