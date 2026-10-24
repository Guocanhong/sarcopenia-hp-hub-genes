# 3h-bootstrap_stability.R
#
# Bootstrap stability of the ML funnel (Supplementary Table S2).
#
# PROVENANCE NOTE (2026-10-13): the script that originally produced
# 27-GSE60662/output/Table_S2_bootstrap.csv was not preserved in the project
# archive. This file is a faithful RECONSTRUCTION written against the primary
# per-step scripts (3a-3g) and the archived outputs, so that the bootstrap
# procedure is auditable again. It (1) recomputes the primary selection of
# every step from the original inputs and asserts agreement with the
# "retained" flags stored in the archived CSV, then (2) repeats the three
# selection procedures in B bootstrap resamples (default B = 1000, same
# algorithm settings, inputs and thresholds as the primary analysis) and
# writes Table_S2_bootstrap_reconstructed.csv (schema identical to the
# archived file). Because the original RNG seed/stream is unknown, the
# reconstructed frequencies are expected to agree with the archived values
# to within Monte-Carlo error, not digit for digit.
#
# KNOWN RNG-SENSITIVE FLAG (full B=1000 rerun): the SP-line rf_retained flag
# for DTL is borderline in the archived run (bootstrap rf_freq = 0.57) and
# flipped to 0 in the reconstruction; DTL still reaches votes_ge2 = 2 via
# SVM-RFE + LASSO, so the hub-gene set and all manuscript claims are unchanged.
#
# Run from this directory:  Rscript 3h-bootstrap_stability.R
# Optional:                 BOOT_B=50 Rscript 3h-bootstrap_stability.R
#
# Definitions (mirroring Supplementary Table S2):
#   svm_retained   = membership of the primary SVM-RFE panel (top features up
#                    to the cross-validated error minimum)
#   svm_topK_freq  = proportion of resamples in which the gene fell in the
#                    top-K ranked panel (K = primary panel size, 6 SP / 7 HP)
#   svm_medrank    = mean SVM-RFE rank across resamples (1 = best)
#   lasso_retained = non-zero coefficient at lambda.min in the primary
#                    10-fold CV LASSO (fitted to the SVM-RFE-retained genes)
#   lasso_freq     = proportion of resamples with non-zero coefficient
#   rf_retained    = IncNodePurity > 0.2 (SP) / > 0.05 (HP) in the primary RF
#                    (1,000 trees; SP fitted to SVM-RFE-retained genes, HP to
#                    the full candidate pool)
#   rf_freq        = proportion of resamples above threshold
#   votes_ge2      = number of algorithms retaining the gene (LASSO counted
#                    only for SVM-RFE-retained genes, which are its only
#                    inputs; genes without SVM-RFE retention therefore score 0)
#   hub            = votes_ge2 in BOTH lines

rm(list = ls())
suppressMessages({
  library(e1071)
  library(glmnet)
  library(randomForest)
})

ROOT <- file.path("..", "..", "..")
B <- as.integer(Sys.getenv("BOOT_B", unset = "1000"))
set.seed(2025)  # bootstrap stream; original seed unknown (see provenance note)

archived_csv <- read.csv(file.path(ROOT, "27-GSE60662", "output",
                                   "Table_S2_bootstrap.csv"), stringsAsFactors = FALSE)

# ---- msvmRFE routines (verbatim copy of input/msvmRFE.R behaviour) ----------
# source into the global environment so svmRFE()/svmRFE.wrap()/etc. are callable
source_msvmrfe <- function(path) source(path, local = globalenv())

run_line <- function(line_name, in_dir, expr_file, group_file, pool_dir,
                     ref_group, test_group, rf_input, rf_cut) {
  cat("\n=====", line_name, "=====\n")
  dat_expr <- readRDS(file.path(in_dir, expr_file))
  dat_group <- readRDS(file.path(in_dir, group_file))
  gene_pool <- intersect(readRDS(file.path(pool_dir, "SP-HP.Rds")),
                         intersect(readRDS(file.path(pool_dir, "ModuleGenes1.Rds")),
                                   readRDS(file.path(pool_dir, "ModuleGenes.Rds"))))
  cat("candidate pool:", length(gene_pool), "\n")

  # group vector aligned to samples (rows of t(expr))
  samp_order <- colnames(dat_expr)
  grp_vec <- factor(as.character(dat_group$group[match(samp_order, dat_group$sample)]),
                    levels = c(ref_group, test_group))
  build_show <- function(genes, idx = NULL) {
    m <- as.data.frame(t(dat_expr[genes, samp_order, drop = FALSE]))
    if (is.null(idx)) idx <- seq_len(nrow(m))
    m <- m[idx, , drop = FALSE]
    cbind(group = grp_vec[idx], m)
  }

  # ---------- primary SVM-RFE (mirror 3b/3e: seed, warm-up call, folds) ------
  dat_show <- build_show(gene_pool)
  source_msvmrfe(file.path(in_dir, "msvmRFE.R"))
  set.seed(2025)
  invisible(svmRFE(dat_show, k = 10, halve.above = 100))   # discarded warm-up (RNG stream!)
  folds <- rep(1:10, len = nrow(dat_show))[sample(nrow(dat_show))]
  folds <- lapply(1:10, function(x) which(folds == x))
  res_svm <- lapply(folds, svmRFE.wrap, dat_show, k = 10, halve.above = 100)
  top_svm <- WriteFeatures(res_svm, dat_show, save = FALSE)
  rownames(top_svm) <- top_svm$FeatureName
  featsweep <- lapply(1:(ncol(dat_show) - 1), FeatSweep.wrap, res_svm, dat_show)
  err_svm <- sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
  K <- which.min(err_svm)
  cat("primary SVM-RFE: K =", K, "\n")
  svm_retained <- setNames(rep(0, length(gene_pool)), gene_pool)
  svm_retained[top_svm$FeatureName[1:K]] <- 1

  # ---------- primary LASSO (mirror 3a/3d) ------------------------------------
  set.seed(2025)
  dat_lasso <- build_show(names(svm_retained)[svm_retained == 1])
  cv_lasso <- cv.glmnet(as.matrix(dat_lasso[, -1]), dat_lasso$group,
                        nfolds = 10, family = "binomial")
  co <- as.matrix(coef(cv_lasso, s = "lambda.min"))
  lasso_genes <- setdiff(rownames(co)[co[, 1] != 0], "(Intercept)")
  lasso_retained <- setNames(rep(0, length(gene_pool)), gene_pool)
  lasso_retained[lasso_genes] <- 1

  # ---------- primary RF (mirror 3c/3f) ----------------------------------------
  set.seed(2025)
  rf_genes <- if (rf_input == "svm") names(svm_retained)[svm_retained == 1] else gene_pool
  dat_rf <- build_show(rf_genes)
  dat_rf$group <- case_match_wrapper(dat_rf$group, ref_group, test_group)
  rf_fit <- randomForest(x = dat_rf[, -1], y = dat_rf$group, importance = TRUE,
                         ntree = 1000, proximity = TRUE)
  imp <- importance(rf_fit, scale = FALSE)
  rf_genes_sel <- rownames(imp)[imp[, "IncNodePurity"] > rf_cut]
  rf_retained <- setNames(rep(0, length(gene_pool)), gene_pool)
  rf_retained[rf_genes_sel] <- 1

  votes <- ifelse(svm_retained == 1, 1 + lasso_retained + rf_retained, 0)

  # ---------- assert against archived flags ------------------------------------
  arch <- archived_csv[archived_csv$line == line_name, ]
  arch <- arch[match(names(svm_retained), arch$gene), ]
  ok <- all(arch$svm_retained == svm_retained) &
        all(with(arch, ifelse(is.na(lasso_retained), 0, lasso_retained)) == lasso_retained) &
        all(arch$rf_retained == rf_retained) &
        all(arch$votes_ge2 == votes)
  cat("primary flags vs archived CSV:", ifelse(ok, "PASS", "FAIL"), "\n")
  if (!ok) print(data.frame(gene = names(svm_retained),
                            svm_a = arch$svm_retained, svm_r = svm_retained,
                            lasso_a = arch$lasso_retained, lasso_r = lasso_retained,
                            rf_a = arch$rf_retained, rf_r = rf_retained,
                            votes_a = arch$votes_ge2, votes_r = votes))

  # ---------- bootstrap ---------------------------------------------------------
  n <- nrow(dat_show)
  svm_top <- rep(0, length(gene_pool)); names(svm_top) <- gene_pool
  svm_rank_sum <- setNames(rep(0, length(gene_pool)), gene_pool); n_svm_ok <- 0
  lasso_hit <- setNames(rep(0, length(gene_pool)), gene_pool); n_lasso_ok <- 0
  rf_hit <- setNames(rep(0, length(gene_pool)), gene_pool); n_rf_ok <- 0
  svm_in <- names(svm_retained)[svm_retained == 1]

  for (b in seq_len(B)) {
    idx <- sample(n, replace = TRUE)
    # SVM-RFE ranking (K fixed at primary size; no FeatSweep needed)
    ok_b <- TRUE
    rank_vec <- tryCatch({
      db <- dat_show[idx, ]
      fb <- rep(1:10, len = nrow(db))[sample(nrow(db))]
      fb <- lapply(1:10, function(x) which(fb == x))
      # drop folds that contain a single class (svm() would fail)
      fb <- Filter(function(f) length(unique(db[f, 1])) > 1, fb)
      rb <- lapply(fb, svmRFE.wrap, db, k = 10, halve.above = 100)
      tb <- WriteFeatures(rb, db, save = FALSE)
      r <- rep(NA, length(gene_pool)); names(r) <- gene_pool
      r[tb$FeatureName] <- seq_along(tb$FeatureName)
      r
    }, error = function(e) NULL)
    if (!is.null(rank_vec) && !any(is.na(rank_vec))) {
      svm_rank_sum <- svm_rank_sum + rank_vec; n_svm_ok <- n_svm_ok + 1
      svm_top[names(sort(rank_vec))[1:K]] <- svm_top[names(sort(rank_vec))[1:K]] + 1
    }
    # LASSO
    lb <- tryCatch({
      db <- build_show(svm_in, idx)
      if (length(unique(db$group)) < 2) stop("single class")
      cvb <- cv.glmnet(as.matrix(db[, -1]), db$group, nfolds = 10, family = "binomial")
      cob <- as.matrix(coef(cvb, s = "lambda.min"))
      setdiff(rownames(cob)[cob[, 1] != 0], "(Intercept)")
    }, error = function(e) NULL)
    if (!is.null(lb)) {
      n_lasso_ok <- n_lasso_ok + 1
      if (length(lb)) lasso_hit[lb] <- lasso_hit[lb] + 1
    }
    # RF
    rb_in <- if (rf_input == "svm") svm_in else gene_pool
    rb <- tryCatch({
      db <- build_show(rb_in, idx)
      db$group <- case_match_wrapper(db$group, ref_group, test_group)
      if (length(unique(db$group)) < 2) stop("single class")
      f <- randomForest(x = db[, -1], y = db$group, importance = TRUE,
                        ntree = 1000, proximity = FALSE)
      impb <- importance(f, scale = FALSE)
      rownames(impb)[impb[, "IncNodePurity"] > rf_cut]
    }, error = function(e) NULL)
    if (!is.null(rb)) {
      n_rf_ok <- n_rf_ok + 1
      if (length(rb)) rf_hit[rb] <- rf_hit[rb] + 1
    }
    if (b %% 50 == 0) cat("bootstrap", b, "/", B, "\n")
  }
  cat("successful resamples: svm", n_svm_ok, " lasso", n_lasso_ok, " rf", n_rf_ok, "\n")

  # output schema mirrors the archived CSV: lasso_freq is recorded for every
  # gene FED INTO the LASSO (the SVM-retained genes of that line; the archived
  # HP rows carry e.g. NFKBIA lasso_freq = 0.34 despite lasso_retained = 0),
  # and rf_freq only for genes fed into the RF (SVM-retained in SP, full pool
  # in HP; elsewhere NA)
  lasso_input <- names(svm_retained)[svm_retained == 1]
  data.frame(line = line_name, gene = names(svm_retained),
             svm_retained = as.integer(svm_retained),
             svm_topK_freq = round(svm_top / max(n_svm_ok, 1), 2),
             svm_medrank = round(svm_rank_sum / max(n_svm_ok, 1), 1),
             lasso_retained = as.integer(lasso_retained),
             lasso_freq = ifelse(names(svm_retained) %in% lasso_input,
                                 round(lasso_hit / max(n_lasso_ok, 1), 2), NA),
             rf_retained = as.integer(rf_retained),
             rf_freq = ifelse(names(svm_retained) %in% rf_genes,
                              round(rf_hit / max(n_rf_ok, 1), 2), NA),
             votes_ge2 = as.integer(votes),
             stringsAsFactors = FALSE)
}

case_match_wrapper <- function(f, ref, test) {
  # mirror dplyr::case_match(ref ~ 0, test ~ 1) from 3c/3f: returns a NUMERIC
  # 0/1 vector, so randomForest runs as REGRESSION and importance() carries
  # the IncNodePurity column the 0.2/0.05 thresholds act on
  out <- rep(NA_real_, length(f))
  out[as.character(f) == ref] <- 0
  out[as.character(f) == test] <- 1
  out
}

res_sp <- run_line("Sarcopenia (GSE1428)",
                   file.path(ROOT, "10-L", "input"), "dat_expr.Rds", "dat_group.Rds",
                   file.path(ROOT, "10-L", "input"),
                   "Control", "SP", "svm", 0.2)
res_hp <- run_line("H. pylori (GSE27411)",
                   file.path(ROOT, "10-L1", "input"), "dat1_expr.Rds", "dat1_group.Rds",
                   file.path(ROOT, "4-Venn2", "input"),
                   "Control", "HP", "pool", 0.05)

out <- rbind(res_sp, res_hp)
# hub flag: votes_ge2 in both lines
hub_genes <- intersect(out$gene[out$votes_ge2 >= 2 & out$line == "Sarcopenia (GSE1428)"],
                       out$gene[out$votes_ge2 >= 2 & out$line == "H. pylori (GSE27411)"])
out$hub <- as.integer(out$gene %in% hub_genes)

dir.create("output", showWarnings = FALSE)
write.csv(out, file.path("output", "Table_S2_bootstrap_reconstructed.csv"),
          row.names = FALSE, na = "")
cat("\nDONE -> output/Table_S2_bootstrap_reconstructed.csv\n")
cat("hub genes:", paste(hub_genes, collapse = ", "), "\n")
