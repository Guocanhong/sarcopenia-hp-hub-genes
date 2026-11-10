# 17_figures.R -- manuscript figures
# Fig 1: forward MR forest (3 Hp exposures x 5 outcomes, IVW RE)  -> Fig1_forward_forest.png
# Fig 4: antigen-spectrum forest (6 antigens x 4 outcomes)        -> Fig4_antigen_forest.png
# Fig 5: NFKBIA -> ALM scatter + forest (rebuild harmonised)      -> Fig5a/Fig5b png
# Fig 6: reverse MR forest (3 traits x 2 Hp outcomes)             -> Fig6_reverse_forest.png
# Fig 2/3 reuse existing line1 scatter/loo PNGs (renamed via shell afterwards).
# English comments only (GBK pitfall). Output to D:/MRout/70-MR (ASCII path).
source(file.path(getwd(), "code", "04_config.R"))
suppressMessages({library(TwoSampleMR); library(data.table); library(ggplot2)})
stopifnot("OPENGWAS_JWT missing" = nzchar(Sys.getenv("OPENGWAS_JWT")))

exp_lab <- c(igg_pos   = "Anti-Hp IgG seropositivity",
             selfreport = "Self-reported Hp infection",
             igg_level = "Anti-Hp IgG titer (ALSPAC children)",
             NFKBIA    = "NFKBIA blood expression",
             C1QA      = "C1QA blood expression",
             CTNNBIP1  = "CTNNBIP1 blood expression")
out_lab <- c(alm = "Appendicular lean mass",
             grip_right = "Right grip strength",
             grip_left  = "Left grip strength",
             walking    = "Walking pace",
             lowgrip60  = "Low grip strength (>=60y)",
             igg_pos    = "Hp seropositivity",
             selfreport = "Self-reported Hp infection")

# ---------- helper: custom forest ----------
forest_df <- function(rows, bonf_keys = character(0)) {
  d <- data.frame(
    exp = rows$exposure_trait, out = rows$outcome_trait,
    b = as.numeric(rows$b), se = as.numeric(rows$se), p = as.numeric(rows$pval),
    stringsAsFactors = FALSE)
  d$exp <- gsub(" \\|\\| id:.*", "", d$exp); d$out <- gsub(" \\|\\| id:.*", "", d$out)
  d$exp <- gsub("blood_cis_eQTL", "blood expression", d$exp)
  d$ci_lo <- d$b - 1.96 * d$se; d$ci_hi <- d$b + 1.96 * d$se
  d$sig <- d$p < 0.05
  d$bonf <- paste(d$exp, d$out) %in% bonf_keys
  d
}
draw_forest <- function(d, title, ncol = 3) {
  d$exp <- factor(d$exp, levels = rev(unique(d$exp)))
  d$out <- factor(d$out, levels = unique(d$out))
  has_bonf <- any(d$bonf)
  p <- ggplot(d, aes(x = b, y = exp)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.18, colour = "grey25") +
    geom_point(aes(shape = sig, fill = sig), size = 2.8, colour = "black") +
    scale_shape_manual(values = c(`TRUE` = 22, `FALSE` = 21),
                       labels = c(`TRUE` = "P < 0.05", `FALSE` = "NS"), name = "") +
    scale_fill_manual(values = c(`TRUE` = "#2C3E50", `FALSE` = "white"), guide = "none") +
    geom_text(aes(label = sprintf("P = %.1e", p)), size = 2.3, vjust = -1.4, colour = "grey20") +
    facet_wrap(~out, ncol = ncol, scales = "free_x") +
    scale_x_continuous(expand = expansion(mult = c(0.18, 0.18))) +
    labs(x = "IVW random-effects estimate (per SD / log-odds, 95% CI)",
         y = NULL, title = title) +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "grey90", colour = NA),
          panel.grid.minor = element_blank(),
          legend.position = "top")
  if (has_bonf) {
    # highlight Bonferroni-passed point in red on top
    db <- d[d$bonf, , drop = FALSE]
    p <- p + geom_point(data = db, aes(x = b, y = exp), shape = 22, size = 3.2,
                        fill = "#C0392B", colour = "black", inherit.aes = FALSE)
  }
  p
}

# ---------- Fig 1: forward line ----------
l1 <- fread(file.path(OUT, "line1_mr_results.csv"))
l1 <- l1[grepl("multiplicative random effects", l1$method, fixed = TRUE)]
# map exposure_trait -> short label via id.exposure
l1[, exp_short := fcase(id.exposure == DATASETS$exposure_hp_igg_pos, "igg_pos",
                        id.exposure == DATASETS$exposure_hp_selfreport, "selfreport",
                        id.exposure == DATASETS$exposure_hp_igg_level, "igg_level",
                        default = id.exposure)]
l1[, out_short := fcase(id.outcome == DATASETS$outcome_alm, "alm",
                        id.outcome == DATASETS$outcome_grip_right, "grip_right",
                        id.outcome == DATASETS$outcome_grip_left, "grip_left",
                        id.outcome == DATASETS$outcome_walk, "walking",
                        id.outcome == DATASETS$outcome_lowgrip, "lowgrip60",
                        default = id.outcome)]
l1$exposure_trait <- exp_lab[l1$exp_short]
l1$outcome_trait  <- out_lab[l1$out_short]
d1 <- forest_df(l1)
g1 <- draw_forest(d1, "Forward MR: H. pylori infection -> sarcopenia-related traits (IVW-RE)", ncol = 3)
ggsave(file.path(OUT, "Fig1_forward_forest.png"), g1, width = 11, height = 7, dpi = 300)

# ---------- Fig 4: antigen line ----------
l1c <- fread(file.path(OUT, "line1c_antigen_mr_results.csv"))
l1c <- l1c[grepl("multiplicative random effects", l1c$method, fixed = TRUE)]
l1c$exposure_trait <- gsub("Helicobacter pylori ", "", l1c$exposure_trait)
l1c$exposure_trait <- gsub(" antibody levels", "", l1c$exposure_trait)
l1c$outcome_trait <- c(alm = "Appendicular lean mass", grip_right = "Right grip strength",
                       lowgrip = "Low grip strength (>=60y)", walking = "Walking pace")[l1c$outcome_trait]
d4 <- forest_df(l1c, bonf_keys = paste("GroEL", "Appendicular lean mass"))
g4 <- draw_forest(d4, "Antigen-spectrum MR: anti-H. pylori antibodies -> sarcopenia-related traits (IVW-RE)", ncol = 2)
ggsave(file.path(OUT, "Fig4_antigen_forest.png"), g4, width = 10, height = 8, dpi = 300)

# ---------- Fig 6: reverse line ----------
l3  <- fread(file.path(OUT, "line3_reverse_mr_results.csv"))
l3b <- fread(file.path(OUT, "line3b_reverse_mr_results.csv"))
l6 <- rbind(l3, l3b, fill = TRUE)
l6 <- l6[grepl("multiplicative random effects", l6$method, fixed = TRUE)]
l6 <- l6[!duplicated(l6[, .(id.exposure, id.outcome, b)])]
l6$exposure_trait <- c(grip_right = "Right grip strength", alm = "Appendicular lean mass",
                       walking = "Walking pace")[l6$exposure_trait]
l6$outcome_trait  <- c(selfreport = "Self-reported Hp infection", igg_pos = "Hp seropositivity")[l6$outcome_trait]
d6 <- forest_df(l6)
g6 <- draw_forest(d6, "Reverse MR: sarcopenia-related traits -> H. pylori infection (IVW-RE)", ncol = 2)
ggsave(file.path(OUT, "Fig6_reverse_forest.png"), g6, width = 10, height = 6, dpi = 300)

# ---------- Fig 5: NFKBIA -> ALM (scatter + forest) ----------
# rebuild harmonised data with the same parameters as 14_nfkbia_heterogeneity.R
har_file <- file.path(OUT, "line2_harmonised_NFKBIA_to_alm.csv")
if (file.exists(har_file)) {
  dat <- fread(har_file)
  cat("NFKBIA harmonised reused from frozen file:", nrow(dat), "\n")
} else {
  dsid <- paste0("eqtl-a-", "ENSG00000100906")
  inst <- extract_instruments(dsid, p1 = FDR_EQTL, clump = TRUE, r2 = CLUMP_R2, kb = CLUMP_KB)
  inst$exposure <- "NFKBIA_blood_cis_eQTL"; inst$id.exposure <- "NFKBIA"
  cat("NFKBIA IVs rebuilt:", nrow(inst), "\n")
  out_dat <- extract_outcome_data(snps = inst$SNP, outcomes = DATASETS$outcome_alm, proxies = TRUE)
  dat <- harmonise_data(inst, out_dat, action = 2)
  cat("harmonised:", nrow(dat), "\n")
  fwrite(as.data.frame(dat), har_file)
}
m <- mr(dat)
print(m[, c("method", "nsnp", "b", "se", "pval")])
g5a <- mr_scatter_plot(m, dat)[[1]] +
  ggtitle("NFKBIA blood expression -> appendicular lean mass (scatter)")
ggsave(file.path(OUT, "Fig5a_NFKBIA_ALM_scatter.png"), g5a, width = 8, height = 6, dpi = 300)
g5b <- mr_forest_plot(mr_singlesnp(dat))[[1]] +
  ggtitle("NFKBIA blood expression -> appendicular lean mass (per-SNP + IVW)")
ggsave(file.path(OUT, "Fig5b_NFKBIA_ALM_forest.png"), g5b, width = 8, height = 7, dpi = 300)

cat("\n=== figures done ===\n")
