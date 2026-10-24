# HP-line donor-level sensitivity analysis (GSE27411)
# 12 biopsies from 6 donors (antrum+corpus each) were treated as independent
# samples in the primary screen. Here we quantify how the HP-line results
# behave when within-donor non-independence is modelled (duplicateCorrelation)
# or removed (donor-level aggregation to 3 vs 3).
suppressMessages({library(limma); library(pROC)})

expr <- readRDS("../1-limma2/input/dat1_expr.Rds")          # 18190 x 12, log2
grp  <- readRDS("../1-limma2/input/dat1_group.Rds")
pool <- read.csv("../2-Venn/output/SP-HP.CSV")$x            # 118 common DEGs
orig <- readRDS("../1-limma2/output/HP-diff_res.Rds")       # biopsy-level limma

# Sample annotation verified against the GSE27411 series matrix
# (titles carry donor IDs: 23A/23C -> donor 23, etc.)
ann <- data.frame(
  gsm   = c("GSM677730","GSM677731","GSM677732","GSM677733","GSM677734","GSM677735",
            "GSM677736","GSM677737","GSM677738","GSM677739","GSM677746","GSM677747"),
  donor = c("D23","D49","D23","D49","D8","D9","D30","D8","D9","D30","D35","D35"),
  site  = c("antrum","antrum","corpus","corpus","antrum","antrum","antrum",
            "corpus","corpus","corpus","antrum","corpus"),
  group = c("Control","Control","Control","Control","HP","HP","HP",
            "HP","HP","HP","Control","Control"),
  stringsAsFactors = FALSE)
stopifnot(setequal(ann$gsm, colnames(expr)), all(ann$group == grp$group[match(ann$gsm, grp$sample)]))
expr <- expr[, ann$gsm]
hubs <- c("C1QA","NFKBIA","DTL","CTNNBIP1")
cat("group x site table:\n"); print(table(ann$group, ann$site))

## Model A: keep all 12 biopsies, block on donor, adjust for site ------------
design <- model.matrix(~ ann$group + ann$site)   # intercept; HP coef = HP vs Control
colnames(design) <- c("Intercept", "HPvsControl", "sitecorpus")
corfit <- duplicateCorrelation(expr, design, block = ann$donor)
cat("consensus within-donor correlation:", round(corfit$consensus, 3), "\n")
fitA <- eBayes(lmFit(expr, design, block = ann$donor, correlation = corfit$consensus))
ttA <- na.omit(topTable(fitA, coef = "HPvsControl", number = Inf))

## Model B: aggregate to donor level (mean of antrum+corpus), 3 vs 3 ----------
agg <- sapply(split(seq_len(ncol(expr)), ann$donor), function(ix) rowMeans(expr[, ix, drop = FALSE]))
dgrp <- ann$group[match(colnames(agg), ann$donor)]
designB <- model.matrix(~ factor(dgrp, levels = c("Control", "HP")))  # intercept
colnames(designB) <- c("Intercept", "HPvsControl")
fitB <- eBayes(lmFit(agg, designB))
ttB <- na.omit(topTable(fitB, coef = "HPvsControl", number = Inf))

## Hub genes under all three views --------------------------------------------
hub_tab <- do.call(rbind, lapply(hubs, function(g) {
  data.frame(gene = g,
    biopsy_logFC = round(orig[g, "logFC"], 2), biopsy_P = signif(orig[g, "P.Value"], 2),
    donorCorr_logFC = round(ttA[g, "logFC"], 2), donorCorr_P = signif(ttA[g, "P.Value"], 2),
    donor3v3_logFC = round(ttB[g, "logFC"], 2), donor3v3_P = signif(ttB[g, "P.Value"], 2),
    donor3v3_AUC = round(as.numeric(pROC::roc(response = dgrp, predictor = agg[g, ],
        levels = c("Control", "HP"), direction = "<", quiet = TRUE)$auc), 3))
}))
print(hub_tab, row.names = FALSE)
write.csv(hub_tab, "output/HP_donor_sensitivity_hub_genes.csv", row.names = FALSE)

## 118-pool survival and direction concordance --------------------------------
pool_tab <- do.call(rbind, lapply(c(`biopsy-level` = "orig", `donor-corrected` = "ttA",
                                    `donor-3v3` = "ttB"), function(nm) {
  tt <- get(nm)
  sig <- rownames(tt)[tt$P.Value < 0.05]
  kept <- pool[pool %in% rownames(tt)]
  surv <- kept[kept %in% sig]
  d_orig <- sign(orig[kept, "logFC"]); d_new <- sign(tt[kept, "logFC"])
  p_hg <- phyper(length(surv) - 1, length(sig), nrow(tt) - length(sig), length(kept),
                 lower.tail = FALSE)
  data.frame(model = nm, n_sig = length(sig), pool_tested = length(kept),
             pool_survived = length(surv),
             fold_enrich = round(length(surv) / length(kept) / (length(sig) / nrow(tt)), 2),
             hypergeom_P = signif(p_hg, 2),
             dir_concordance = round(mean(d_orig == d_new), 2))
}))
print(pool_tab, row.names = FALSE)
write.csv(pool_tab, "output/HP_donor_sensitivity_pool118.csv", row.names = FALSE)
cat("[OK] done\n")
