# 1-validation.R  GSE167186/GPL20301 bulk 肌少症验证（75 例）
# 四基因 C1QA / NFKBIA / DTL / CTNNBIP1
#   1) Sarcopenia vs Old Healthy（年龄匹配主对比）组间比较 + 小提琴图
#   2) 单基因 ROC（Sarcopenia vs Old Healthy；另存 vs 全部对照）
#   3) 表型相关性：握力 / SPPB / TUG / 6min walk（老年人群 Sarco+Old）
# 运行方式：在本模块目录内 Rscript 1-validation.R（相对路径，规避中文路径问题）

suppressMessages({
  library(data.table); library(ggplot2); library(pROC); library(patchwork)
})
dir.create("output", showWarnings = FALSE)
set.seed(2026)

HUB <- c("C1QA", "NFKBIA", "DTL", "CTNNBIP1")

## ---------- 1. 读入表达矩阵 ----------
cnt <- fread("raw/bulk_counts.csv.gz", header = TRUE)
mat <- as.matrix(cnt[, -1]); rownames(mat) <- as.character(cnt[[1]])
cat("counts matrix:", nrow(mat), "genes x", ncol(mat), "samples\n")

## ---------- 2. 解析 series matrix 表型 ----------
con <- gzfile("raw/GSE167186-GPL20301_series_matrix.txt.gz")
lines <- readLines(con)
close(con)
pick <- function(pfx) {
  ln <- lines[grepl("^!Sample_characteristics_ch1", lines)]
  rows <- lapply(ln, function(l) {
    v <- strsplit(l, "\t", fixed = TRUE)[[1]][-1]
    v <- gsub('^"|"$', "", v)
  })
  out <- rep(NA_character_, length(rows[[1]]))
  for (r in rows) {
    hit <- grep(pfx, r)
    if (length(hit)) {
      val <- sub(pfx, "", r[hit])
      val <- trimws(val)
      out[!is.na(out) & out != "" & out != val & FALSE] <- out[FALSE]  # no-op
      have <- is.na(out) | out == ""
      out[have] <- val[have]
    }
  }
  out
}
gsms <- gsub('^"|"$', "", strsplit(lines[grepl("^!Sample_geo_accession", lines)][1], "\t")[[1]][-1])
titles <- gsub('^"|"$', "", strsplit(lines[grepl("^!Sample_title", lines)][1], "\t")[[1]][-1])

meta <- data.frame(
  gsm   = gsms,
  title = titles,
  group = pick("^group: "),
  age   = suppressWarnings(as.numeric(pick("^age: "))),
  grip  = suppressWarnings(as.numeric(pick("grip strength \\(kg\\): "))),
  sppb  = suppressWarnings(as.numeric(pick("sppb: "))),
  tug   = suppressWarnings(as.numeric(pick("time up and go \\(s\\): "))),
  walk  = suppressWarnings(as.numeric(pick("6 min walk test \\(m/s\\): ")))
)
cat("group counts:\n"); print(table(meta$group))

## ---------- 3. CPM 归一化 ----------
lib <- colSums(mat)
cpm <- sweep(mat, 2, lib, "/") * 1e6
lcpm <- log2(cpm + 1)
expr <- lcpm[HUB, , drop = FALSE]
df <- data.frame(title = colnames(expr), t(expr), check.names = FALSE)
meta <- merge(meta, df, by = "title", sort = FALSE)

ana <- meta[!meta$group %in% "UNCLASSIFIED" & !is.na(meta$group), ]
ana$group <- factor(ana$group, levels = c("Young Healthy", "Old Healthy", "Sarcopenia"))
cat("analyzed n:\n"); print(table(ana$group))

gcol <- c("Young Healthy" = "#00A087", "Old Healthy" = "#4DBBD5", "Sarcopenia" = "#E64B35")

## ---------- 4. 组间比较（Sarcopenia vs Old Healthy 为主对比） ----------
stat <- do.call(rbind, lapply(HUB, function(g) {
  so <- ana$group %in% c("Old Healthy", "Sarcopenia")
  x1 <- ana[[g]][ana$group == "Sarcopenia"]; x0 <- ana[[g]][ana$group == "Old Healthy"]
  w <- wilcox.test(x1, x0)
  sy <- ana[[g]][ana$group == "Young Healthy"]
  wy <- wilcox.test(x1, sy)
  data.frame(gene = g,
             mean_Sarco = mean(x1), mean_Old = mean(x0), mean_Young = mean(sy),
             log2FC_Sarco_vs_Old = mean(x1) - mean(x0),
             p_Sarco_vs_Old = w$p.value, p_Sarco_vs_Young = wy$p.value)
}))
stat$p_Sarco_vs_Old  <- p.adjust(stat$p_Sarco_vs_Old,  "BH")
stat$p_Sarco_vs_Young <- p.adjust(stat$p_Sarco_vs_Young, "BH")
fwrite(stat, "output/group_comparison_GPL20301.csv")
cat("\n=== group comparison (BH-adjusted) ===\n"); print(stat, digits = 3)

## 小提琴图（4 基因，3 组）
plots <- lapply(HUB, function(g) {
  ggplot(ana, aes(x = group, y = .data[[g]], fill = group)) +
    geom_violin(alpha = .8, color = "grey30", linewidth = .3) +
    geom_boxplot(width = .12, fill = "white", outlier.size = .6, linewidth = .3) +
    scale_fill_manual(values = gcol) +
    labs(title = g, x = NULL, y = expression(log[2](CPM + 1))) +
    theme_bw(base_size = 11) + theme(legend.position = "none")
})
p <- wrap_plots(plots, nrow = 1) + plot_annotation("GPL20301 (GSE167186 bulk): hub gene expression")
ggsave("output/Violin_4genes_GPL20301.pdf", p, width = 11, height = 3.2)

## ---------- 5. ROC：Sarcopenia vs Old Healthy ----------
roc_tab <- do.call(rbind, lapply(HUB, function(g) {
  so <- ana[ana$group %in% c("Old Healthy", "Sarcopenia"), ]
  r <- roc(response = so$group == "Sarcopenia", predictor = so[[g]], quiet = TRUE)
  ci <- as.numeric(ci.auc(r))
  data.frame(gene = g, AUC = as.numeric(auc(r)), AUC_low = ci[1], AUC_high = ci[3],
             n_sarco = sum(so$group == "Sarcopenia"), n_old = sum(so$group == "Old Healthy"))
}))
fwrite(roc_tab, "output/ROC_GPL20301_Sarco_vs_Old.csv")
cat("\n=== ROC Sarcopenia vs Old Healthy ===\n"); print(roc_tab, digits = 3)

for (g in HUB) {
  so <- ana[ana$group %in% c("Old Healthy", "Sarcopenia"), ]
  r <- roc(response = so$group == "Sarcopenia", predictor = so[[g]], quiet = TRUE)
  auc.v <- auc(r)
  p <- ggroc(r, color = "#E64B35", linewidth = .8) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey50") +
    annotate("text", x = .25, y = .15,
             label = sprintf("AUC = %.3f (95%%CI %.3f-%.3f)", auc.v, roc_tab$AUC_low[roc_tab$gene == g], roc_tab$AUC_high[roc_tab$gene == g]),
             size = 3.4) +
    coord_equal() +
    labs(title = sprintf("%s: Sarcopenia vs Old Healthy (n=%d/%d)", g, roc_tab$n_sarco[roc_tab$gene == g], roc_tab$n_old[roc_tab$gene == g]),
         x = "1 - Specificity", y = "Sensitivity") +
    theme_bw(base_size = 11)
  ggsave(sprintf("output/ROC_%s_GPL20301.pdf", g), p, width = 4.6, height = 4.6)
}

## ROC vs 全部对照（Sarco vs Old+Young）
roc_all <- do.call(rbind, lapply(HUB, function(g) {
  r <- roc(response = ana$group == "Sarcopenia", predictor = ana[[g]], quiet = TRUE)
  ci <- as.numeric(ci.auc(r))
  data.frame(gene = g, AUC = as.numeric(auc(r)), AUC_low = ci[1], AUC_high = ci[3])
}))
fwrite(roc_all, "output/ROC_GPL20301_Sarco_vs_AllControls.csv")

## ROC：老年（Sarco+Old）vs Young —— 对应 GSE25941 的 young-vs-old 设计
oldvyoung <- do.call(rbind, lapply(HUB, function(g) {
  o <- ana[ana$group %in% c("Old Healthy", "Sarcopenia", "Young Healthy"), ]
  r <- roc(response = o$group == "Young Healthy", predictor = o[[g]], quiet = TRUE)
  ci <- as.numeric(ci.auc(r))
  data.frame(gene = g, AUC = as.numeric(auc(r)), AUC_low = ci[1], AUC_high = ci[3])
}))
fwrite(oldvyoung, "output/ROC_GPL20301_Older_vs_Young.csv")
cat("\n=== ROC Older(Sarco+Old) vs Young ===\n"); print(oldvyoung, digits = 3)

for (g in HUB) {
  o <- ana[ana$group %in% c("Old Healthy", "Sarcopenia", "Young Healthy"), ]
  r <- roc(response = o$group == "Young Healthy", predictor = o[[g]], quiet = TRUE)
  auc.v <- auc(r); ci <- as.numeric(ci.auc(r))
  p <- ggroc(r, color = "#E64B35", linewidth = .8) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey50") +
    annotate("text", x = .25, y = .15,
             label = sprintf("AUC = %.3f (95%%CI %.3f-%.3f)", auc.v, ci[1], ci[3]), size = 3.4) +
    coord_equal() +
    labs(title = sprintf("%s: Older (Sarcopenia+Old) vs Young (n=%d/%d)", g,
                         sum(o$group != "Young Healthy"), sum(o$group == "Young Healthy")),
         x = "1 - Specificity", y = "Sensitivity") +
    theme_bw(base_size = 11)
  ggsave(sprintf("output/ROC_%s_GPL20301_OlderVsYoung.pdf", g), p, width = 4.6, height = 4.6)
}

## ---------- 6. 表型相关性（老年人群：Sarco + Old） ----------
old <- ana[ana$group %in% c("Old Healthy", "Sarcopenia"), ]
pheno <- c(grip = "Grip strength (kg)", sppb = "SPPB", tug = "TUG (s)", walk = "6-min walk speed (m/s)")
cor_tab <- do.call(rbind, lapply(HUB, function(g) do.call(rbind, lapply(names(pheno), function(p) {
  ok <- !is.na(old[[g]]) & !is.na(old[[p]])
  ct <- suppressWarnings(cor.test(old[[g]][ok], old[[p]][ok], method = "spearman"))
  data.frame(gene = g, phenotype = pheno[p], n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}))))
cor_tab$n <- as.integer(cor_tab$n)
cor_tab$p_adj <- p.adjust(cor_tab$p, "BH")
fwrite(cor_tab, "output/phenotype_correlation_GPL20301.csv")
cat("\n=== phenotype correlation (older adults, Spearman) ===\n"); print(cor_tab, digits = 3)

## C1QA × 四表型散点面板（其余基因按需改）
sc <- lapply(names(pheno), function(p) {
  ok <- !is.na(old$C1QA) & !is.na(old[[p]])
  ggplot(old[ok, ], aes(x = .data[[p]], y = C1QA)) +
    geom_point(size = 1.6, color = "#4DBBD5", alpha = .85) +
    geom_smooth(method = "lm", se = TRUE, color = "#E64B35", fill = "#E64B35", alpha = .15, linewidth = .5) +
    labs(title = pheno[p], x = NULL, y = expression(C1QA~log[2](CPM+1))) +
    theme_bw(base_size = 10) + theme(plot.title = element_text(size = 10))
})
p2 <- wrap_plots(sc, nrow = 1) + plot_annotation("C1QA vs physical performance (older adults, GPL20301)")
ggsave("output/Scatter_C1QA_phenotype_GPL20301.pdf", p2, width = 11, height = 3)

saveRDS(ana, "output/GPL20301_meta_expr.Rds")
cat("\nDONE. outputs in ./output\n")
