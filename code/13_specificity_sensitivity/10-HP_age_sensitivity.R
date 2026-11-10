# 10-HP_age_sensitivity.R
# R2-M2 (C4) HP-line age-confounding sensitivity analysis
# GSE5081: donor ages/sex recovered from GEO series_matrix characteristics_ch1
# Goal: test whether HP+ vs HP- gastric signal can be explained by donor age
# Output: 60-Review/Medicine_Submission_Package/HP_age_sensitivity_GSE5081.csv (+ console summary)

# run from project root; no setwd (Chinese path breaks Rscript locale)
options(stringsAsFactors = FALSE)

sm <- gzfile("8-HP/input/GSE5081_series_matrix.txt.gz", "rt")
lines <- readLines(sm, n = 200)
close(sm)

get_row <- function(key) {
  ln <- grep(paste0("^!", key), lines, value = TRUE)[1]
  strsplit(ln, "\t")[[1]][-1] |> gsub('"', '', x = _)
}
gsm     <- get_row("Sample_geo_accession")
char_rows <- lines[grepl("^!Sample_characteristics_ch1", lines)]
char_list <- lapply(char_rows, function(l) gsub('"', '', strsplit(l, "\t")[[1]][-1]))
names(char_list) <- sapply(char_list, function(x) sub(":.*", "", x[1]))
age_v    <- as.numeric(sub("Age: ", "", char_list[["Age"]]))
sex_v    <- sub("Gender: ", "", char_list[["Gender"]])
tissue_v <- sub("Tissue: ", "", char_list[["Tissue"]])

meta <- data.frame(gsm = gsm, age = age_v, sex = sex_v, tissue = tissue_v)
stopifnot(nrow(meta) == 32)

grp <- readRDS("8-HP/input/dat_group.Rds")
expr <- readRDS("8-HP/input/dat_expr.Rds")
grp_v <- grp$group
names(grp_v) <- grp$sample

meta$group <- grp_v[match(meta$gsm, names(grp_v))]
stopifnot(!any(is.na(meta$group)))
cat("== GSE5081 sample metadata (n=32) ==\n"); print(meta)

# donor pairing check: erosive vs non-erosive same age within HP status
chk <- aggregate(age ~ group, data = meta[grepl("erosive", meta$tissue), ], FUN = function(x) paste(x, collapse = ","))
cat("\nAge by tissue-site group:\n"); print(chk)

# donor level: one row per (group, age) pair collapsed
donors <- data.frame(
  donor_id = 1:16,
  group    = c(rep("Control", 8), rep("HP", 8)),
  age      = c(meta$age[meta$tissue == "HP- erosive gastritis, site of erosion"],
               meta$age[meta$tissue == "HP+ erosive gastritis, site of erosion"]),
  sex      = c(meta$sex[meta$tissue == "HP- erosive gastritis, site of erosion"],
               meta$sex[meta$tissue == "HP+ erosive gastritis, site of erosion"])
)
cat("\n== Donor-level age/sex (8 HP- vs 8 HP+) ==\n"); print(donors)
wt_age <- wilcox.test(age ~ group, data = donors)
cat(sprintf("Age: HP- mean %.1f (SD %.1f), HP+ mean %.1f (SD %.1f), Wilcoxon P = %.3f\n",
    mean(donors$age[donors$group == "Control"]), sd(donors$age[donors$group == "Control"]),
    mean(donors$age[donors$group == "HP"]),      sd(donors$age[donors$group == "HP"]),
    wt_age$p.value))
ft <- fisher.test(table(donors$group, donors$sex))
cat(sprintf("Sex: HP- %dF/%dM, HP+ %dF/%dM, Fisher P = %.3f\n",
    sum(donors$sex == "female" & donors$group == "Control"), sum(donors$sex == "male" & donors$group == "Control"),
    sum(donors$sex == "female" & donors$group == "HP"),      sum(donors$sex == "male" & donors$group == "HP"),
    ft$p.value))

# hub genes: verify against Key_gene.Rds
kg <- readRDS("8-HP/input/Key_gene.Rds")
cat("\nKey_gene object:\n"); print(kg)
hub <- c("C1QA", "NFKBIA", "DTL", "CTNNBIP1")
cat("rownames class:", class(expr[1:5, 1]), " sample probe rows:", sum(rownames(expr) %in% hub), "\n")

get_expr <- function(g) as.numeric(expr[rownames(expr) == g, meta$gsm])

out <- data.frame()
for (g in hub) {
  x <- get_expr(g)
  if (length(x) == 0 || all(is.na(x))) { cat(g, ": NOT FOUND in expr\n"); next }
  # biopsy level (n=32), matches Table 3
  p_nom <- wilcox.test(x ~ meta$group)$p.value
  auc_nom <- suppressMessages(pROC::roc(meta$group, x, quiet = TRUE,
                    levels = c("Control", "HP"), direction = "auto"))$auc
  sp_age <- suppressWarnings(cor(x, meta$age, method = "spearman"))
  p_agecorr <- suppressWarnings(cor.test(x, meta$age, method = "spearman"))$p.value
  # biopsy-level logistic: gene + age
  df32 <- data.frame(y = as.numeric(meta$group == "HP"), gene = scale(x), age = scale(meta$age))
  m_gene <- glm(y ~ gene, data = df32, family = binomial)
  m_age  <- glm(y ~ age,  data = df32, family = binomial)
  m_both <- glm(y ~ gene + age, data = df32, family = binomial)
  lrt <- anova(m_age, m_both, test = "LRT")
  auc_gene <- as.numeric(pROC::roc(df32$y, fitted(m_gene), quiet = TRUE)$auc)
  auc_age  <- as.numeric(pROC::roc(df32$y, fitted(m_age),  quiet = TRUE)$auc)
  auc_both <- as.numeric(pROC::roc(df32$y, fitted(m_both), quiet = TRUE)$auc)
  p_gene_adj <- summary(m_both)$coefficients["gene", "Pr(>|z|)"]
  p_age_adj  <- summary(m_both)$coefficients["age", "Pr(>|z|)"]
  # donor level (n=16): mean of erosive + non-erosive biopsies per donor
  # meta rows: 1-8 HP- erosive, 9-16 HP+ erosive, 17-24 HP+ non-erosive, 25-32 HP- non-erosive
  xdonor <- c((x[1:8] + x[25:32]) / 2, (x[9:16] + x[17:24]) / 2)
  df16 <- data.frame(y = as.numeric(donors$group == "HP"), gene = scale(xdonor), age = scale(donors$age))
  m16 <- glm(y ~ gene + age, data = df16, family = binomial)
  m16a <- glm(y ~ age, data = df16, family = binomial)
  p16 <- tryCatch(summary(m16)$coefficients["gene", "Pr(>|z|)"], error = function(e) NA)
  lrt16 <- tryCatch(anova(m16a, m16, test = "LRT")$`Pr(>Chi)`[2], error = function(e) NA)
  p16_nom <- wilcox.test(xdonor ~ donors$group)$p.value
  out <- rbind(out, data.frame(
    Gene = g,
    AUC_biopsy_nominal = round(auc_nom, 3),
    P_wilcox_biopsy = signif(p_nom, 3),
    Spearman_vs_age = round(sp_age, 3),
    P_spearman_age = signif(p_agecorr, 3),
    AUC_age_only = round(auc_age, 3),
    AUC_gene_only = round(auc_gene, 3),
    AUC_gene_plus_age = round(auc_both, 3),
    P_gene_age_adjusted_biopsy = signif(p_gene_adj, 3),
    P_LRT_gene_beyond_age = signif(lrt$`Pr(>Chi)`[2], 3),
    P_wilcox_donor = signif(p16_nom, 3),
    P_gene_age_adjusted_donor = signif(p16, 3)
  ))
}
cat("\n== Hub-gene age sensitivity summary (GSE5081) ==\n"); print(out)
dir.create("60-Review/Medicine_Submission_Package", showWarnings = FALSE)
write.csv(out, "60-Review/Medicine_Submission_Package/HP_age_sensitivity_GSE5081.csv", row.names = FALSE)
write.csv(donors, "60-Review/Medicine_Submission_Package/GSE5081_donor_age_sex.csv", row.names = FALSE)
cat("\nSaved: HP_age_sensitivity_GSE5081.csv, GSE5081_donor_age_sex.csv\n")
