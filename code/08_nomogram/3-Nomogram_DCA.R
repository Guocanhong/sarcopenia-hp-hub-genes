# =====================================================================
# Nomogram + 校准曲线 + DCA
# 数据: 9-Lo SP 数据集 (GSE1428, Control 10 / SP 12)
# 模型: 多因素 Logistic (NFKBIA + CTNNBIP1)，参考 9-Lo VIF 筛选结果
# 输出: Fig10A-Nomogram.pdf / Fig10B-Calibration.pdf / Fig10C-DCA.pdf
# =====================================================================
rm(list = ls())
library(rms)
library(ggplot2)

dat_expr  <- readRDS('input/dat_expr.Rds')
dat_group <- readRDS('input/dat_group.Rds')

genes <- c('C1QA', 'NFKBIA', 'DTL', 'CTNNBIP1')
stopifnot(all(genes %in% rownames(dat_expr)))

# 整理建模数据 ---------------------------------------------------------------
dat <- as.data.frame(t(dat_expr[genes, ]))
dat$group <- factor(dat_group$group, levels = c('Control', 'SP'))
# ggDCA 要求数值型结局 (0=Control, 1=SP)
dat$y <- as.numeric(dat$group == 'SP')
stopifnot(!any(is.na(dat)), table(dat$group)[['Control']] == 10)

dd <- datadist(dat)
options(datadist = 'dd')

# 主模型：多因素保留的 2 基因 (NFKBIA + CTNNBIP1)
# 注：4 基因全模型存在完全分离，lrm 拟合失败，故不纳入
fit2 <- lrm(y ~ NFKBIA + CTNNBIP1, data = dat, x = TRUE, y = TRUE)
print(fit2)

# ---- Fig10A 列线图 ----------------------------------------------------------
pdf(file = 'output/Fig10A-Nomogram.pdf', width = 8, height = 6, onefile = FALSE)
nom <- nomogram(
  fit2,
  fun = plogis,
  funlabel = 'Risk of Sarcopenia',
  fun.at = c(0.05, 0.25, 0.5, 0.75, 0.95),
  lp = FALSE
)
plot(nom, xfrac = 0.45, cex.axis = 0.75, cex.var = 0.9)
title(main = 'Sarcopenia Risk Nomogram (NFKBIA + CTNNBIP1)')
dev.off()

# ---- Fig10B 校准曲线 (Bootstrap, B=1000) ------------------------------------
pdf(file = 'output/Fig10B-Calibration.pdf', width = 6, height = 6, onefile = FALSE)
cal <- calibrate(fit2, method = 'boot', B = 1000)
plot(
  cal,
  xlab = 'Predicted Probability of Sarcopenia',
  ylab = 'Observed Probability of Sarcopenia',
  main = 'Calibration Curve (Bootstrap B = 1000)',
  subtitles = FALSE
)
dev.off()

# 校准指标由 plot(cal) 自带输出（Mean absolute error）

# ---- Fig10C 决策曲线分析 DCA ------------------------------------------------
# 手动计算净获益率 (Net Benefit)，不依赖 ggDCA/rmda
pred <- as.numeric(predict(fit2, type = 'fitted'))
y    <- dat$y
pts  <- seq(0.01, 0.99, 0.01)
nb_model <- sapply(pts, function(pt) {
  tp <- mean(pred >= pt & y == 1)
  fp <- mean(pred >= pt & y == 0)
  tp - fp * pt / (1 - pt)
})
nb_all  <- sapply(pts, function(pt) mean(y == 1) - mean(y == 0) * pt / (1 - pt))
nb_none <- rep(0, length(pts))

dat_dca <- rbind(
  data.frame(threshold = pts, NB = nb_model, Model = '2-gene model (NFKBIA+CTNNBIP1)'),
  data.frame(threshold = pts, NB = nb_all,   Model = 'Treat All'),
  data.frame(threshold = pts, NB = nb_none,  Model = 'Treat None')
)
p_dca <- ggplot(dat_dca, aes(threshold, NB, color = Model)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c('#E64B35', '#3C5488', 'grey50')) +
  coord_cartesian(xlim = c(0, 0.8), ylim = c(-0.05, 0.6)) +
  theme_bw() +
  theme(legend.position = c(0.7, 0.75), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(
    title = 'Decision Curve Analysis',
    x = 'Threshold Probability',
    y = 'Net Benefit'
  )
ggsave(file = 'output/Fig10C-DCA.pdf', p_dca, width = 7, height = 6)

# ---- 模型参数保存 ------------------------------------------------------------
saveRDS(fit2, file = 'output/lrm_2gene_SP.Rds')
write.csv(
  data.frame(
    Gene = names(coef(fit2)),
    Coefficient = as.numeric(coef(fit2)),
    OR = exp(as.numeric(coef(fit2)))
  ),
  file = 'output/nomogram_coef_SP.csv', row.names = FALSE
)

# 模型内部表现
cat('\n2基因模型  AUC(Dxy推导) =', round((fit2$stats['Dxy'] / 2 + 0.5), 3), '\n')
cat('完成：Fig10A-Nomogram.pdf / Fig10B-Calibration.pdf / Fig10C-DCA.pdf\n')
