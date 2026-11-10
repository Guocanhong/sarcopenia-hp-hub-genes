# 06b_plots.R -- make scatter/forest/LOO plots from harmonised CSVs (igg_pos pairs)
source(file.path(getwd(), "code", "04_config.R"))
suppressMessages({library(TwoSampleMR); library(data.table); library(ggplot2)})
files <- list.files(OUT, pattern = "^line1_harmonised_igg_pos_.*\\.csv$", full.names = TRUE)
for (f in files) {
  nm <- sub("^line1_harmonised_igg_pos_(.*)\\.csv$", "\\1", basename(f))
  dat <- as.data.frame(fread(f))
  tryCatch({
    p1 <- mr_scatter_plot(mr(dat), dat)
    ggsave(file.path(OUT, paste0("line1_scatter_igg_pos_", nm, ".png")), p1[[1]],
           width = 9, height = 7, dpi = 300)
    p4 <- mr_forest_plot(mr_singlesnp(dat))
    ggsave(file.path(OUT, paste0("line1_forest_igg_pos_", nm, ".png")), p4[[1]],
           width = 9, height = 7, dpi = 300)
    p3 <- mr_leaveoneout_plot(mr_leaveoneout(dat))
    ggsave(file.path(OUT, paste0("line1_loo_igg_pos_", nm, ".png")), p3[[1]],
           width = 9, height = 7, dpi = 300)
    cat(nm, ": plots saved\n")
  }, error = function(e) cat(nm, "plot fail:", conditionMessage(e), "\n"))
}
