# res_scan_20261012.R  Resolution sweep verification for Supp Note 2 / Methods 2.10
# Loads the archived probe_seurat.rds (harmony + neighbours already computed)
# and re-runs Louvain clustering at resolutions 0.2-0.6. Read-only; rds not re-saved.
suppressMessages({library(Seurat)})
obj <- readRDS("probe_seurat.rds")
cat("cells:", ncol(obj), "\n")
for (r in c(0.2, 0.3, 0.4, 0.5, 0.6)) {
  o <- FindClusters(obj, resolution = r, verbose = FALSE)
  cat("resolution", r, "-> clusters:", nlevels(Idents(o)), "\n")
}
cat("RES SCAN DONE\n")
