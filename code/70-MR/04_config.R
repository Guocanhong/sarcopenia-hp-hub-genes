# 04_config.R -- shared config for all MR scripts (ASCII only!)
# OpenGWAS JWT is stored in data/token.txt (from api.opengwas.io profile)
MR_DIR  <- getwd()
DATA    <- file.path(MR_DIR, "data")
# data.table::fwrite cannot create files under a non-ASCII (Chinese) path on this
# Windows box, so all outputs go to a pure-ASCII directory instead.
OUT     <- "D:/MRout/70-MR"
dir.create(DATA, showWarnings = FALSE); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

TOKEN_FILE <- file.path(DATA, "token.txt")
OPENGWAS_JWT <- if (file.exists(TOKEN_FILE)) trimws(readLines(TOKEN_FILE, warn = FALSE)[1]) else ""
if (nzchar(OPENGWAS_JWT)) Sys.setenv(OPENGWAS_JWT = OPENGWAS_JWT)

# Four hub genes -- versioned GENCODE ids verified via GTEx /api/v2/reference/gene
HUB_GENES <- list(
  C1QA     = "ENSG00000173372.16",
  NFKBIA   = "ENSG00000100906.10",
  DTL      = "ENSG00000143476.17",
  CTNNBIP1 = "ENSG00000178585.14"
)

# OpenGWAS dataset ids (verified via gwasinfo + literature)
DATASETS <- list(
  # exposures (Hp)
  exposure_hp_selfreport = "ukb-b-531",            # H. pylori self-report, N=462,933
  exposure_hp_igg_pos    = "ebi-a-GCST90006910",   # anti-Hp IgG seropositivity, N=8,735
  exposure_hp_igg_level  = "ieu-b-4905",           # anti-Hp IgG levels, N=4,683
  # outcomes (sarcopenia traits)
  outcome_grip_right  = "ukb-b-10215",             # right hand grip, N=461,089
  outcome_grip_left   = "ukb-b-7478",              # left hand grip, N=461,026
  outcome_alm         = "ebi-a-GCST90000025",      # ALM, N=450,243
  outcome_walk        = "ukb-b-4711",              # usual walking pace, N=459,915
  outcome_lowgrip     = "ebi-a-GCST90007526"       # low grip strength 60+ (EWGSOP), binary
)

# Analysis parameters
P_EXP     <- 5e-8    # IV threshold; relax to 5e-6 for Hp if too few SNPs
CLUMP_R2  <- 0.001
CLUMP_KB  <- 10000
FDR_EQTL  <- 5e-8    # eQTL IV threshold
