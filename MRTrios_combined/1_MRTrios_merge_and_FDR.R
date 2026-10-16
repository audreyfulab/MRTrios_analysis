#### MRTrios: Step 1 (merge output shards) + Step 2 (BH FDR & Storey q-value)
#### One config table + two functions, looped over all cancer types.

library(data.table)
if (!requireNamespace("qvalue", quietly = TRUE)) BiocManager::install("qvalue")
library(qvalue)
library(MRGN)

# ════════════════════════════════════════════════════════════════
# Settings
# ════════════════════════════════════════════════════════════════
env       <- "local"          # "ondemand" or "local"
run_merge <- TRUE                # Step 1
run_fdr   <- TRUE                # Step 2
alpha     <- 0.05

ondemand_root <- "/wsu/home/hb/hb68/hb6890/fulab/MRTrios"
local_root    <- "/Users/lianzuo/LZ/ResearchProject/Fulab"

# OnDemand: everything lives under MRTrios/
# Local:    each cancer has its own project folder (MRTrios_BRCA, MRTrios_LIHC, ...)
proj_dir <- function(project) {
  if (env == "local") file.path(local_root, project) else ondemand_root
}

# ── One entry per dataset ──────────────────────────────────
# tag is used for both the shard pattern (trio_<tag>_results_part*.txt)
# and the FDR output name (<tag>_trio_Model_results_ALL_with_BH_fdr_qval_byLZ.txt)
cancers <- list(
  posER_BRCA = list(project     = "MRTrios_BRCA", tag = "posER_BRCA",
                    part_dir    = "Output_posER_BRCA_part_Oct09_2026",
                    out_dir     = "Output_posER_BRCA_Oct09_2026",
                    merged_file = "trio_posER_BRCA_results_ALL_v2.txt"),
  negER_BRCA = list(project     = "MRTrios_BRCA", tag = "negER_BRCA",
                    part_dir    = "Output_negER_BRCA_part_Oct09_2026",
                    out_dir     = "Output_negER_BRCA_Oct09_2026",
                    merged_file = "trio_negER_BRCA_results_ALL_v2.txt"),
  LIHC       = list(project     = "MRTrios_LIHC", tag = "LIHC",
                    part_dir    = "Output_LIHC_part",
                    out_dir     = "Output_LIHC",
                    merged_file = "trio_results_LIHC_ALL.txt"),
  BLCA       = list(project     = "MRTrios_BLCA", tag = "BLCA",
                    part_dir    = "Output_BLCA_part",
                    out_dir     = "Output_BLCA",
                    merged_file = "trio_results_BLCA_ALL.txt")
)

# ════════════════════════════════════════════════════════════════
# Step 1: merge shards
# ════════════════════════════════════════════════════════════════
merge_shards <- function(cfg) {
  dir_in  <- file.path(proj_dir(cfg$project), cfg$part_dir)
  dir_out <- file.path(proj_dir(cfg$project), cfg$out_dir)
  dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
  
  shard_files <- sort(list.files(dir_in,
                                 pattern = sprintf("trio_%s_results_part.*\\.txt$", cfg$tag),
                                 full.names = TRUE))
  if (length(shard_files) == 0) stop("No shard files found in ", dir_in)
  cat(sprintf("[%s] Found %d shard files. Merging...\n", cfg$tag, length(shard_files)))
  
  res <- rbindlist(lapply(shard_files, fread), fill = TRUE)
  cat(sprintf("[%s] Total rows: %d | Total cols: %d\n", cfg$tag, nrow(res), ncol(res)))
  
  out_file <- file.path(dir_out, cfg$merged_file)
  fwrite(res, file = out_file, sep = "\t")
  cat(sprintf("[%s] Saved: %s\n", cfg$tag, out_file))
  print(table(res$Confounders))
  invisible(res)
}

# ════════════════════════════════════════════════════════════════
# Step 2: FDR adjustment
# ════════════════════════════════════════════════════════════════
# p-values adjusted together in pairs (same as before)
p_pairs   <- list(c("b11", "b21"), c("b12", "b22"), c("V1:T1", "V1:T2"))
# order of indicators passed to MRGN::class.vec
ind_order <- c("b11", "b12", "b21", "b22", "V1:T1", "V1:T2")

# method = "BH" -> columns *_BH_fdr ; method = "qvalue" -> columns *_qval
add_adjusted <- function(dt, method) {
  suffix <- if (method == "BH") "BH_fdr" else "qval"
  n <- nrow(dt)
  
  # adjusted p-values: pb11_BH_fdr, pV1:T1_qval, ...
  for (pr in p_pairs) {
    p   <- c(dt[[paste0("p", pr[1])]], dt[[paste0("p", pr[2])]])
    adj <- if (method == "BH") p.adjust(p, method = "BH") else qvalue(p = p)$qvalues
    dt[, (paste0("p", pr[1], "_", suffix)) := adj[1:n]]
    dt[, (paste0("p", pr[2], "_", suffix)) := adj[(n + 1):(2 * n)]]
  }
  
  # binary indicators: b11_BH_fdr, V1:T1_qval, ...
  for (b in ind_order) {
    dt[, (paste0(b, "_", suffix)) := as.integer(get(paste0("p", b, "_", suffix)) < alpha)]
  }
  
  # re-inferred model: Inferred.Model.BH_fdr / Inferred.Model.qval
  ind_mat <- as.matrix(dt[, paste0(ind_order, "_", suffix), with = FALSE])
  dt[, (paste0("Inferred.Model.", suffix)) := apply(ind_mat, 1, class.vec)]
  invisible(dt)
}

fdr_adjust <- function(cfg, dt = NULL) {
  dir_out <- file.path(proj_dir(cfg$project), cfg$out_dir)
  if (is.null(dt)) dt <- fread(file.path(dir_out, cfg$merged_file), sep = "\t")
  
  add_adjusted(dt, "BH")
  add_adjusted(dt, "qvalue")
  
  cat(sprintf("\n===== %s =====\n", cfg$tag))
  cat("Original:\n");   print(table(dt$Inferred.Model))
  cat("\nBH (FDR):\n"); print(table(dt$Inferred.Model.BH_fdr))
  cat("\nqvalue:\n");   print(table(dt$Inferred.Model.qval))
  cat("\nConfounders:\n"); print(table(dt$Confounders))
  
  out_file <- file.path(dir_out, sprintf("%s_trio_Model_results_ALL_with_BH_fdr_qval_v2.txt", cfg$tag))
  fwrite(dt, file = out_file, sep = "\t")
  cat(sprintf("[%s] Saved: %s\n", cfg$tag, out_file))
  invisible(dt)
}

# ════════════════════════════════════════════════════════════════
# Run (subset names(cancers) to run only some, e.g. c("LIHC", "BLCA"))
# ════════════════════════════════════════════════════════════════
results <- list()
for (nm in names(cancers)) {
  cfg <- cancers[[nm]]
  dt  <- if (run_merge) merge_shards(cfg) else NULL
  if (run_fdr) dt <- fdr_adjust(cfg, dt)
  results[[nm]] <- dt
}

# ════════════════════════════════════════════════════════════════
# Output print
# ════════════════════════════════════════════════════════════════
# [posER_BRCA] Found 600 shard files. Merging...
# [posER_BRCA] Total rows: 292525 | Total cols: 24
# [posER_BRCA] Saved: /Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_posER_BRCA_Oct09_2026/trio_posER_BRCA_results_ALL_v2.txt
# 
# age|race|U_exp|U_meth       age|race|U_meth 
# 292372                   153 
# 
# ===== posER_BRCA =====
#   Original:
#   
#   M0.1  M0.2  M1.1  M1.2  M2.1  M2.2    M3    M4 Other 
# 83091 19476 34315  7255  8006  2086 44927 16363 77006 
# 
# BH (FDR):
#   
#   M0.1  M0.2  M1.1  M1.2  M2.1  M2.2    M3    M4 Other 
# 74262 20750 36960  8890 10386  2999 53335 22960 61983 
# 
# qvalue:
#   
#   M0.1  M0.2  M1.1  M1.2  M2.1  M2.2    M3    M4 Other 
# 53316 19946 39576 11416 14632  4851 68700 43089 36999 
# 
# Confounders:
#   
#   age|race|U_exp|U_meth       age|race|U_meth 
# 292372                   153 
# [posER_BRCA] Saved: /Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_posER_BRCA_Oct09_2026/posER_BRCA_trio_Model_results_ALL_with_BH_fdr_qval_v2.txt
# [negER_BRCA] Found 600 shard files. Merging...
# [negER_BRCA] Total rows: 292525 | Total cols: 24
# [negER_BRCA] Saved: /Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_negER_BRCA_Oct09_2026/trio_negER_BRCA_results_ALL_v2.txt
# 
# age|race|U_exp age|race|U_exp|U_meth       age|race|U_meth 
# 127                290235                  2163 
# 
# ===== negER_BRCA =====
#   Original:
#   
#   M0.1   M0.2   M1.1   M1.2   M2.1   M2.2     M3     M4  Other 
# 118069  13697  13714   1980   1678    413  16707    881 125386 
# 
# BH (FDR):
#   
#   M0.1   M0.2   M1.1   M1.2   M2.1   M2.2     M3     M4  Other 
# 122553  16137  11486   1920   1922    571  22953    971 114012 
# 
# qvalue:
#   
#   M0.1   M0.2   M1.1   M1.2   M2.1   M2.2     M3     M4  Other 
# 118033  20066  15805   3157   3433   1014  37345   2855  90817 
# 
# Confounders:
#   
#   age|race|U_exp age|race|U_exp|U_meth       age|race|U_meth 
# 127                290235                  2163 
# [negER_BRCA] Saved: /Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_negER_BRCA_Oct09_2026/negER_BRCA_trio_Model_results_ALL_with_BH_fdr_qval_v2.txt
# [LIHC] Found 620 shard files. Merging...
# [LIHC] Total rows: 308028 | Total cols: 24
# [LIHC] Saved: /Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_LIHC/Output_LIHC/trio_results_LIHC_ALL.txt
# 
# age|sex|race        age|sex|race|U_exp age|sex|race|U_exp|U_meth       age|sex|race|U_meth 
# 13                      2416                    302983                      2616 
# 
# ===== LIHC =====
#   Original:
#   
#   M0.1   M0.2   M1.1   M1.2   M2.1   M2.2     M3     M4  Other 
# 101753  13913  29166   6389   4704   1907  38847   5716 105633 
# 
# BH (FDR):
#   
#   M0.1  M0.2  M1.1  M1.2  M2.1  M2.2    M3    M4 Other 
# 94678 15987 32408  8250  6247  2467 47698  9009 91284 
# 
# qvalue:
#   
#   M0.1  M0.2  M1.1  M1.2  M2.1  M2.2    M3    M4 Other 
# 76978 18105 38363 11978 10320  4345 64124 19423 64392 
# 
# Confounders:
#   
#   age|sex|race        age|sex|race|U_exp age|sex|race|U_exp|U_meth       age|sex|race|U_meth 
# 13                      2416                    302983                      2616 
# [LIHC] Saved: /Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_LIHC/Output_LIHC/LIHC_trio_Model_results_ALL_with_BH_fdr_qval_v2.txt
# [BLCA] Found 620 shard files. Merging...
# [BLCA] Total rows: 308405 | Total cols: 24
# [BLCA] Saved: /Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BLCA/Output_BLCA/trio_results_BLCA_ALL.txt
# 
# age|sex|race        age|sex|race|U_exp age|sex|race|U_exp|U_meth       age|sex|race|U_meth 
# 4                      1605                    305106                      1690 
# 
# ===== BLCA =====
#   Original:
#   
#   M0.1  M0.2  M1.1  M1.2  M2.1  M2.2    M3    M4 Other 
# 94588 15541 37606  8697  9987  2486 44780 11042 83678 
# 
# BH (FDR):
#   
#   M0.1  M0.2  M1.1  M1.2  M2.1  M2.2    M3    M4 Other 
# 83391 16751 41195 10539 13062  4080 54285 16994 68108 
# 
# qvalue:
#   
#   M0.1  M0.2  M1.1  M1.2  M2.1  M2.2    M3    M4 Other 
# 58414 16083 43900 13325 19364  7615 72616 37270 39818 
# 
# Confounders:
#   
#   age|sex|race        age|sex|race|U_exp age|sex|race|U_exp|U_meth       age|sex|race|U_meth 
# 4                      1605                    305106                      1690 
# [BLCA] Saved: /Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BLCA/Output_BLCA/BLCA_trio_Model_results_ALL_with_BH_fdr_qval_v2.txt