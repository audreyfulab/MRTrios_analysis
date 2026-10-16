library(data.table)
# ============================================================
# Per-trio summary statistics for posER BRCA, negER BRCA, LIHC, BLCA
#   - mean of gene expression / methylation
#   - pairwise correlations: CNA-Exp, CNA-Meth, Exp-Meth
# ============================================================

# ============================================================
# Settings
# ============================================================
env    <- "local"            # "local" or "ondemand"
method <- "pearson"          # or "spearman"

ondemand_root <- "/wsu/home/hb/hb68/hb6890/fulab/MRTrios"
local_root    <- "/Users/lianzuo/LZ/ResearchProject/Fulab"

# Trio file (shared by all cancers unless a cancer sets its own `trio_file` below)
trio_file <- if (env == "local") {
  file.path(local_root, "MRTrios_analysis/raw_Data_Methyl/trio.final.protein.coding.txt")
} else {
  file.path(ondemand_root, "raw_Data_Methyl/trio.final.protein.coding.txt")
}

# Where the combined summary is written
summary_dir <- if (env == "local") file.path(local_root, "MRTrios_analysis") else ondemand_root

# ── One entry per dataset ──────────────────────────────────
# rds_local / rds_ondemand: folder holding the *_data_filter_list.rds
# (stats file is written to the same folder)
cancers <- list(
  posER_BRCA = list(project = "MRTrios_BRCA",
                    rds_local = "process_data_BRCA", rds_ondemand = "Output_posER_BRCA",
                    rds_file  = "posER_BRCA_data_filter_list.rds"),
  negER_BRCA = list(project = "MRTrios_BRCA",
                    rds_local = "process_data_BRCA", rds_ondemand = "Output_posER_BRCA",
                    rds_file  = "negER_BRCA_data_filter_list.rds"),
  LIHC       = list(project = "MRTrios_LIHC",
                    rds_local = "process_data_LIHC", rds_ondemand = "Output_LIHC",
                    rds_file  = "LIHC_data_filter_list.rds"),
  BLCA       = list(project = "MRTrios_BLCA",
                    rds_local = "process_data_BLCA", rds_ondemand = "Output_BLCA",
                    rds_file  = "BLCA_data_filter_list.rds")
)

rds_dir <- function(cfg) {
  if (env == "local") file.path(local_root, cfg$project, cfg$rds_local)
  else                file.path(ondemand_root, cfg$rds_ondemand)
}

# ============================================================
# Helper functions (unchanged)
# ============================================================

# Build a trio-aligned numeric matrix:
# row i of the output = the row used by trio i (NA row if not found)
to_trio_matrix <- function(df, id_col, trio_ids, sample_ids) {
  idx <- match(as.character(trio_ids), as.character(df[[id_col]]))
  m <- as.matrix(df[idx, sample_ids, drop = FALSE])
  storage.mode(m) <- "double"
  m
}

# Fast row-wise correlation with pairwise NA removal
# (equivalent to cor(x, y, use = "complete.obs") for each row)
row_cor <- function(X, Y, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  ok <- !is.na(X) & !is.na(Y)
  X[!ok] <- NA; Y[!ok] <- NA
  if (method == "spearman") {
    X <- t(apply(X, 1, rank, na.last = "keep"))
    Y <- t(apply(Y, 1, rank, na.last = "keep"))
  }
  n  <- rowSums(ok)
  Xc <- X - rowMeans(X, na.rm = TRUE)
  Yc <- Y - rowMeans(Y, na.rm = TRUE)
  num <- rowSums(Xc * Yc, na.rm = TRUE)
  den <- sqrt(rowSums(Xc^2, na.rm = TRUE) * rowSums(Yc^2, na.rm = TRUE))
  r <- num / den
  r[!is.finite(r) | n < 3] <- NA   # zero variance (e.g. constant CNA) or too few samples
  list(r = r, n = n)
}

compute_trio_stats <- function(d, trios, method = "pearson") {
  ids <- d$unique_id
  
  C <- to_trio_matrix(d$cna_filter,  "cna.row",  trios$cna.row,  ids)
  E <- to_trio_matrix(d$exp_filter,  "gene.row", trios$gene.row, ids)
  M <- to_trio_matrix(d$meth_filter, "meth.row", trios$meth.row, ids)
  
  n_miss <- c(CNA  = sum(rowSums(!is.na(C)) == 0),
              Exp  = sum(rowSums(!is.na(E)) == 0),
              Meth = sum(rowSums(!is.na(M)) == 0))
  message("  Trios with no data -> ",
          paste(names(n_miss), n_miss, sep = ": ", collapse = ", "))
  
  ce <- row_cor(C, E, method)
  cm <- row_cor(C, M, method)
  em <- row_cor(E, M, method)
  
  mean_exp  <- rowMeans(E, na.rm = TRUE); mean_exp[is.nan(mean_exp)]   <- NA
  mean_meth <- rowMeans(M, na.rm = TRUE); mean_meth[is.nan(mean_meth)] <- NA
  
  data.table(
    trio_id      = seq_len(nrow(trios)),
    trios,
    mean_exp     = mean_exp,
    #    sd_exp       = apply(E, 1, sd, na.rm = TRUE),
    mean_meth    = mean_meth,
    #    sd_meth      = apply(M, 1, sd, na.rm = TRUE),
    cor_cna_exp  = ce$r,
    cor_cna_meth = cm$r,
    cor_exp_meth = em$r
  )
}

# ============================================================
# Run for every cancer and save
# (subset names(cancers) to run only some, e.g. c("LIHC", "BLCA"))
# ============================================================
trio_cache <- list()
stats_list <- list()

for (nm in names(cancers)) {
  cfg <- cancers[[nm]]
  rds_path <- file.path(rds_dir(cfg), cfg$rds_file)
  if (!file.exists(rds_path)) {
    warning("[", nm, "] skipped: file not found ", rds_path)
    next
  }
  message("[", nm, "] computing stats...")
  
  tf <- if (!is.null(cfg$trio_file)) cfg$trio_file else trio_file
  if (is.null(trio_cache[[tf]])) trio_cache[[tf]] <- fread(tf)
  
  d <- readRDS(rds_path)
  stats <- compute_trio_stats(d, trio_cache[[tf]], method)
  rm(d); gc()
  
  out_file <- file.path(rds_dir(cfg), paste0(nm, "_trio_stats.txt"))
  fwrite(stats, out_file, sep = "\t")
  message("  Saved: ", out_file)
  stats_list[[nm]] <- stats
}

# ============================================================
# Side-by-side summary across all cancers
# ============================================================
all_stats <- rbindlist(stats_list, idcol = "cancer")

summary_all <- all_stats[, .(
  n_trios             = .N,
  mean_of_mean_exp    = mean(mean_exp,       na.rm = TRUE),
  mean_of_mean_meth   = mean(mean_meth,      na.rm = TRUE),
  median_cor_cna_exp  = median(cor_cna_exp,  na.rm = TRUE),
  median_cor_cna_meth = median(cor_cna_meth, na.rm = TRUE),
  median_cor_exp_meth = median(cor_exp_meth, na.rm = TRUE),
  n_NA_cor_cna_exp    = sum(is.na(cor_cna_exp)),
  n_NA_cor_cna_meth   = sum(is.na(cor_cna_meth)),
  n_NA_cor_exp_meth   = sum(is.na(cor_exp_meth))
), by = cancer]

print(summary_all)
fwrite(summary_all, file.path(summary_dir, "ALL_cancer_trio_stats_summary.txt"), sep = "\t")
