library(data.table)
library(tidyverse)
# ============================================================
# Per-trio summary statistics for posER / negER BRCA
#   - mean & SD of gene expression
#   - mean & SD of methylation
#   - pairwise correlations: CNA-Exp, CNA-Meth, Exp-Meth
# ============================================================


# ============================================================
# Step 1: Load trios + saved preprocessed lists
# ============================================================
trios <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/raw_Data_Methyl/trio.final.protein.coding.txt")

setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/process_data_BRCA")
posER_BRCA_filter <- readRDS("posER_BRCA_data_filter_list.rds")
negER_BRCA_filter <- readRDS("negER_BRCA_data_filter_list.rds")

# ============================================================
# Step 2: Helper functions
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

# ============================================================
# Step 3: Main function - stats for every trio
# ============================================================
compute_trio_stats <- function(d, trios, method = "pearson") {
  ids <- d$unique_id
  
  C <- to_trio_matrix(d$cna_filter,  "cna.row",  trios$cna.row,  ids)
  E <- to_trio_matrix(d$exp_filter,  "gene.row", trios$gene.row, ids)
  M <- to_trio_matrix(d$meth_filter, "meth.row", trios$meth.row, ids)
   
  # Sanity check: trios whose rows were not found (all-NA)
  n_miss <- c(CNA  = sum(rowSums(!is.na(C)) == 0),
              Exp  = sum(rowSums(!is.na(E)) == 0),
              Meth = sum(rowSums(!is.na(M)) == 0))
  message("Trios with no data -> ",
          paste(names(n_miss), n_miss, sep = ": ", collapse = ", "))
  
  ce <- row_cor(C, E, method)
  cm <- row_cor(C, M, method)
  em <- row_cor(E, M, method)
  
  mean_exp  <- rowMeans(E, na.rm = TRUE); mean_exp[is.nan(mean_exp)]   <- NA
  mean_meth <- rowMeans(M, na.rm = TRUE); mean_meth[is.nan(mean_meth)] <- NA
  
  data.table(
    trio_id      = seq_len(nrow(trios)),
    trios,                                   # keeps cna.row / gene.row / meth.row etc.
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
# Step 4: Run for both subgroups and save
# ============================================================
posER_stats <- compute_trio_stats(posER_BRCA_filter, trios)
negER_stats <- compute_trio_stats(negER_BRCA_filter, trios)

fwrite(posER_stats, "posER_BRCA_trio_stats.txt", sep = "\t")
fwrite(negER_stats, "negER_BRCA_trio_stats.txt", sep = "\t")

# ============================================================
# Step 5: Side-by-side summary (posER vs negER)
# ============================================================
all_stats <- rbindlist(list(posER = posER_stats, negER = negER_stats),
                       idcol = "ER")

summary_ER <- all_stats[, .(
  n_trios             = .N,
  mean_of_mean_exp    = mean(mean_exp,      na.rm = TRUE),
  mean_of_mean_meth   = mean(mean_meth,     na.rm = TRUE),
  median_cor_cna_exp  = median(cor_cna_exp,  na.rm = TRUE),
  median_cor_cna_meth = median(cor_cna_meth, na.rm = TRUE),
  median_cor_exp_meth = median(cor_exp_meth, na.rm = TRUE),
  n_NA_cor_cna_exp    = sum(is.na(cor_cna_exp)),
  n_NA_cor_cna_meth   = sum(is.na(cor_cna_meth)),
  n_NA_cor_exp_meth   = sum(is.na(cor_exp_meth))
), by = ER]

print(summary_ER)
fwrite(summary_ER, "BRCA_ER_trio_stats_summary.txt", sep = "\t") 
