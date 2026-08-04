########################################################################################
###################      FDR    ########################################################
########################################################################################

########################################################################################
### ADD BH p.ajust
########################################################################################


Model_posER_BRCA <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/trio_posER_BRCA_results_ALL.txt", sep = "\t")

## ── Pos model FDR ──────────────────────────────────────────────────
output_dir = "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/"

# ── Install if needed ──────────────────────────────────────
if (!requireNamespace("qvalue", quietly = TRUE)) {
  BiocManager::install("qvalue")
}
library(qvalue)

n <- nrow(Model_posER_BRCA)

## ════════════════════════════════════════════════════════════════
## BH (FDR) — Benjamini-Hochberg via p.adjust
## ════════════════════════════════════════════════════════════════

# pb11 & pb21 together
combined_1 <- c(Model_posER_BRCA$pb11, Model_posER_BRCA$pb21)
fdr_1      <- p.adjust(combined_1, method = "BH")
Model_posER_BRCA$pb11_BH_fdr <- fdr_1[1:n]
Model_posER_BRCA$pb21_BH_fdr <- fdr_1[(n + 1):(2 * n)]

# pb12 & pb22 together
combined_2 <- c(Model_posER_BRCA$pb12, Model_posER_BRCA$pb22)
fdr_2      <- p.adjust(combined_2, method = "BH")
Model_posER_BRCA$pb12_BH_fdr <- fdr_2[1:n]
Model_posER_BRCA$pb22_BH_fdr <- fdr_2[(n + 1):(2 * n)]

# pV1:T1 & pV1:T2 together
combined_3 <- c(Model_posER_BRCA$`pV1:T1`, Model_posER_BRCA$`pV1:T2`)
fdr_3      <- p.adjust(combined_3, method = "BH")
Model_posER_BRCA$`pV1:T1_BH_fdr` <- fdr_3[1:n]
Model_posER_BRCA$`pV1:T2_BH_fdr` <- fdr_3[(n + 1):(2 * n)]

alpha <- 0.05
Model_posER_BRCA$b11_BH_fdr     <- as.integer(Model_posER_BRCA$pb11_BH_fdr     < alpha)
Model_posER_BRCA$b12_BH_fdr     <- as.integer(Model_posER_BRCA$pb12_BH_fdr     < alpha)
Model_posER_BRCA$b21_BH_fdr     <- as.integer(Model_posER_BRCA$pb21_BH_fdr     < alpha)
Model_posER_BRCA$b22_BH_fdr     <- as.integer(Model_posER_BRCA$pb22_BH_fdr     < alpha)
Model_posER_BRCA$`V1:T1_BH_fdr` <- as.integer(Model_posER_BRCA$`pV1:T1_BH_fdr` < alpha)
Model_posER_BRCA$`V1:T2_BH_fdr` <- as.integer(Model_posER_BRCA$`pV1:T2_BH_fdr` < alpha)

library(MRGN)
Model_posER_BRCA$Inferred.Model.BH_fdr <- apply(Model_posER_BRCA, 1, function(row) {
  vec <- as.numeric(c(
    row["b11_BH_fdr"],
    row["b12_BH_fdr"],
    row["b21_BH_fdr"],
    row["b22_BH_fdr"],
    row["V1:T1_BH_fdr"],
    row["V1:T2_BH_fdr"]
  ))
  class.vec(vec)
})

## ════════════════════════════════════════════════════════════════
## qvalue (Storey)
## ════════════════════════════════════════════════════════════════

# pb11 & pb21 together
qobj_1 <- qvalue(p = combined_1)
Model_posER_BRCA$pb11_qval <- qobj_1$qvalues[1:n]
Model_posER_BRCA$pb21_qval <- qobj_1$qvalues[(n + 1):(2 * n)]

# pb12 & pb22 together
qobj_2 <- qvalue(p = combined_2)
Model_posER_BRCA$pb12_qval <- qobj_2$qvalues[1:n]
Model_posER_BRCA$pb22_qval <- qobj_2$qvalues[(n + 1):(2 * n)]

# pV1:T1 & pV1:T2 together
qobj_3 <- qvalue(p = combined_3)
Model_posER_BRCA$`pV1:T1_qval` <- qobj_3$qvalues[1:n]
Model_posER_BRCA$`pV1:T2_qval` <- qobj_3$qvalues[(n + 1):(2 * n)]

# ── Re-infer binary indicators ─────────────────────────────
Model_posER_BRCA$b11_qval    <- as.integer(Model_posER_BRCA$pb11_qval    < alpha)
Model_posER_BRCA$b12_qval    <- as.integer(Model_posER_BRCA$pb12_qval    < alpha)
Model_posER_BRCA$b21_qval    <- as.integer(Model_posER_BRCA$pb21_qval    < alpha)
Model_posER_BRCA$b22_qval    <- as.integer(Model_posER_BRCA$pb22_qval    < alpha)
Model_posER_BRCA$`V1:T1_qval` <- as.integer(Model_posER_BRCA$`pV1:T1_qval` < alpha)
Model_posER_BRCA$`V1:T2_qval` <- as.integer(Model_posER_BRCA$`pV1:T2_qval` < alpha)

# ── Re-infer model ─────────────────────────────────────────
Model_posER_BRCA$Inferred.Model.qval <- apply(Model_posER_BRCA, 1, function(row) {
  vec <- as.numeric(c(
    row["b11_qval"],
    row["b12_qval"],
    row["b21_qval"],
    row["b22_qval"],
    row["V1:T1_qval"],
    row["V1:T2_qval"]
  ))
  class.vec(vec)
})

# ── Compare all three methods ──────────────────────────────
cat("Original:\n");    print(table(Model_posER_BRCA$Inferred.Model))
cat("\nBH (FDR):\n");  print(table(Model_posER_BRCA$Inferred.Model.BH_fdr))
cat("\nqvalue:\n");    print(table(Model_posER_BRCA$Inferred.Model.qval))

# ── Save ───────────────────────────────────────────────────
fwrite(Model_posER_BRCA,
       file = file.path(output_dir, "posER_BRCA_trio_Model_results_ALL_with_BH_fdr_qval_byLZ.txt"),
       sep  = "\t")

########################################################################################
