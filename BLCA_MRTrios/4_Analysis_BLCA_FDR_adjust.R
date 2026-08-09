########################################################################################
###################      FDR    ########################################################
########################################################################################

BLCA_Model <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BLCA/Output_BLCA/trio_results_BLCA_ALL.txt",sep ="\t")
## ── BLCA model FDR ──────────────────────────────────────────────────
output_dir ="/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BLCA/Output_BLCA/"

# ── Install if needed ──────────────────────────────────────
if (!requireNamespace("qvalue", quietly = TRUE)) {
  BiocManager::install("qvalue")
}
library(qvalue)

n <- nrow(BLCA_Model)

## ════════════════════════════════════════════════════════════════
## BH (FDR) — Benjamini-Hochberg via p.adjust
## ════════════════════════════════════════════════════════════════

# pb11 & pb21 together
combined_1 <- c(BLCA_Model$pb11, BLCA_Model$pb21)
fdr_1      <- p.adjust(combined_1, method = "BH")
BLCA_Model$pb11_BH_fdr <- fdr_1[1:n]
BLCA_Model$pb21_BH_fdr <- fdr_1[(n + 1):(2 * n)]

# pb12 & pb22 together
combined_2 <- c(BLCA_Model$pb12, BLCA_Model$pb22)
fdr_2      <- p.adjust(combined_2, method = "BH")
BLCA_Model$pb12_BH_fdr <- fdr_2[1:n]
BLCA_Model$pb22_BH_fdr <- fdr_2[(n + 1):(2 * n)]

# pV1:T1 & pV1:T2 together
combined_3 <- c(BLCA_Model$`pV1:T1`, BLCA_Model$`pV1:T2`)
fdr_3      <- p.adjust(combined_3, method = "BH")
BLCA_Model$`pV1:T1_BH_fdr` <- fdr_3[1:n]
BLCA_Model$`pV1:T2_BH_fdr` <- fdr_3[(n + 1):(2 * n)]

alpha <- 0.05
BLCA_Model$b11_BH_fdr     <- as.integer(BLCA_Model$pb11_BH_fdr     < alpha)
BLCA_Model$b12_BH_fdr     <- as.integer(BLCA_Model$pb12_BH_fdr     < alpha)
BLCA_Model$b21_BH_fdr     <- as.integer(BLCA_Model$pb21_BH_fdr     < alpha)
BLCA_Model$b22_BH_fdr     <- as.integer(BLCA_Model$pb22_BH_fdr     < alpha)
BLCA_Model$`V1:T1_BH_fdr` <- as.integer(BLCA_Model$`pV1:T1_BH_fdr` < alpha)
BLCA_Model$`V1:T2_BH_fdr` <- as.integer(BLCA_Model$`pV1:T2_BH_fdr` < alpha)

library(MRGN)
BLCA_Model$Inferred.Model.BH_fdr <- apply(BLCA_Model, 1, function(row) {
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
BLCA_Model$pb11_qval <- qobj_1$qvalues[1:n]
BLCA_Model$pb21_qval <- qobj_1$qvalues[(n + 1):(2 * n)]

# pb12 & pb22 together
qobj_2 <- qvalue(p = combined_2)
BLCA_Model$pb12_qval <- qobj_2$qvalues[1:n]
BLCA_Model$pb22_qval <- qobj_2$qvalues[(n + 1):(2 * n)]

# pV1:T1 & pV1:T2 together
qobj_3 <- qvalue(p = combined_3)
BLCA_Model$`pV1:T1_qval` <- qobj_3$qvalues[1:n]
BLCA_Model$`pV1:T2_qval` <- qobj_3$qvalues[(n + 1):(2 * n)]

# ── Re-infer binary indicators ─────────────────────────────
BLCA_Model$b11_qval    <- as.integer(BLCA_Model$pb11_qval    < alpha)
BLCA_Model$b12_qval    <- as.integer(BLCA_Model$pb12_qval    < alpha)
BLCA_Model$b21_qval    <- as.integer(BLCA_Model$pb21_qval    < alpha)
BLCA_Model$b22_qval    <- as.integer(BLCA_Model$pb22_qval    < alpha)
BLCA_Model$`V1:T1_qval` <- as.integer(BLCA_Model$`pV1:T1_qval` < alpha)
BLCA_Model$`V1:T2_qval` <- as.integer(BLCA_Model$`pV1:T2_qval` < alpha)

# ── Re-infer model ─────────────────────────────────────────
BLCA_Model$Inferred.Model.qval <- apply(BLCA_Model, 1, function(row) {
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
cat("Original:\n");    print(table(BLCA_Model$Inferred.Model))
cat("\nBH (FDR):\n");  print(table(BLCA_Model$Inferred.Model.BH_fdr))
cat("\nqvalue:\n");    print(table(BLCA_Model$Inferred.Model.qval))

# ── Save ───────────────────────────────────────────────────
fwrite(BLCA_Model,
       file = file.path(output_dir, "BLCA_trio_Model_results_ALL_with_BH_fdr_qval_byLZ.txt"),
       sep  = "\t")

########################################################################################
