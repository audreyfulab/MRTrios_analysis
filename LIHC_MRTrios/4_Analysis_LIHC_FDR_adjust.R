########################################################################################
###################      FDR    ########################################################
########################################################################################

LIHC_Model <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_LIHC/Output_LIHC/trio_results_LIHC_ALL.txt",sep ="\t")
## ── LIHC model FDR ──────────────────────────────────────────────────
output_dir ="/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_LIHC/Output_LIHC/"

# ── Install if needed ──────────────────────────────────────
if (!requireNamespace("qvalue", quietly = TRUE)) {
  BiocManager::install("qvalue")
}
library(qvalue)

n <- nrow(LIHC_Model)

## ════════════════════════════════════════════════════════════════
## BH (FDR) — Benjamini-Hochberg via p.adjust
## ════════════════════════════════════════════════════════════════

# pb11 & pb21 together
combined_1 <- c(LIHC_Model$pb11, LIHC_Model$pb21)
fdr_1      <- p.adjust(combined_1, method = "BH")
LIHC_Model$pb11_BH_fdr <- fdr_1[1:n]
LIHC_Model$pb21_BH_fdr <- fdr_1[(n + 1):(2 * n)]

# pb12 & pb22 together
combined_2 <- c(LIHC_Model$pb12, LIHC_Model$pb22)
fdr_2      <- p.adjust(combined_2, method = "BH")
LIHC_Model$pb12_BH_fdr <- fdr_2[1:n]
LIHC_Model$pb22_BH_fdr <- fdr_2[(n + 1):(2 * n)]

# pV1:T1 & pV1:T2 together
combined_3 <- c(LIHC_Model$`pV1:T1`, LIHC_Model$`pV1:T2`)
fdr_3      <- p.adjust(combined_3, method = "BH")
LIHC_Model$`pV1:T1_BH_fdr` <- fdr_3[1:n]
LIHC_Model$`pV1:T2_BH_fdr` <- fdr_3[(n + 1):(2 * n)]

alpha <- 0.05
LIHC_Model$b11_BH_fdr     <- as.integer(LIHC_Model$pb11_BH_fdr     < alpha)
LIHC_Model$b12_BH_fdr     <- as.integer(LIHC_Model$pb12_BH_fdr     < alpha)
LIHC_Model$b21_BH_fdr     <- as.integer(LIHC_Model$pb21_BH_fdr     < alpha)
LIHC_Model$b22_BH_fdr     <- as.integer(LIHC_Model$pb22_BH_fdr     < alpha)
LIHC_Model$`V1:T1_BH_fdr` <- as.integer(LIHC_Model$`pV1:T1_BH_fdr` < alpha)
LIHC_Model$`V1:T2_BH_fdr` <- as.integer(LIHC_Model$`pV1:T2_BH_fdr` < alpha)

library(MRGN)
LIHC_Model$Inferred.Model.BH_fdr <- apply(LIHC_Model, 1, function(row) {
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
LIHC_Model$pb11_qval <- qobj_1$qvalues[1:n]
LIHC_Model$pb21_qval <- qobj_1$qvalues[(n + 1):(2 * n)]

# pb12 & pb22 together
qobj_2 <- qvalue(p = combined_2)
LIHC_Model$pb12_qval <- qobj_2$qvalues[1:n]
LIHC_Model$pb22_qval <- qobj_2$qvalues[(n + 1):(2 * n)]

# pV1:T1 & pV1:T2 together
qobj_3 <- qvalue(p = combined_3)
LIHC_Model$`pV1:T1_qval` <- qobj_3$qvalues[1:n]
LIHC_Model$`pV1:T2_qval` <- qobj_3$qvalues[(n + 1):(2 * n)]

# ── Re-infer binary indicators ─────────────────────────────
LIHC_Model$b11_qval    <- as.integer(LIHC_Model$pb11_qval    < alpha)
LIHC_Model$b12_qval    <- as.integer(LIHC_Model$pb12_qval    < alpha)
LIHC_Model$b21_qval    <- as.integer(LIHC_Model$pb21_qval    < alpha)
LIHC_Model$b22_qval    <- as.integer(LIHC_Model$pb22_qval    < alpha)
LIHC_Model$`V1:T1_qval` <- as.integer(LIHC_Model$`pV1:T1_qval` < alpha)
LIHC_Model$`V1:T2_qval` <- as.integer(LIHC_Model$`pV1:T2_qval` < alpha)

# ── Re-infer model ─────────────────────────────────────────
LIHC_Model$Inferred.Model.qval <- apply(LIHC_Model, 1, function(row) {
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
cat("Original:\n");    print(table(LIHC_Model$Inferred.Model))
cat("\nBH (FDR):\n");  print(table(LIHC_Model$Inferred.Model.BH_fdr))
cat("\nqvalue:\n");    print(table(LIHC_Model$Inferred.Model.qval))

# ── Save ───────────────────────────────────────────────────
fwrite(LIHC_Model,
       file = file.path(output_dir, "LIHC_trio_Model_results_ALL_with_BH_fdr_qval_byLZ.txt"),
       sep  = "\t")

# ── Check Confounders categoriers ─────────────────────────────────────────────────── 
table(LIHC_Model$Confounders)


########################################################################################
# ─────────────────────────────────────────────────────────────────────────────────────────  
# LIHC_Model saved:
# /Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_LIHC/Output_LIHC/LIHC_trio_Model_results_ALL_with_BH_fdr_qval_byLZ.txt

# Compare all three methods:
# Original with pvalue:
#  M0.1   M0.2   M1.1   M1.2   M2.1   M2.2     M3     M4  Other 
# 101753  13913  29166   6389   4704   1907  38847   5716 105633

# BH (FDR):
#  M0.1   M0.2   M1.1   M1.2   M2.1   M2.2     M3     M4  Other 
# 94678   15987  32408   8250  6247    2467  47698   9009 91284 

# qvalue:
#  M0.1   M0.2   M1.1   M1.2   M2.1   M2.2     M3     M4  Other  
# 76978   18105  38363  11978  10320   4345  64124   19423 64392

# table(LIHC_Model$Confounders)

# age|sex|race        age|sex|race|U_exp     age|sex|race|U_exp|U_meth       age|sex|race|U_meth
#    13                      2416                    302983                         2616
 
# ─────────────────────────────────────────────────────────────────────────────────────────





# ─────────────────────────────────────────────────────────────────────────────────────────  
