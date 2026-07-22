########################################################################################
###################      FDR    ########################################################
########################################################################################

BLCA_Model <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/Output_BLCA/trio_results_BLCA_ALL.txt",sep ="\t")

## ── Pos model FDR ──────────────────────────────────────────────────
output_dir ="/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_BLCA/"

# ── Install if needed ──────────────────────────────────────
if (!requireNamespace("qvalue", quietly = TRUE)) {
  BiocManager::install("qvalue")
}
library(qvalue)

# ── Basic usage ────────────────────────────────────────────
qobj <- qvalue(p = BLCA_Model$pb11)
BLCA_Model$pb11_qval <- qobj$qvalues

# ── Apply to all p-values ──────────────────────────────────
# pb11 & pb21 together
combined_1       <- c(BLCA_Model$pb11, BLCA_Model$pb21)
qobj_1           <- qvalue(p = combined_1)
n                <- nrow(BLCA_Model)
BLCA_Model$pb11_qval <- qobj_1$qvalues[1:n]
BLCA_Model$pb21_qval <- qobj_1$qvalues[(n+1):(2*n)]

# pb12 & pb22 together
combined_2       <- c(BLCA_Model$pb12, BLCA_Model$pb22)
qobj_2           <- qvalue(p = combined_2)
BLCA_Model$pb12_qval <- qobj_2$qvalues[1:n]
BLCA_Model$pb22_qval <- qobj_2$qvalues[(n+1):(2*n)]

# pV1:T1 & pV1:T2 together
combined_3       <- c(BLCA_Model$`pV1:T1`, BLCA_Model$`pV1:T2`)
qobj_3           <- qvalue(p = combined_3)
BLCA_Model$`pV1:T1_qval` <- qobj_3$qvalues[1:n]
BLCA_Model$`pV1:T2_qval` <- qobj_3$qvalues[(n+1):(2*n)]

# ── Re-infer binary indicators ─────────────────────────────
alpha <- 0.05
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
cat("\nBH (FDR):\n");  print(table(BLCA_Model$Inferred.Model.fdr))
cat("\nqvalue:\n");    print(table(BLCA_Model$Inferred.Model.qval))

# ── Save ───────────────────────────────────────────────────
fwrite(BLCA_Model,
       file = file.path(output_dir, "LZ_trio_results_BLCA_ALL_qval.txt"),
       sep  = "\t")

########################################################################################
