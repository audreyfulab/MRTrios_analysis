########################################################################################
###################      FDR    ########################################################
########################################################################################

LIHC_Model <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/Output_LIHC/trio_results_LIHC_ALL.txt",sep ="\t")

## ── LIHC model FDR ──────────────────────────────────────────────────
output_dir ="/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_LIHC/"

# ── Install if needed ──────────────────────────────────────
if (!requireNamespace("qvalue", quietly = TRUE)) {
  BiocManager::install("qvalue")
}
library(qvalue)

# ── Basic usage ────────────────────────────────────────────
qobj <- qvalue(p = LIHC_Model$pb11)
LIHC_Model$pb11_qval <- qobj$qvalues

# ── Apply to all p-values ──────────────────────────────────
# pb11 & pb21 together
combined_1       <- c(LIHC_Model$pb11, LIHC_Model$pb21)
qobj_1           <- qvalue(p = combined_1)
n                <- nrow(LIHC_Model)
LIHC_Model$pb11_qval <- qobj_1$qvalues[1:n]
LIHC_Model$pb21_qval <- qobj_1$qvalues[(n+1):(2*n)]

# pb12 & pb22 together
combined_2       <- c(LIHC_Model$pb12, LIHC_Model$pb22)
qobj_2           <- qvalue(p = combined_2)
LIHC_Model$pb12_qval <- qobj_2$qvalues[1:n]
LIHC_Model$pb22_qval <- qobj_2$qvalues[(n+1):(2*n)]

# pV1:T1 & pV1:T2 together
combined_3       <- c(LIHC_Model$`pV1:T1`, LIHC_Model$`pV1:T2`)
qobj_3           <- qvalue(p = combined_3)
LIHC_Model$`pV1:T1_qval` <- qobj_3$qvalues[1:n]
LIHC_Model$`pV1:T2_qval` <- qobj_3$qvalues[(n+1):(2*n)]

# ── Re-infer binary indicators ─────────────────────────────
alpha <- 0.05
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
cat("\nBH (FDR):\n");  print(table(LIHC_Model$Inferred.Model.fdr))
cat("\nqvalue:\n");    print(table(LIHC_Model$Inferred.Model.qval))

# ── Save ───────────────────────────────────────────────────
fwrite(LIHC_Model,
       file = file.path(output_dir, "LZ_trio_results_LIHC_ALL_qval.txt"),
       sep  = "\t")

########################################################################################
