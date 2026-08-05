# ========================================================================================================
# Load posER_BRCA trios Model result
# ========================================================================================================

## Without FDR
Model_posER_BRCA <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/trio_posER_BRCA_results_ALL.txt")

## With FDR
Model_posER_BRCA <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/posER_BRCA_trio_Model_results_ALL_with_BH_fdr_qval_byLZ.txt")

# ========================================================================================================
# Figure out which trios Failed to infer models (C,E,M data Missing from raw data)
# ========================================================================================================

## Failed: [WARN] trio 277040: C/E/M did not match exactly one row each (cna=0, exp=1, meth=1); skipping trio.
## Succeed: trio 277041: running infer.trio() with confounders: age, race, U_exp, U_meth

range(Model_posER_BRCA$trio_row)
[1]      1 295492
nrow(Model_posER_BRCA)
[1] 292469

missing_trios <- setdiff(1:295492, Model_posER_BRCA$trio_row)
length(missing_trios)   # how many are missing
[1] 3023

head(missing_trios)     # a peek at which ones
[1] 9219 9220 9221 9222 9223 9224

missing_sorted <- sort(missing_trios)

# group into runs of consecutive integers
run_id <- cumsum(c(1, diff(missing_sorted) != 1))
runs <- split(missing_sorted, run_id)

run_summary <- data.frame(
     start  = sapply(runs, min),
     end    = sapply(runs, max),
     length = sapply(runs, length)
 )

run_summary[order(-run_summary$length), ]

     start    end length
8    29322  29495    174
31   82271  82409    139
7    28216  28351    136
