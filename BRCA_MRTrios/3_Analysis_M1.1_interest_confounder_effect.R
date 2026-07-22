
# ==============================================================================
# Check m1.1_interest examples
# (confounder effect)
# ==============================================================================
library(tidyverse)
library(data.table)
setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026")
m1.1_interest <- fread("m1.1_interest.txt")

## Drop some cols,remove dup rows
m1.1_filter <- m1.1_interest %>% select(-c(Group,b11,b12,b21,b22,`V1:T1`,`V1:T2`,pb11,pb12,pb21,pb22,`pV1:T1`,`pV1:T2`,Inferred.Model)) %>% distinct()
### Head 6 examples
head(m1.1_filter %>% select(c(trio_row, meth.row, cna.row, gene.row,)))
# trio_row meth.row cna.row gene.row
# <int>    <int>   <int>    <int>
#   1:   180593   240726    1349    13218
# 2:   151855    29318    3010     5201
# 3:     8201   270287    9918    14391
# 4:   123253   390588   21411     2094
# 5:   123249    62115   21411     2094
# 6:   135220   440845   13618    17112

# ========================================================================================

setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026")
pos_meth <- readRDS("pos_meth.rds")
pos_cna <- readRDS("pos_cna.rds")
pos_exp <- readRDS("pos_exp.rds")

trio_180593=m1.1_filter$trio_row[1]
meth_240726=m1.1_filter$meth.row[1]
gene_13218=m1.1_filter$gene.row[1]
cna_1349=m1.1_filter$cna.row[1]


M_240726 <- pos_meth %>% filter(meth.row==240726)
C_1349 <- pos_cna %>% filter(cna.row==1349)
E_13218  <- pos_exp %>% filter(gene.row==13218 )

df_trio_180593_noU <- cbind(
  C  = t(C_1349[,-1]) ,
  E  = t(E_13218[,-1]),
  M  = t(M_240726[,-1]))

df_trio_180593_noU <- as.data.frame(df_trio_180593_noU)

library(MRGN)
### check the inferred model without confounder variables
infer_model_trio_180593_noU =infer.trio(df_trio_180593_noU)
infer_model_trio_180593_noU$Inferred.Model ## "Other"

### compare with the inferred model with confounder  variables
Model_pos=readRDS("LZ_trio_results_pos_ALL_qval.rds")
Model_pos %>% filter(trio_row==180593)
## with confounder the inferred model is "Other", but it changed to M1.1 after using qValue

# ========================================================================================

check_trio <- function(i, m1.1_filter, pos_meth, pos_cna, pos_exp, Model_pos) {
  meth_r <- m1.1_filter$meth.row[i]
  gene_r <- m1.1_filter$gene.row[i]
  cna_r  <- m1.1_filter$cna.row[i]
  trio_r <- m1.1_filter$trio_row[i]
  
  M <- pos_meth %>% filter(meth.row == meth_r)
  C <- pos_cna  %>% filter(cna.row  == cna_r)
  E <- pos_exp  %>% filter(gene.row == gene_r)
  
  df <- as.data.frame(cbind(
    C = t(C[,-1]),
    E = t(E[,-1]),
    M = t(M[,-1])
  ))
  
  fit <- infer.trio(df)
  
  data.frame(
    trio_row       = trio_r,
    meth.row       = meth_r,
    gene.row       = gene_r,
    cna.row        = cna_r,
    model_noU      = fit$Inferred.Model,
    model_withU    = Model_pos$Inferred.Model[match(trio_r, Model_pos$trio_row)],
    model_withU_qval  = Model_pos$Inferred.Model.qval[match(trio_r, Model_pos$trio_row)],
    stringsAsFactors = FALSE
  )
}

# Run over all trios (or a subset, e.g. 1:20)
# results <- do.call(rbind, lapply(seq_len(nrow(m1.1_filter)), function(i)
#   check_trio(i, m1.1_filter, pos_meth, pos_cna, pos_exp, Model_pos)))

results <- do.call(rbind, lapply(1:20, function(i)
  check_trio(i, m1.1_filter, pos_meth, pos_cna, pos_exp, Model_pos)))

results

# Compare where the no-confounder model disagrees with the confounder model
results %>% filter(model_noU != model_withU)

