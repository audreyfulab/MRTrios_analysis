
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

# ============================ OUTPUT ======================================================
infer_model_trio_180593_noU
  b11 b12 b21 b22 V1:T1 V1:T2       pb11         pb12       pb21         pb22    pV1:T1   pV1:T2
1   0   1   0   1     0     0 0.08157669 1.197528e-99 0.02926805 1.197528e-99 0.8519377 0.185065
  Minor.freq Inferred.Model
1  0.3262411          Other

### compare with the inferred model with confounder  variables
Model_pos=readRDS("LZ_trio_results_pos_ALL_qval.rds")
## Model_pos %>% filter(trio_row==180593)
## with confounder the inferred model is "Other", but it changed to M1.1 after using qValue

Model_pos %>% filter(trio_row==180593)
   trio_row meth.row cna.row gene.row Total.PC.Count   b11   b12   b21   b22 V1:T1 V1:T2       pb11
      <int>    <int>   <int>    <int>          <int> <int> <int> <int> <int> <int> <int>      <num>
1:   180593   240726    1349    13218             32     0     1     0     1     0     0 0.05548527
          pb12      pb21        pb22    pV1:T1   pV1:T2 Minor.freq Inferred.Model  pb11_qval   pb12_qval
         <num>     <num>       <num>     <num>    <num>      <num>         <char>      <num>       <num>
1: 0.000999001 0.1511994 0.000999001 0.8519377 0.185065 0.04787234          Other 0.02389754 0.001671113
    pb21_qval   pb22_qval pV1:T1_qval pV1:T2_qval b11_qval b12_qval b21_qval b22_qval V1:T1_qval
        <num>       <num>       <num>       <num>    <int>    <int>    <int>    <int>      <int>
1: 0.05719164 0.001671113   0.2337852   0.0669408        1        1        0        1          0
   V1:T2_qval Inferred.Model.qval
        <int>              <char>
1:          0                M1.1

# ========================================================================================

### Write a Function to compare more trios' result:
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

# ============================ OUTPUT ======================================================
  trio_row meth.row gene.row cna.row model_noU model_withU model_withU_qval
1    180593   240726    13218    1349     Other       Other             M1.1
2    151855    29318     5201    3010     Other       Other             M1.1
3      8201   270287    14391    9918     Other       Other             M1.1
4    123253   390588     2094   21411     Other       Other             M1.1
5    123249    62115     2094   21411     Other       Other             M1.1
6    135220   440845    17112   13618      M2.1        M1.1             M1.1
7    123250   149374     2094   21411     Other       Other             M1.1
8    270299   382681     5103     615     Other        M1.1             M1.1
9    112467   151299    12984    3281     Other       Other             M1.1
10     8658   472810    14391    9918      M1.2       Other             M1.1
11   277295   156032     8255   16790     Other        M1.1             M1.1
12   200834    26265    13940    5540      M2.1       Other             M1.1
13   233022    32849    17018   23492     Other        M1.1             M1.1
14    56542    82372      400    3179      M1.2        M1.1             M1.1
15     8117   230419    14391    9918     Other       Other             M1.1
16     8088   217199    14391    9918     Other       Other             M1.1
17    56539    68672      400    3179     Other        M1.1             M1.1
18   123252   219458     2094   21411     Other       Other             M1.1
19   270296   296528     5103     615     Other        M1.1             M1.1
20     7700    42674    14391    9918     Other       Other             M1.1                                 
# ==================================================================================
                                 
# Compare where the no-confounder model disagrees with the confounder model
results %>% filter(model_noU != model_withU)
                                 
# === OUTPUT ========================================================================
  trio_row meth.row gene.row cna.row model_noU model_withU model_withU_qval
1   135220   440845    17112   13618      M2.1        M1.1             M1.1
2   270299   382681     5103     615     Other        M1.1             M1.1
3     8658   472810    14391    9918      M1.2       Other             M1.1
4   277295   156032     8255   16790     Other        M1.1             M1.1
5   200834    26265    13940    5540      M2.1       Other             M1.1
6   233022    32849    17018   23492     Other        M1.1             M1.1
7    56542    82372      400    3179      M1.2        M1.1             M1.1
8    56539    68672      400    3179     Other        M1.1             M1.1
9   270296   296528     5103     615     Other        M1.1             M1.1                                 
# ==================================================================================
