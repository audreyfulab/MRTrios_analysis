# ============================================================
# Step 0: Load required libraries
# ============================================================
# library (devtools)
# devtools::install_github("audreyqyfu/MRTrios")
# devtools::install_github("audreyfulab/MRGN")

library(MRTrios)
library(data.table)
library(tidyverse)
library(MRGN)
# ============================================================
# Step 1: Load data Files
# ============================================================

# ----- set the directory -----
 setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/raw_Data_Methyl") ### Run Locally
#setwd("/wsu/home/hb/hb68/hb6890/fulab/MRTrios_analysis") ### Run on Server

# ----- load Main data files -----
cna <- fread('data_CNA.txt',data.table = F)
gene.exp <- fread('data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt',data.table = F)
TCGA.meth <- fread("split.names.TCGA.meth.logit.txt", data.table=F)
trios <- fread("trio.final.protein.coding.txt")

# ================ Load neg ER+ patient data ================
clinical.pos  <- fread("names.pos.patient2.txt", header = FALSE) # Clinical ER+ patients ID

pc.gene.pos <- fread("PCA.gene.exp.posER.txt") ## gene.exp PC scores
pc.meth.pos <- fread("PCA.meth.posER.txt")  ## meth PC scores

### Note:meth,gene.exp table : # V1,col1, col2 , the col2 is the final pcscore index.
gene.table.pos <- fread("gene.exp.posER.table.txt") ## gene.exp PC raw and corrected index
meth.table.pos <- fread("meth.posER.table.txt")# meth indices values

gene_indice.pos  <- readRDS("gene_exp_posER_sig_pcs_indices.rds") # gene.exp indices values
meth_indice.pos  <- readRDS("meth_posER_sig_pcs_indices.rds") ## meth PC raw and corrected index

# length 787: meth, gene.exp
# length 777: Finding common individuals among the 3 datasets:cna,gene.exp, meth
# length 572: meth, gene.exp,clinical.pos
# lenth  564: cna,meth, gene.exp,clinical.pos
model_pos <- read.delim("model.trio.MRGN.all.posER.reclassify2.txt", header = TRUE)

# ===================  Load neg ER- patient data ================

clinical.neg  <- fread("names.neg.patient2.txt", header = FALSE) # ER- patients ID

pc.gene.neg <- fread("PCA.gene.exp.negER.txt")
pc.meth.neg <- fread("PCA.meth.negER.txt")

### meth,gene.exp table : # V1,col1, col2 , the col2 is the final pcscore index.
meth.table.neg <- fread("meth.negER.table.txt")
gene.table.neg <- fread("gene.exp.negER.table.txt")

gene_indice.neg <- readRDS("gene_exp_negER_sig_pcs_indices.rds")
meth_indice.neg <- readRDS("meth_negER_sig_pcs_indices.rds")

model_neg <- read.delim("model.trio.MRGN.all.negER.reclassify2.txt", header = TRUE)

# ===================   Usage  ===================  
setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/Methyl_Analysis_LZ/Data_Result_Methyl_LZ")

# final_pos_M1.1 <- fread("final_pos_M1.1.txt",sep = "\t")

# cor_M1.1=final_pos_M1.1  %>% 
#   filter(abs(cor_meth_exp) > 0.5 & abs(cor_cna_exp) > 0.5) %>%
#   arrange(desc(abs(cor_meth_exp)), desc(abs(cor_cna_exp)))

final_pos_M1.2<- fread("final_pos_M1.2.txt",sep = "\t")

cor_M1.2=final_pos_M1.2  %>% 
  filter(abs(cor_meth_exp) > 0.5 & abs(cor_meth_cna) > 0.3) %>%
  arrange(desc(abs(cor_meth_exp)), desc(abs(cor_meth_cna)))

corM1.2_models_11 <- get_trio_models(cor_M1.2$pos_trios_row, 
                                     trios, cna, gene.exp, TCGA.meth,
                                     clinical = clinical.pos,
                                     gene_indice = gene_indice.pos,
                                     meth_indice = meth_indice.pos,
                                     pc.gene = pc.gene.pos,
                                     pc.meth = pc.meth.pos,
                                     gene.table = gene.table.pos,
                                     meth.table = meth.table.pos)

print(corM1.2_models_11)
### check corM1.2
print(cor_M1.2)

#### =============================================================================

### Check M3 model:




#### =============================================================================

### Test to infer head10 pos models
pos_models_head10 <- get_trio_models(1:10, 
                                     trios, cna, gene.exp, TCGA.meth,
                                     clinical = clinical.pos,
                                     gene_indice = gene_indice.pos,
                                     meth_indice = meth_indice.pos,
                                     pc.gene = pc.gene.pos,
                                     pc.meth = pc.meth.pos,
                                     gene.table = gene.table.pos,
                                     meth.table = meth.table.pos)
print(pos_models_head10) ## new inferred
old_model_pos_head10 <- model_pos %>% arrange(trio_row) %>% slice(1:10) ## old inferred_Model3

### write new function: get_model,set is.CNA = TRUE,nperms=100
pos_models_head10_v2 <- get_trio_models_V2(1:10, 
                                     trios, cna, gene.exp, TCGA.meth,
                                     clinical = clinical.pos,
                                     gene_indice = gene_indice.pos,
                                     meth_indice = meth_indice.pos,
                                     pc.gene = pc.gene.pos,
                                     pc.meth = pc.meth.pos,
                                     gene.table = gene.table.pos,
                                     meth.table = meth.table.pos)

### compare 3 inferred version:
print(pos_models_head10) ## new inferred with is.CNA = TRUE,nperms=1000
print(pos_models_head10_v2) ## new inferred with is.CNA = TRUE,nperms=100
print(old_model_pos_head10)  ## old inferred without set is.CNA 





### Run to infer the whole pos/neg models =========
### set trio_row = 1:nrow(trios)

# ── run head 100 ────────────────────────────────────────────────────────────────────
chunk_size <- 20
chunks <- split(101:10000, ceiling((101:10000) / chunk_size))
# chunks = list(1:20, 21:40, 41:60, 61:80, 81:100)

results_list <- lapply(chunks, function(idx) {
  message("Running rows ", idx[1], " to ", tail(idx, 1))
  get_trio_models(idx,
                  trios, cna, gene.exp, TCGA.meth,
                  clinical,
                  gene_indice,
                  meth_indice,
                  pc.gene,
                  pc.meth,
                  gene.table,
                  meth.table)
})

neg_models <- do.call(rbind, results_list)


saveRDS(neg_models, "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/Infer_negER_models_head101_10000.rds")


# ============================================================
# Step 3: MRTrios function : get_trio_models()
# ============================================================
# ----- load Main data files -----
cna <- fread('data_CNA.txt',data.table = F)
gene.exp <- fread('data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt',data.table = F)
TCGA.meth <- fread("split.names.TCGA.meth.logit.txt", data.table=F)
trios <- fread("trio.final.protein.coding.txt")

# ================ Load neg ER+ patient data ================
clinical  <- fread("names.pos.patient2.txt", header = FALSE) # Clinical ER+ patients ID

pc.gene <- fread("PCA.gene.exp.posER.txt") ## gene.exp PC scores
pc.meth <- fread("PCA.meth.posER.txt")  ## meth PC scores

### Note:meth,gene.exp table : # V1,col1, col2 , the col2 is the final pcscore index.
gene.table <- fread("gene.exp.posER.table.txt") ## gene.exp PC raw and corrected index
meth.table <- fread("meth.posER.table.txt")# meth indices values

gene_indice  <- readRDS("gene_exp_posER_sig_pcs_indices.rds") # gene.exp indices values
meth_indice  <- readRDS("meth_posER_sig_pcs_indices.rds") ## meth PC raw and corrected index

get_trio_models <- function(trio_row,
                            trios,
                            cna,
                            gene.exp,
                            TCGA.meth,
                            clinical,
                            gene_indice,
                            meth_indice,
                            pc.gene,
                            pc.meth,
                            gene.table,
                            meth.table) {
  
  assign("isTrue", isTRUE, envir = .GlobalEnv)
  
  # ── Step 1: Find common individuals across all 4 datasets ──────────────────
  
  com.ind <- intersect(
    intersect(colnames(gene.exp)[3:ncol(gene.exp)],
              colnames(TCGA.meth)[5:ncol(TCGA.meth)]),
    colnames(cna)[3:ncol(cna)])
  
  unique_id <- intersect(com.ind, clinical$V1)
  
  # ── Step 2: Extract clinical covariates ────────────────────────────────────
  
  clinical_filter <- clinical %>%
    arrange(match(V1, unique_id)) %>%
    filter(V1 %in% unique_id)
  
  age  <- clinical_filter$V2
  race <- clinical_filter$V3
  
  # ── Step 3: Filter datasets to common individuals ──────────────────────────
  
  cna_filter <- cna %>%
    mutate(cna.row = rownames(cna)) %>%
    filter(cna.row %in% unique(trios$cna.row)) %>%
    select(cna.row, all_of(unique_id))
  
  exp_filter <- gene.exp %>%
    mutate(gene.row = rownames(gene.exp)) %>%
    filter(gene.row %in% unique(trios$gene.row)) %>%
    select(gene.row, all_of(unique_id))
  
  meth_filter <- TCGA.meth %>%
    mutate(meth.row = rownames(TCGA.meth)) %>%
    filter(meth.row %in% unique(trios$meth.row)) %>%
    select(meth.row, all_of(unique_id))
  
  # ── Step 4: Align PCA score matrices to common individuals ─────────────────
  
  pca_E_final <- pc.gene %>%
    filter(V1 %in% unique_id) %>%
    arrange(match(V1, unique_id)) %>%
    select(-V1)
  
  pca_E_final <- pc.gene %>%
    filter(V1 %in% unique_id) %>%
    arrange(match(V1, unique_id)) %>%
    select(-V1)
  
  pca_M_final <- pc.meth %>%
    filter(V1 %in% unique_id) %>%
    arrange(match(V1, unique_id)) %>%
    select(-V1)
  
  # ── Step 5: Loop over trio indices and infer models ────────────────────────
  
  do.call(rbind, lapply(trio_row, function(i) {
    
    trio     <- trios[i, ]
    meth_row <- trio$meth.row
    cna_row  <- trio$cna.row
    gene_row <- trio$gene.row
    
    # Extract relevant rows from each dataset
    M <- meth_filter %>% filter(meth.row == meth_row)
    C <- cna_filter  %>% filter(cna.row  == cna_row)
    E <- exp_filter  %>% filter(gene.row == gene_row)
    
    # Retrieve significant PCs for expression
    gene_pc_index <- gene.table %>% filter(col1 == gene_row) %>% pull(col2)
    U_exp <- pca_E_final %>% select(all_of(as.numeric(gene_indice[gene_pc_index][[1]])))
    
    # gene_pc_index <- gene.table %>% filter(col1 == gene_row) %>% pull(col2)
    # gene_indices <- gene_indice[gene_pc_index][[1]]
    # gene_indices_num <- as.numeric(gene_indices)
    # U_exp <- pca_E_final %>% select(all_of(gene_indices_num))
    # 
    # Retrieve significant PCs for methylation
    meth_pc_index <- meth.table %>% filter(col1 == meth_row) %>% pull(col2)
    U_meth <- pca_M_final %>% select(all_of(as.numeric(meth_indice[meth_pc_index][[1]])))
    
    Total.PC.Count <- sum(ncol(U_exp), ncol(U_meth), na.rm = TRUE)
    
    # Build regression data frame
    df_trio <- as.data.frame(cbind(
      C = t(C[, -1]),
      E = t(E[, -1]),
      M = t(M[, -1]),
      age,
      race,
      U_exp,
      U_meth
    ))
    
    colnames(df_trio)[1:3] <- c("C", "E", "M")
    df_trio$race <- as.factor(df_trio$race)
    
    # Infer causal model ## default  alpha=0.01, nperms=1000
    model <- as.data.frame(infer.trio(df_trio,is.CNA = TRUE))
    
    # Return one row of results
    cbind(
      data.frame(
        trio_row       = i,
        meth.row       = meth_row,
        cna.row        = cna_row,
        gene.row       = gene_row,
        Total.PC.Count = Total.PC.Count
      ),
      model
    )
  }))
}



#### ================= ER pos ===================
trio_pos_ALL <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/trio_results_pos_ALL.txt")
a1=trio_pos_ALL %>% filter(is.na(cna.row)|is.na(gene.row))
dim(trio_pos_ALL)## 294913 
dim(a1)## 2417
trio_pos_model <- anti_join(trio_pos_ALL,a1)### remove (cna.row==NA) inferred results
dim(trio_pos_model) ## 292496 

table(trio_pos_model$Inferred.Model)
dim(model_pos) ### 292465
table(model_pos$Inferred.Model)
table(model_pos$Inferred.Model3)


trio_dif <- setdiff(trio_pos_model$trio_row,model_pos$trio_row)
length(trio_dif)
df1=trio_pos_model %>% filter(trio_row %in% trio_dif)
trios[df1$trio_row,]

aa1=pos_exp_filter %>% filter(gene.row==9100)
aa2=pos_meth_filter %>% filter(meth.row %in% df1[21:24]$meth.row)

hist(c(trio_pos_model$pb11,trio_pos_model$pb21))
hist(c(trio_pos_model$`pV1:T1`,trio_pos_model$`pV1:T2`))

trio_pos_model %>% filter(is.na())
# NAs per column
colSums(is.na(trio_pos_model)) ###292496
LZ_pos_model <- trio_pos_model %>% filter(!is.na(`V1:T1`)). ###292453
fwrite(LZ_pos_model,"/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/LZ_pos_model.txt",sep ="\t")

## Rows in LZ_pos_model but NOT in model_pos:
setdiff(LZ_pos_model$trio_row, model_pos$trio_row)

## Rows in model_pos but NOT in LZ_pos_model:
setdiff(model_pos$trio_row, LZ_pos_model$trio_row)

missing=setdiff(model_pos$trio_row, LZ_pos_model$trio_row)

### missing trios:  291020 291021 291022 291023 291189 291190 291191 291192 291193 291194 
###                 294601 295282 295283 295388 259513 259514

trios[missing,]


df_trio_291020 <- get_trio_data(291020, trios, pos_meth_filter, pos_cna_filter,
                                 pos_exp_filter, pca_E_final, pca_M_final,
                                 gene_indice.pos, meth_indice.pos,gene.table.pos,meth.table.pos,
                                 pos_age, pos_race)
library(MRGN)
model_trio_291020_with_U <- infer.trio(df_trio_291020,is.CNA = TRUE,compute.nominal=F)
## the result is the same as using get_trio_models( ) from MRTrios

mis_meth_c=meth.table.pos %>% filter(col1 %in% trios[missing,]$meth.row)
mis_gene_c=gene.table.pos %>% filter(col1 %in% trios[missing,]$gene.row)

mis_meth=pos_meth_filter %>% filter(meth.row %in% mis_meth_c$col2)
mis_gene=pos_exp_filter %>% filter(gene.row %in% mis_gene_c$col2)

pos_exp_filter[12330,]
colSums(is.na(mis_meth))

df_a1=c(291020:291023,291189:291194,294601,295282 ,295283 ,295388 ,259513 ,259514)
results_df_a1 <- get_trio_models(
  trio_row       = df_a1,
  trios          = trios,
  cna            = cna,
  gene.exp       = gene.exp,
  TCGA.meth      = TCGA.meth,
  clinical   = clinical.pos,
  gene_indice = gene_indice.pos,
  meth_indice = meth_indice.pos,
  pc.gene    = pc.gene.pos,
  pc.meth    = pc.meth.pos,
  gene.table = gene.table.pos,
  meth.table = meth.table.pos
)
fwrite(results_df_a1,"/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/missing_16trios_pos_model.txt",sep ="\t")

LZ_pos_model <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/LZ_pos_model.txt")
missing_16trios_pos_model <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/missing_16trios_pos_model.txt")
LZ_final_pos_model= rbind(LZ_pos_model,missing_16trios_pos_model) %>% arrange(trio_row)
fwrite(LZ_final_pos_model,"/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/LZ_pos_model_final_fix.txt",sep ="\t")


### there are # NA positions (non-null but all NA) in gene_indice.pos
na_pos = which(sapply(gene_indice.pos, function(x) !is.null(x) && all(is.na(x))))
na_pos
# [1]    15  1624  4268  4806  5737  6159  6911  9271  9274  9278  9750 10051 10122 11631 12005 12323 12330 12363
# [19] 12466 12472 12564 12607 12652 12658 13921 16508 16509 16531

which(sapply(meth_indice.pos, function(x) !is.null(x) && all(is.na(x))))
# [1] 396024 396060
which(sapply(gene_indice.pos, function(x) !is.null(x) && all(is.na(x))))
# [1]    15  1624  4268  4806  5737  6159  6911  9271  9274  9278  9750 10051 10122 11631 12005 12323 12330 12363
# [19] 12466 12472 12564 12607 12652 12658 13921 16508 16509 16531

########## example: trio==291020  POS======================
########## it not infered in LZ code, but they are infered by previous model:
model_pos %>% filter(trio_row==291020)
trios[291020,]
gene.table.pos %>% filter(col1==12378) %>% pull(col2)
# 12323
meth.table.pos %>% filter(col1==137485) %>% pull(col2)
# 112570
gene_indice.pos[[12323]] ## there is null values in gene_indice.pos[[12323]]

meth_indice.pos[[112570]]

#### ================= ER neg ===================
trio_results_neg_ALL <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/trio_results_neg_ALL.txt")
neg_mis=trio_results_neg_ALL %>% filter(is.na(cna.row)|is.na(gene.row))###2189

trio_neg_model <- anti_join(trio_results_neg_ALL,neg_mis)### remove (cna.row==NA) inferred results
dim(trio_neg_model) ## 290235  ##previous:292331

LZ_neg_model <- trio_neg_model %>% filter(!is.na(`V1:T1`)) ###290235
fwrite(LZ_neg_model,"/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/LZ_neg_model.txt",sep ="\t")

a1 <- which(sapply(meth_indice.neg, function(x) !is.null(x) && all(is.na(x))))### 266
a2 <- which(sapply(gene_indice.neg, function(x) !is.null(x) && all(is.na(x))))## 525
missing_meth_row <- meth.table.neg %>% filter(col2 %in% a1) %>% pull(col1)
missing_gene_row <- gene.table.neg %>% filter(col2 %in% a2) %>% pull(col1)
## 2296
missing_neg_trios <- trios %>% mutate(trio_row=1:nrow(trios)) %>% filter(meth.row %in% missing_meth_row|gene.row %in% missing_gene_row)
missing_neg_trios_final <- missing_neg_trios %>% filter(!is.na(cna.row))## 2101
fwrite(missing_neg_trios_final,"/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/missing_2101trios_neg.txt",sep ="\t")



##2101
neg_miss_2=setdiff(model_neg$trio_row, trio_results_neg_ALL$trio_row)

length(setdiff(neg_miss_2,missing_neg_trios_final$trio_row)) ##0 

sum(is.na(gene.table.neg$col2))
sum(is.na(meth.table.neg$col2))

LZ_neg_model <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/LZ_neg_model.txt")
neg_miss_2=setdiff(model_neg$trio_row, LZ_neg_model$trio_row)### 2101
neg_trio_miss <- trios %>% slice(neg_miss_2)
results_neg_miss_2 <- get_trio_models(
  trio_row       = neg_miss_2,
  trios          = trios,
  cna            = cna,
  gene.exp       = gene.exp,
  TCGA.meth      = TCGA.meth,
  clinical   = clinical.pos,
  gene_indice = gene_indice.pos,
  meth_indice = meth_indice.pos,
  pc.gene    = pc.gene.pos,
  pc.meth    = pc.meth.pos,
  gene.table = gene.table.pos,
  meth.table = meth.table.pos
)
# fwrite(results_neg_miss_2,"/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/missing_2101trios_neg_model.txt",sep ="\t")
# trio_neg_miss <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/missing_2101trios_neg_model.txt")
missing_2101trios_neg <- fread("/wsu/home/hb/hb68/hb6890/fulab/MRTrios/raw_Data_Methyl/missing_2101trios_neg.txt")
neg_2101miss_head500 <- get_trio_models(
  trio_row       = missing_2101trios_neg$trio_row[1:500],
  trios          = trios,
  cna            = cna,
  gene.exp       = gene.exp,
  TCGA.meth      = TCGA.meth,
  clinical   = clinical.pos,
  gene_indice = gene_indice.pos,
  meth_indice = meth_indice.pos,
  pc.gene    = pc.gene.pos,
  pc.meth    = pc.meth.pos,
  gene.table = gene.table.pos,
  meth.table = meth.table.pos
)
fwrite(neg_2101miss_head500,"/wsu/home/hb/hb68/hb6890/fulab/MRTrios/Output_neg/head500_missing_2101trios_neg_model.txt",sep ="\t"))

# neg_missed_models <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/trio_results_neg_missed.txt")
neg_missed_models_v2 <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/trio_results_neg_missed_v2.txt")
## combine the inferred + missed model result 
df0= neg_missed_models_v2 %>% select(-trio_row)
df_trio= trios %>% mutate(trio_row=1:nrow(trios))

df1=df_trio %>% filter(meth.row %in% df0$meth.row & cna.row %in% df0$cna.row & gene.row %in% df0$gene.row)
dim(df1)
identical(df0$meth.row,df1$meth.row)
neg_missed_models_v2$trio_row=df1$trio_row

LZ_final_neg_model <- rbind(LZ_neg_model,neg_missed_models_v2)
fwrite(LZ_final_neg_model,"/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/LZ_final_neg_model.txt",sep ="\t")


##


########## example: trio==274100 neg======================
########## it not infered in LZ code, but they are infered by previous model:

model_neg %>% filter(trio_row==274100)
trios[274100,]
 gene.table.neg %>% filter(col1==5930) %>% pull(col2)
 # 5889
 meth.table.neg %>% filter(col1==160435) %>% pull(col2)
# 131306
gene_indice.neg[[5889]] ## there is null values in gene_indice.neg[[5889]]

meth_indice.neg[[131306]]

########################################################################################
###################      FDR    ########################################################
########################################################################################

NegModel <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/LZ_final_neg_model.txt",sep ="\t")
PosModel <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/LZ_pos_model_final_fix.txt",sep ="\t")

## ── Pos model FDR ──────────────────────────────────────────────────
output_dir ="/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/"

# ── Install if needed ──────────────────────────────────────
if (!requireNamespace("qvalue", quietly = TRUE)) {
  BiocManager::install("qvalue")
}
library(qvalue)

# ── Basic usage ────────────────────────────────────────────
qobj <- qvalue(p = LZ_final_pos_model$pb11)
LZ_final_pos_model$pb11_qval <- qobj$qvalues

# ── Apply to all p-values ──────────────────────────────────
# pb11 & pb21 together
combined_1       <- c(LZ_final_pos_model$pb11, LZ_final_pos_model$pb21)
qobj_1           <- qvalue(p = combined_1)
n                <- nrow(LZ_final_pos_model)
LZ_final_pos_model$pb11_qval <- qobj_1$qvalues[1:n]
LZ_final_pos_model$pb21_qval <- qobj_1$qvalues[(n+1):(2*n)]

# pb12 & pb22 together
combined_2       <- c(LZ_final_pos_model$pb12, LZ_final_pos_model$pb22)
qobj_2           <- qvalue(p = combined_2)
LZ_final_pos_model$pb12_qval <- qobj_2$qvalues[1:n]
LZ_final_pos_model$pb22_qval <- qobj_2$qvalues[(n+1):(2*n)]

# pV1:T1 & pV1:T2 together
combined_3       <- c(LZ_final_pos_model$`pV1:T1`, LZ_final_pos_model$`pV1:T2`)
qobj_3           <- qvalue(p = combined_3)
LZ_final_pos_model$`pV1:T1_qval` <- qobj_3$qvalues[1:n]
LZ_final_pos_model$`pV1:T2_qval` <- qobj_3$qvalues[(n+1):(2*n)]

# ── Re-infer binary indicators ─────────────────────────────
alpha <- 0.05
LZ_final_pos_model$b11_qval    <- as.integer(LZ_final_pos_model$pb11_qval    < alpha)
LZ_final_pos_model$b12_qval    <- as.integer(LZ_final_pos_model$pb12_qval    < alpha)
LZ_final_pos_model$b21_qval    <- as.integer(LZ_final_pos_model$pb21_qval    < alpha)
LZ_final_pos_model$b22_qval    <- as.integer(LZ_final_pos_model$pb22_qval    < alpha)
LZ_final_pos_model$`V1:T1_qval` <- as.integer(LZ_final_pos_model$`pV1:T1_qval` < alpha)
LZ_final_pos_model$`V1:T2_qval` <- as.integer(LZ_final_pos_model$`pV1:T2_qval` < alpha)

# ── Re-infer model ─────────────────────────────────────────
LZ_final_pos_model$Inferred.Model.qval <- apply(LZ_final_pos_model, 1, function(row) {
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
cat("Original:\n");    print(table(LZ_final_pos_model$Inferred.Model))
cat("\nBH (FDR):\n");  print(table(LZ_final_pos_model$Inferred.Model.fdr))
cat("\nqvalue:\n");    print(table(LZ_final_pos_model$Inferred.Model.qval))

# ── Save ───────────────────────────────────────────────────
fwrite(LZ_final_pos_model,
       file = file.path(output_dir, "LZ_trio_results_pos_ALL_qval.txt"),
       sep  = "\t")

########################################################################################

## ── Neg model FDR ──────────────────────────────────────────────────
output_dir ="/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/"

# ── Install if needed ──────────────────────────────────────
if (!requireNamespace("qvalue", quietly = TRUE)) {
  BiocManager::install("qvalue")
}
library(qvalue)

# ── Basic usage ────────────────────────────────────────────
qobj <- qvalue(p = LZ_final_neg_model$pb11)
LZ_final_neg_model$pb11_qval <- qobj$qvalues

# ── Apply to all p-values ──────────────────────────────────
# pb11 & pb21 together
combined_1       <- c(LZ_final_neg_model$pb11, LZ_final_neg_model$pb21)
qobj_1           <- qvalue(p = combined_1)
n                <- nrow(LZ_final_neg_model)
LZ_final_neg_model$pb11_qval <- qobj_1$qvalues[1:n]
LZ_final_neg_model$pb21_qval <- qobj_1$qvalues[(n+1):(2*n)]

# pb12 & pb22 together
combined_2       <- c(LZ_final_neg_model$pb12, LZ_final_neg_model$pb22)
qobj_2           <- qvalue(p = combined_2)
LZ_final_neg_model$pb12_qval <- qobj_2$qvalues[1:n]
LZ_final_neg_model$pb22_qval <- qobj_2$qvalues[(n+1):(2*n)]

# pV1:T1 & pV1:T2 together
combined_3       <- c(LZ_final_neg_model$`pV1:T1`, LZ_final_neg_model$`pV1:T2`)
qobj_3           <- qvalue(p = combined_3)
LZ_final_neg_model$`pV1:T1_qval` <- qobj_3$qvalues[1:n]
LZ_final_neg_model$`pV1:T2_qval` <- qobj_3$qvalues[(n+1):(2*n)]

# ── Re-infer binary indicators ─────────────────────────────
alpha <- 0.05
LZ_final_neg_model$b11_qval    <- as.integer(LZ_final_neg_model$pb11_qval    < alpha)
LZ_final_neg_model$b12_qval    <- as.integer(LZ_final_neg_model$pb12_qval    < alpha)
LZ_final_neg_model$b21_qval    <- as.integer(LZ_final_neg_model$pb21_qval    < alpha)
LZ_final_neg_model$b22_qval    <- as.integer(LZ_final_neg_model$pb22_qval    < alpha)
LZ_final_neg_model$`V1:T1_qval` <- as.integer(LZ_final_neg_model$`pV1:T1_qval` < alpha)
LZ_final_neg_model$`V1:T2_qval` <- as.integer(LZ_final_neg_model$`pV1:T2_qval` < alpha)

# ── Re-infer model ─────────────────────────────────────────
LZ_final_neg_model$Inferred.Model.qval <- apply(LZ_final_neg_model, 1, function(row) {
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
cat("Original:\n");    print(table(LZ_final_neg_model$Inferred.Model))
cat("\nBH (FDR):\n");  print(table(LZ_final_neg_model$Inferred.Model.fdr))
cat("\nqvalue:\n");    print(table(LZ_final_neg_model$Inferred.Model.qval))

# ── Save ───────────────────────────────────────────────────
fwrite(LZ_final_neg_model,
       file = file.path(output_dir, "LZ_trio_results_neg_ALL_qval.txt"),
       sep  = "\t")


### ================== Read final results ====================
output_dir ="/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/"
neg_model_qval <- fread(file = file.path(output_dir, "LZ_trio_results_neg_ALL_qval.txt"))
pos_model_qval <- fread(file = file.path(output_dir, "LZ_trio_results_pos_ALL_qval.txt"))
