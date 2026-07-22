# ============================================================
# Step 0: Load required libraries
# ============================================================
# library (devtools)
# install_github("audreyqyfu/MRTrios")
library(MRTrios)
library(data.table)
library(tidyverse)
library(MRGN)
# ============================================================
# Step 1: Load data Files
# ============================================================

# ----- set the directory -----
setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/raw_Data_Methyl")

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

models_head100 <- do.call(rbind, results_list)


saveRDS(models_head100, "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026/Infer_models_head100.rds")

