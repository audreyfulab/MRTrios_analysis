
# ============================================================
# Step 0: Load required libraries
# ============================================================
library(data.table)
library(tidyverse)
library(MRGN)
# ============================================================
# Step 1: Load data Files
# ============================================================

# ----- set the directory -----
setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/raw_Data_Methyl")

# ----- load data files -----
gene_indice.pos  <- readRDS("gene_exp_posER_sig_pcs_indices.rds") # gene.exp indices values
meth_indice.pos  <- readRDS("meth_posER_sig_pcs_indices.rds") ## meth PC raw and corrected index
pc.gene.pos <- fread("PCA.gene.exp.posER.txt")
pc.meth.pos <- fread("PCA.meth.posER.txt")

### Note:meth,gene.exp table : # V1,col1, col2 , the col2 is the final pcscore index.
gene.table.pos <- fread("gene.exp.posER.table.txt") ## gene.exp PC raw and corrected index
meth.table.pos <- fread("meth.posER.table.txt")# meth indices values

cna <- fread('data_CNA.txt',data.table = F)
gene.exp <- fread('data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt',data.table = F)
TCGA.meth <- fread("split.names.TCGA.meth.logit.txt", data.table=F)
trios <- fread("trio.final.protein.coding.txt")
clinical.pos  <- fread("names.pos.patient2.txt", header = FALSE) # ER+ patients ID

# length 787: meth, gene.exp
# length 777: Finding common individuals among the 3 datasets:cna,gene.exp, meth
# length 572: meth, gene.exp,clinical.pos
# lenth  564: cna,meth, gene.exp,clinical.pos

# ============================================================
# Step 2: process data
# ============================================================

## ----- extract unique IDs --------------
## com.ind from 3 datasets:cna,gene.exp, meth
com.ind = intersect(
  intersect(colnames(gene.exp)[3:ncol(gene.exp)],
            colnames(TCGA.meth)[5:ncol(TCGA.meth)]),
  colnames(cna)[3:ncol(cna)]) # length 777

## com_id of 4 datasets:cna,gene.exp, meth,clinical 
id_pos_filter <- intersect(com.ind,clinical.pos$V1) # length 564

clinical.pos_filter <- clinical.pos  %>% arrange(match(V1, id_pos_filter)) %>% filter(V1 %in% id_pos_filter)

## ----- extract AGE, RACE --------------
## extract Pos age, race column as confounder 
pos_age <- clinical.pos_filter$V2
pos_race <- clinical.pos_filter$V3


# ------- Filter cna, gene.exp, meth ------------
pos_cna_filter <-  cna %>% mutate(cna.row= rownames(cna)) %>% 
  filter(cna.row %in% unique(trios$cna.row)) %>% 
  select(cna.row,all_of(id_pos_filter))

pos_exp_filter <-  gene.exp %>% mutate(gene.row= rownames(gene.exp)) %>% 
  filter(gene.row %in% unique(trios$gene.row)) %>% 
  select(gene.row,all_of(id_pos_filter))

pos_meth_filter <-  TCGA.meth %>% mutate(meth.row= rownames(TCGA.meth)) %>% 
  filter(meth.row %in% unique(trios$meth.row)) %>% 
  select(meth.row,all_of(id_pos_filter))

# ------- Filter meth_pca and gene_pca of ERpos ------------
pca_E_sorted =pc.gene.pos %>% filter(V1 %in% id_pos_filter) %>% arrange(match(V1, id_pos_filter))
pca_M_sorted =pc.meth.pos %>% filter(V1 %in% id_pos_filter) %>% arrange(match(V1, id_pos_filter))

# ------- keep only pca value, remove the id col ------------
pca_E_final <- pca_E_sorted[,-1]
pca_M_final <- pca_M_sorted[,-1]

# ============================================================
# Step 3: write function : get_trio_data()
# ============================================================
get_trio_data <- function(trio_row, 
                          trios, 
                          pos_meth_filter, 
                          pos_cna_filter, 
                          pos_exp_filter,
                          pca_E_final,
                          pca_M_final,
                          gene_indice.pos,
                          meth_indice.pos,
                          gene.table.pos,
                          meth.table.pos,
                          pos_age,
                          pos_race) {
  
  # --- Extract row indices ---
  trio        <- trios[trio_row, ]
  meth_row    <- trio$meth.row
  cna_row     <- trio$cna.row
  gene_row    <- trio$gene.row
  
  # --- Filter data rows ---
  M <- pos_meth_filter %>% filter(meth.row == meth_row)
  C <- pos_cna_filter  %>% filter(cna.row  == cna_row)
  E <- pos_exp_filter  %>% filter(gene.row == gene_row)
  
  # Retrieve significant PCs for expression
  gene_pc_index <- gene.table.pos %>% filter(col1 == gene_row) %>% pull(col2)
  U_exp <- pca_E_final %>% select(all_of(as.numeric(gene_indice.pos[gene_pc_index][[1]])))
  
  # Retrieve significant PCs for methylation
  meth_pc_index <- meth.table.pos %>% filter(col1 == meth_row) %>% pull(col2)
  U_meth <- pca_M_final %>% select(all_of(as.numeric(meth_indice.pos[meth_pc_index][[1]])))
  
  Total.PC.Count <- sum(ncol(U_exp), ncol(U_meth), na.rm = TRUE)
  
  # --- Build regression dataframe ---
  df <- cbind(
    C      = t(C[, -1]),
    E      = t(E[, -1]),
    M     = t(M[, -1]),
    pos_age,
    pos_race,
    U_exp,
    U_meth
  )
  
  df <- as.data.frame(df)
  colnames(df)[1:3] <- c("C", "E", "M")
  df$pos_race <- as.factor(df$pos_race)
  
  return(df)
}

# ============================================================
# Step 4: Find a specific trio's model: Use MRGN
# ============================================================
# --- a trio example ----
df_trio_144897 <- get_trio_data(144897, trios, pos_meth_filter, pos_cna_filter,
                                pos_exp_filter, pca_E_final, pca_M_final,
                                gene_indice.pos, meth_indice.pos,gene.table.pos,meth.table.pos,
                                pos_age, pos_race)

# some error running MGRN, need Define isTrue as it's missing from your environment
assign("isTrue", isTRUE, envir = .GlobalEnv)

# --- Use MRGN to find the trio belong to which model ----
### Find the model of the trio

model_trio_144897_with_U <- infer.trio(df_trio_144897,is.CNA = TRUE,compute.nominal=TRUE)

model_trio_144897_without_U  <- infer.trio(df_trio_144897[1:3],is.CNA = TRUE,compute.nominal=TRUE)

### with Confounder, Inferred to M1.2
print(model_trio_144897_with_U)

### without Confounder, Inferred to M4
print(model_trio_144897_without_U)







