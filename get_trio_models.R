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

cna <- fread('data_CNA.txt',data.table = F)
gene.exp <- fread('data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt',data.table = F)
TCGA.meth <- fread("split.names.TCGA.meth.logit.txt", data.table=F)
trios <- fread("trio.final.protein.coding.txt")

gene_pos_indice <- readRDS("gene_exp_posER_sig_pcs_indices.rds")
meth_pos_indice <- readRDS("meth_posER_sig_pcs_indices.rds")
pc.gene.pos <- fread("PCA.gene.exp.posER.txt")
pc.meth.pos <- fread("PCA.meth.posER.txt")

### meth,gene.exp table : # V1,col1, col2 , the col2 is the final pcscore index.
meth.table.pos <- fread("meth.posER.table.txt")
gene.table.pos <- fread("gene.exp.posER.table.txt")


clinical.pos  <- fread("names.pos.patient2.txt", header = FALSE) # ER+ patients ID

# length 787: meth, gene.exp
# length 777: Finding common individuals among the 3 datasets:cna,gene.exp, meth
# length 572: meth, gene.exp,clinical.pos
# lenth  564: cna,meth, gene.exp,clinical.pos


# ============================================================
# Step 2: write function : get_trio_models()
# ============================================================
get_trio_models <- function(trio_row,
                            trios, cna, gene.exp, TCGA.meth,
                            clinical.pos,
                            gene_pos_indice,
                            meth_pos_indice,
                            pc.gene.pos,
                            pc.meth.pos,
                            gene.table.pos,
                            meth.table.pos) {
  
  assign("isTrue", isTRUE, envir = .GlobalEnv)
  
  # ============================================================
  # Step 2.1: process data ONCE outside the loop
  # ============================================================
  
  # --- extract unique IDs ---
  com.ind <- intersect(
    intersect(colnames(gene.exp)[3:ncol(gene.exp)],
              colnames(TCGA.meth)[5:ncol(TCGA.meth)]),
    colnames(cna)[3:ncol(cna)])
  
  unique_id <- intersect(com.ind, clinical.pos$V1)
  
  # --- extract AGE, RACE as confounder ---
  clinical.pos_filter <- clinical.pos %>%
    arrange(match(V1, unique_id)) %>%
    filter(V1 %in% unique_id)
  
  pos_age  <- clinical.pos_filter$V2
  pos_race <- clinical.pos_filter$V3
  
  # --- Filter cna, gene.exp, meth ---
  pos_cna_filter <- cna %>%
    mutate(cna.row = rownames(cna)) %>%
    filter(cna.row %in% unique(trios$cna.row)) %>%
    select(cna.row, all_of(unique_id))
  
  pos_exp_filter <- gene.exp %>%
    mutate(gene.row = rownames(gene.exp)) %>%
    filter(gene.row %in% unique(trios$gene.row)) %>%
    select(gene.row, all_of(unique_id))
  
  pos_meth_filter <- TCGA.meth %>%
    mutate(meth.row = rownames(TCGA.meth)) %>%
    filter(meth.row %in% unique(trios$meth.row)) %>%
    select(meth.row, all_of(unique_id))
  
  # --- Filter and sort PCA matrices ---
  pca_E_sorted <- pc.gene.pos %>%
    filter(V1 %in% unique_id) %>%
    arrange(match(V1, unique_id))
  
  pca_M_sorted <- pc.meth.pos %>%
    filter(V1 %in% unique_id) %>%
    arrange(match(V1, unique_id))
  
  # --- Remove ID column ---
  pca_E_final <- pca_E_sorted[, -1]
  pca_M_final <- pca_M_sorted[, -1]
  
  # ============================================================
  # Step 2.2: loop over trio indices
  # ============================================================
  
  do.call(rbind, lapply(trio_row, function(i) {
    
    # BUG 1 FIX: use i, not trio_row
    trio     <- trios[i, ]
    meth_row <- trio$meth.row
    cna_row  <- trio$cna.row
    gene_row <- trio$gene.row
    
    # --- Filter data rows ---
    M <- pos_meth_filter %>% filter(meth.row == meth_row)
    C <- pos_cna_filter  %>% filter(cna.row  == cna_row)
    E <- pos_exp_filter  %>% filter(gene.row == gene_row)
    
    # --- PCA covariates for expression ---
    gene_pc_index <- gene.table.pos %>% filter(col1 == gene_row) %>% pull(col2)
    U_E   <- gene_pos_indice[gene_pc_index]
    U_exp <- pca_E_final %>% select(all_of(as.numeric(U_E[[1]])))
    
    # --- PCA covariates for methylation ---
    meth_pc_index <- meth.table.pos %>% filter(col1 == meth_row) %>% pull(col2)
    U_M    <- meth_pos_indice[meth_pc_index]
    U_meth <- pca_M_final %>% select(all_of(as.numeric(U_M[[1]])))
    
    # --- Build regression dataframe ---
    df_trio <- cbind(
      C = t(C[, -1]),
      E = t(E[, -1]),
      M = t(M[, -1]),
      pos_age,
      pos_race,
      U_exp,
      U_meth
    )
    
    # BUG 2 FIX: use df_trio consistently
    df_trio <- as.data.frame(df_trio)
    colnames(df_trio)[1:3] <- c("C", "E", "M")
    df_trio$pos_race <- as.factor(df_trio$pos_race)
    
    Total.PC.Count <- sum(ncol(U_meth), ncol(U_exp), na.rm = TRUE)
    
    # --- Infer trio model ---
    model <- as.data.frame(infer.trio(df_trio))
    
    # BUG 3 FIX: save cbind result and return it
    result_row <- cbind(
      data.frame(trio_row = i,
                 meth.row = meth_row,
                 cna.row  = cna_row,
                 gene.row = gene_row,
                 Total.PC.Count = Total.PC.Count),
      model
    )
    
    return(result_row)
  }))
}


# ==============================================
# Step 3: Model Application and Validation
# ==============================================
corM1.2_models_11 <- get_trio_models(cor_M1.2$pos_trios_row, 
                                     trios, cna, gene.exp, TCGA.meth,
                                     clinical.pos,
                                     gene_pos_indice,
                                     meth_pos_indice,
                                     pc.gene.pos,
                                     pc.meth.pos,
                                     gene.table.pos,
                                     meth.table.pos)

# Original case (Pos_M1.2 model)

# Any other set of trios
corM1.1_models <- get_trio_models(cor_M1.1$pos_trios_row, 
                                     trios, cna, gene.exp, TCGA.meth,
                                     clinical.pos,
                                     gene_pos_indice,
                                     meth_pos_indice,
                                     pc.gene.pos,
                                     pc.meth.pos,
                                     gene.table.pos,
                                     meth.table.pos)



### === check all the Pos_models of M1.2 ====

final_pos_M1.2_models <- get_trio_models(final_pos_M1.2$pos_trios_row,
                                     trios, cna, gene.exp, TCGA.meth,
                                     clinical.pos,
                                     gene_pos_indice,
                                     meth_pos_indice,
                                     pc.gene.pos,
                                     pc.meth.pos,
                                     gene.table.pos,
                                     meth.table.pos) 
