# ============================================================
# Step 1: Load required libraries
# ============================================================
library(data.table)
library(dplyr)
library(tidyverse)
library(MRTrios)

# ============================================================
# Step 2: Load data Files
# ============================================================
setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/raw_Data_Methyl")

humanmeth         <- read.csv("GPL13534_HumanMethylation450_15017482_v.1.1 2.csv", skip = 7, header = TRUE)
biomart           <- read.delim("ensembl37_genes_p13_biomart.txt", header = TRUE)

setwd("/wsu/home/hb/hb68/hb6890/fulab/GDCdata/TCGA-BLCA/raw_data")
setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/GDCdata/TCGA-BLCA/raw_data")

## Truncate a TCGA barcode down to the patient-level ID (first 3 "-" fields)
truncate_barcode <- function(x) {
  vapply(strsplit(x, "-"), function(parts) paste(parts[1:3], collapse = "-"), character(1))
}

prep_dataset <- function(df, id_col, sample_col_start) {
  df <- as.data.frame(df)
  sample_idx <- sample_col_start:ncol(df)
  colnames(df)[sample_idx] <- truncate_barcode(colnames(df)[sample_idx])
  
  n_dup <- sum(duplicated(names(df)))
  if (n_dup > 0) {
    message(sprintf("%s: %d duplicate sample columns after barcode truncation -- disambiguating with make.unique()",
                    id_col, n_dup))
    names(df) <- make.unique(names(df))
  }
  
  df[[id_col]] <- rownames(df)
  df
}

## ---- Methylation -----------------------------------------------------------
BLCA.meth <- readRDS("split.names.BLCA.meth.logit.rds")                                  # [421414, 444]
BLCA.meth <- prep_dataset(BLCA.meth, id_col = "meth.row", sample_col_start = 5)

## ---- Gene expression --------------------------------------------------------
BLCA.gene <- fread("blca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt")
BLCA.gene <- prep_dataset(BLCA.gene, id_col = "gene.row", sample_col_start = 3)

## ---- CNA ---------------------------------------------------------------
BLCA.cna <- fread("blca_tcga_pan_can_atlas_2018/data_cna.txt")
BLCA.cna <- prep_dataset(BLCA.cna, id_col = "cna.row", sample_col_start = 3)

## ---- Trios ---------------------------------------------------------------
trios <- as.data.frame(fread("trio.final.protein.coding.txt"))

## ---- Clinical ---------------------------------------------------------
clinical.BLCA <- fread("split_new_data_clinical_patient.txt")
clinical <- clinical.BLCA %>% select(1, 5, 6, 26)
colnames(clinical) <- c("person_id", "age", "sex", "race")

cna         <- BLCA.cna
gene.exp    <- BLCA.gene
TCGA.meth   <- BLCA.meth
#trios


## load new inferred FDR models
Model_BLCA <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BLCA/Output_BLCA/BLCA_trio_Model_results_ALL_with_BH_fdr_qval_byLZ.txt", sep = "\t")

# ============================================================
# Step 3a: Preprocess once
# ============================================================
BLCA.ind <- clinical[,1] 

# Use the defined function to get the M0 and M1 models
# add Inferred.Model3=Inferred.Model.BH_fdr, use for extractHumanProbeInfo function:

Model_BLCA$Inferred.Model3=Model_BLCA$Inferred.Model.BH_fdr
Model_BLCA$Inferred.Model3=Model_BLCA$Inferred.Model.BH_fdr

# ============================================================
# Step 3b: extractHumanMethProbeInfo()
# ============================================================
M0.1_BLCA = extractHumanMethProbeInfo(Model_BLCA, "M0.1", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, BLCA.ind)
M0.2_BLCA = extractHumanMethProbeInfo(Model_BLCA, "M0.2", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, BLCA.ind)

M1.1_BLCA = extractHumanMethProbeInfo(Model_BLCA, "M1.1", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, BLCA.ind)
M1.2_BLCA = extractHumanMethProbeInfo(Model_BLCA, "M1.2", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, BLCA.ind)


# ============================================================
# Step 4: Save data
# ============================================================
setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BLCA/Output_BLCA/")
write.table(M0.1_BLCA,"M0.1_BLCA_extract.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(M0.2_BLCA,"M0.2_BLCA_extract.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(M1.1_BLCA,"M1.1_BLCA_extract.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(M1.2_BLCA,"M1.2_BLCA_extract.txt",sep = "\t",quote = F,row.names = F,col.names = T)



         
