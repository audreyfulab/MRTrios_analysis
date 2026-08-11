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

cna               <- fread('data_CNA.txt',data.table = F)
gene.exp          <- fread('data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt',data.table = F) 
TCGA.meth         <- readRDS("split.names.TCGA.meth.logit.rds") 
trios             <- fread("trio.final.protein.coding.txt")
humanmeth         <- read.csv("GPL13534_HumanMethylation450_15017482_v.1.1 2.csv", skip = 7, header = TRUE)
biomart           <- read.delim("ensembl37_genes_p13_biomart.txt", header = TRUE)
clinical.neg      <- fread("names.neg.patient2.txt", header = FALSE)
clinical.pos      <- fread("names.pos.patient2.txt", header = FALSE)
# --- load new inferred FDR models ------
Model_posER_BRCA  <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_posER_BRCA/trio_posER_BRCA_results_ALL.txt", sep = "\t")
Model_posER_BRCA  <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_posER_BRCA/trio_posER_BRCA_results_ALL.txt", sep = "\t")

# ============================================================
# Step 3a: Preprocess once
# ============================================================
neg.ind <- clinical.neg[,1]
pos.ind <- clinical.pos[,1]

# Use the defined function to get the M0 and M1 models
# add Inferred.Model3=Inferred.Model.BH_fdr, use for extractHumanProbeInfo function:

Model_posER_BRCA$Inferred.Model3=Model_posER_BRCA$Inferred.Model.BH_fdr
Model_negER_BRCA$Inferred.Model3=Model_negER_BRCA$Inferred.Model.BH_fdr

# ============================================================
# Step 3b: extractHumanMethProbeInfo()
# ============================================================
M0.1_pos = extractHumanMethProbeInfo(Model_posER_BRCA, "M0.1", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, pos.ind)
M0.2_pos = extractHumanMethProbeInfo(Model_posER_BRCA, "M0.2", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, pos.ind)

M1.1_pos = extractHumanMethProbeInfo(Model_posER_BRCA, "M1.1", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, pos.ind)
M1.2_pos = extractHumanMethProbeInfo(Model_posER_BRCA, "M1.2", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, pos.ind)

M0.1_neg = extractHumanMethProbeInfo(Model_negER_BRCA, "M0.1", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, neg.ind)
M0.2_neg = extractHumanMethProbeInfo(Model_negER_BRCA, "M0.2", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, neg.ind)

M1.1_neg = extractHumanMethProbeInfo(Model_negER_BRCA, "M1.1", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, neg.ind)
M1.2_neg = extractHumanMethProbeInfo(Model_negER_BRCA, "M1.2", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, neg.ind)

# ============================================================
# Step 4: Save data
# ============================================================
setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_posER_BRCA/")
write.table(M0.1_pos,"M0.1_pos_extract.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(M0.2_pos,"M0.2_pos_extract.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(M1.1_pos,"M1.1_pos_extract.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(M1.2_pos,"M1.2_pos_extract.txt",sep = "\t",quote = F,row.names = F,col.names = T)

setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_negER_BRCA/")
write.table(M0.1_neg,"M0.1_neg_extract.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(M0.2_neg,"M0.2_neg_extract.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(M1.1_neg,"M1.1_neg_extract.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(M1.2_neg,"M1.2_neg_extract.txt",sep = "\t",quote = F,row.names = F,col.names = T)

