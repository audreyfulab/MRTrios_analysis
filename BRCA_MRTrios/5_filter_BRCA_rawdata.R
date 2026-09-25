# ============================================================
# Step 0: Read command line arguments
# ============================================================
args    <- commandArgs(trailingOnly = TRUE)
part_id <- as.integer(args[1])  # which part (1 to 10)
n_parts <- as.integer(args[2])  # total parts (10)

cat(sprintf("Running part %d of %d\n", part_id, n_parts))

# ============================================================
# Step 1: Load required libraries
# ============================================================
.libPaths(c('/wsu/home/hb/hb68/hb6890/fulab/MRGNgeneral',
            '/wsu/home/hb/hb68/hb6890/fulab/MRGN',
            .libPaths()))
library(data.table)
library(tidyverse)
library(MRGN)


# ============================================================
# Step 2: Load data Files
# ============================================================
setwd("/wsu/home/hb/hb68/hb6890/fulab/MRTrios/raw_Data_Methyl") ## On Grid
setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/raw_Data_Methyl") # local

cna         <- fread('data_CNA.txt', data.table = F)
gene.exp    <- fread('data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt', data.table = F)
TCGA.meth   <- fread("split.names.TCGA.meth.logit.txt", data.table = F)
trios       <- fread("trio.final.protein.coding.txt")
clinical    <- fread("names.pos.patient2.txt", header = FALSE)
pc.gene     <- fread("PCA.gene.exp.posER.txt")
pc.meth     <- fread("PCA.meth.posER.txt")
gene.table  <- fread("gene.exp.posER.table.txt")
meth.table  <- fread("meth.posER.table.txt")
gene_indice <- readRDS("gene_exp_posER_sig_pcs_indices.rds")
meth_indice <- readRDS("meth_posER_sig_pcs_indices.rds")

# ── Normalize empty/NULL PC-index list entries to NA ────────
# gene_indice / meth_indice are lists keyed by PC-table row. An
# entry can come back as NULL, character(0)/integer(0), or "" —
# collapse all of those forms to a single NA so every downstream
# check only has to test for NA.
normalize_indice <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA)
  x[x == ""] <- NA
  x
}
gene_indice <- lapply(gene_indice, normalize_indice)
meth_indice <- lapply(meth_indice, normalize_indice)


# ============================================================
# Step 3a: Preprocess once
# ============================================================
preprocess_trio_data <- function(trios, cna, gene.exp, TCGA.meth, clinical,
                                 gene_indice, meth_indice, pc.gene, pc.meth,
                                 gene.table, meth.table) {
  
  com.ind <- intersect(
    intersect(colnames(gene.exp)[3:ncol(gene.exp)],
              colnames(TCGA.meth)[5:ncol(TCGA.meth)]),
    colnames(cna)[3:ncol(cna)])
  
  unique_id <- intersect(com.ind, clinical$V1)
  
  clinical_filter <- clinical %>%
    arrange(match(V1, unique_id)) %>%
    filter(V1 %in% unique_id)
  
  age  <- clinical_filter$V2
  race <- clinical_filter$V3
  
  # ── Treat empty/blank/placeholder values as NA ─────────────
  # fread() already turns truly empty numeric cells into NA, but
  # values stored as blank strings ("", " ") or literal "NA"/"NULL"
  # text survive as character data. Normalize all of these to NA
  # before they're used as confounders.
  age  <- trimws(as.character(age))
  race <- trimws(as.character(race))
  age[age   %in% c("", "NA", "NULL", "N/A")]  <- NA
  race[race %in% c("", "NA", "NULL", "N/A")]  <- NA
  age <- suppressWarnings(as.numeric(age))
  
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
  
  pca_E_final <- pc.gene %>%
    filter(V1 %in% unique_id) %>%
    arrange(match(V1, unique_id)) %>%
    select(-V1)
  
  pca_M_final <- pc.meth %>%
    filter(V1 %in% unique_id) %>%
    arrange(match(V1, unique_id)) %>%
    select(-V1)
    
  list(
    unique_id   = unique_id,
    age         = age,
    race        = race,
    cna_filter  = cna_filter,
    exp_filter  = exp_filter,
    meth_filter = meth_filter,
    pca_E_final = pca_E_final,
    pca_M_final = pca_M_final,
    gene.table  = gene.table,
    meth.table  = meth.table,
    gene_indice = gene_indice,
    meth_indice = meth_indice
  )
}

# ============================================================
# Step 3b: Preprocess once
# ============================================================
cat("Preprocessing data...\n")
preprocessed <- preprocess_trio_data(
  trios, cna, gene.exp, TCGA.meth, clinical,
  gene_indice, meth_indice, pc.gene, pc.meth,
  gene.table, meth.table
)
cat("Preprocessing done. Total trios:", nrow(trios), "\n")

# ============================================================
# Step 4:  Save Preprocess Data ----- posER_BRCA
# ============================================================

setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/process_data_BRCA")
saveRDS(preprocessed,"posER_BRCA_data_filter_list.rds")

posER_BRCA_filter <- readRDS("posER_BRCA_data_filter_list.rds") #load it


## negER

# ============================================================
# Step 2: Load data Files
# ============================================================
setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/raw_Data_Methyl")

#cna         <- fread('data_CNA.txt', data.table = F)
#gene.exp    <- fread('data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt', data.table = F)
#TCGA.meth   <- fread("split.names.TCGA.meth.logit.txt", data.table = F)
#trios       <- fread("trio.final.protein.coding.txt")
clinical    <- fread("names.neg.patient2.txt", header = FALSE)
pc.gene     <- fread("PCA.gene.exp.negER.txt")
pc.meth     <- fread("PCA.meth.negER.txt")
gene.table  <- fread("gene.exp.negER.table.txt")
meth.table  <- fread("meth.negER.table.txt")
gene_indice <- readRDS("gene_exp_negER_sig_pcs_indices.rds")
meth_indice <- readRDS("meth_negER_sig_pcs_indices.rds")

gene_indice <- lapply(gene_indice, normalize_indice)
meth_indice <- lapply(meth_indice, normalize_indice)

# ============================================================
# Step 4:  Save Preprocess Data ----- negER_BRCA
# ============================================================
cat("Preprocessing data...\n")
preprocessed <- preprocess_trio_data(
  trios, cna, gene.exp, TCGA.meth, clinical,
  gene_indice, meth_indice, pc.gene, pc.meth,
  gene.table, meth.table
)
cat("Preprocessing done. Total trios:", nrow(trios), "\n")

setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/process_data_BRCA")
saveRDS(preprocessed,"negER_BRCA_data_filter_list.rds")

negER_BRCA_filter <- readRDS("negER_BRCA_data_filter_list.rds") #load it

