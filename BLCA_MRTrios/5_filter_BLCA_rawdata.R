# ============================================================
# Step 2: Load data Files
# ============================================================

setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/GDCdata/TCGA-BLCA/")

library(data.table)
library(dplyr)
library(tibble)

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

## ---- PCA score matrices ------------------------------------------------------
pc.meth <- read.table("PCA.meth.txt", header = TRUE)
rownames(pc.meth) <- truncate_barcode(rownames(pc.meth))

pc.gene <- read.table("PCA.gene.exp.txt", header = TRUE)
rownames(pc.gene) <- truncate_barcode(rownames(pc.gene))

## ---- Indices ---------------------------------------------------------------
meth.table   <- fread("meth.table.txt", drop = 1)
gene.table   <- fread("gene.exp.table.txt", drop = 1)
meth_indice  <- readRDS("meth_sig_pcs_indices.rds")
gene_indice  <- readRDS("gene_exp_sig_pcs_indices.rds")

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


cna         <- BLCA.cna
gene.exp    <- BLCA.gene
TCGA.meth   <- BLCA.meth

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
  
  unique_id <- intersect(com.ind, clinical$person_id)
  
  clinical_filter <- clinical %>%
    arrange(match(person_id, unique_id)) %>%
    filter(person_id %in% unique_id)
  
  age  <- clinical_filter$age
  sex  <- clinical_filter$sex
  race <- clinical_filter$race
  
  # ── Treat empty/blank/placeholder values as NA ─────────────
  # fread() already turns truly empty numeric cells into NA, but
  # values stored as blank strings ("", " ") or literal "NA"/"NULL"
  # text survive as character data. Normalize all of these to NA
  # before they're used as confounders.
  age  <- trimws(as.character(age))
  sex  <- trimws(as.character(sex))
  race <- trimws(as.character(race))
  age[age   %in% c("", "NA", "NULL", "N/A")]  <- NA
  sex[sex   %in% c("", "NA", "NULL", "N/A")]  <- NA
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
  
  names(TCGA.meth) <- make.unique(names(TCGA.meth))
  TCGA.meth$meth.row <- rownames(TCGA.meth)
  
  meth_filter <- TCGA.meth %>%
    filter(meth.row %in% unique(trios$meth.row)) %>%
    select(meth.row, all_of(unique_id))
  
  pca_E_final <-  pc.gene %>% rownames_to_column(var = "person_id") %>%
    filter(person_id %in% unique_id) %>%
    arrange(match(person_id, unique_id)) %>%
    select(-person_id)
  
  pca_M_final <- pc.meth %>% rownames_to_column(var = "person_id") %>%
    filter(person_id %in% unique_id) %>%
    arrange(match(person_id, unique_id)) %>%
    select(-person_id)
  
  
  list(
    unique_id   = unique_id,
    age         = age,
    sex         = sex,
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
# Step 4: Preprocess once
# ============================================================
cat("Preprocessing data...\n")
preprocessed <- preprocess_trio_data(
  trios, cna, gene.exp, TCGA.meth, clinical,
  gene_indice, meth_indice, pc.gene, pc.meth,
  gene.table, meth.table
)
cat("Preprocessing done. Total trios:", nrow(trios), "\n")

setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BLCA/process_data_BLCA")
saveRDS(preprocessed,"BLCA_data_filter_list.rds")

BLCA_filter <- readRDS("BLCA_data_filter_list.rds") #load it
