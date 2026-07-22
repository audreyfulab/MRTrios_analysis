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

setwd("/wsu/home/hb/hb68/hb6890/fulab/GDCdata/TCGA-BLCA/raw_data")

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
  sex <- clinical_filter$sex
  race <- clinical_filter$race
  
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
# Step 3b: Infer models using preprocessed data
# ============================================================
 get_trio_models <- function(trio_row, trios, preprocessed) {
  
  assign("isTrue", isTRUE, envir = .GlobalEnv)
  list2env(preprocessed, envir = environment())
  
  do.call(rbind, lapply(trio_row, function(i) {
    
    trio     <- trios[i, ]
    meth_row <- trio$meth.row
    cna_row  <- trio$cna.row
    gene_row <- trio$gene.row
    
    # ── Check NA in row identifiers first ─────────────────
    if (is.na(gene_row)) {
      cat(sprintf("  [WARN] trio %d: gene_row is NA, skipping.\n", i))
      return(NULL)
    }
    if (is.na(cna_row)) {
      cat(sprintf("  [WARN] trio %d: cna_row is NA, skipping.\n", i))
      return(NULL)
    }
    if (is.na(meth_row)) {
      cat(sprintf("  [WARN] trio %d: meth_row is NA, skipping.\n", i))
      return(NULL)
    }
    
    # ── Extract relevant rows from each dataset ────────────
    M <- meth_filter %>% filter(meth.row == meth_row)
    C <- cna_filter  %>% filter(cna.row  == cna_row)
    E <- exp_filter  %>% filter(gene.row == gene_row)
    
    # ── Get gene PCA confounders (optional) ───────────────
    gene_pc_index <- gene.table %>% filter(col1 == gene_row) %>% pull(col2)
    
    if (length(gene_pc_index) == 0) {
      # No significant PCs → run without gene confounders
      cat(sprintf("  [INFO] trio %d: gene_row '%s' has no significant PCs.\n", i, gene_row))
      U_exp <- NULL
    } else {
      gene_indices <- gene_indice[gene_pc_index][[1]]
      if (is.null(gene_indices) || all(is.na(gene_indices))) {
        # NULL/NA → run without gene confounders
        cat(sprintf("  [INFO] trio %d: gene_indice[%d] is NULL/NA, no gene confounders.\n",
                    i, gene_pc_index))
        U_exp <- NULL
      } else {
        gene_indices_num <- as.numeric(gene_indices)
        if (any(gene_indices_num > ncol(pca_E_final))) {
          cat(sprintf("  [WARN] trio %d: gene index out of range (max col = %d), skipping.\n",
                      i, ncol(pca_E_final)))
          return(NULL)  # ← real data error
        }
        U_exp <- pca_E_final %>% select(all_of(gene_indices_num))
      }
    }
    
    # ── Get meth PCA confounders (optional) ───────────────
    meth_pc_index <- meth.table %>% filter(col1 == meth_row) %>% pull(col2)
    
    if (length(meth_pc_index) == 0) {
      # No significant PCs → run without meth confounders
      cat(sprintf("  [INFO] trio %d: meth_row '%s' has no significant PCs.\n", i, meth_row))
      U_meth <- NULL
    } else {
      meth_indices <- meth_indice[meth_pc_index][[1]]
      if (is.null(meth_indices) || all(is.na(meth_indices))) {
        # NULL/NA → run without meth confounders
        cat(sprintf("  [INFO] trio %d: meth_indice[%d] is NULL/NA, no meth confounders.\n",
                    i, meth_pc_index))
        U_meth <- NULL
      } else {
        meth_indices_num <- as.numeric(meth_indices)
        if (any(meth_indices_num > ncol(pca_M_final))) {
          cat(sprintf("  [WARN] trio %d: meth index out of range (max col = %d), skipping.\n",
                      i, ncol(pca_M_final)))
          return(NULL)  # ← real data error
        }
        U_meth <- pca_M_final %>% select(all_of(meth_indices_num))
      }
    }
    
    # ── Count total PCs used ───────────────────────────────
    Total.PC.Count <- sum(
      if (!is.null(U_exp))  ncol(U_exp)  else 0,
      if (!is.null(U_meth)) ncol(U_meth) else 0
    )
    
    # ── Build regression dataframe ─────────────────────────
    base_df <- data.frame(
      C    = as.numeric(t(C[, -1])),
      E    = as.numeric(t(E[, -1])),
      M    = as.numeric(t(M[, -1])),
      age  = age,
      sex  = sex,
      race = race
    )
    if (!is.null(U_exp))  base_df <- cbind(base_df, U_exp)
    if (!is.null(U_meth)) base_df <- cbind(base_df, U_meth)
    
    df_trio       <- as.data.frame(base_df)
    df_trio$race  <- as.factor(df_trio$race)
    df_trio$sex  <- as.factor(df_trio$sex)
    
    # ── Infer causal model ─────────────────────────────────
    model <- as.data.frame(infer.trio(df_trio, is.CNA = TRUE))
    
    # ── Return one row of results ──────────────────────────
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

# ============================================================
# Step 5: Compute this part's trio range
# ============================================================
total_trios  <- nrow(trios)
already_done <- 0                         # ← add this
part_size    <- ceiling((total_trios - already_done) / n_parts)
part_start   <- already_done + (part_id - 1) * part_size + 1
part_end     <- min(already_done + part_id * part_size, total_trios)

cat(sprintf("Part %d: trios %d to %d (%d trios)\n",
            part_id, part_start, part_end, part_end - part_start + 1))

# ============================================================
# Step 6: Run sharded batches within this part
# ============================================================
output_dir    <- "/wsu/home/hb/hb68/hb6890/fulab/MRTrios/Output_BLCA_part"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

batch_size    <- 500
batch_starts  <- seq(part_start, part_end, by = batch_size)
total_batches <- length(batch_starts)

cat(sprintf("Batch size: %d | Total batches in this part: %d\n",
            batch_size, total_batches))

for (b in seq_along(batch_starts)) {
  
  start    <- batch_starts[b]
  end      <- min(start + batch_size - 1, part_end)
  out_file <- file.path(output_dir,
                        sprintf("trio_results_part%02d_shard_%04d.txt", part_id, b))
  
  if (file.exists(out_file)) {
    cat(sprintf("[Part %02d | Batch %d / %d] Already exists, skipping.\n",
                part_id, b, total_batches))
    next
  }
  
  cat(sprintf("[Part %02d | Batch %d / %d] Trios %d to %d ...\n",
              part_id, b, total_batches, start, end))
  t_start <- proc.time()
  
  batch_result <- tryCatch(
    get_trio_models(start:end, trios, preprocessed),
    error = function(e) {
      cat(sprintf("  -> ERROR in batch %d (rows %d-%d): %s\n",
                  b, start, end, e$message))
      NULL
    }
  )
  
  # ── NEW - replace with this ────────────────────────────────
  if (!is.null(batch_result)) {
    fwrite(batch_result, file = out_file, sep = "\t")
    elapsed <- round((proc.time() - t_start)["elapsed"], 1)
    cat(sprintf("  -> Saved: %s (%d rows) [%.1f sec]\n",
                basename(out_file), nrow(batch_result), elapsed))
  } else {
    cat(sprintf("  -> SKIPPED: batch %d produced no results (all trios warned/errored)\n", b))
  }
  
} # end for loop


cat(sprintf("\nPart %d complete. Shards saved to: %s\n", part_id, output_dir))