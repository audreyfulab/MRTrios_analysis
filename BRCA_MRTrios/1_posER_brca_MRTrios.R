
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
setwd("/wsu/home/hb/hb68/hb6890/fulab/MRTrios/raw_Data_Methyl")

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
# Step 3b: Infer models using preprocessed data
# ============================================================
# Confounders (age, race, gene PCs, meth PCs) are now OPTIONAL per
# trio: instead of skipping a trio when a confounder is missing/NA,
# we drop just that confounder and still run infer.trio() on
# whatever is left (at minimum C, E, M).

get_trio_models <- function(trio_row, trios, preprocessed) {
  
  assign("isTrue", isTRUE, envir = .GlobalEnv)
  list2env(preprocessed, envir = environment())
  
  do.call(rbind, lapply(trio_row, function(i) {
    
    trio     <- trios[i, ]
    meth_row <- trio$meth.row
    cna_row  <- trio$cna.row
    gene_row <- trio$gene.row
    
    M <- meth_filter %>% filter(meth.row == meth_row)
    C <- cna_filter  %>% filter(cna.row  == cna_row)
    E <- exp_filter  %>% filter(gene.row == gene_row)
    
    # ── Required trio identifiers: cna_row/gene_row/meth_row must
    # each resolve to exactly one row of real C/E/M data. This is
    # NOT a confounder and cannot be made optional — if any of these
    # is NA or fails to match (e.g. NA == NA never matches in R),
    # there is no trio data to run infer.trio() on, so skip it.
    if (nrow(C) != 1 || nrow(E) != 1 || nrow(M) != 1) {
      cat(sprintf("  [WARN] trio %d: C/E/M did not match exactly one row each (cna=%d, exp=%d, meth=%d); skipping trio.\n",
                  i, nrow(C), nrow(E), nrow(M)))
      return(NULL)
    }
    
    # ── gene PC confounder (U_exp) ──────────────────────────
    U_exp <- NULL
    gene_pc_index <- gene.table %>% filter(col1 == gene_row) %>% pull(col2)
    
    if (length(gene_pc_index) == 0) {
      cat(sprintf("  [INFO] trio %d: gene_row '%s' not found in gene.table; running without gene PCs.\n",
                  i, gene_row))
    } else {
      gene_indices <- gene_indice[gene_pc_index][[1]]
      
      if (all(is.na(gene_indices))) {
        cat(sprintf("  [INFO] trio %d: gene_indice[%d] is NA/empty; running without gene PCs.\n",
                    i, gene_pc_index))
      } else {
        gene_indices_num <- suppressWarnings(as.numeric(gene_indices))
        gene_indices_num <- gene_indices_num[!is.na(gene_indices_num)]
        
        if (length(gene_indices_num) == 0) {
          cat(sprintf("  [INFO] trio %d: gene_indice[%d] has no valid numeric index; running without gene PCs.\n",
                      i, gene_pc_index))
        } else if (any(gene_indices_num > ncol(pca_E_final))) {
          cat(sprintf("  [WARN] trio %d: gene index %s out of range (max col = %d); running without gene PCs.\n",
                      i, paste(gene_indices_num, collapse=","), ncol(pca_E_final)))
        } else {
          U_exp <- pca_E_final %>% select(all_of(gene_indices_num))
        }
      }
    }
    
    # ── meth PC confounder (U_meth) ─────────────────────────
    U_meth <- NULL
    meth_pc_index <- meth.table %>% filter(col1 == meth_row) %>% pull(col2)
    
    if (length(meth_pc_index) == 0) {
      cat(sprintf("  [INFO] trio %d: meth_row '%s' not found in meth.table; running without meth PCs.\n",
                  i, meth_row))
    } else {
      meth_indices <- meth_indice[meth_pc_index][[1]]
      
      if (all(is.na(meth_indices))) {
        cat(sprintf("  [INFO] trio %d: meth_indice[%d] is NA/empty; running without meth PCs.\n",
                    i, meth_pc_index))
      } else {
        meth_indices_num <- suppressWarnings(as.numeric(meth_indices))
        meth_indices_num <- meth_indices_num[!is.na(meth_indices_num)]
        
        if (length(meth_indices_num) == 0) {
          cat(sprintf("  [INFO] trio %d: meth_indice[%d] has no valid numeric index; running without meth PCs.\n",
                      i, meth_pc_index))
        } else if (any(meth_indices_num > ncol(pca_M_final))) {
          cat(sprintf("  [WARN] trio %d: meth index %s out of range (max col = %d); running without meth PCs.\n",
                      i, paste(meth_indices_num, collapse=","), ncol(pca_M_final)))
        } else {
          U_meth <- pca_M_final %>% select(all_of(meth_indices_num))
        }
      }
    }
    
    Total.PC.Count <- sum(ncol(U_exp), ncol(U_meth), na.rm = TRUE)
    
    # ── Build the confounder set, dropping any that are NA/missing ──
    df_trio <- as.data.frame(cbind(
      C = t(C[, -1]),
      E = t(E[, -1]),
      M = t(M[, -1])
    ))
    colnames(df_trio)[1:3] <- c("C", "E", "M")
    
    used_confounders <- character(0)
    
    if (!all(is.na(age))) {
      df_trio$age <- age
      used_confounders <- c(used_confounders, "age")
    }
    if (!all(is.na(race))) {
      df_trio$race <- as.factor(race)
      used_confounders <- c(used_confounders, "race")
    }
    if (!is.null(U_exp)) {
      df_trio <- cbind(df_trio, U_exp)
      used_confounders <- c(used_confounders, "U_exp")
    }
    if (!is.null(U_meth)) {
      df_trio <- cbind(df_trio, U_meth)
      used_confounders <- c(used_confounders, "U_meth")
    }
    
    if (length(used_confounders) == 0) {
      cat(sprintf("  [INFO] trio %d: no confounders available; running infer.trio() on C/E/M only.\n", i))
    } else {
      cat(sprintf("  [INFO] trio %d: running infer.trio() with confounders: %s\n",
                  i, paste(used_confounders, collapse = ", ")))
    }
    
    model <- as.data.frame(infer.trio(df_trio, is.CNA = TRUE, compute.nominal = FALSE, use.perm = TRUE))
    
    cbind(
      data.frame(
        trio_row       = i,
        meth.row       = meth_row,
        cna.row        = cna_row,
        gene.row       = gene_row,
        Total.PC.Count = Total.PC.Count,
        Confounders    = paste(used_confounders, collapse = "|")
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
output_dir    <- "/wsu/home/hb/hb68/hb6890/fulab/MRTrios/Output_posER_BRCA_part"
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
                        sprintf("trio_posER_BRCA_results_part%02d_shard_%04d.txt", part_id, b))
  
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
