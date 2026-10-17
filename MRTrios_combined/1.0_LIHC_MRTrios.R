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

setwd("/wsu/home/hb/hb68/hb6890/fulab/GDCdata/TCGA-LIHC/raw_data")

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
LIHC.meth <- readRDS("split.names.LIHC.meth.logit.rds")                                  
LIHC.meth <- prep_dataset(LIHC.meth, id_col = "meth.row", sample_col_start = 5)

## ---- Gene expression --------------------------------------------------------
LIHC.gene <- fread("lihc_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt")
LIHC.gene <- prep_dataset(LIHC.gene, id_col = "gene.row", sample_col_start = 3)

## ---- CNA ---------------------------------------------------------------
LIHC.cna <- fread("lihc_tcga_pan_can_atlas_2018/data_cna.txt")
LIHC.cna <- prep_dataset(LIHC.cna, id_col = "cna.row", sample_col_start = 3)

## ---- Trios ---------------------------------------------------------------
trios <- as.data.frame(fread("trio.final.final.protein.coding.txt"))

## ---- Clinical ---------------------------------------------------------
clinical.LIHC <- fread("split_new_data_clinical_patient.txt")
clinical <- clinical.LIHC %>% select(1, 5, 6, 26)
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


cna         <- LIHC.cna
gene.exp    <- LIHC.gene
TCGA.meth   <- LIHC.meth

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
# Step 3b: Infer models using preprocessed data
# ============================================================
# Confounders (age, sex, race, gene PCs, meth PCs) are OPTIONAL per
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
    if (!all(is.na(sex))) {
      df_trio$sex <- as.factor(sex)
      used_confounders <- c(used_confounders, "sex")
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

    # ── Guard against confounder factors collapsing to 1 level ──
    # infer.trio() does its own internal NA/row filtering before
    # fitting the model. If that filtering happens to leave only one
    # category of a factor confounder (sex or race) among the
    # retained samples for THIS trio's C/E/M data, R's model-fitting
    # code throws "contrasts can be applied only to factors with 2 or
    # more levels". This is a per-trio data issue, not a bug — but it
    # must not be allowed to propagate, because a single such trio
    # would otherwise abort the whole do.call(rbind, lapply(...))
    # batch and take every other (successful) trio in the batch down
    # with it. (This is what silently dropped the posER_BRCA
    # part10/shard0057 output — trio 294363's race confounder
    # collapsed to 1 level and the resulting error killed the entire
    # 500-trio batch.)
    for (fac in c("sex", "race")) {
      if (fac %in% used_confounders && length(unique(na.omit(df_trio[[fac]]))) < 2) {
        cat(sprintf("  [WARN] trio %d: %s has < 2 levels after filtering; dropping %s from confounders.\n",
                    i, fac, fac))
        df_trio[[fac]] <- NULL
        used_confounders <- setdiff(used_confounders, fac)
      }
    }

    if (length(used_confounders) == 0) {
      cat(sprintf("  [INFO] trio %d: no confounders available; running infer.trio() on C/E/M only.\n", i))
    } else {
      cat(sprintf("  [INFO] trio %d: running infer.trio() with confounders: %s\n",
                  i, paste(used_confounders, collapse = ", ")))
    }

    # ── Per-trio error isolation ─────────────────────────────
    # infer.trio() can still fail for reasons not caught above (e.g.
    # rank-deficient design after MRGN's own internal filtering, or
    # other data idiosyncrasies for this specific trio). Catch the
    # error HERE, per trio, and return NULL for just this one row —
    # NOT at the batch level — so one bad trio never nukes the other
    # ~499 successful trios in the same shard/batch.
    model <- tryCatch(
      as.data.frame(infer.trio(df_trio, is.CNA = TRUE, compute.nominal = FALSE,
                                use.perm = TRUE)),
      error = function(e) {
        cat(sprintf("  [WARN] trio %d: infer.trio() failed (%s); skipping trio.\n", i, e$message))
        NULL
      }
    )
    if (is.null(model)) return(NULL)

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
output_dir    <- "/wsu/home/hb/hb68/hb6890/fulab/MRTrios/Output_LIHC_part"
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
                        sprintf("trio_LIHC_results_part%02d_shard_%04d.txt", part_id, b))
  
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
