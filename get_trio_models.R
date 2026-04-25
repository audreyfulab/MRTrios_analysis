#' Infer Causal Trio Models with Confounding Variables
#'
#' @description
#' For each trio (methylation probe, CNA, gene expression), this function:
#' \itemize{
#'   \item Aligns individuals across all datasets
#'   \item Extracts significant PCA covariates for methylation and expression
#'   \item Builds a regression-ready data frame
#'   \item Calls \code{infer.trio()} from the MRGN package to infer the causal model
#' }
#'
#' @param trio_row Integer or integer vector. Row index/indices into \code{trios}
#'   to process. E.g., \code{1}, \code{1:100}, or \code{c(1, 5, 20)}.
#' @param trios A data frame with at least three columns:
#'   \code{meth.row}, \code{cna.row}, \code{gene.row} — row indices
#'   into the methylation, CNA, and gene expression matrices respectively.
#' @param cna A data frame of copy number alteration data;
#'   genes in rows, individuals in columns.
#' @param gene.exp A data frame of gene expression data (Z-scores);
#'   genes in rows, individuals in columns.
#' @param TCGA.meth A data frame of DNA methylation data (logit-transformed);
#'   probes in rows, individuals in columns.
#' @param clinical.pos A data frame of clinical data for ER+ patients with columns:
#'   \code{V1} = patient ID, \code{V2} = age, \code{V3} = race.
#' @param gene_pos_indice A list of significant PC indices for gene expression,
#'   as returned by \code{findPCs()} or loaded from \code{.rds}.
#' @param meth_pos_indice A list of significant PC indices for methylation,
#'   as returned by \code{findPCs()} or loaded from \code{.rds}.
#' @param pc.gene.pos A data frame of PC scores for gene expression;
#'   first column is patient ID (\code{V1}), remaining columns are PC scores.
#' @param pc.meth.pos A data frame of PC scores for methylation;
#'   first column is patient ID (\code{V1}), remaining columns are PC scores.
#' @param gene.table.pos A data frame mapping gene expression row indices to
#'   their corresponding index in \code{gene_pos_indice}.
#'   Columns: \code{col1} = original row index, \code{col2} = list index.
#' @param meth.table.pos A data frame mapping methylation row indices to
#'   their corresponding index in \code{meth_pos_indice}.
#'   Columns: \code{col1} = original row index, \code{col2} = list index.
#'
#' @return A data frame with one row per trio and the following columns:
#' \describe{
#'   \item{trio_row}{Row index of the trio in \code{trios}}
#'   \item{meth.row}{Methylation probe row index}
#'   \item{cna.row}{CNA row index}
#'   \item{gene.row}{Gene expression row index}
#'   \item{Total.PC.Count}{Total number of PCs used as covariates}
#'   \item{b11}{Binary: 1 if conditional test T1 ~ V | \{T2, U\} rejects null}
#'   \item{b12}{Binary: outcome for conditional test T1 ~ T2 | \{V, U\}}
#'   \item{b21}{Binary: outcome for conditional test T2 ~ V | \{T1, U\}}
#'   \item{b22}{Binary: outcome for conditional test T2 ~ T1 | \{V, U\}}
#'   \item{pb11}{P-value for conditional test T1 ~ V | \{T2, U\}}
#'   \item{pb12}{P-value for conditional test T1 ~ T2 | \{V, U\}}
#'   \item{pb21}{P-value for conditional test T2 ~ V | \{T1, U\}}
#'   \item{pb22}{P-value for conditional test T2 ~ T1 | \{V, U\}}
#'   \item{Inferred.Model}{String indicating the causal model type}
#' }
#'
#' @details
#' The function performs the following steps internally:
#' \enumerate{
#'   \item Finds common individuals across CNA, gene expression, methylation,
#'         and clinical datasets.
#'   \item Filters and aligns each dataset to the common individuals.
#'   \item For each trio index \code{i} in \code{trio_row}:
#'     \itemize{
#'       \item Extracts the relevant rows from each filtered dataset.
#'       \item Retrieves significant PCs for methylation and expression.
#'       \item Builds a complete-case regression data frame including
#'             age, race, and PC covariates.
#'       \item Calls \code{\link[MRGN]{infer.trio}} to infer the causal model.
#'     }
#' }
#'
#' @seealso
#' \code{\link[MRGN]{infer.trio}} for the underlying causal inference method.
#'
#' @examples
#' \dontrun{
#' # Run on first 10 trios
#' results <- get_trio_models(
#'   trio_row       = 1:10,
#'   trios          = trios,
#'   cna            = cna,
#'   gene.exp       = gene.exp,
#'   TCGA.meth      = TCGA.meth,
#'   clinical.pos   = clinical.pos,
#'   gene_pos_indice = gene_pos_indice,
#'   meth_pos_indice = meth_pos_indice,
#'   pc.gene.pos    = pc.gene.pos,
#'   pc.meth.pos    = pc.meth.pos,
#'   gene.table.pos = gene.table.pos,
#'   meth.table.pos = meth.table.pos
#' )
#'
#' # View results
#' head(results)
#'
#' # Run on a single trio
#' result_one <- get_trio_models(trio_row = 5, ...)
#'
#' # Save results
#' saveRDS(results, "trio_models_results.rds")
#' }
#'
#' @importFrom dplyr filter select pull arrange match all_of mutate
#' @importFrom MRGN infer.trio
#' @export

get_trio_models <- function(trio_row,
                            trios,
                            cna,
                            gene.exp,
                            TCGA.meth,
                            clinical.pos,
                            gene_pos_indice,
                            meth_pos_indice,
                            pc.gene.pos,
                            pc.meth.pos,
                            gene.table.pos,
                            meth.table.pos) {

  assign("isTrue", isTRUE, envir = .GlobalEnv)

  # ── Step 1: Find common individuals across all 4 datasets ──────────────────

  com.ind <- intersect(
    intersect(colnames(gene.exp)[3:ncol(gene.exp)],
              colnames(TCGA.meth)[5:ncol(TCGA.meth)]),
    colnames(cna)[3:ncol(cna)])

  unique_id <- intersect(com.ind, clinical.pos$V1)

  # ── Step 2: Extract clinical covariates ────────────────────────────────────

  clinical.pos_filter <- clinical.pos %>%
    arrange(match(V1, unique_id)) %>%
    filter(V1 %in% unique_id)

  pos_age  <- clinical.pos_filter$V2
  pos_race <- clinical.pos_filter$V3

  # ── Step 3: Filter datasets to common individuals ──────────────────────────

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

  # ── Step 4: Align PCA score matrices to common individuals ─────────────────

  pca_E_final <- pc.gene.pos %>%
    filter(V1 %in% unique_id) %>%
    arrange(match(V1, unique_id)) %>%
    select(-V1)

  pca_M_final <- pc.meth.pos %>%
    filter(V1 %in% unique_id) %>%
    arrange(match(V1, unique_id)) %>%
    select(-V1)

  # ── Step 5: Loop over trio indices and infer models ────────────────────────

  do.call(rbind, lapply(trio_row, function(i) {

    trio     <- trios[i, ]
    meth_row <- trio$meth.row
    cna_row  <- trio$cna.row
    gene_row <- trio$gene.row

    # Extract relevant rows from each dataset
    M <- pos_meth_filter %>% filter(meth.row == meth_row)
    C <- pos_cna_filter  %>% filter(cna.row  == cna_row)
    E <- pos_exp_filter  %>% filter(gene.row == gene_row)

    # Retrieve significant PCs for expression
    gene_pc_index <- gene.table.pos %>% filter(col1 == gene_row) %>% pull(col2)
    U_exp <- pca_E_final %>% select(all_of(as.numeric(gene_pos_indice[gene_pc_index][[1]])))

    # Retrieve significant PCs for methylation
    meth_pc_index <- meth.table.pos %>% filter(col1 == meth_row) %>% pull(col2)
    U_meth <- pca_M_final %>% select(all_of(as.numeric(meth_pos_indice[meth_pc_index][[1]])))

    Total.PC.Count <- sum(ncol(U_exp), ncol(U_meth), na.rm = TRUE)

    # Build regression data frame
    df_trio <- as.data.frame(cbind(
      C = t(C[, -1]),
      E = t(E[, -1]),
      M = t(M[, -1]),
      pos_age,
      pos_race,
      U_exp,
      U_meth
    ))

    colnames(df_trio)[1:3] <- c("C", "E", "M")
    df_trio$pos_race <- as.factor(df_trio$pos_race)

    # Infer causal model
    model <- as.data.frame(infer.trio(df_trio))

    # Return one row of results
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
