# =======================================================================

# ========================================================================
# Merge trio stats with MRGN model results for all 3 cancer type datasets:
#   posER BRCA, negER BRCA, BLCA, LIHC
# ========================================================================

library(data.table)

base_dir <- "/Users/lianzuo/LZ/ResearchProject/Fulab"

# ============================================================
# Step 1: File paths for each dataset
# ============================================================
datasets <- list(
  posER_BRCA = list(
    stats = file.path(base_dir, "MRTrios_BRCA/process_data_BRCA/posER_BRCA_trio_stats.txt"),
    model = file.path(base_dir, "MRTrios_BRCA/Output_posER_BRCA/posER_BRCA_trio_Model_results_ALL_with_BH_fdr_qval_byLZ.txt"),
    out   = file.path(base_dir, "MRTrios_BRCA/Output_posER_BRCA/posER_BRCA_trio_Model_results_with_stats.txt")
  ),
  negER_BRCA = list(
    stats = file.path(base_dir, "MRTrios_BRCA/process_data_BRCA/negER_BRCA_trio_stats.txt"),
    model = file.path(base_dir, "MRTrios_BRCA/Output_negER_BRCA/negER_BRCA_trio_Model_results_ALL_with_BH_fdr_qval_byLZ.txt"),
    out   = file.path(base_dir, "MRTrios_BRCA/Output_negER_BRCA/negER_BRCA_trio_Model_results_with_stats.txt")
  ),
  BLCA = list(
    stats = file.path(base_dir, "MRTrios_BLCA/process_data_BLCA/BLCA_trio_stats.txt"),
    model = file.path(base_dir, "MRTrios_BLCA/Output_BLCA/BLCA_trio_Model_results_ALL_with_BH_fdr_qval_byLZ.txt"),
    out   = file.path(base_dir, "MRTrios_BLCA/Output_BLCA/BLCA_trio_Model_results_with_stats.txt")
  ),
  LIHC = list(
    stats = file.path(base_dir, "MRTrios_LIHC/process_data_LIHC/LIHC_trio_stats.txt"),
    model = file.path(base_dir, "MRTrios_LIHC/Output_LIHC/LIHC_trio_Model_results_ALL_with_BH_fdr_qval_byLZ.txt"),
    out   = file.path(base_dir, "MRTrios_LIHC/Output_LIHC/LIHC_trio_Model_results_with_stats.txt")
  )
)

front_cols <- c("trio_row", "Gene.name", "meth.row", "cna.row", "gene.row",
                "Inferred.Model.BH_fdr",
                "mean_exp", "mean_meth",
                "cor_cna_exp", "cor_cna_meth", "cor_exp_meth")

# ============================================================
# Step 2: Merge function for one dataset
# ============================================================
merge_stats_model <- function(stats, model, name) {
  merged <- merge(
    model, stats,
    by.x  = c("trio_row", "meth.row", "cna.row", "gene.row"),
    by.y  = c("trio_id",  "meth.row", "cna.row", "gene.row"),
    all.x = TRUE
  )
  
  # Checks
  n_missing_model <- length(setdiff(stats$trio_id, model$trio_row))
  check <- data.table(
    Dataset             = name,
    stats_rows          = nrow(stats),
    model_rows          = nrow(model),
    merged_rows         = nrow(merged),
    unmatched_in_merge  = merged[is.na(mean_exp) & is.na(mean_meth), .N],
    trios_without_model = n_missing_model,
    dup_trio_row_model  = model[duplicated(trio_row), .N]
  )
  
  setcolorder(merged, intersect(front_cols, names(merged)))
  setorder(merged, trio_row)
  list(merged = merged, check = check)
}

# ============================================================
# Step 3: Run for every dataset and save each merged table
# ============================================================
merged_list <- list()
check_list  <- list()

for (nm in names(datasets)) {
  cat("Processing", nm, "...\n")
  p     <- datasets[[nm]]
  stats <- fread(p$stats)
  model <- fread(p$model)
  
  res <- merge_stats_model(stats, model, nm)
  merged_list[[nm]] <- res$merged
  check_list[[nm]]  <- res$check
  
  fwrite(res$merged, p$out, sep = "\t")
  cat("  saved:", p$out, "\n")
}

# ============================================================
# Step 4: Merge checks for all datasets
# ============================================================
merge_checks <- rbindlist(check_list)
print(merge_checks)

# ============================================================
# Step 5: Combine all datasets into one long table
# ============================================================
all_merged <- rbindlist(merged_list, idcol = "Dataset", fill = TRUE)

# ============================================================
# Step 6: Summaries across datasets
# ============================================================

# 6a. Overall stats per dataset
summary_dataset <- all_merged[, .(
  n_trios             = .N,
  mean_of_mean_exp    = mean(mean_exp,      na.rm = TRUE),
  mean_of_mean_meth   = mean(mean_meth,     na.rm = TRUE),
  median_cor_cna_exp  = median(cor_cna_exp,  na.rm = TRUE),
  median_cor_cna_meth = median(cor_cna_meth, na.rm = TRUE),
  median_cor_exp_meth = median(cor_exp_meth, na.rm = TRUE)
), by = Dataset]
print(summary_dataset)

# 6b. Same stats broken down by inferred model
summary_model <- all_merged[, .(
  n_trios             = .N,
  mean_of_mean_exp    = mean(mean_exp,      na.rm = TRUE),
  mean_of_mean_meth   = mean(mean_meth,     na.rm = TRUE),
  median_cor_cna_exp  = median(cor_cna_exp,  na.rm = TRUE),
  median_cor_cna_meth = median(cor_cna_meth, na.rm = TRUE),
  median_cor_exp_meth = median(cor_exp_meth, na.rm = TRUE)
), by = .(Dataset, Inferred.Model.BH_fdr)]
setorder(summary_model, Dataset, Inferred.Model.BH_fdr)
print(summary_model)

# 6c. Model counts: one row per model, one column per dataset
model_counts <- dcast(all_merged, Inferred.Model.BH_fdr ~ Dataset,
                      value.var = "trio_row", fun.aggregate = length)
print(model_counts)

# ============================================================
# Step 7: Save combined outputs
# ============================================================
combined_dir <- file.path(base_dir, "MRTrios_combined")
dir.create(combined_dir, showWarnings = FALSE)

fwrite(all_merged,      file.path(combined_dir, "ALL_trio_Model_results_with_stats.txt"), sep = "\t")
fwrite(merge_checks,    file.path(combined_dir, "ALL_merge_checks.txt"),                   sep = "\t")
fwrite(summary_dataset, file.path(combined_dir, "ALL_summary_by_dataset.txt"),             sep = "\t")
fwrite(summary_model,   file.path(combined_dir, "ALL_summary_by_dataset_model.txt"),       sep = "\t")
fwrite(model_counts,    file.path(combined_dir, "ALL_model_counts.txt"),                   sep = "\t")
