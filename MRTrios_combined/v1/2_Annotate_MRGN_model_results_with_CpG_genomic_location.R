# ============================================================
# Annotate trio results with CpG genomic location
# (UCSC_RefGene_Group), collapsed into 3 groups: TSS / Body / UTR
# for posER BRCA, negER BRCA, BLCA, LIHC
# ============================================================

library(data.table)

base_dir <- "/Users/lianzuo/LZ/ResearchProject/Fulab"

# ============================================================
# Step 0: Settings
# ============================================================

# Annotation file (450K manifest) -- read once, used by every cancer
humanmeth_file <- file.path(base_dir, "MRTrios_analysis/raw_Data_Methyl",
                            "GPL13534_HumanMethylation450_15017482_v.1.1 2.csv")

# Column in each methylation file that holds the CpG probe ID
probe_col <- "Row.names"

# TRUE  = one row per (trio, Group): a CpG annotated "TSS200;TSS1500"
#         counts once as TSS, not twice
# FALSE = one row per raw region entry (same as the original LIHC script)
dedup_groups <- TRUE

# ---- Per-cancer source files (trios + methylation) ----------
# A cancer's trio numbering and meth.row indices refer to ITS OWN trio
# and methylation files, so the probe lookup is built per cancer.
# posER and negER BRCA share the same trios and meth file.
sources <- list(
  BRCA = list(
    trios = file.path(base_dir, "MRTrios_analysis/raw_Data_Methyl/trio.final.protein.coding.txt"),
    meth  = file.path(base_dir, "MRTrios_analysis/raw_Data_Methyl/split.names.TCGA.meth.logit.rds")
  ),
  BLCA = list(   # <-- CHECK these two BLCA paths
    trios = file.path(base_dir, "GDCdata/TCGA-BLCA/trio.final.protein.coding.txt"),
    meth  = file.path(base_dir, "GDCdata/TCGA-BLCA/split.names.BLCA.meth.logit.rds")
  ),
  LIHC = list(
    trios = file.path(base_dir, "GDCdata/TCGA-LIHC/trio.final.final.protein.coding.txt"),
    meth  = file.path(base_dir, "GDCdata/TCGA-LIHC/split.names.LIHC.meth.logit.rds")
  )
)

# ---- Results to annotate ------------------------------------
# Uses the model + stats merged files from the previous script.
# To annotate the model-only files instead, point `input` at
# *_trio_Model_results_ALL_with_BH_fdr_qval_byLZ.txt
datasets <- list(
  posER_BRCA = list(
    source = "BRCA",
    input  = file.path(base_dir, "MRTrios_BRCA/Output_posER_BRCA/posER_BRCA_trio_Model_results_with_stats.txt"),
    out    = file.path(base_dir, "MRTrios_BRCA/Output_posER_BRCA/LZ_trio_results_posER_BRCA_with_stats_location.rds")
  ),
  negER_BRCA = list(
    source = "BRCA",
    input  = file.path(base_dir, "MRTrios_BRCA/Output_negER_BRCA/negER_BRCA_trio_Model_results_with_stats.txt"),
    out    = file.path(base_dir, "MRTrios_BRCA/Output_negER_BRCA/LZ_trio_results_negER_BRCA_with_stats_location.rds")
  ),
  BLCA = list(
    source = "BLCA",
    input  = file.path(base_dir, "MRTrios_BLCA/Output_BLCA/BLCA_trio_Model_results_with_stats.txt"),
    out    = file.path(base_dir, "MRTrios_BLCA/Output_BLCA/LZ_trio_results_BLCA_with_stats_location.rds")
  ),
  LIHC = list(
    source = "LIHC",
    input  = file.path(base_dir, "MRTrios_LIHC/Output_LIHC/LIHC_trio_Model_results_with_stats.txt"),
    out    = file.path(base_dir, "MRTrios_LIHC/Output_LIHC/LZ_trio_results_LIHC_with_stats_location.rds")
  )
)

# ---- Check all input files exist before doing any work -----
all_inputs <- c(humanmeth_file,
                unlist(lapply(sources, unlist)),
                sapply(datasets, `[[`, "input"))
missing_files <- all_inputs[!file.exists(all_inputs)]
if (length(missing_files) > 0) {
  stop("These input files were not found:\n  ",
       paste(missing_files, collapse = "\n  "))
}

# ============================================================
# Step 1: Load the 450K annotation once
# ============================================================
humanmeth <- read.csv(humanmeth_file, skip = 7, header = TRUE)
annot <- as.data.table(humanmeth)[, .(Name, UCSC_RefGene_Group)]
annot <- annot[!is.na(Name) & Name != ""]
rm(humanmeth); invisible(gc())
cat(sprintf("Annotation: %d probes\n", nrow(annot)))

# ============================================================
# Step 2: Helper functions
# ============================================================

# Read only the probe-ID column of a methylation file
get_probe_names <- function(file, col = probe_col) {
  if (grepl("\\.rds$", file, ignore.case = TRUE)) {
    m <- readRDS(file)
    ids <- as.character(m[[col]])
    rm(m); invisible(gc())
  } else {
    ids <- as.character(fread(file, select = col)[[1]])
  }
  ids
}

# Build lookup: trio_row -> meth.row -> probe Name -> Group
build_location_lookup <- function(trios_file, meth_file, annot, dedup = TRUE) {
  trios  <- fread(trios_file)
  probes <- get_probe_names(meth_file)
  
  stopifnot(
    "trios$meth.row exceeds number of rows in methylation file" =
      max(trios$meth.row, na.rm = TRUE) <= length(probes)
  )
  
  lk <- data.table(trio_row = seq_len(nrow(trios)),
                   meth.row = trios$meth.row,
                   Name     = probes[trios$meth.row])
  stopifnot("some meth.row values map to a missing probe Name" = !anyNA(lk$Name))
  
  # Attach raw UCSC_RefGene_Group
  a <- annot[Name %in% lk$Name]
  stopifnot("annotation has duplicate Names" = !anyDuplicated(a$Name))
  lk <- merge(lk, a, by = "Name", all.x = TRUE)
  cat(sprintf("  trios without annotation match: %d of %d\n",
              lk[is.na(UCSC_RefGene_Group), .N], nrow(lk)))
  
  # Expand semicolon-separated regions (blank / NA -> one NA region)
  parts <- strsplit(as.character(lk$UCSC_RefGene_Group), ";", fixed = TRUE)
  parts[lengths(parts) == 0] <- NA_character_
  reps <- lengths(parts)
  
  long <- data.table(
    trio_row           = rep(lk$trio_row, reps),
    meth.row           = rep(lk$meth.row, reps),
    Name               = rep(lk$Name,     reps),
    UCSC_RefGene_Group = rep(lk$UCSC_RefGene_Group, reps),   # full original string
    Region             = unlist(parts)
  )
  
  # Collapse 6 categories into 3
  long[, Group := fcase(
    Region %in% c("TSS200", "TSS1500"), "TSS",
    Region %in% c("Body", "1stExon"),   "Body",
    Region %in% c("3'UTR", "5'UTR"),    "5'/3'UTR",
    default = NA_character_              # blank / intergenic
  )]
  
  if (dedup) {
    long <- unique(long, by = c("trio_row", "Group"))
    long[, Region := NULL]
  }
  setorder(long, trio_row)
  long
}

# ============================================================
# Step 3: Build one lookup per cancer source
# ============================================================
lookups <- list()
for (s in names(sources)) {
  cat("Building location lookup for", s, "...\n")
  lookups[[s]] <- build_location_lookup(sources[[s]]$trios, sources[[s]]$meth,
                                        annot, dedup = dedup_groups)
}

# ============================================================
# Step 4: Annotate each dataset's results
# ============================================================
annotated_list <- list()
check_list     <- list()

for (nm in names(datasets)) {
  cat("Annotating", nm, "...\n")
  p   <- datasets[[nm]]
  res <- fread(p$input)
  lk  <- lookups[[p$source]]
  
  out <- merge(res, lk, by = c("trio_row", "meth.row"),
               all.x = TRUE, allow.cartesian = TRUE)
  setorder(out, trio_row)
  
  check_list[[nm]] <- data.table(
    Dataset          = nm,
    input_rows       = nrow(res),
    input_trios      = uniqueN(res$trio_row),
    output_rows      = nrow(out),
    trios_no_probe   = out[is.na(Name), uniqueN(trio_row)],
    trios_no_group   = out[, all(is.na(Group)), by = trio_row][V1 == TRUE, .N]
  )
  
  saveRDS(out, p$out, compress = "xz")
  cat("  saved:", p$out, "\n")
  annotated_list[[nm]] <- out
}

# ============================================================
# Step 5: Checks and summaries across datasets
# ============================================================
annot_checks <- rbindlist(check_list)
print(annot_checks)

all_annot <- rbindlist(annotated_list, idcol = "Dataset", fill = TRUE)

# Trio counts per Group x inferred model, one column per dataset
group_model_counts <- dcast(
  all_annot[!is.na(Group)],
  Group + Inferred.Model.BH_fdr ~ Dataset,
  value.var = "trio_row", fun.aggregate = uniqueN
)
print(group_model_counts)

# Mean stats per Group (uses the stats columns from the merge step)
group_stats <- all_annot[!is.na(Group), .(
  n_trios             = uniqueN(trio_row),
  mean_of_mean_exp    = mean(mean_exp,      na.rm = TRUE),
  mean_of_mean_meth   = mean(mean_meth,     na.rm = TRUE),
  median_cor_cna_exp  = median(cor_cna_exp,  na.rm = TRUE),
  median_cor_cna_meth = median(cor_cna_meth, na.rm = TRUE),
  median_cor_exp_meth = median(cor_exp_meth, na.rm = TRUE)
), by = .(Dataset, Group)]
setorder(group_stats, Dataset, Group)
print(group_stats)

combined_dir <- file.path(base_dir, "MRTrios_combined")
dir.create(combined_dir, showWarnings = FALSE)
saveRDS(all_annot, file.path(combined_dir, "ALL_trio_results_with_stats_location.rds"), compress = "xz")
fwrite(annot_checks,       file.path(combined_dir, "ALL_location_checks.txt"),       sep = "\t")
fwrite(group_model_counts, file.path(combined_dir, "ALL_group_model_counts.txt"),    sep = "\t")
fwrite(group_stats,        file.path(combined_dir, "ALL_group_stats.txt"),           sep = "\t")
