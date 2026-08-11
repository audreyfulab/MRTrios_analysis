# ============================================================
# Annotate posER-BRCA trio results with CpG genomic location
# (UCSC_RefGene_Group), collapsed into 3 groups: TSS / Body / UTR.
# ============================================================

library(data.table)
library(dplyr)
library(tidyr)

# ============================================================
# Step 0: Paths
# ============================================================
raw_dir    <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/raw_Data_Methyl"
out_dir    <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_posER_BRCA"
model_file <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_posER_BRCA/posER_BRCA_trio_Model_results_ALL_with_BH_fdr_qval_byLZ.txt"

negER_model_file <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_negER_BRCA/negER_BRCA_trio_Model_results_ALL_with_BH_fdr_qval_byLZ.txt"


# NOTE: filename below says "v.1.1" but the original script's comments
# describe switching to a "v1.2" manifest to get a full Name match
# (0 mismatches instead of 29,462). Worth double-checking this is
# actually the v1.2 file and the "v.1.1" in the name is just a stale
# label — the assertion in Step 2 will catch it either way if it's
# actually the older manifest with unmatched Names.
humanmeth_file <- "GPL13534_HumanMethylation450_15017482_v.1.1 2.csv"

# ============================================================
# Step 1: Load data
# ============================================================
humanmeth <- fread(file.path(raw_dir, humanmeth_file), skip = 7, data.table = FALSE)
cat(sprintf("humanmeth: %d rows, %d cols\n", nrow(humanmeth), ncol(humanmeth)))

trios     <- fread(file.path(raw_dir, "trio.final.protein.coding.txt"), data.table = FALSE)
TCGA.meth <- readRDS(file.path(raw_dir, "split.names.TCGA.meth.logit.rds"))

# ============================================================
# Step 2: Map each trio's meth.row to its CpG probe Name
# ============================================================
# trios$meth.row indexes into TCGA.meth; subsetting must return
# exactly one row per trio, in trio order, with no unmatched (NA)
# rows — this is what lets us attach `Name` by position below.
meth_trio <- TCGA.meth[trios$meth.row, ]

stopifnot(
  "meth_trio row count doesn't match trios row count" =
    nrow(meth_trio) == nrow(trios),
  "some trios$meth.row values did not match TCGA.meth rows (NA Row.names)" =
    !anyNA(meth_trio$Row.names)
)

trios_add_Name <- trios
trios_add_Name$Name      <- meth_trio$Row.names
trios_add_Name$trios_row <- seq_len(nrow(trios_add_Name))

# ============================================================
# Step 3: Join CpG probe Name -> UCSC_RefGene_Group location
# ============================================================
info <- humanmeth %>%
  select(Name, UCSC_RefGene_Group) %>%
  filter(Name %in% trios_add_Name$Name)

stopifnot(
  "humanmeth has duplicate Name entries after filtering to trio Names — merge below assumes 1 row per Name" =
    !anyNA(info$Name) && !any(duplicated(info$Name))
)

n_unmatched <- sum(is.na(match(trios_add_Name$Name, info$Name)))
cat(sprintf("Trio meth Names not found in humanmeth annotation: %d of %d\n",
            n_unmatched, nrow(trios_add_Name)))

merge_trio_loc <- merge(trios_add_Name, info, by = "Name", all = FALSE) %>%
  arrange(trios_row)

cat(sprintf("merge_trio_loc: %d rows (expected %d if humanmeth fully covers trios)\n",
            nrow(merge_trio_loc), nrow(trios_add_Name)))

# Quick summary of the raw (pre-split) location categories present
print(table(unlist(strsplit(as.character(merge_trio_loc$UCSC_RefGene_Group), ";"))))

# ============================================================
# Step 4: Save the location lookup table
# ============================================================
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
loc_file <- file.path(out_dir, "final_protein_coding_location_byLZ.txt")

fwrite(merge_trio_loc, file = loc_file, sep = "\t",
       quote = FALSE, row.names = FALSE, col.names = TRUE)

# merge_trio_loc is already in memory, so re-use it directly below
# instead of re-reading the file we just wrote. The file itself is
# still saved to disk for reuse in later sessions/scripts.
loc <- merge_trio_loc

# ============================================================
# Step 5: Merge trio model results with CpG location
# ============================================================
Model_posER_BRCA <- fread(model_file, sep = "\t", data.table = FALSE)

PosModel <- Model_posER_BRCA
PosModel$Inferred.Model3 <- PosModel$Inferred.Model.BH_fdr

n_model_rows <- nrow(PosModel)

merge_pos_df <- left_join(
  PosModel, loc,
  by = c("trio_row" = "trios_row", "cna.row", "gene.row", "meth.row")
)

if (nrow(merge_pos_df) != n_model_rows) {
  warning(sprintf(
    "left_join changed row count: %d -> %d (check for duplicate keys in `loc`)",
    n_model_rows, nrow(merge_pos_df)
  ))
}

# One CpG can belong to multiple gene regions (semicolon-separated);
# expand so each region gets its own row.
pos_df <- merge_pos_df %>%
  separate_rows(UCSC_RefGene_Group, sep = ";")

# ============================================================
# Step 6: Collapse the 6 UCSC_RefGene_Group categories into 3
# ============================================================
pos_df_2 <- pos_df %>%
  mutate(Group = case_when(
    UCSC_RefGene_Group %in% c("TSS200", "TSS1500") ~ "TSS",
    UCSC_RefGene_Group %in% c("Body", "1stExon")    ~ "Body",
    UCSC_RefGene_Group %in% c("3'UTR", "5'UTR")     ~ "5'/3'UTR",
    TRUE ~ NA_character_   # includes blank/IGR (intergenic) entries
  ))

# ============================================================
# Step 7: Save final annotated results
# ============================================================
# final_file <- file.path(out_dir, "LZ_trio_results_pos_ALL_with_BH_fdr_qval_with_location.txt")
# fwrite(pos_df_2, final_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

final_file <- file.path(out_dir, "LZ_trio_results_pos_ALL_with_BH_fdr_qval_with_location.rds")
saveRDS(pos_df_2, final_file, compress = "xz")

cat(sprintf("Done. Final file: %s (%d rows)\n", final_file, nrow(pos_df_2)))

## Final file: /Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_posER_BRCA/LZ_trio_results_pos_ALL_with_BH_fdr_qval_with_location.txt (675071 rows)


# ============================================================
# Step 5.2: Merge trio model results with CpG location
# ============================================================
Model_negER_BRCA <- fread(negER_model_file, sep = "\t", data.table = FALSE)

negModel <- Model_negER_BRCA
negModel$Inferred.Model3 <- negModel$Inferred.Model.BH_fdr

n_model_rows <- nrow(negModel)

merge_neg_df <- left_join(
  negModel, loc,
  by = c("trio_row" = "trios_row", "cna.row", "gene.row", "meth.row")
)

if (nrow(merge_neg_df) != n_model_rows) {
  warning(sprintf(
    "left_join changed row count: %d -> %d (check for duplicate keys in `loc`)",
    n_model_rows, nrow(merge_neg_df)
  ))
}

# One CpG can belong to multiple gene regions (semicolon-separated);
# expand so each region gets its own row.
neg_df <- merge_neg_df %>%
  separate_rows(UCSC_RefGene_Group, sep = ";")

# ============================================================
# Step 6: Collapse the 6 UCSC_RefGene_Group categories into 3
# ============================================================
neg_df_2 <- neg_df %>%
  mutate(Group = case_when(
    UCSC_RefGene_Group %in% c("TSS200", "TSS1500") ~ "TSS",
    UCSC_RefGene_Group %in% c("Body", "1stExon")    ~ "Body",
    UCSC_RefGene_Group %in% c("3'UTR", "5'UTR")     ~ "5'/3'UTR",
    TRUE ~ NA_character_   # includes blank/IGR (intergenic) entries
  ))

# ============================================================
# Step 7: Save final annotated results
# ============================================================
out_dir    <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_negER_BRCA"

final_file <- file.path(out_dir, "LZ_trio_results_neg_ALL_with_BH_fdr_qval_with_location.txt")
fwrite(neg_df_2, final_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

final_file <- file.path(out_dir, "LZ_trio_results_neg_ALL_with_BH_fdr_qval_with_location.rds")
saveRDS(neg_df_2, final_file, compress = "xz")

cat(sprintf("Done. Final file: %s (%d rows)\n", final_file, nrow(neg_df_2)))

