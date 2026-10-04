# ============================================================
# Add human methylation probe + biomart gene information to all_annot
# (vectorized version of MRTrios::extractHumanMethProbeInfo)
# for posER BRCA, negER BRCA, BLCA, LIHC
#
# New columns:
#   IlmnID, UCSC_RefGene_Name, Genome_Build, CHR, MAPINFO,
#   UCSC_CpG_Islands_Name, Relation_to_UCSC_CpG_Island, Trio_row,
#   Biomart_GeneName, <all biomart columns>, diff_cpG_mapinfo,
#   diff_mapinfo_geneStart, gene_length, Methyl_mean, Methyl_sd
# ============================================================

library(data.table)

base_dir     <- "/Users/lianzuo/LZ/ResearchProject/Fulab"
raw_dir      <- file.path(base_dir, "MRTrios_analysis/raw_Data_Methyl")
combined_dir <- file.path(base_dir, "MRTrios_combined")

# ============================================================
# Step 0: Settings
# ============================================================

# Output of annotate_location_all_cancers.R
all_annot_file <- file.path(combined_dir, "ALL_trio_results_with_stats_location.rds")

humanmeth_file <- file.path(raw_dir, "GPL13534_HumanMethylation450_15017482_v.1.1 2.csv")
biomart_file   <- file.path(raw_dir, "ensembl37_genes_p13_biomart.txt")

# Which probe genes to keep:
#   "trio_gene"     = only the gene of that trio (recommended; one gene row per trio)
#   "any_trio_gene" = any gene that appears in any trio
#                     (what the original function did:
#                      filter(Biomart_GeneName %in% trios$Gene.name))
gene_filter <- "trio_gene"

# TRUE  = keep one Ensembl entry per trio gene (the one containing, or nearest to, the probe)
# FALSE = keep every Ensembl entry on the probe's chromosome
one_biomart_per_gene <- TRUE

# Preprocessed lists (meth_filter + unique_id) used for Methyl_mean / Methyl_sd.
# These hold exactly the samples used for each subgroup.
# If a file is missing, Methyl_mean falls back to mean_meth and Methyl_sd is NA.
filter_lists <- list(
  posER_BRCA = file.path(base_dir, "MRTrios_BRCA/process_data_BRCA/posER_BRCA_data_filter_list.rds"),
  negER_BRCA = file.path(base_dir, "MRTrios_BRCA/process_data_BRCA/negER_BRCA_data_filter_list.rds"),
  BLCA       = file.path(base_dir, "MRTrios_BLCA/process_data_BLCA/BLCA_data_filter_list.rds"),   # <-- CHECK
  LIHC       = file.path(base_dir, "MRTrios_LIHC/process_data_LIHC/LIHC_data_filter_list.rds")    # <-- CHECK
)

stopifnot(file.exists(all_annot_file), file.exists(humanmeth_file), file.exists(biomart_file))

# ============================================================
# Step 1: Load data
# ============================================================
all_annot <- readRDS(all_annot_file)
cat(sprintf("all_annot: %d rows, %d cols\n", nrow(all_annot), ncol(all_annot)))

humanmeth <- read.csv(humanmeth_file, skip = 7, header = TRUE)
hm <- as.data.table(humanmeth)[, .(Name, IlmnID, UCSC_RefGene_Name, Genome_Build, CHR,
                                   MAPINFO, UCSC_CpG_Islands_Name, Relation_to_UCSC_CpG_Island)]
hm <- hm[!is.na(Name) & Name != ""]
hm <- unique(hm, by = "Name")
rm(humanmeth); invisible(gc())

biomart <- unique(as.data.table(read.delim(biomart_file, header = TRUE)))
stopifnot(all(c("Gene.name", "Gene.start..bp.", "Gene.end..bp.") %in% names(biomart)))
cat(sprintf("humanmeth: %d probes | biomart: %d rows\n", nrow(hm), nrow(biomart)))

# ============================================================
# Step 2: Probe + gene info for one dataset
# ============================================================
build_probe_info <- function(tri, hm, biomart, gene_filter = "trio_gene") {
  # tri: one row per trio with trio_row, meth.row, Name (probe), Gene.name (trio gene)
  
  # -- humanmeth columns for each trio's probe --
  pi <- merge(tri, hm, by = "Name", all.x = TRUE)
  
  # -- distance from probe to CpG island midpoint --
  # island name looks like "chr1:2004858-2005346"
  isl <- as.character(pi$UCSC_CpG_Islands_Name)
  has_isl <- !is.na(isl) & grepl(":\\d+-\\d+$", isl)
  isl_start <- rep(NA_real_, nrow(pi)); isl_end <- isl_start
  isl_start[has_isl] <- as.numeric(sub("^.*:(\\d+)-(\\d+)$", "\\1", isl[has_isl]))
  isl_end[has_isl]   <- as.numeric(sub("^.*:(\\d+)-(\\d+)$", "\\2", isl[has_isl]))
  pi[, diff_cpG_mapinfo := (isl_start + isl_end) / 2 - as.numeric(MAPINFO)]
  
  # -- one row per unique gene in UCSC_RefGene_Name --
  genes <- lapply(strsplit(as.character(pi$UCSC_RefGene_Name), ";", fixed = TRUE), unique)
  genes[lengths(genes) == 0] <- NA_character_
  long <- pi[rep(seq_len(.N), lengths(genes))]
  long[, Biomart_GeneName := unlist(genes)]
  long <- long[!is.na(Biomart_GeneName) & Biomart_GeneName != ""]
  
  # -- keep only relevant genes --
  if (gene_filter == "trio_gene") {
    long <- long[Biomart_GeneName == Gene.name]
  } else {
    long <- long[Biomart_GeneName %in% unique(tri$Gene.name)]
  }
  
  # -- join biomart (inner: only genes found in biomart, as in the original) --
  bm <- copy(biomart)
  bm[, Biomart_GeneName := Gene.name]
  long <- merge(long, bm, by = "Biomart_GeneName", allow.cartesian = TRUE,
                suffixes = c("", ".biomart"))
  
  # Keep only the gene entry on the probe's own chromosome; this drops
  # copies on patch / alternate-haplotype scaffolds (e.g. HG991_PATCH)
  long <- long[as.character(Chromosome.scaffold.name) == as.character(CHR)]
  
  # If a gene name still has several Ensembl entries on the same chromosome,
  # optionally keep the one whose gene body contains the probe, or else the nearest
  if (one_biomart_per_gene) {
    pos <- as.numeric(long$MAPINFO)
    long[, dist_to_gene := pmax(Gene.start..bp. - pos, pos - Gene.end..bp., 0)]
    setorder(long, trio_row, Biomart_GeneName, dist_to_gene)
    long <- unique(long, by = c("trio_row", "Biomart_GeneName"))
    long[, dist_to_gene := NULL]
  }
  
  long[, diff_mapinfo_geneStart := as.numeric(MAPINFO) - Gene.start..bp.]
  long[, gene_length            := Gene.end..bp. - Gene.start..bp.]
  
  # Gene.name.biomart is identical to Biomart_GeneName; keep the trio Gene.name only
  long[, c("Gene.name", "Gene.name.biomart") := NULL]
  long[, Trio_row := trio_row]
  long
}

# ============================================================
# Step 3: Methylation mean / SD per meth.row for one dataset
# ============================================================
methyl_mean_sd <- function(filter_file) {
  fl <- readRDS(filter_file)
  M  <- as.matrix(fl$meth_filter[, fl$unique_id])
  storage.mode(M) <- "double"
  mu <- rowMeans(M, na.rm = TRUE)
  n  <- rowSums(!is.na(M))
  s  <- sqrt(rowSums((M - mu)^2, na.rm = TRUE) / (n - 1))
  mu[is.nan(mu)] <- NA; s[n < 2] <- NA
  data.table(meth.row    = as.integer(fl$meth_filter$meth.row),
             Methyl_mean = mu,
             Methyl_sd   = s)
}

# ============================================================
# Step 4: Run for every dataset
# ============================================================
bm_cols  <- setdiff(names(biomart), "Gene.name")
new_cols <- c("IlmnID", "UCSC_RefGene_Name", "Genome_Build", "CHR", "MAPINFO",
              "UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island", "Trio_row",
              "Biomart_GeneName", bm_cols,
              "diff_cpG_mapinfo", "diff_mapinfo_geneStart", "gene_length",
              "Methyl_mean", "Methyl_sd")

out_list   <- list()
check_list <- list()

for (nm in unique(all_annot$Dataset)) {
  cat("Processing", nm, "...\n")
  ds  <- all_annot[Dataset == nm]
  tri <- unique(ds[, .(trio_row, meth.row, Name, Gene.name)])
  
  info <- build_probe_info(tri, hm, biomart, gene_filter)
  
  # Methyl_mean / Methyl_sd
  if (!is.null(filter_lists[[nm]]) && file.exists(filter_lists[[nm]])) {
    info <- merge(info, methyl_mean_sd(filter_lists[[nm]]), by = "meth.row", all.x = TRUE)
  } else {
    message("  filter list not found for ", nm,
            " -> Methyl_mean = mean_meth, Methyl_sd = NA")
    mm <- unique(ds[, .(trio_row, Methyl_mean = mean_meth)])
    info <- merge(info, mm, by = "trio_row", all.x = TRUE)
    info[, Methyl_sd := NA_real_]
  }
  
  # Attach to all_annot rows (left join keeps trios without probe/gene info)
  info <- info[, c("trio_row", "meth.row", new_cols), with = FALSE]
  out  <- merge(ds, info, by = c("trio_row", "meth.row"),
                all.x = TRUE, allow.cartesian = TRUE)
  setcolorder(out, c(names(ds), new_cols))
  setorder(out, trio_row)
  
  check_list[[nm]] <- data.table(
    Dataset               = nm,
    rows_before           = nrow(ds),
    rows_after            = nrow(out),
    trios                 = uniqueN(ds$trio_row),
    trios_with_probe_info = out[!is.na(IlmnID), uniqueN(trio_row)],
    trios_with_biomart    = out[!is.na(Biomart_GeneName), uniqueN(trio_row)],
    trios_multi_biomart   = info[, .N, by = trio_row][N > 1, .N],
    max_diff_Methyl_mean_vs_mean_meth =
      out[, suppressWarnings(max(abs(Methyl_mean - mean_meth), na.rm = TRUE))]
  )
  out_list[[nm]] <- out
}

probe_checks <- rbindlist(check_list)
print(probe_checks)

# ============================================================
# Step 5: Combine and save
# ============================================================
all_annot_probe <- rbindlist(out_list, fill = TRUE)

saveRDS(all_annot_probe,
        file.path(combined_dir, "ALL_trio_results_with_stats_location_probeinfo.rds"),
        compress = "xz")
fwrite(probe_checks, file.path(combined_dir, "ALL_probeinfo_checks.txt"), sep = "\t")

# Optional: one file per dataset
for (nm in names(out_list)) {
  saveRDS(out_list[[nm]],
          file.path(combined_dir, paste0(nm, "_trio_results_with_stats_location_probeinfo.rds")),
          compress = "xz")
}
