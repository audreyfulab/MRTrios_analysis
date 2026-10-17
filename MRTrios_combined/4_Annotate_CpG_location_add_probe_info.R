# ============================================================
# Annotate trio results with Illumina 450K probe information
# and Ensembl (GRCh37) gene information
#
# Output: one row per trio, with
#   - probe ID and coordinates (TCGA.meth)
#   - probe annotation (Illumina 450K manifest)
#   - gene region / accession for the trio's own gene
#   - Ensembl gene info (biomart, GRCh37)
#   - diff_cpG_mapinfo, diff_mapinfo_geneStart, gene_length
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)   # dplyr, tidyr, stringr
})

# ------------------------------------------------------------
# Step 0: Paths and settings
# ------------------------------------------------------------
base_dir     <- "/Users/lianzuo/LZ/ResearchProject/Fulab"
raw_dir      <- file.path(base_dir, "MRTrios_analysis/raw_Data_Methyl")
combined_dir <- file.path(base_dir, "MRTrios_combined")

all_annot_file <- file.path(combined_dir, "ALL_trio_Model_results_with_stats_v2.rds")
meth_file      <- file.path(raw_dir, "split.names.TCGA.meth.logit.rds")
humanmeth_file <- file.path(raw_dir, "GPL13534_HumanMethylation450_15017482_v.1.1 2.csv")
biomart_file   <- file.path(raw_dir, "ensembl37_genes_p13_biomart.txt")

# Dataset(s) whose meth.row refers to rows of the TCGA.meth file above
datasets_use <- "posER_BRCA"

# TRUE  = keep only trios whose gene is listed in the probe's UCSC_RefGene_Name
# FALSE = keep all trios; gene-region columns are NA when Illumina does not
#         annotate the probe to that gene
require_gene_match <- TRUE

# Trio gene names that are named differently in the GRCh37 biomart export
gene_alias <- c(LINS1  = "LINS",
                LVRN   = "AQPEP",
                PHF16  = "JADE3",
                GGTA1  = "GGTA1P",
                PLSCR3 = "TMEM256-PLSCR3")

# Genes that exist only on a patch scaffold in biomart
# (GRCh37 patch rows are given in chromosome coordinates, so they can be used)
patch_genes <- c("SEBOX")

primary_chr <- c(1:22, "X", "Y", "MT")

out_file <- file.path(combined_dir,
                      paste0("trio_probe_gene_annot_",
                             paste(datasets_use, collapse = "_"), ".rds"))
out_file_2 <- file.path(combined_dir,
                      paste0("trio_probe_gene_annot_loc_split_",
                             paste(datasets_use, collapse = "_"), ".rds"))
stopifnot(file.exists(all_annot_file), file.exists(meth_file),
          file.exists(humanmeth_file), file.exists(biomart_file))

# ------------------------------------------------------------
# Step 1: Load data
# ------------------------------------------------------------
all_annot <- readRDS(all_annot_file) %>%
  as_tibble() %>%
  filter(Dataset %in% datasets_use)
stopifnot(all(c("Dataset", "trio_row", "meth.row", "Gene.name") %in% names(all_annot)))
cat(sprintf("all_annot (%s): %d rows\n", paste(datasets_use, collapse = ", "), nrow(all_annot)))

TCGA.meth <- readRDS(meth_file)
humanmeth <- read.csv(humanmeth_file, skip = 7, header = TRUE, stringsAsFactors = FALSE)
biomart   <- read.delim(biomart_file, header = TRUE, stringsAsFactors = FALSE) %>% distinct()

# ------------------------------------------------------------
# Step 2: Methylation probe key (meth.row -> probe ID, position)
# ------------------------------------------------------------
meth_key <- tibble(
  meth.row           = as.integer(rownames(TCGA.meth)),
  Row.names          = as.character(TCGA.meth$Row.names),
  Chromosome         = sub("^chr", "", as.character(TCGA.meth$Chromosome)),
  Genomic_Coordinate = as.integer(TCGA.meth$Genomic_Coordinate)
)
stopifnot(!anyNA(meth_key$meth.row), !anyDuplicated(meth_key$meth.row))

rm(TCGA.meth); invisible(gc())   # only the four key columns are needed

# ------------------------------------------------------------
# Step 3: Illumina 450K manifest
# ------------------------------------------------------------
gene_cols <- c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")

humanmeth_probe <- humanmeth %>%
  select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
         Forward_Sequence, Genome_Build, CHR, MAPINFO, SourceSeq, Strand,
         all_of(gene_cols),
         UCSC_CpG_Islands_Name, Relation_to_UCSC_CpG_Island) %>%
  filter(IlmnID == Name) %>%            # drops the [Controls] section (851 rows)
  distinct() %>%
  mutate(CHR = as.character(CHR), MAPINFO = as.integer(MAPINFO))
stopifnot(!anyDuplicated(humanmeth_probe$Name))

rm(humanmeth); invisible(gc())

# 3a. Probe-level annotation: one row per probe
probe_annot <- humanmeth_probe %>% select(-all_of(gene_cols))

# 3b. Gene-level annotation: one row per probe + gene
#     (a gene listed several times for one probe = several transcripts;
#      their accessions / regions are collapsed with ";")
probe_gene <- humanmeth_probe %>%
  select(Name, all_of(gene_cols)) %>%
  separate_rows(all_of(gene_cols), sep = ";") %>%
  filter(UCSC_RefGene_Name != "") %>%
  group_by(Name, UCSC_RefGene_Name) %>%
  summarise(UCSC_RefGene_Accession = paste(unique(UCSC_RefGene_Accession), collapse = ";"),
            UCSC_RefGene_Group     = paste(unique(UCSC_RefGene_Group),     collapse = ";"),
            .groups = "drop")

# ------------------------------------------------------------
# Step 4: biomart, one row per gene name
#         (primary chromosomes first; patch rows only for patch_genes)
# ------------------------------------------------------------
biomart_gene <- biomart %>%
  filter(Chromosome.scaffold.name %in% primary_chr | Gene.name %in% patch_genes) %>%
  arrange(!Chromosome.scaffold.name %in% primary_chr) %>%
  distinct(Gene.name, .keep_all = TRUE)

# ------------------------------------------------------------
# Step 5: Combine
# ------------------------------------------------------------
gene_join <- if (require_gene_match) inner_join else left_join

df_comb_all <- all_annot %>%
  left_join(meth_key, by = "meth.row", relationship = "many-to-one") %>%
  left_join(probe_annot, by = c("Row.names" = "Name"), relationship = "many-to-one") %>%
  gene_join(probe_gene,
            by = c("Row.names" = "Name", "Gene.name" = "UCSC_RefGene_Name"),
            relationship = "many-to-one") %>%
  mutate(biomart_name = coalesce(unname(gene_alias[Gene.name]), Gene.name)) %>%
  left_join(biomart_gene, by = c("biomart_name" = "Gene.name"),
            relationship = "many-to-one")

# ------------------------------------------------------------
# Step 6: Distance / length columns
# ------------------------------------------------------------
cpg_pos <- str_match(df_comb_all$UCSC_CpG_Islands_Name, ":(\\d+)-(\\d+)")

df_comb_all <- df_comb_all %>%
  mutate(
    diff_cpG_mapinfo       = (as.numeric(cpg_pos[, 2]) + as.numeric(cpg_pos[, 3])) / 2 - MAPINFO,
    diff_mapinfo_geneStart = MAPINFO - as.numeric(Gene.start..bp.),
    gene_length            = as.numeric(Gene.end..bp.) - as.numeric(Gene.start..bp.)
  )
df_comb_all_2 <- df_comb_all %>%
  separate_rows(UCSC_RefGene_Group,
                sep = ";")

df_comb_all_final <- df_comb_all_2 %>% select(-UCSC_RefGene_Accession) %>% distinct()
# ------------------------------------------------------------
# Step 7: Checks and save
# ------------------------------------------------------------
cat(sprintf("Trios in:  %d\n", nrow(all_annot)))
cat(sprintf("Trios out: %d (%s)\n", nrow(df_comb_all),
            if (require_gene_match) "gene matched to probe annotation" else "all trios kept"))

if (anyDuplicated(df_comb_all[c("Dataset", "trio_row")]) > 0)
  warning("Some trios appear more than once in the output")

cat("Missing values:\n")
print(colSums(is.na(df_comb_all[c("MAPINFO", "diff_cpG_mapinfo",
                                  "diff_mapinfo_geneStart", "gene_length")])))

cat("Trio genes without a biomart match:\n")
print(df_comb_all %>% filter(is.na(gene_length)) %>% count(Gene.name, sort = TRUE))

saveRDS(df_comb_all, out_file)
cat("Saved:", out_file, "\n")

saveRDS(df_comb_all_2, out_file_2)
cat("Saved:", out_file_2, "\n")

# ================================================================================================
# Draft Version:
# ================================================================================================
# Add human methylation probe + biomart gene information to all_annot
# (vectorized version of MRTrios::extractHumanMethProbeInfo)

library(data.table)
library(tidyverse)

setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/raw_Data_Methyl")
TCGA.meth   <- readRDS("split.names.TCGA.meth.logit.rds")
trios       <- fread("trio.final.protein.coding.txt")
base_dir     <- "/Users/lianzuo/LZ/ResearchProject/Fulab"
raw_dir      <- file.path(base_dir, "MRTrios_analysis/raw_Data_Methyl")
combined_dir <- file.path(base_dir, "MRTrios_combined")
all_annot_file <- file.path(combined_dir, "ALL_trio_Model_results_with_stats_v2.rds")
humanmeth_file <- file.path(raw_dir, "GPL13534_HumanMethylation450_15017482_v.1.1 2.csv")
biomart_file   <- file.path(raw_dir, "ensembl37_genes_p13_biomart.txt")


all_annot <- readRDS(all_annot_file)
##  biomart: one row per gene name, primary chromosomes only (drops *_PATCH scaffolds)
biomart <- unique(as.data.table(read.delim(biomart_file, header = TRUE)))
biomart_1 <- biomart %>%
  filter(Chromosome.scaffold.name %in% c(1:22, "X", "Y", "MT")) %>%
  distinct(Gene.name, .keep_all = TRUE)

humanmeth <- read.csv(humanmeth_file, skip = 7, header = TRUE)
humanmeth_filter <- humanmeth %>% select(c(1:6,10:14,17,22:26)) %>% distinct()
humanmeth_filter2 <- humanmeth_filter %>%
  filter(IlmnID == Name)
## IlmnID != Name:  851 rows 

meth_filter=TCGA.meth  %>%  select(1:4)%>% mutate(meth.row=rownames(TCGA.meth))
meth_key <- meth_filter %>%
  mutate(meth.row           = as.integer(meth.row),
         Row.names          = as.character(Row.names),
         Chromosome         = sub("^chr", "", as.character(Chromosome)),
         Genomic_Coordinate = as.integer(Genomic_Coordinate)) %>%
  select(meth.row, Row.names, Chromosome, Genomic_Coordinate)

comb1=meth_key %>% left_join(humanmeth_filter2,by = c("Row.names" ="Name" ))
all_annot_posER_BRCA=all_annot %>% filter(Dataset=="posER_BRCA")
## 4. Combine

comb2 <- all_annot_posER_BRCA %>%
  as_tibble() %>%
  left_join(comb1, by =  c("meth.row"="meth.row"))
comb2_long <- comb2 %>%
  separate_rows(UCSC_RefGene_Name,
                UCSC_RefGene_Accession,
                UCSC_RefGene_Group,
                sep = ";")
comb_final=comb2_long %>% filter(Gene.name==UCSC_RefGene_Name)

df_comb <- comb_final %>% left_join(biomart_1,by="Gene.name")

df_comb_2 <- df_comb %>%
  mutate(
    # CpG island "chr16:53467000-53468000" -> start and end
    cpg_start = as.numeric(str_match(UCSC_CpG_Islands_Name, ":(\\d+)-(\\d+)")[, 2]),
    cpg_end   = as.numeric(str_match(UCSC_CpG_Islands_Name, ":(\\d+)-(\\d+)")[, 3]),
    
    diff_cpG_mapinfo       = (cpg_start + cpg_end) / 2 - MAPINFO,
    diff_mapinfo_geneStart = MAPINFO - as.numeric(Gene.start..bp.),
    gene_length            = as.numeric(Gene.end..bp.) - as.numeric(Gene.start..bp.)
  ) %>%
  select(-cpg_start, -cpg_end)

##511,966 rows
df_comb_3 <- df_comb_2 %>% select(-UCSC_RefGene_Accession) %>% distinct()
# ================================================================================================
