# ============================================================
# Annotate trio results with Illumina 450K probe information
# and Ensembl (GRCh37) gene information -- ALL datasets
#
# Datasets: posER_BRCA, negER_BRCA, LIHC, BLCA
#   - Illumina 450K manifest and biomart are shared (loaded once)
#   - meth.row is dataset-specific: each dataset's meth.row refers to
#     rows of THAT cancer's TCGA.meth file, so the probe key is built
#     per dataset and joined by (Dataset, meth.row)
#
# Outputs:
#   out_file   : one row per trio, with
#                - probe ID and coordinates (TCGA.meth)
#                - probe annotation (Illumina 450K manifest)
#                - gene region / accession for the trio's own gene
#                - Ensembl gene info (biomart, GRCh37)
#                - diff_cpG_mapinfo, diff_mapinfo_geneStart, gene_length
#                - position_status flags, n_accession
#   out_file_2 : one row per trio x Location (see Step 7), where
#                - trios with ONE region            -> 1 row
#                - multi_region, 1 accession        -> split, 1 row per region
#                  (5'UTR + 1stExon or 1stExon + 3'UTR of one transcript)
#                - multi_region, >1 accession       -> REMOVED
#                  (different transcripts disagree on the region)
#
# Overall flow:
#   trio (Dataset, meth.row, Gene.name)
#     --[Dataset + meth.row]-->  probe ID (cg...) from that cancer's TCGA.meth
#     --[probe ID]----------->   probe annotation from the 450K manifest
#     --[probe ID + gene]---->   gene region (TSS200, Body, ...) for the trio gene
#     --[gene name]---------->   gene start / end / Ensembl ID from biomart
# ============================================================

## Load packages quietly.
##   data.table : kept for compatibility (TCGA.meth may be a data.table)
##   tidyverse  : dplyr (joins, mutate), tidyr (separate_rows, pivot_wider),
##                stringr (str_match / str_split), purrr (map_chr, walk)
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)   # dplyr, tidyr, stringr, purrr
})

# ------------------------------------------------------------
# Step 0: Paths and settings
# ------------------------------------------------------------
## All user-editable settings live here so the rest of the script
## never needs to be touched when files move or datasets change.

## Project folders
base_dir     <- "/Users/lianzuo/LZ/ResearchProject/Fulab"
raw_dir      <- file.path(base_dir, "MRTrios_analysis/raw_Data_Methyl")   # raw inputs
combined_dir <- file.path(base_dir, "MRTrios_combined")                   # results

## Inputs shared by every dataset:
##   all_annot_file : trio model results for ALL datasets (has a Dataset column)
##   humanmeth_file : Illumina HumanMethylation450 manifest (GPL13534);
##                    probe -> chromosome, position, genes, CpG island
##   biomart_file   : Ensembl GRCh37 gene table exported from biomart;
##                    gene name -> Ensembl ID, chromosome, start, end
all_annot_file <- file.path(combined_dir, "ALL_trio_Model_results_with_stats_v2.rds")
humanmeth_file <- file.path(raw_dir, "GPL13534_HumanMethylation450_15017482_v.1.1 2.csv")
biomart_file   <- file.path(raw_dir, "ensembl37_genes_p13_biomart.txt")

## Dataset-specific input: the methylation matrix used to build each
## dataset's trios. This matters because meth.row in all_annot is only a
## ROW NUMBER -- row 1000 in the BRCA file is a different probe from row
## 1000 in the LIHC file. Each dataset must therefore be matched to the
## exact file its meth.row values were taken from.
##
## >>> EDIT: the TCGA.meth file whose ROW NUMBERS each dataset's meth.row refers to.
##     posER_BRCA and negER_BRCA can point to the same file if both were run on
##     the same BRCA methylation matrix (it is only read once).
meth_files <- c(
  posER_BRCA = file.path(raw_dir, "split.names.TCGA.meth.logit.rds"),
  negER_BRCA = file.path(raw_dir, "split.names.TCGA.meth.logit.rds"),        # EDIT if different
  LIHC       = file.path("/Users/lianzuo/LZ/ResearchProject/Fulab/GDCdata/TCGA-LIHC/split.names.LIHC.meth.logit.rds"),   # EDIT
  BLCA       = file.path("/Users/lianzuo/LZ/ResearchProject/Fulab/GDCdata/TCGA-BLCA/split.names.BLCA.meth.logit.rds")    # EDIT
)

## Which datasets to process. By default all four (the names of meth_files).
## To rerun only some, give a subset, e.g. c("LIHC", "BLCA").
datasets_use <- names(meth_files)   # or a subset, e.g. c("LIHC", "BLCA")

## How strictly to match the trio gene to the probe:
## TRUE  = keep only trios whose gene is listed in the probe's UCSC_RefGene_Name
##         (i.e. Illumina says the probe lies in/near that gene)
## FALSE = keep all trios; gene-region columns are NA when Illumina does not
##         annotate the probe to that gene
require_gene_match <- TRUE

## Trio gene names (current HGNC symbols) that are named differently in the
## older GRCh37 biomart export. Without this, these genes would get NA for
## gene start/end. Format: c(trio_name = "biomart_name").
gene_alias <- c(LINS1  = "LINS",
                LVRN   = "AQPEP",
                PHF16  = "JADE3",
                GGTA1  = "GGTA1P",
                PLSCR3 = "TMEM256-PLSCR3")
## (if LIHC/BLCA trios bring new unmatched genes, add them here -- see the
##  "Trio genes without a biomart match" printout in Step 8)

## Genes that exist only on a patch scaffold in biomart (e.g. HG..._PATCH)
## are now found automatically in Step 4 (patch_genes).
## Ensembl stores each GRCh37 patch as a full-length alternative chromosome,
## so patch rows use whole-chromosome coordinates and can be used for these
## genes instead of being dropped. They are validated against the probe
## position in Step 8 (patch_check).
## For posER_BRCA this found: SEBOX (HG883_PATCH), PECAM1 (HG183_PATCH),
## BAHCC1 (HG271_PATCH) -- all chr17.

## The normal chromosomes; anything else in biomart is a patch / alt scaffold.
primary_chr <- c(1:22, "X", "Y", "MT")

## Output file names. If every dataset is processed the tag is "ALL";
## otherwise the dataset names are pasted together, e.g. "LIHC_BLCA".
##   out_file   : one row per trio (main result)
##   out_file_2 : one row per trio x Location (see Step 7)
tag <- if (setequal(datasets_use, names(meth_files))) "ALL" else paste(datasets_use, collapse = "_")
out_file   <- file.path(combined_dir, paste0("trio_probe_gene_annot_", tag, ".rds"))
#out_file_2 <- file.path(combined_dir, paste0("trio_probe_gene_annot_loc_split_", tag, ".rds"))
out_file_4 <- file.path(combined_dir, paste0("trio_probe_gene_annot_loc_split_allCancer", tag, ".rds"))
## Fail early (before any slow reading) if a setting or file is wrong:
##   1. every requested dataset must have an entry in meth_files
##   2. the shared inputs must exist
##   3. every dataset's methylation file must exist (names the missing ones)
stopifnot(all(datasets_use %in% names(meth_files)))
stopifnot(file.exists(all_annot_file), file.exists(humanmeth_file), file.exists(biomart_file))
missing_meth <- meth_files[datasets_use][!file.exists(meth_files[datasets_use])]
if (length(missing_meth) > 0)
  stop("Methylation file(s) not found:\n", paste(names(missing_meth), missing_meth, sep = ": ", collapse = "\n"))

# ------------------------------------------------------------
# Step 1: Load trio results + shared annotation
# ------------------------------------------------------------
## all_annot: one row per trio. Keep only the datasets being processed.
## Required columns:
##   Dataset   : which cancer / subtype the trio came from
##   trio_row  : trio ID within the dataset
##   meth.row  : row number of the methylation probe in that dataset's TCGA.meth
##   Gene.name : the gene in the trio
all_annot <- readRDS(all_annot_file) %>%
  as_tibble() %>%
  filter(Dataset %in% datasets_use)
stopifnot(all(c("Dataset", "trio_row", "meth.row", "Gene.name") %in% names(all_annot)))

## Print the number of trios per dataset so you can compare with Step 8.
cat("all_annot rows per dataset:\n"); print(count(all_annot, Dataset))

## Illumina manifest: the first 7 lines are a header block ([Heading] etc.),
## so skip them; line 8 holds the column names.
humanmeth <- read.csv(humanmeth_file, skip = 7, header = TRUE, stringsAsFactors = FALSE)

## biomart export: tab-delimited. distinct() removes exact duplicate rows
## (biomart repeats a gene once per transcript/attribute when exported).
biomart   <- read.delim(biomart_file, header = TRUE, stringsAsFactors = FALSE) %>% distinct()

# ------------------------------------------------------------
# Step 2: Methylation probe key per dataset
#         (Dataset, meth.row) -> probe ID, position
# ------------------------------------------------------------
## The TCGA.meth files are large (probes x samples). Only four columns are
## needed to turn a row number into a probe, so this function reads one
## file, keeps those columns, and frees the memory straight away.
##
##   meth.row           : the row name (= row number used in all_annot)
##   Row.names          : the probe ID, e.g. "cg00000029"
##   Chromosome         : "chr16" -> "16" (to match the manifest's CHR style)
##   Genomic_Coordinate : probe position in the TCGA file
read_meth_key <- function(f) {
  cat("Reading", f, "\n")
  m <- readRDS(f)
  key <- tibble(
    meth.row           = as.integer(rownames(m)),
    Row.names          = as.character(m$Row.names),
    Chromosome         = sub("^chr", "", as.character(m$Chromosome)),
    Genomic_Coordinate = as.integer(m$Genomic_Coordinate)
  )
  rm(m); invisible(gc())   # only the four key columns are needed
  ## Every row must have a valid, unique row number, otherwise the
  ## join in Step 5 could attach the wrong probe or duplicate trios.
  stopifnot(!anyNA(key$meth.row), !anyDuplicated(key$meth.row))
  key
}

## Read each DISTINCT file only once. posER_BRCA and negER_BRCA usually share
## the BRCA file, so it is read once and reused for both.
## keys_by_file is a list named by file path: list("<path>" = key tibble, ...)
uniq_files   <- unique(unname(meth_files[datasets_use]))
keys_by_file <- setNames(lapply(uniq_files, read_meth_key), uniq_files)

## Build one long key table with a Dataset column:
## for each dataset, look up its file's key and label it with the dataset
## name, then stack them. A shared file appears once per dataset using it.
## Result columns: Dataset, meth.row, Row.names, Chromosome, Genomic_Coordinate
meth_key <- bind_rows(lapply(datasets_use, function(d)
  keys_by_file[[ meth_files[[d]] ]] %>% mutate(Dataset = d, .before = 1)))
rm(keys_by_file); invisible(gc())

# ------------------------------------------------------------
# Step 3: Illumina 450K manifest (shared)
# ------------------------------------------------------------
## The three gene columns are semicolon-separated lists that line up by
## position, e.g. for one probe:
##   UCSC_RefGene_Name      = "TP53;TP53;WRAP53"
##   UCSC_RefGene_Accession = "NM_000546;NM_001126112;NR_..."
##   UCSC_RefGene_Group     = "TSS1500;TSS1500;Body"
## They are handled separately in 3b.
gene_cols <- c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")

## Keep the useful manifest columns:
##   IlmnID / Name            : probe ID
##   Address*/AlleleA/B Seq   : probe sequences and bead addresses
##   Forward_Sequence, SourceSeq, Strand, Genome_Build
##   CHR, MAPINFO             : chromosome and position (GRCh37)
##   gene columns (above)
##   UCSC_CpG_Islands_Name    : nearest CpG island, e.g. "chr16:53467000-53468000"
##   Relation_to_UCSC_CpG_Island : Island / N_Shore / S_Shelf / ... / blank (open sea)
## Then:
##   IlmnID == Name  -> keeps real probes; the [Controls] section at the end
##                      of the file has different values in these columns
##   distinct()      -> drop repeated rows
##   CHR as text (so "X"/"Y" and numbers share one type), MAPINFO as integer
humanmeth_probe <- humanmeth %>%
  select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
         Forward_Sequence, Genome_Build, CHR, MAPINFO, SourceSeq, Strand,
         all_of(gene_cols),
         UCSC_CpG_Islands_Name, Relation_to_UCSC_CpG_Island) %>%
  filter(IlmnID == Name) %>%            # drops the [Controls] section (851 rows)
  distinct() %>%
  mutate(CHR = as.character(CHR), MAPINFO = as.integer(MAPINFO))

## Each probe ID must appear once, so the probe join in Step 5 is one-to-one.
stopifnot(!anyDuplicated(humanmeth_probe$Name))
rm(humanmeth); invisible(gc())

## 3a. Probe-level annotation: one row per probe.
##     Everything except the gene columns (position, sequence, CpG island).
probe_annot <- humanmeth_probe %>% select(-all_of(gene_cols))

## 3b. Gene-level annotation: one row per probe + gene.
##     separate_rows() splits the three gene columns together, so the k-th
##     name stays paired with the k-th accession and k-th region:
##         cg01  TP53    NM_000546     TSS1500
##         cg01  TP53    NM_001126112  TSS1500
##         cg01  WRAP53  NR_...        Body
##     Probes with no gene (intergenic, empty name) are dropped.
##     A gene listed several times for one probe = several transcripts;
##     group by probe + gene and collapse their accessions / regions with ";"
##     (unique() avoids "TSS1500;TSS1500"):
##         cg01  TP53    NM_000546;NM_001126112  TSS1500
##         cg01  WRAP53  NR_...                  Body
##     This gives exactly one row per (probe, gene), which is what the
##     trio join (probe + trio gene) needs.
probe_gene <- humanmeth_probe %>%
  select(Name, all_of(gene_cols)) %>%
  separate_rows(all_of(gene_cols), sep = ";") %>%
  filter(UCSC_RefGene_Name != "") %>%
  group_by(Name, UCSC_RefGene_Name) %>%
  summarise(UCSC_RefGene_Accession = paste(unique(UCSC_RefGene_Accession), collapse = ";"),
            UCSC_RefGene_Group     = paste(unique(UCSC_RefGene_Group),     collapse = ";"),
            .groups = "drop")

# ------------------------------------------------------------
# Step 4: biomart, one row per gene name (shared)
# ------------------------------------------------------------
## A gene name can appear on several rows in biomart (primary chromosome
## plus patch/alt scaffolds). To get exactly one row per gene:
##   1. keep primary-chromosome rows, plus patch rows only for patch_genes
##   2. sort primary-chromosome rows first
##      (!x %in% primary_chr is FALSE for primary -> FALSE sorts before TRUE)
##   3. distinct(Gene.name, .keep_all = TRUE) keeps the FIRST row per gene,
##      i.e. the primary-chromosome row when one exists
## Columns used later: Gene.stable.ID, Gene.start..bp., Gene.end..bp., ...

## 4a. Find patch-only trio genes automatically.
##     A gene is "patch-only" if it never appears on a primary chromosome.
##     Genes that are on a patch AND a primary chromosome (e.g. TMEM97 on
##     HG883_PATCH) are not included -- they keep the primary coordinates.
trio_biomart_names <- unique(coalesce(unname(gene_alias[all_annot$Gene.name]),
                                      all_annot$Gene.name))

on_primary <- biomart %>%
  filter(Chromosome.scaffold.name %in% primary_chr) %>%
  pull(Gene.name) %>% unique()

patch_only_trio <- biomart %>%
  filter(!Chromosome.scaffold.name %in% primary_chr,
         !Gene.name %in% on_primary,
         Gene.name %in% trio_biomart_names) %>%
  distinct(Gene.name, Gene.stable.ID, Chromosome.scaffold.name,
           Gene.start..bp., Gene.end..bp.)

cat("Patch-only trio genes:\n"); print(patch_only_trio)

# Patch-only trio genes:
#   Gene.name  Gene.stable.ID Chromosome.scaffold.name Gene.start..bp. Gene.end..bp.
# 1    PECAM1 ENSG00000261371              HG183_PATCH        62399863      62491136
# 2    BAHCC1 ENSG00000266074              HG271_PATCH        79369275      79440126
# 3     SEBOX ENSG00000263225              HG883_PATCH        26691290      26692265

## A gene on several patches would get several rows; distinct() below would
## then keep an arbitrary one, so warn if that happens.
multi_patch <- count(patch_only_trio, Gene.name) %>% filter(n > 1)
if (nrow(multi_patch) > 0) {
  warning("Some trio genes are on more than one patch -- check which to use")
  print(multi_patch)
}

patch_genes <- unique(patch_only_trio$Gene.name)

## 4b. One row per gene name
biomart_gene <- biomart %>%
  filter(Chromosome.scaffold.name %in% primary_chr | Gene.name %in% patch_genes) %>%
  arrange(!Chromosome.scaffold.name %in% primary_chr) %>%
  distinct(Gene.name, .keep_all = TRUE)

# ------------------------------------------------------------
# Step 5: Combine  (meth key joined by Dataset + meth.row)
# ------------------------------------------------------------
## Pick the gene join based on the setting in Step 0:
##   inner_join -> trios whose gene is not annotated to the probe are dropped
##   left_join  -> they are kept with NA gene-region columns
gene_join <- if (require_gene_match) inner_join else left_join

## Chain of joins (each adds columns to the trio table):
##   1. meth_key    by Dataset + meth.row -> probe ID and TCGA position.
##                  Joining on BOTH columns is what makes this multi-cancer
##                  safe: a LIHC trio can only find a LIHC probe row.
##   2. probe_annot by probe ID          -> manifest position, CpG island, ...
##   3. probe_gene  by probe ID + trio gene -> the region(s) / accessions for
##                  THIS trio's gene only (not other genes at the same probe)
##   4. biomart_name: use the alias if the gene has one, otherwise the
##                  original name (coalesce takes the first non-NA)
##   5. biomart_gene by biomart_name     -> Ensembl ID, gene start / end
##
## relationship = "many-to-one" makes dplyr throw an error if the right-hand
## table has duplicate keys, so a join can never silently duplicate trios.
df_comb_all <- all_annot %>%
  left_join(meth_key, by = c("Dataset", "meth.row"), relationship = "many-to-one") %>%
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
## Pull the CpG island start and end out of names like
## "chr16:53467000-53468000". str_match returns a matrix:
##   column 1 = whole match, column 2 = start, column 3 = end
## (NA when the probe has no CpG island, i.e. open sea)
cpg_pos <- str_match(df_comb_all$UCSC_CpG_Islands_Name, ":(\\d+)-(\\d+)")

## New columns (all in bp, GRCh37):
##   diff_cpG_mapinfo       = CpG island midpoint - probe position
##                            (+ : island is downstream of the probe on the
##                             reference strand, - : upstream; NA if no island)
##   diff_mapinfo_geneStart = probe position - gene start
##                            NOTE: uses biomart "Gene start", which is the
##                            lower coordinate. For minus-strand genes the TSS
##                            is actually at Gene end, so read this with the
##                            gene strand in mind.
##   gene_length            = gene end - gene start
df_comb_all <- df_comb_all %>%
  mutate(
    diff_cpG_mapinfo       = (as.numeric(cpg_pos[, 2]) + as.numeric(cpg_pos[, 3])) / 2 - MAPINFO,
    diff_mapinfo_geneStart = MAPINFO - as.numeric(Gene.start..bp.),
    gene_length            = as.numeric(Gene.end..bp.) - as.numeric(Gene.start..bp.)
  )

# ------------------------------------------------------------
# Step 7: Mark multi-position probes (within each dataset)
#         and build the Location table
# ------------------------------------------------------------
## Some probes are ambiguous: they sit in several genes, or in several
## regions of one gene (different transcripts). These flags let you filter
## or stratify them later.

## Recognised gene regions (also used to order region combinations and
## the Location counts in Step 8).
prio <- c("TSS200", "TSS1500", "5'UTR", "1stExon", "Body", "3'UTR")

df_comb_all <- df_comb_all %>%
  ## (1) Several regions inside one cell, e.g. "TSS1500;Body":
  ##     number of regions = number of ";" + 1 (0 when empty / NA).
  ##     n_accession counts the RefSeq transcripts the same way:
  ##       1  -> one transcript (several regions = overlapping labels,
  ##             e.g. 5'UTR;1stExon of the same first exon)
  ##       >1 -> several transcripts (several regions = transcripts disagree)
  mutate(
    n_refgene_group     = if_else(UCSC_RefGene_Group == "" | is.na(UCSC_RefGene_Group),
                                  0L, str_count(UCSC_RefGene_Group, ";") + 1L),
    multi_refgene_group = n_refgene_group > 1,
    n_accession         = if_else(UCSC_RefGene_Accession == "" | is.na(UCSC_RefGene_Accession),
                                  0L, str_count(UCSC_RefGene_Accession, ";") + 1L)
  ) %>%
  ## (2) The same probe used by several trios in the same dataset.
  ##     Grouping by Dataset AND probe so that a probe shared by BRCA and LIHC
  ##     is not counted as "duplicated" just because it is in two cancers.
  ##       n_rows_probe : trios in this dataset using this probe
  ##       n_genes      : different Ensembl genes those trios point to
  group_by(Dataset, IlmnID) %>%          # per dataset, so BRCA/LIHC/BLCA don't mix
  mutate(
    n_rows_probe = n(),
    n_genes      = n_distinct(Gene.stable.ID, na.rm = TRUE),
    multi_row    = n_rows_probe > 1,
    multi_gene   = n_genes > 1
  ) %>%
  ungroup() %>%
  mutate(
    ## (3) One overall label. case_when checks in order, so a probe gets
    ##     the most "serious" label that applies:
    ##       multi_gene     : probe linked to >1 gene in this dataset
    ##       multi_region   : one gene, but >1 region
    ##       duplicated_row : same probe + gene in several trios
    ##       unique         : none of the above
    position_status = case_when(
      multi_gene          ~ "multi_gene",
      multi_refgene_group ~ "multi_region",
      multi_row           ~ "duplicated_row",
      TRUE                ~ "unique"
    )
  )

## (4) Location: one row per trio x region, used for region-level analysis.
##     a. Remove trios whose regions come from >1 transcript
##        (n_refgene_group > 1 AND n_accession > 1). Their location is
##        genuinely ambiguous, e.g. TSS200 for one transcript, Body for another.
##     b. Keep trios with one region (Location = that region).
##     c. Keep multi-region trios with ONE transcript (n_accession == 1),
##        i.e. 5'UTR;1stExon or 1stExon;3'UTR, and split them into one row
##        per region:
##          "5'UTR;1stExon"  ->  Location = "5'UTR"   (row 1)
##                               Location = "1stExon" (row 2)
##     The rule is applied to every trio regardless of position_status, so
##     multi_gene trios follow the same rule for their own gene's regions.
##     UCSC_RefGene_Group keeps the original string; Location is the split value.

a1=df_comb_all %>%  filter(!(n_refgene_group >1 & n_accession > 1))

table(a1$UCSC_RefGene_Group)
# 
# 1stExon 1stExon;3'UTR 1stExon;5'UTR         3'UTR 3'UTR;1stExon         5'UTR 5'UTR;1stExon 
# 34559          1047         27787         54472          1083         89907         22592 
# 5'UTR;Body          Body    Body;5'UTR       TSS1500        TSS200 
# 2        491688             8        209198        165027

a1 %>% filter(UCSC_RefGene_Group %in% c("5'UTR;Body","Body;5'UTR")) %>% select(trio_row,Row.names,Inferred.Model.BH_fdr,UCSC_RefGene_Accession,UCSC_RefGene_Group,biomart_name,MAPINFO) %>% arrange(MAPINFO)
# # A tibble: 10 × 7
# trio_row Row.names Inferred.Model.BH_fdr UCSC_RefGene_Accession UCSC_RefGene_Group biomart_name MAPINFO
# <int> <chr>     <chr>                 <chr>                  <chr>              <chr>          <int>
#   1   174586 cg206763… Other                 NM_001127212           Body;5'UTR         GAGE2A        4.92e7
#  2   289023 cg206763… Other                 NM_001472              5'UTR;Body         GAGE2C        4.92e7
# 3   174586 cg206763… M0.1                  NM_001127212           Body;5'UTR         GAGE2A        4.92e7
#  4   289023 cg206763… M0.1                  NM_001472              5'UTR;Body         GAGE2C        4.92e7
# 5   174582 cg132811… Other                 NM_001127212           Body;5'UTR         GAGE2A        4.92e7
#  6   289021 cg132811… Other                 NM_001472              Body;5'UTR         GAGE2C        4.92e7
# 7   174582 cg132811… Other                 NM_001127212           Body;5'UTR         GAGE2A        4.92e7
#  8   289021 cg132811… Other                 NM_001472              Body;5'UTR         GAGE2C        4.92e7
# 9   174585 cg205030… Other                 NM_001127212           Body;5'UTR         GAGE2A        4.92e7
# 10   174585 cg205030… M0.1                  NM_001127212           Body;5'UTR         GAGE2A        4.92e7

## 1) n_accession >1 & n_refgene_group =1 ==> keep but not split two rows (NM_022763;NM_001135095  TSS1500)
## 2) n_accession >1 & n_refgene_group > 1 ==> remove (NM_022763;NM_001135095.  TSS200;TSS1500)
## 3.1) n_accession =1 & n_refgene_group > 1 ==> loc=5'UTR or 3'UTR;1stExon , keep and split rows( NM_005611.  5'UTR;1stExon)
## 3.2) n_accession =1 & n_refgene_group > 1 ==> "5'UTR;Body","Body;5'UTR" ,remove (NM_001127212,   Body;5'UTR)
df_comb_all_2 <- df_comb_all %>%
  filter(!(n_refgene_group >1 & n_accession > 1)) %>%
  mutate(Location = UCSC_RefGene_Group, .after = UCSC_RefGene_Group) %>%
  separate_rows(Location, sep = ";") %>%
  mutate(Location = na_if(Location, ""))

# -------------------------------
## Multi-region combinations allowed to be split (single transcript only):
## the first exon overlaps the 5'UTR or the 3'UTR, so both labels are true
allowed_combo <- c("1stExon;5'UTR", "1stExon;3'UTR")    # sorted form

df_comb_all_3 <- df_comb_all %>%
  mutate(region_combo = map_chr(str_split(UCSC_RefGene_Group, ";"),
                                ~ paste(sort(unique(.x)), collapse = ";")),.after = UCSC_RefGene_Group)

df_comb_all_4 <- df_comb_all_3 %>%
  filter(
    n_refgene_group <= 1 |                                   # 1)  one region (any n_accession): keep
      (n_accession == 1 & region_combo %in% allowed_combo)   # 3.1) 1 accession, 5'UTR/3'UTR + 1stExon: keep
  ) %>%                                                      # 2) and 3.2) everything else: removed
  mutate(Location = UCSC_RefGene_Group, .after = UCSC_RefGene_Group) %>%
  separate_rows(Location, sep = ";") %>%                    # only 3.1) rows are split
  mutate(Location = na_if(Location, ""))




# df_comb_all_2 %>% filter(Dataset=="posER_BRCA"&biomart_name=="RBL2") %>% select(trio_row,Row.names,Inferred.Model.BH_fdr,Location,biomart_name,MAPINFO) %>% arrange(MAPINFO)
# # A tibble: 13 × 6
# trio_row Row.names  Inferred.Model.BH_fdr Location biomart_name  MAPINFO
# <int> <chr>      <chr>                 <chr>    <chr>           <int>
#   1        1 cg00000029 M4                    TSS1500  RBL2         53468112
# 2        7 cg02851944 M4                    TSS1500  RBL2         53468125
# 3        5 cg02021817 M4                    TSS200   RBL2         53468150
# 4        9 cg03943270 M0.1                  TSS200   RBL2         53468205
# 5        8 cg02976952 M3                    5'UTR    RBL2         53468350
#  6        8 cg02976952 M3                    1stExon  RBL2         53468350
#  7       12 cg13426503 M0.1                  Body     RBL2         53468712
#  8       10 cg07396495 M0.1                  Body     RBL2         53468716
#  9        3 cg01016459 M0.1                  Body     RBL2         53468960
# 10       11 cg12190057 M3                    Body     RBL2         53469041
# 11        2 cg00011284 M0.1                  Body     RBL2         53469344
# 12        6 cg02119764 M0.1                  Body     RBL2         53469525
# 13        4 cg01783195 M4                    3'UTR    RBL2         53524545

## Same as the Location table without the accession column, de-duplicated.
## (Kept in memory for downstream use; not saved.)
df_comb_all_final <- df_comb_all_2 %>% select(-UCSC_RefGene_Accession) %>% distinct()

# ------------------------------------------------------------
# Step 8: Checks (per dataset) and save
# ------------------------------------------------------------
## (a) How many trios went in vs. came out, per dataset.
##     With require_gene_match = TRUE, "out" < "in" is expected: trios whose
##     gene is not annotated to their probe by Illumina are dropped.
cat("\nTrios in / out per dataset:\n")
print(full_join(count(all_annot, Dataset, name = "trios_in"),
                count(df_comb_all, Dataset, name = "trios_out"), by = "Dataset"))
cat(if (require_gene_match) "(gene matched to probe annotation)\n" else "(all trios kept)\n")

## (b) meth.row values that have no row in the dataset's TCGA.meth file.
##     These would show NA probe IDs; usually means meth_files points a
##     dataset to the wrong file. (With inner_join they are already dropped
##     in Step 5, so this mainly matters when require_gene_match = FALSE.)
no_probe <- df_comb_all %>% filter(is.na(Row.names)) %>% count(Dataset)
if (nrow(no_probe) > 0) {
  warning("Some meth.row values were not found in the dataset's methylation file")
  print(no_probe)
}

## (c) Sanity check on the file mapping: the position stored in TCGA.meth
##     should equal the manifest MAPINFO for the same probe. A large count
##     for one dataset means its meth.row values were matched to the wrong
##     methylation file.
cat("\nCoordinate mismatches (Genomic_Coordinate != MAPINFO) per dataset:\n")
print(df_comb_all %>% group_by(Dataset) %>%
        summarise(n_mismatch = sum(Genomic_Coordinate != MAPINFO, na.rm = TRUE)))

## (d) Each trio should appear once. Duplicates would mean a join key
##     wasn't unique (should not happen thanks to many-to-one).
if (anyDuplicated(df_comb_all[c("Dataset", "trio_row")]) > 0)
  warning("Some trios appear more than once in the output")

## (e) Missing values in the key numeric columns, per dataset:
##     diff_cpG_mapinfo NA = probe not near a CpG island (expected, open sea)
##     gene_length NA      = trio gene not found in biomart (see (f))
cat("\nMissing values per dataset:\n")
print(df_comb_all %>% group_by(Dataset) %>%
        summarise(across(c(MAPINFO, diff_cpG_mapinfo, diff_mapinfo_geneStart, gene_length),
                         ~ sum(is.na(.x)))))

## (f) Trio genes with no biomart row. For each one, look up its GRCh37 name
##     (old symbol) and add it to gene_alias in Step 0. (Patch-only genes are
##     now picked up automatically in Step 4.)
cat("\nTrio genes without a biomart match:\n")
print(df_comb_all %>% filter(is.na(gene_length)) %>% count(Dataset, Gene.name, sort = TRUE))

## (g) Validate patch-gene coordinates: the probe should be on the same
##     chromosome and inside the gene or within 5 kb of it (the probe was
##     matched to this gene by Illumina, so it should be in/near the gene).
##     near_gene FALSE -> patch coordinates don't fit; set them to NA or take
##     them from another GRCh37 source.
patch_check <- df_comb_all %>%
  filter(biomart_name %in% patch_genes) %>%
  mutate(near_gene = MAPINFO >= as.numeric(Gene.start..bp.) - 5000 &
           MAPINFO <= as.numeric(Gene.end..bp.)   + 5000) %>%
  select(Dataset, trio_row, Gene.name, CHR, MAPINFO,
         Chromosome.scaffold.name, Gene.start..bp., Gene.end..bp., near_gene)

cat("\nPatch-gene coordinate check:\n")
print(patch_check)
if (any(!patch_check$near_gene, na.rm = TRUE))
  warning("Some patch-gene coordinates are far from the probe -- check these genes")

## (h) Distribution of the position flags from Step 7, per dataset.
cat("\nposition_status per dataset:\n")
print(count(df_comb_all, Dataset, position_status))

## (i) Location table: how many trios were kept / split / removed, per dataset.
##       single_region   : one region, 1 row each
##       split_1acc      : multi-region, 1 accession, split into several rows
##       removed_gt1acc  : multi-region, >1 accession, dropped
cat("\nLocation table: trios kept / split / removed per dataset:\n")
print(df_comb_all %>%
        mutate(loc_rule = case_when(
          multi_refgene_group & n_accession > 1 ~ "removed_gt1acc",
          multi_refgene_group                   ~ "split_1acc",
          TRUE                                  ~ "single_region")) %>%
        count(Dataset, loc_rule) %>%
        pivot_wider(names_from = loc_rule, values_from = n, values_fill = 0))

## (j) The split trios should only be 5'UTR + 1stExon or 1stExon + 3'UTR.
##     Any other combination is printed here so you can decide on it.
cat("\nRegion combinations among split (1-accession) trios:\n")
print(df_comb_all %>%
        filter(multi_refgene_group, n_accession == 1) %>%
        mutate(region_combo = map_chr(str_split(UCSC_RefGene_Group, ";"),
                                      ~ paste(prio[prio %in% .x], collapse = " + "))) %>%
        count(Dataset, region_combo))

## (k) Rows per Location, per dataset (a split trio counts once in each region).
cat("\nRows per Location per dataset:\n")
print(df_comb_all_2 %>%
        count(Dataset, Location = factor(Location, levels = prio)) %>%
        pivot_wider(names_from = Dataset, values_from = n, values_fill = 0))

## Save both tables (all datasets together; filter by Dataset when reading).
saveRDS(df_comb_all, out_file)
cat("Saved:", out_file, "\n")

# saveRDS(df_comb_all_2, out_file_2)
# cat("Saved:", out_file_2, "\n")

saveRDS(df_comb_all_4, out_file_4)
cat("Saved:", out_file_4, "\n")

## Optional: also write one file per dataset.
## split() breaks the table into a list by Dataset; walk() saves each piece.
# walk(split(df_comb_all, df_comb_all$Dataset), function(d)
#   saveRDS(d, file.path(combined_dir, paste0("trio_probe_gene_annot_", d$Dataset[1], ".rds"))))
# ================================================================================================
