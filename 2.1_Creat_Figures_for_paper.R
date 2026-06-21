# ==============================================================================
# MRTrios Methylation Analysis
# ==============================================================================

# ── Libraries (loaded once) ───────────────────────────────────────────────────
library(data.table)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)

# ── Paths (defined once) ──────────────────────────────────────────────────────
DATA_DIR <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026"
RAW_DIR  <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/raw_Data_Methyl"
FIG_DIR  <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Figure_2026"

# ── Load model results (loaded once) ─────────────────────────────────────────
model_pos     <- readRDS(file.path(DATA_DIR, "LZ_trio_results_pos_ALL_qval.rds"))
model_neg     <- readRDS(file.path(DATA_DIR, "LZ_trio_results_neg_ALL_qval.rds"))
model_pos_loc <- readRDS(file.path(DATA_DIR, "LZ_trio_results_pos_ALL_qval_with_location.rds"))
model_neg_loc <- readRDS(file.path(DATA_DIR, "LZ_trio_results_neg_ALL_qval_with_location.rds"))

# ── Load extracted model sub-tables (loaded once) ────────────────────────────
# BUG FIX: original used M1.2_neg for pie_dat_M1.2_pos — corrected to M1.2_pos
model_names <- c("M0.1", "M0.2", "M1.1", "M1.2")
for (nm in model_names) {
  assign(paste0(nm, "_pos"), fread(file.path(DATA_DIR, paste0(nm, "_pos_extract.txt")), sep = "\t"))
  assign(paste0(nm, "_neg"), fread(file.path(DATA_DIR, paste0(nm, "_neg_extract.txt")), sep = "\t"))
}
# Objects now in environment: M0.1_pos, M0.1_neg, M0.2_pos, ... M1.2_neg

# ── Colours ───────────────────────────────────────────────────────────────────
COL_POS <- "#f781bf"
COL_NEG <- "#377eb8"

# ── Helper: histogram for a single model/variable ────────────────────────────
make_hist <- function(df, col, transform = NULL, fill, title,
                      xlab, bins = 35, binwidth = NULL,
                      xlim = NULL, ylim = NULL) {
  vals <- df[[col]]
  vals <- vals[!is.na(vals)]
  if (!is.null(transform)) vals <- transform(vals)
  p <- ggplot(data.frame(x = vals), aes(x = x)) +
    (if (is.null(binwidth))
      geom_histogram(aes(y = ..density..), bins = bins,
                     fill = fill, color = "black", alpha = 0.7)
     else
       geom_histogram(aes(y = ..density..), binwidth = binwidth,
                      fill = fill, color = "black", alpha = 0.7)) +
    labs(title = title, x = xlab, y = "Density") +
    theme_bw() +
    theme(plot.title = element_text(size = 10, hjust = 0.5, face = "plain"))
  if (!is.null(xlim)) p <- p + xlim(xlim)
  if (!is.null(ylim)) p <- p + ylim(ylim)
  p
}

# ── Helper: build 2×2 patchwork + combined ER+/ER- grid ──────────────────────
make_2x2 <- function(p_M0.1, p_M1.1, p_M0.2, p_M1.2, er_label) {
  (p_M0.1 | p_M1.1) / (p_M0.2 | p_M1.2) +
    plot_annotation(title = er_label) &
    theme(plot.title = element_text(hjust = 0.5, size = 12))
}

make_combined <- function(layout_pos, layout_neg) {
  plot_grid(layout_pos, layout_neg, labels = "AUTO", label_size = 15, ncol = 1)
}

# ── Helper: build all 8 histograms for one variable across model types ────────
build_plots <- function(col, transform = NULL, fill_pos = COL_POS, fill_neg = COL_NEG,
                        xlab, bins = 35, binwidth = NULL,
                        xlim = NULL, ylim = NULL,
                        pre_filter = NULL,   # optional function(df) → df
                        ylim_override = list()) {
  result <- list()
  for (er in c("pos", "neg")) {
    fill <- if (er == "pos") fill_pos else fill_neg
    for (nm in model_names) {
      key <- paste0(nm, "_", er)
      df  <- get(key)
      if (!is.null(pre_filter)) df <- pre_filter(df)
      yl <- if (key %in% names(ylim_override)) ylim_override[[key]] else ylim
      result[[key]] <- make_hist(df, col, transform = transform, fill = fill,
                                 title = nm, xlab = xlab, bins = bins,
                                 binwidth = binwidth, xlim = xlim, ylim = yl)
    }
  }
  result
}

# ==============================================================================
# Figure A: log10(distance) — CpG probe to nearby CpG island
# ==============================================================================

cpg_plots <- build_plots(
  col       = "diff_cpG_mapinfo",
  transform = function(x) log10(abs(x[x != 0])),
  pre_filter = function(df) df %>% filter(!is.na(diff_cpG_mapinfo), diff_cpG_mapinfo != 0),
  xlab = "log10 (distance)", bins = 40,
  xlim = c(-2, 4), ylim = c(0, 1)
)

layout_cpg_pos <- make_2x2(cpg_plots$M0.1_pos, cpg_plots$M1.1_pos,
                           cpg_plots$M0.2_pos, cpg_plots$M1.2_pos, "ER+")
layout_cpg_neg <- make_2x2(cpg_plots$M0.1_neg, cpg_plots$M1.1_neg,
                           cpg_plots$M0.2_neg, cpg_plots$M1.2_neg, "ER-")
make_combined(layout_cpg_pos, layout_cpg_neg)
# ggsave(file.path(FIG_DIR, "Figure6_hist_log10_distance_between_nearby_CpG_island_and_CpG_in_methy_probe.pdf"), width = 6, height = 8)

# ==============================================================================
# Figure B: log10(distance) — CpG probe to gene start position
# ==============================================================================

gs_plots <- build_plots(
  col       = "diff_mapinfo_geneStart",
  transform = function(x) log10(abs(x[x != 0])),
  pre_filter = function(df) df %>% filter(!is.na(diff_mapinfo_geneStart), diff_mapinfo_geneStart != 0),
  xlab = "log10 (distance)", bins = 35,
  xlim = c(-2.5, 7.5), ylim = c(0, 0.9)
)

layout_gs_pos <- make_2x2(gs_plots$M0.1_pos, gs_plots$M1.1_pos,
                          gs_plots$M0.2_pos, gs_plots$M1.2_pos, "ER+")
layout_gs_neg <- make_2x2(gs_plots$M0.1_neg, gs_plots$M1.1_neg,
                          gs_plots$M0.2_neg, gs_plots$M1.2_neg, "ER-")
make_combined(layout_gs_pos, layout_gs_neg)
ggsave(file.path(FIG_DIR, "Figure8_hist_log10_distance_between_CpG_in_methy_probe_and_startposition_of_a_gene.pdf"), width = 6, height = 8)

# ==============================================================================
# Figure C: log10(gene length)
# ==============================================================================

gl_plots <- build_plots(
  col       = "gene_length",
  transform = log10,
  pre_filter = function(df) df %>% select(gene_length) %>% distinct() %>%
    filter(!is.na(gene_length), gene_length > 0),
  xlab = "log10 (length)", bins = 30,
  xlim = c(1, 7), ylim = c(0, 0.8)
)

layout_gl_pos <- make_2x2(gl_plots$M0.1_pos, gl_plots$M1.1_pos,
                          gl_plots$M0.2_pos, gl_plots$M1.2_pos, "ER+")
layout_gl_neg <- make_2x2(gl_plots$M0.1_neg, gl_plots$M1.1_neg,
                          gl_plots$M0.2_neg, gl_plots$M1.2_neg, "ER-")
make_combined(layout_gl_pos, layout_gl_neg)
ggsave(file.path(FIG_DIR, "Figure10_log10_length_between_startposition_and_endposition_of_a_gene.pdf"), width = 6, height = 8)

# # ==============================================================================
# # Figure D: CpG island relation — pie charts
# # ==============================================================================
# 
# PIE_COLS   <- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f')
# 
# make_pie_data <- function(df) {
#   df %>%
#     count(Relation_to_UCSC_CpG_Island) %>%
#     mutate(percentage = n / sum(n) * 100,
#            Relation_to_UCSC_CpG_Island = if_else(
#              Relation_to_UCSC_CpG_Island == "" | is.na(Relation_to_UCSC_CpG_Island),
#              "No Info", Relation_to_UCSC_CpG_Island))
# }
# 
# draw_pie <- function(pie_df, title) {
#   pie(pie_df$n,
#       labels = paste0(pie_df$Relation_to_UCSC_CpG_Island,
#                       " (", round(pie_df$percentage, 1), "%)"),
#       col = PIE_COLS, main = title, cex = 0.75)
# }
# 
# pie_data <- lapply(setNames(
#   paste0(rep(model_names, each = 2), "_", c("pos", "neg")),
#   paste0(rep(model_names, each = 2), "_", c("pos", "neg"))
# ), function(key) make_pie_data(get(key)))
# 
# par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
# for (nm in model_names) draw_pie(pie_data[[paste0(nm, "_pos")]], paste(nm, "pos"))
# # Figure9.1: 5×7 landscape
# 
# par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
# for (nm in model_names) draw_pie(pie_data[[paste0(nm, "_neg")]], paste(nm, "neg"))
# # Figure9.2: 5×7 landscape

# ==============================================================================
# Figure D: CpG island relation — pie charts (ER+ and ER- combined, 2×4 grid)
# ==============================================================================

PIE_COLS <- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f')

make_pie_data <- function(df) {
  df %>%
    count(Relation_to_UCSC_CpG_Island) %>%
    mutate(percentage = n / sum(n) * 100,
           Relation_to_UCSC_CpG_Island = if_else(
             Relation_to_UCSC_CpG_Island == "" | is.na(Relation_to_UCSC_CpG_Island),
             "No Info", Relation_to_UCSC_CpG_Island))
}

draw_pie <- function(pie_df, title) {
  pie(pie_df$n,
      labels = paste0(pie_df$Relation_to_UCSC_CpG_Island,
                      " (", round(pie_df$percentage, 1), "%)"),
      col = PIE_COLS, main = title, cex = 1.1, radius = 0.65)
}

# draw_pie <- function(pie_df, title, cex = 1.1) {
#   pie(pie_df$n,
#       labels = paste0(pie_df$Relation_to_UCSC_CpG_Island,
#                       " (", round(pie_df$percentage, 1), "%)"),
#       col = PIE_COLS, main = title, cex = cex, radius = 0.65)
# }

pie_data <- lapply(setNames(
  paste0(rep(model_names, each = 2), "_", c("pos", "neg")),
  paste0(rep(model_names, each = 2), "_", c("pos", "neg"))
), function(key) make_pie_data(get(key)))

# Combined figure: ER+ top row, ER- bottom row
# Save as Figure9_pie_CpG_island.pdf, 10×7 landscape
par(mfrow = c(2, 4), mar = c(2, 2, 2, 2))
for (nm in model_names) draw_pie(pie_data[[paste0(nm, "_pos")]], paste(nm, "ER+"))
for (nm in model_names) draw_pie(pie_data[[paste0(nm, "_neg")]], paste(nm, "ER-"))

ggsave(file.path(FIG_DIR, "Figure9_pie_Distribution_relation_CpG_island_five_sections.pdf"), width = 10, height = 7)

# ==============================================================================
# Figure E: Mean methylation distribution
# ==============================================================================

mm_plots <- build_plots(
  col  = "Methyl_mean",
  xlab = "Mean methylation", bins = 30,
  xlim = c(-5, 5), ylim = c(0, 0.45)
)

layout_mm_pos <- make_2x2(mm_plots$M0.1_pos, mm_plots$M1.1_pos,
                          mm_plots$M0.2_pos, mm_plots$M1.2_pos, "ER+")
layout_mm_neg <- make_2x2(mm_plots$M0.1_neg, mm_plots$M1.1_neg,
                          mm_plots$M0.2_neg, mm_plots$M1.2_neg, "ER-")
make_combined(layout_mm_pos, layout_mm_neg)
ggsave(file.path(FIG_DIR, "Figure11_hist_meanMethylation_distribution.pdf"), width = 6, height = 8)

# ==============================================================================
# Figure F: GC content distribution
# NOTE: M0.2_neg and M1.2_neg have a wider ylim (0.25) than the others (0.1)
# ==============================================================================

gc_plots <- build_plots(
  col      = "Gene...GC.content",
  xlab     = "GC content", binwidth = 2.5,
  xlim     = c(28, 75), ylim = c(0, 0.1),
  ylim_override = list(M0.2_neg = c(0, 0.25), M1.2_neg = c(0, 0.25))
)

layout_gc_pos <- make_2x2(gc_plots$M0.1_pos, gc_plots$M1.1_pos,
                          gc_plots$M0.2_pos, gc_plots$M1.2_pos, "ER+")
layout_gc_neg <- make_2x2(gc_plots$M0.1_neg, gc_plots$M1.1_neg,
                          gc_plots$M0.2_neg, gc_plots$M1.2_neg, "ER-")
make_combined(layout_gc_pos, layout_gc_neg)
ggsave(file.path(FIG_DIR, "Figure13_GC_content_distribution.pdf"), width = 6, height = 8)

# ==============================================================================
# Summary statistics
# ==============================================================================

log10_summary <- function(df, col, min_val = -2) {
  x <- df[[col]]
  x <- x[!is.na(x) & x > 0]
  summary(log10(x)[log10(x) > min_val])
}

cat("\n── log10(|diff_cpG_mapinfo|) ──\n")
for (key in paste0(rep(model_names, each=2), "_", c("pos","neg")))
{ cat(key, ":\n"); print(log10_summary(get(key), "diff_cpG_mapinfo")); cat("\n") }

cat("\n── log10(|diff_mapinfo_geneStart|) ──\n")
for (key in paste0(rep(model_names, each=2), "_", c("pos","neg")))
{ cat(key, ":\n"); print(log10_summary(get(key), "diff_mapinfo_geneStart")); cat("\n") }

cat("\n── log10(gene_length) ──\n")
for (key in paste0(rep(model_names, each=2), "_", c("pos","neg")))
{ cat(key, ":\n"); print(log10_summary(get(key), "gene_length")); cat("\n") }

cat("\n── Methyl_mean ──\n")
for (key in paste0(rep(model_names, each=2), "_", c("pos","neg")))
{ cat(key, ":\n"); print(summary(get(key)$Methyl_mean)); cat("\n") }

cat("\n── GC content ──\n")
for (key in paste0(rep(model_names, each=2), "_", c("pos","neg")))
{ cat(key, ":\n"); print(summary(get(key)$Gene...GC.content)); cat("\n") }

# ==============================================================================
# Location bar plots (UCSC_RefGene_Group)
# ==============================================================================

LOC_COLS <- c("TSS" = "#377eb8", "Body" = "#4daf4a", "5'/3'UTR" = "#e7298a")

recode_location <- function(df) {
  df %>%
    separate_rows(UCSC_RefGene_Group, sep = ";") %>%
    mutate(Group = case_when(
      UCSC_RefGene_Group %in% c("TSS200", "TSS1500") ~ "TSS",
      UCSC_RefGene_Group %in% c("Body",   "1stExon") ~ "Body",
      UCSC_RefGene_Group %in% c("3'UTR",  "5'UTR")   ~ "5'/3'UTR",
      TRUE ~ NA_character_
    ))
}

# Build location table for TCGA probe → gene mapping (run once, then save)
humanmeth <- read.csv(file.path(RAW_DIR, "GPL13534_HumanMethylation450_15017482_v.1.1 2.csv"),
                      skip = 7, header = TRUE)
trios      <- fread(file.path(RAW_DIR, "trio.final.protein.coding.txt"), data.table = FALSE)
TCGA.meth  <- readRDS(file.path(RAW_DIR, "split.names.TCGA.meth.logit.rds"))

meth_trio      <- TCGA.meth[trios$meth.row, ]
trios_add_Name <- cbind(trios, Name = meth_trio$Row.names)
trios_add_Name$trios_row <- seq_len(nrow(trios_add_Name))

info <- humanmeth %>% select(Name, UCSC_RefGene_Group) %>%
  filter(Name %in% trios_add_Name$Name)
merge_trio_loc <- merge(trios_add_Name, info, by = "Name", all = FALSE) %>%
  arrange(trios_row)

write.table(merge_trio_loc,
            file.path(DATA_DIR, "final_protein_coding_location_LZ_03_19_25.csv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

loc <- fread(file.path(DATA_DIR, "final_protein_coding_location_LZ_03_19_25.csv"),
             data.table = FALSE)

make_loc_barplot <- function(model_file, loc, er_label, y_max) {
  df <- fread(file.path(DATA_DIR, model_file), header = TRUE) %>%
    left_join(loc, by = c("trio_row" = "trios_row")) %>%
    recode_location()
  counts <- df %>% group_by(Inferred.Model3, Group) %>%
    summarise(Count = n(), .groups = "drop")
  ggplot(counts, aes(x = Inferred.Model3, y = Count, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_y_continuous(limits = c(0, y_max),
                       breaks = seq(0, y_max, 50000),
                       labels = seq(0, y_max / 10000, 5)) +
    scale_fill_manual(name = "Location", values = LOC_COLS) +
    labs(title = paste0("(", er_label, ") Location by Model Type"),
         x = "Model Type", y = "Count (thousand)") +
    theme_bw() +
    theme(legend.position = "top", plot.title = element_text(size = 14, hjust = 0.5))
}

make_loc_barplot("LZ_trio_results_pos_ALL_qval.txt", loc, "ER+", 180000)
make_loc_barplot("LZ_trio_results_neg_ALL_qval.txt", loc, "ER-", 150000)

# ==============================================================================
# Load raw data for correlation analysis (loaded once)
# ==============================================================================

cna      <- fread(file.path(RAW_DIR, "data_CNA.txt"), data.table = FALSE)
gene.exp <- fread(file.path(RAW_DIR, "data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt"),
                  data.table = FALSE)
TCGA.meth <- readRDS(file.path(RAW_DIR, "split.names.TCGA.meth.logit.rds"))
trios     <- fread(file.path(RAW_DIR, "trio.final.protein.coding.txt"))
biomart   <- read.delim(file.path(RAW_DIR, "ensembl37_genes_p13_biomart.txt"), header = TRUE)

id_pos_patient <- fread(file.path(RAW_DIR, "names.pos.patient2.txt"), header = FALSE)$V1
id_neg_patient <- fread(file.path(RAW_DIR, "names.neg.patient2.txt"), header = FALSE)$V1

# Common patient IDs across all three assays
E_M_C_Id <- colnames(gene.exp)[3:ncol(gene.exp)] %>%
  intersect(colnames(TCGA.meth)[5:ncol(TCGA.meth)]) %>%
  intersect(colnames(cna)[3:ncol(cna)])

E_M_C_id_pos <- intersect(E_M_C_Id, id_pos_patient)
E_M_C_id_neg <- intersect(E_M_C_Id, id_neg_patient)

# Subset matrices to common patient columns, indexed by trios
# BUG FIX: original code duplicated this entire block verbatim with slightly
# different variable names (E_M_C_id_pos_patient vs E_M_C_patient_pos_Id etc.)
meth_join <- TCGA.meth %>% left_join(humanmeth, by = c("Row.names" = "Name"))

cna_in_trios      <- cna      %>% mutate(cna.row  = as.numeric(rownames(.))) %>%
  filter(cna.row  %in% unique(trios$cna.row))
gene.exp_in_trios <- gene.exp %>% mutate(gene.row = as.numeric(rownames(.))) %>%
  filter(gene.row %in% unique(trios$gene.row))
meth_in_trios     <- meth_join %>% mutate(meth.row = as.numeric(rownames(.))) %>%
  filter(meth.row %in% unique(trios$meth.row))

pos_cna  <- cna_in_trios      %>% select(cna.row,  all_of(E_M_C_id_pos))
pos_meth <- meth_in_trios     %>% select(meth.row, all_of(E_M_C_id_pos))
pos_exp  <- gene.exp_in_trios %>% select(gene.row, all_of(E_M_C_id_pos))

neg_cna  <- cna_in_trios      %>% select(cna.row,  all_of(E_M_C_id_neg))
neg_meth <- meth_in_trios     %>% select(meth.row, all_of(E_M_C_id_neg))
neg_exp  <- gene.exp_in_trios %>% select(gene.row, all_of(E_M_C_id_neg))

saveRDS(pos_cna,  file.path(DATA_DIR, "pos_cna.rds"))
saveRDS(pos_meth, file.path(DATA_DIR, "pos_meth.rds"))
saveRDS(pos_exp,  file.path(DATA_DIR, "pos_exp.rds"))
saveRDS(neg_cna,  file.path(DATA_DIR, "neg_cna.rds"))
saveRDS(neg_meth, file.path(DATA_DIR, "neg_meth.rds"))
saveRDS(neg_exp,  file.path(DATA_DIR, "neg_exp.rds"))

pos_cna  <- readRDS(file.path(DATA_DIR, "pos_cna.rds"))
pos_meth <- readRDS(file.path(DATA_DIR, "pos_meth.rds"))
pos_exp <- readRDS(file.path(DATA_DIR, "pos_exp.rds"))
neg_cna <- readRDS(file.path(DATA_DIR, "neg_cna.rds"))
neg_meth <- readRDS(file.path(DATA_DIR, "neg_meth.rds"))
neg_exp <- readRDS(file.path(DATA_DIR, "neg_exp.rds"))

# ==============================================================================
# Process ER+ and ER- trios, compute E-M correlations
# ==============================================================================

process_er_group <- function(model, trios, TCGA.meth, gene.exp,
                             id_patients, E_M_C_id,
                             meth_n_meta = 5, gene_n_meta = 3) {
  idx_M <- match(E_M_C_id, colnames(TCGA.meth))
  idx_E <- match(E_M_C_id, colnames(gene.exp))
  n_pat <- length(E_M_C_id)
  
  sub_trios <- trios %>% mutate(trios_row = as.integer(rownames(trios))) %>%
    slice(model$trio_row)
  
  meth_sub <- TCGA.meth %>%
    mutate(Methy_row = as.integer(rownames(TCGA.meth))) %>%
    slice(sub_trios$meth.row) %>%
    select(c(1:4, 900), all_of(idx_M))
  meth_sub$Mean_Methy <- rowMeans(meth_sub[, 6:ncol(meth_sub)], na.rm = TRUE)
  meth_sub <- meth_sub %>% select(1:5, n_pat + 6, 6:(n_pat + 5))
  
  gene_sub <- gene.exp %>%
    mutate(Gene_row = as.integer(rownames(gene.exp))) %>%
    slice(sub_trios$gene.row) %>%
    select(c(1:2, 1103), all_of(idx_E))
  gene_sub$Mean_Gene <- rowMeans(gene_sub[, 4:ncol(gene_sub)], na.rm = TRUE)
  gene_sub <- gene_sub %>% select(1:3, n_pat + 4, 4:(n_pat + 3))
  
  final <- cbind(gene_sub[, 1:4], meth_sub[, 1:6], model)
  
  # Compute E-M correlation and p-value per trio row
  gene_vals <- gene_sub[, 5:(n_pat + 4)]
  meth_vals <- meth_sub[, 7:(n_pat + 6)]
  corr_df <- map_dfr(seq_len(nrow(gene_vals)), function(i) {
    res <- cor.test(as.numeric(gene_vals[i, ]), as.numeric(meth_vals[i, ]))
    data.frame(trio_row = model$trio_row[i],
               corr     = res$estimate,
               p_val    = res$p.value)
  })
  
  final %>% left_join(corr_df, by = "trio_row")
}

final_model_pos <- process_er_group(model_pos, trios, TCGA.meth, gene.exp,
                                    id_pos_patient, E_M_C_id_pos)
final_model_neg <- process_er_group(model_neg, trios, TCGA.meth, gene.exp,
                                    id_neg_patient, E_M_C_id_neg)

# Rename correlation columns to match downstream usage
final_model_pos <- final_model_pos %>% rename(corr_pos = corr, p_val_pos = p_val)
final_model_neg <- final_model_neg %>% rename(corr_neg = corr, p_val_neg = p_val)

colnames(final_model_pos)[1] <- "Gene.name"
colnames(final_model_neg)[1] <- "Gene.name"

saveRDS(final_model_pos, file.path(DATA_DIR, "final_model_pos_with_correlation.rds"))
saveRDS(final_model_neg, file.path(DATA_DIR, "final_model_neg_with_correlation.rds"))


final_model_pos <- readRDS(file.path(DATA_DIR, "final_model_pos_with_correlation.rds"))
final_model_neg <- readRDS(file.path(DATA_DIR, "final_model_neg_with_correlation.rds"))

final_model_combined_pos <- readRDS(file.path(DATA_DIR, "final_model_combined_pos.rds"))
final_model_combined_neg <- readRDS(file.path(DATA_DIR, "final_model_combined_neg.rds"))


# E-M correlation boxplots
for (er in list(list(df = final_model_pos, col = "corr_pos", label = "ER+"),
                list(df = final_model_neg, col = "corr_neg", label = "ER-"))) {
  print(
    ggplot(er$df, aes(x = Inferred.Model.qval, y = .data[[er$col]],
                      fill = Inferred.Model.qval)) +
      geom_boxplot() +
      labs(title = paste0(er$label, ": E-M correlation by model type"),
           x = "Inferred model types", y = "Correlation") +
      theme_bw() + theme(legend.position = "none")
  )
}

# ==============================================================================
# Figures K/L: Mean gene expression and mean methylation histograms
# (from final_model_pos / final_model_neg, filtered by model type)
# ==============================================================================

make_final_model_plots <- function(final_df, col, fill, xlim, ylim = NULL, bins = 40) {
  lapply(setNames(model_names, model_names), function(nm) {
    sub <- final_df %>%
      filter(Inferred.Model.qval == nm) %>%
      select(Inferred.Model.qval, Gene_Symbol, all_of(col)) %>%
      distinct()
    make_hist(sub, col, fill = fill, title = nm,
              xlab = col, bins = bins, xlim = xlim, ylim = ylim)
  })
}

# Mean gene expression
ge_pos_plots <- make_final_model_plots(final_model_pos, "Mean_Gene", COL_POS, c(-2, 2))
ge_neg_plots <- make_final_model_plots(final_model_neg, "Mean_Gene", COL_NEG, c(-2, 2), c(0, 2))

make_combined(
  make_2x2(ge_pos_plots$M0.1, ge_pos_plots$M1.1, ge_pos_plots$M0.2, ge_pos_plots$M1.2, "ER+"),
  make_2x2(ge_neg_plots$M0.1, ge_neg_plots$M1.1, ge_neg_plots$M0.2, ge_neg_plots$M1.2, "ER-")
)
# ggsave(file.path(FIG_DIR, "hist_GeneExpMean.pdf"), width = 6, height = 8)

# Mean methylation
mm_pos_plots <- make_final_model_plots(final_model_pos, "Mean_Methy", COL_POS, c(-5, 5), c(0, 0.4))
mm_neg_plots <- make_final_model_plots(final_model_neg, "Mean_Methy", COL_NEG, c(-5, 5), c(0, 0.4))

make_combined(
  make_2x2(mm_pos_plots$M0.1, mm_pos_plots$M1.1, mm_pos_plots$M0.2, mm_pos_plots$M1.2, "ER+"),
  make_2x2(mm_neg_plots$M0.1, mm_neg_plots$M1.1, mm_neg_plots$M0.2, mm_neg_plots$M1.2, "ER-")
)
# ggsave(file.path(FIG_DIR, "hist_MethylMean.pdf"), width = 6, height = 8)

# ==============================================================================
# M1, M2, M4 correlation boxplots by genomic location
# ==============================================================================

combined_pos <- inner_join(final_model_pos, model_pos_loc %>%
                             select(trio_row, meth.row, cna.row, gene.row, Inferred.Model3, Name,
                                    UCSC_RefGene_Group, Group),
                           by = c("trio_row", "meth.row", "cna.row", "gene.row"))

combined_neg <- inner_join(final_model_neg, model_neg_loc %>%
                             select(trio_row, meth.row, cna.row, gene.row, Inferred.Model3, Name,
                                    UCSC_RefGene_Group, Group),
                           by = c("trio_row", "meth.row", "cna.row", "gene.row"))

saveRDS(combined_pos, file.path(DATA_DIR, "final_model_combined_pos.rds"))
saveRDS(combined_neg, file.path(DATA_DIR, "final_model_combined_neg.rds"))

LOC_MODELS <- c("M1.1", "M1.2", "M2.1", "M2.2", "M4")

make_loc_corr_boxplots <- function(combined, corr_col, fill, er_label) {
  plots <- lapply(setNames(LOC_MODELS, LOC_MODELS), function(m) {
    ggplot(combined %>% filter(Inferred.Model3 == m),
           aes(x = Group, y = .data[[corr_col]])) +
      geom_boxplot(fill = fill, color = "black") +
      labs(title = m, x = "Location", y = "Correlation") +
      ylim(-0.9, 0.9) + theme_bw()
  })
  wrap_plots(plots, nrow = 1) +
    plot_annotation(title = er_label) &
    theme(plot.title = element_text(hjust = 0.5, size = 12))
}

make_loc_corr_boxplots(combined_pos, "corr_pos", COL_POS, "ER+")
# ggsave(file.path(FIG_DIR, "pos_M1_M2_M4_correlation_boxplot.pdf"), width = 10, height = 3.5)

make_loc_corr_boxplots(combined_neg, "corr_neg", COL_NEG, "ER-")
# ggsave(file.path(FIG_DIR, "neg_M1_M2_M4_correlation_boxplot.pdf"), width = 10, height = 3.5)

# ==============================================================================
# Trio-level pairwise correlations (M~E, C~E, M~C) for M1.1, M1.2, M3
# BUG FIX: original cor_results_posM3 used pos_M1.1_final row indices — fixed.
# ==============================================================================

# Pre-index matrices for fast lookup (avoids filter() inside loop)
to_mat <- function(df, id_col) {
  m <- as.matrix(df[, -which(names(df) == id_col)])
  rownames(m) <- as.character(df[[id_col]])
  m
}
pm <- to_mat(pos_meth, "meth.row")
pc <- to_mat(pos_cna,  "cna.row")
pe <- to_mat(pos_exp,  "gene.row")

compute_trio_cors <- function(model_sub, pm, pc, pe) {
  t(vapply(seq_len(nrow(model_sub)), function(i) {
    m <- pm[as.character(model_sub$meth.row[i]), ]
    c <- pc[as.character(model_sub$cna.row[i]),  ]
    e <- pe[as.character(model_sub$gene.row[i]), ]
    c(cor_meth_exp = cor(m, e, use = "complete.obs"),
      cor_cna_exp  = cor(c, e, use = "complete.obs"),
      cor_meth_cna = cor(m, c, use = "complete.obs"))
  }, numeric(3)))
}

for (mt in c("M1.1", "M1.2", "M3")) {
  sub    <- final_model_pos %>% filter(Inferred.Model.qval == mt) %>% arrange(corr_pos)
  cors   <- compute_trio_cors(sub, pm, pc, pe)
  result <- cbind(sub, cors)
  fwrite(result, file.path(DATA_DIR, paste0("pos_", mt, "_all.txt")), sep = "\t")
  assign(paste0("pos_", gsub("\\.", "", mt), "_all"), result)
}

# M1.1_all   <- fread(file.path(DATA_DIR, "pos_M1.1_all.txt"))

loc_cols <- model_pos_loc %>% 
  select(trio_row, UCSC_RefGene_Group, Group)

for (mt in c("M1.1", "M1.2", "M3")) {
  fname <- paste0("pos_", gsub("\\.", "", mt), "_all")   # pos_M11_all, pos_M12_all, pos_M3_all
  df    <- get(fname)
  
  df_with_loc <- df %>%
    left_join(loc_cols, by = "trio_row")
  
  assign(fname, df_with_loc)
  fwrite(df_with_loc, file.path(DATA_DIR, paste0("pos_", mt, "_all.txt")), sep = "\t")
}

# Example filtering on the saved results
a1 <- pos_M11_all %>%
  filter(abs(cor_meth_exp) > 0.5 & abs(cor_cna_exp) > 0.5) %>%
  arrange(desc(abs(cor_meth_exp)), desc(abs(cor_cna_exp)))

a2 <- pos_M12_all %>%
  filter(abs(cor_meth_exp) > 0.5 & abs(cor_meth_cna) > 0.3) %>%
  arrange(desc(abs(cor_meth_exp)), desc(abs(cor_meth_cna)))

a_M3 <- pos_M3_all %>%
  filter(abs(cor_meth_exp) > 0.5 & abs(cor_meth_cna) > 0.3)



# ==============================================================================
# Find M1.1 trios where cor(CNA, EXP) < cor(METH, EXP)
# ==============================================================================

# Step 1: Pull M1.1 trios from pos_M11_all (already has all three correlations)
# m1.1_interest <- pos_M11_all %>%
#   filter(abs(cor_cna_exp) < abs(cor_meth_exp)) %>%
#   arrange(desc(abs(cor_meth_exp) - abs(cor_cna_exp)))

m1.1_interest <- pos_M11_all %>%
  filter(abs(cor_cna_exp) < abs(cor_meth_exp)) %>%
  arrange(desc(abs(cor_meth_exp) - abs(cor_cna_exp))) %>% select(-UCSC_RefGene_Group) %>% distinct() # biggest gap first

cat("Number of M1.1 trios where |cor(CNA,EXP)| < |cor(METH,EXP)|:", nrow(m1.1_interest), "\n")

# ==============================================================================
# Correlation scatter plots for a selected trio
# (change i to look at different trios)
# ==============================================================================

plot_trio_cors <- function(i, model_sub, pm, pc, pe) {
  m        <- pm[as.character(model_sub$meth.row[i]), ]
  cna_vals <- pc[as.character(model_sub$cna.row[i]),  ]
  e        <- pe[as.character(model_sub$gene.row[i]), ]
  
  gene_name <- model_sub$Hugo_Symbol[i]
  loc       <- model_sub$Group[i]
  
  r_me <- round(cor(m,        e, use = "complete.obs"), 3)
  r_ce <- round(cor(cna_vals, e, use = "complete.obs"), 3)
  r_mc <- round(cor(m, cna_vals, use = "complete.obs"), 3)
  
  df <- data.frame(meth = m, cna = cna_vals, exp = e)
  
  p1 <- ggplot(df, aes(x = meth, y = exp)) +
    geom_point(alpha = 0.5, color = COL_POS) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    labs(title = paste0("METH vs EXP  (r = ", r_me, ")"),
         x = "Methylation", y = "Gene expression") +
    theme_bw()
  
  p2 <- ggplot(df, aes(x = factor(cna), y = exp)) +
    geom_boxplot(fill = COL_NEG, color = "black", alpha = 0.7, outlier.size = 0.8) +
    labs(title = paste0("CNA vs EXP  (r = ", r_ce, ")"),
         x = "CNA", y = "Gene expression") +
    theme_bw()
  
  p3 <- ggplot(df, aes(x = factor(cna), y = meth)) +
    geom_boxplot(fill = "#4daf4a", color = "black", alpha = 0.7, outlier.size = 0.8) +
    labs(title = paste0("METH vs CNA  (r = ", r_mc, ")"),
         x = "CNA", y = "Methylation") +
    theme_bw()
  
  wrap_plots(p1, p2, p3, nrow = 1) +
    plot_annotation(
      title    = paste0("Gene: ", gene_name, " | Location: ", loc, " | Model: M1.1"),
      subtitle = paste0("trio_row = ", model_sub$trio_row[i])
    ) & theme(plot.title    = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5))
}

# Plot the top trio (largest gap between cor_meth_exp and cor_cna_exp)
plot_trio_cors(1, m1.1_interest, pm, pc, pe)
# ggsave(file.path(FIG_DIR, "M1.1_trio_cor_top1.pdf"), width = 10, height = 4)

# Plot multiple trios at once (e.g. top 6)
library(cowplot)
top_plots <- lapply(1:6, function(i) plot_trio_cors(i, m1.1_interest, pm, pc, pe))
plot_grid(plotlist = top_plots, ncol = 1)

top_plots <- lapply(1:20, function(i) plot_trio_cors(i, m1.1_interest, pm, pc, pe))
plot_grid(plotlist = top_plots, ncol = 1)

# ggsave(file.path(FIG_DIR, "M1.1_trio_cor_top6.pdf"), width = 10, height = 24)
## 6interset_trios_pos_M1.1.pdf ,width = 18, height = 10
