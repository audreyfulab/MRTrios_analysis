# ===============================================================================================
## Figures:
# 1. Correlation landscape by genomic location across three cancers.
# - boxplots of correlations between E and M at difference locations: TSS1500 → TSS200 → 5'UTR → 1st exon → Body → 3'UTR.
# 2. Causal-model distribution by genomic location across three cancers.
# 3. E–M correlation does not determine causal direction, BRCA high-confidence trios.
# 4. Direction, sign, and magnitude of causal effects, using the regression coefficients.
# - for M1 and M2 trios:
#   - For M->E trios, we know that the M-E correlations can be positive and negative, 
#     but those correlations are marginal ones. Need to look at the distribution of the effect sizes, 
#     which account for confounding variables. Does that depend on location: TSS vs gene body vs UTR?
#     The classical model is M represses E in the promoter region.
# - For E->M trios, marginal cors are also positive and negative. Need to look at the effect sizes. 
# - Does effect magnitude differ by genomic location or causal direction?
#   5. Biological characterization of mediation trios, restricted to the most informative features.

# ===============================================================================================

##############################################################
# Figure 1. Correlation landscape by genomic location
# Figure 2. Causal-model distribution by genomic location
#
# Input: the Location table made by trio_probe_gene_annot_all_datasets.R
#        (out_file_2 = trio_probe_gene_annot_loc_split_ALL.rds)
#   - one row per trio x Location
#   - multi-region trios with n_refgene_group >1 & >1 accession were already removed there
#     (different transcripts disagree on the region -> ambiguous)
#   - n_accession >1 & n_refgene_group =1  keep as the raw rows (NM_022763;NM_001135095  TSS1500)
#   - 5'UTR;1stExon and 1stExon;3'UTR trios (1 accession) were split into
#     one row per region, so they are counted in BOTH locations
#   - The locations are in gene order (TSS1500 -> TSS200 -> 5'UTR -> 1st exon -> Body -> 3'UTR)
#
# Output (written to MRTrios_combined/figures):
#   Fig1_correlation_by_location_6regions.pdf
#   Fig2_model_distribution_by_location_6regions.pdf
#   TableS1_n_pairs_by_cancer_model_location_6regions.csv
#   + statistics printed to the console
#
# Overall flow:
#   Location table
#     --> keep trios with a correlation and one of the 6 locations
#     --> relabel cancers / locations, fix their order (factors)
#     --> Figure 1: E-M correlation per location (A boxplots, B medians,
#                   C boxplots split by inferred model)
#     --> Figure 2: share of each inferred model per location
#                   (A stacked bars, B dodged bars)
#     --> statistics for the results text
##############################################################







## Packages:
##   data.table : fast grouping (dat[, ..., by = ...]), dcast, CJ, fwrite
##   ggplot2    : all plots
##   patchwork  : combine several ggplots into one figure with | and /
library(data.table)
library(ggplot2)
library(patchwork)

## Read the Location table. It was saved as a tibble by the annotation
## script; setDT() turns it into a data.table in place (no copy), so the
## data.table syntax below works.
combined_dir <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_combined"
loc_df <- readRDS(file.path(combined_dir, "trio_probe_gene_annot_loc_split_allCancerALL.rds"))
setDT(loc_df)

## All figures and tables go to MRTrios_combined/figures.
## dir.create() makes the folder if needed; showWarnings = FALSE keeps it
## quiet when the folder already exists.
fig_dir <- file.path(combined_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)
setwd(fig_dir)

## Names of the two analysis columns. They come from all_annot (the trio
## model results) and were carried through the annotation script:
##   cor_col   : per-trio Pearson correlation between expression (E) and
##               methylation (M) across samples
##   model_col : inferred causal model per trio (M0.1, M1.1, ..., Other)
## If your columns are named differently, change them here only; the rest
## of the script uses these variables.
cor_col   <- "cor_exp_meth"
model_col <- "Inferred.Model.BH_fdr"

## Stop early with a clear error if any needed column is missing,
## instead of failing halfway through with a confusing message.
stopifnot(all(c("Dataset", "Location", "IlmnID", "Gene.name",
                cor_col, model_col) %in% names(loc_df)))

# ======================= 1. Settings =======================
## Order of the cancer panels / boxes. These are the DISPLAY names; the
## raw names (negER_BRCA, posER_BRCA) are converted in section 2.
## "\u2212" is a true minus sign, which looks better than a hyphen in "ER-".
cancer_levels <- c("BLCA", "LIHC", "BRCA ER\u2212", "BRCA ER+")

## The six Illumina gene regions, in gene order (upstream -> downstream):
##   TSS1500 : 200-1500 bp upstream of the transcription start site
##   TSS200  : 0-200 bp upstream of the TSS
##   5'UTR   : 5' untranslated region
##   1stExon : first exon
##   Body    : rest of the gene body (introns and later exons)
##   3'UTR   : 3' untranslated region
## loc_levels  = the values as they appear in the Location column
## loc_labels  = the names shown on the x-axis (same order)
## Ordering by gene structure makes the medians in Fig 1B read as a trend
## along the gene.
loc_levels <- c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR")
loc_labels <- c("TSS1500", "TSS200", "5'UTR", "1st exon", "Gene body", "3'UTR")

## Colors.
##   cancer_cols : same as the earlier Figure 1, so cancers match across figures
##   loc_cols    : one colour per location (only used in Fig 2B); promoter
##                 regions in orange/yellow, UTR/exon in greens, body in blue,
##                 3'UTR in pink. setNames() attaches the location labels as
##                 names, which scale_fill_manual uses to match colours.
##   model_cols  : same palette as before; each model family shares a hue
##                 (M0 blues, M1 greens, M2 reds, M3/M4 oranges, Other grey)
## Every name must match the data exactly, otherwise ggplot shows grey.
cancer_cols <- c("BLCA" = "#66C2A5", "LIHC" = "#FC8D62",
                 "BRCA ER\u2212" = "#8DA0CB", "BRCA ER+" = "#E78AC3")
loc_cols    <- setNames(c("#E69F00", "#F0C35A", "#009E73", "#7FCDBB", "#56B4E9", "#CC79A7"),
                        loc_labels)
model_cols  <- c(M0.1 = "#A6CEE3", M0.2 = "#1F78B4", M1.1 = "#B2DF8A", M1.2 = "#33A02C",
                 M2.1 = "#FB9A99", M2.2 = "#E31A1C", M3 = "#FDBF6F", M4 = "#FF7F00",
                 Other = "grey75")

## Small tuning values:
##   small_cut : in Fig 1C, boxes based on fewer CpG-gene pairs than this
##               get a "*" (sample-size warning, NOT a significance mark)
##   dodge_w   : dodge width in Fig 1C; the boxes and the "*" must use the
##               same width, or the stars drift away from their boxes
##   label_cut : in Fig 2A, % labels are printed only inside segments that
##               are at least this large (8%), so small slivers stay clean
small_cut <- 1000
dodge_w   <- 0.85
label_cut <- 0.08

## Shared theme for every panel: white background, light-grey facet strips
## with bold titles, no minor grid and no vertical grid lines (the x-axis
## is categorical, so vertical lines add nothing), black x-axis text.
base_theme <- theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(colour = "black"))

# ======================= 2. Data =======================
## Build the plotting table:
##   - keep rows with a correlation value (NA = correlation could not be
##     computed, e.g. constant expression) and one of the six locations
##   - keep only the columns the figures need
##   - get() looks up the column whose NAME is stored in cor_col / model_col,
##     and the result is renamed to cor_exp_meth / Model
dat <- loc_df[!is.na(get(cor_col)) & Location %in% loc_levels,
              .(Dataset, Gene.name, IlmnID, Location,
                cor_exp_meth = get(cor_col),
                Model        = get(model_col))]

## Convert raw dataset names to display names. fcase() is data.table's
## vectorised if/else; default = Dataset keeps BLCA and LIHC as they are.
## Then drop anything not in cancer_levels (in case other datasets exist).
dat[, Dataset := fcase(Dataset == "negER_BRCA", "BRCA ER\u2212",
                       Dataset == "posER_BRCA", "BRCA ER+",
                       default = Dataset)]
dat <- dat[Dataset %in% cancer_levels]

## Turn Dataset and Location into factors. Factor levels set the ORDER of
## panels, boxes and legend entries (otherwise ggplot sorts alphabetically,
## e.g. 1stExon before TSS1500). For Location, levels = the raw values and
## labels = the display names, so "Body" is shown as "Gene body".
dat[, Dataset  := factor(Dataset, levels = cancer_levels)]
dat[, Location := factor(Location, levels = loc_levels, labels = loc_labels)]

## Model table for Fig 1C and Fig 2:
##   - model_levels = models in the palette that actually occur in the data,
##     kept in palette order (intersect keeps the order of its first argument)
##   - dat_m keeps only trios with a recognised model call (drops NA / other
##     labels not in model_cols)
##   - Model becomes a factor so the legend follows M0.1, M0.2, ..., Other
model_levels <- intersect(names(model_cols), unique(dat$Model))
dat_m <- dat[Model %in% model_levels]
dat_m[, Model := factor(Model, levels = model_levels)]

## Checks before plotting:
##   - every cancer x location cell should have a count and there should be
##     no NA column (an NA column means a label did not match loc_levels)
##   - the number of trios per model; a missing model here means its name
##     in the data differs from the name in model_cols
print(table(dat$Dataset, dat$Location, useNA = "ifany"))
print(table(dat_m$Model, useNA = "ifany"))


##############################################################
# Figure 1. Correlation landscape by genomic location
##############################################################

# ---- 1A: E-M correlation by location, cancers side by side ----
## One box per cancer at each location (dodged side by side).
##   - dashed line at 0 separates negative (expected for promoter
##     methylation silencing expression) from positive correlations
##   - outlier.shape = NA hides the outlier points: with ~300k trios they
##     would cover the plot. The whiskers still show the 1.5 x IQR range.
##   - coord_cartesian() ZOOMS the y-axis to -0.6..0.65 without removing
##     data, so medians and boxes are still computed from all values
##     (scale_y_continuous(limits = ...) would drop data and change boxes)
p1A <- ggplot(dat, aes(x = Location, y = cor_exp_meth, fill = Dataset)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8),
               width = 0.7, linewidth = 0.4) +
  scale_fill_manual(values = cancer_cols, name = "Cancer") +
  scale_y_continuous(breaks = round(seq(-0.6, 0.6, 0.2), 1)) +
  coord_cartesian(ylim = c(-0.6, 0.65)) +
  labs(x = NULL, y = "E\u2013M correlation") +
  base_theme +
  theme(legend.position = "top")

# ---- 1B: medians per location, one line per cancer (trend along the gene) ----
## Summarise each cancer x location to its median (q1/q3 are kept for the
## table / text if needed). Joining the medians with a line shows how the
## correlation changes from upstream of the TSS into the gene body.
## The y-axis is NOT fixed here, so small differences between locations
## are visible (panel A's wide axis hides them).
## The legend is hidden because panel A already shows the cancer colours.
med_df <- dat[, .(median_r = median(cor_exp_meth),
                  q1 = quantile(cor_exp_meth, 0.25),
                  q3 = quantile(cor_exp_meth, 0.75)),
              by = .(Dataset, Location)]

p1B <- ggplot(med_df, aes(x = Location, y = median_r,
                          colour = Dataset, group = Dataset)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = cancer_cols, name = "Cancer") +
  labs(x = NULL, y = "Median E\u2013M correlation") +
  base_theme +
  theme(legend.position = "none")

# ---- 1C: E-M correlation by location x inferred model, one row per cancer ----
## n_m     : number of CpG-gene pairs behind each box (also used for Table S1)
## small_n : boxes below small_cut, which get a "*"
## dodge_m : shared dodge for boxes and stars. preserve = "single" keeps
##           every box the same width even when a model is missing at a
##           location, so boxes of the same model line up across locations.
n_m     <- dat_m[, .N, by = .(Dataset, Location, Model)]
small_n <- n_m[N < small_cut]
dodge_m <- position_dodge(width = dodge_w, preserve = "single")

## Boxes are filled by model and dodged within each location; one facet
## row per cancer (ncol = 1) gives each panel the full width it needs for
## 6 locations x up to 9 models.
## The "*" is drawn at y = 0.95 above its box. group = Model tells the
## dodge which slot the star belongs to; inherit.aes = FALSE stops it from
## inheriting the y = cor_exp_meth mapping of the main plot.
## ylim goes slightly above 1 so the stars are not cut off.
p1C <- ggplot(dat_m, aes(x = Location, y = cor_exp_meth, fill = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_boxplot(position = dodge_m, width = 0.8, outlier.shape = NA, linewidth = 0.25) +
  geom_text(data = small_n,
            aes(x = Location, y = 0.95, label = "*", group = Model),
            position = dodge_m, size = 3.2, inherit.aes = FALSE) +
  facet_wrap(~ Dataset, ncol = 1) +
  scale_fill_manual(values = model_cols, name = "Model") +
  scale_y_continuous(breaks = seq(-1, 1, 0.5)) +
  coord_cartesian(ylim = c(-1, 1.02)) +
  labs(x = NULL, y = "E\u2013M correlation",
       caption = paste0("* n < ", format(small_cut, big.mark = ","),
                        " CpG\u2013gene pairs; full counts in Table S1")) +
  guides(fill = guide_legend(nrow = 1)) +
  base_theme +
  theme(legend.position = "bottom",
        plot.caption = element_text(size = 8.5, colour = "grey40", hjust = 0))

# ---- combine & save ----
## Layout with patchwork:
##   top row    : A | B, with A about twice as wide as B (widths = c(2.2, 1))
##   bottom row : C, about 2.6 times as tall as the top row
## tag_levels = "A" labels the panels A, B, C automatically.
## cairo_pdf embeds fonts and renders the Unicode minus / en dash correctly.
fig1 <- ((p1A | p1B) + plot_layout(widths = c(2.2, 1))) / p1C +
  plot_layout(heights = c(1, 2.6)) +
  plot_annotation(
    title = "Correlation landscape by genomic location across cancers",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 14)))

fig1
ggsave("Fig1_correlation_by_location_6regions.pdf", fig1,
       width = 14, height = 15, device = cairo_pdf)

# ---- Figure 1 statistics and Table S1 ----
## (a) One summary cell per cancer x location:
##       median correlation (% of trios with a negative correlation, n)
##     dcast() puts cancers in rows and locations in columns, so the table
##     reads like panel A. Use these numbers in the results text.
print(dcast(dat[, .(val = sprintf("%.3f (%.0f%% neg, n=%s)",
                                  median(cor_exp_meth),
                                  100 * mean(cor_exp_meth < 0),
                                  format(.N, big.mark = ","))),
                by = .(Dataset, Location)],
            Dataset ~ Location, value.var = "val"))

## (b) Does the correlation differ between locations? One test per cancer.
##     Kruskal-Wallis = rank-based ANOVA (no normality assumption).
##     With hundreds of thousands of trios, p will be tiny almost always,
##     so also report the effect size epsilon^2 = H / (n - 1):
##       ~0.01 small, ~0.08 medium, >= 0.26 large (common rule of thumb).
##     Note: split trios appear in two locations, so the rows are not fully
##     independent; treat p as descriptive.
print(dat[, .(kw_p = kruskal.test(cor_exp_meth ~ Location)$p.value,
              eps2 = {h <- kruskal.test(cor_exp_meth ~ Location)$statistic
              round(as.numeric(h) / (.N - 1), 3)}),   # effect size epsilon^2
          by = Dataset])

## (c) Table S1: number of CpG-gene pairs behind every box in panel C.
##     One row per cancer x model, one column per location; fill = 0 for
##     model/location combinations that do not occur.
fwrite(dcast(n_m, Dataset + Model ~ Location, value.var = "N", fill = 0),
       "TableS1_n_pairs_by_cancer_model_location_6regions.csv")


##############################################################
# Figure 2. Causal-model distribution by genomic location
##############################################################

# ---- proportions on a complete grid (zeros where a model is absent) ----
## 1. Count trios per cancer x location x model.
## 2. CJ() ("cross join") builds every possible cancer x location x model
##    combination; joining the counts onto it adds rows for combinations
##    with no trios. Without this, a missing model would make bars and
##    labels shift position in panel B.
## 3. Missing counts (NA after the join) become 0.
## 4. total = all trios at that cancer x location; prop = share of each model.
##    So the models in one bar of panel A add up to 100%.
## 5. lab = % text for panel A, blank for segments below label_cut. Labels
##    are made for ALL rows (some empty) so the text stacks in exactly the
##    same order as the bar segments.
prop_df <- dat_m[, .N, by = .(Dataset, Location, Model)]
prop_df <- prop_df[CJ(Dataset  = factor(cancer_levels, levels = cancer_levels),
                      Location = factor(loc_labels,    levels = loc_labels),
                      Model    = factor(model_levels,  levels = model_levels)),
                   on = .(Dataset, Location, Model)]
prop_df[is.na(N), N := 0L]
prop_df[, total := sum(N), by = .(Dataset, Location)]
prop_df[, prop  := N / total]
prop_df[, lab   := fifelse(prop >= label_cut, scales::percent(prop, accuracy = 1), "")]

## Total pairs per bar, shown above each bar in thousands (e.g. "42k"),
## so readers can see how many trios each 100% bar represents.
tot_df <- unique(prop_df[, .(Dataset, Location, total)])
tot_df[, tot_lab := paste0(round(total / 1000), "k")]

# ---- 2A: 100% stacked bars, all models ----
## One bar per location, segments = models, one panel per cancer.
##   - position_stack(vjust = 0.5) puts each % label in the middle of its
##     segment
##   - the totals are drawn at y = 1.04, just above the bars; clip = "off"
##     and ylim up to 1.07 make room for them without a y-axis tick at 104%
##   - expand = c(0, 0.02) makes the bars start right at 0 (no gap)
##   - x labels are angled because six location names do not fit side by
##     side in a quarter-width panel
p2A <- ggplot(prop_df, aes(x = Location, y = prop, fill = Model)) +
  geom_col(width = 0.75, colour = "grey30", linewidth = 0.2) +
  geom_text(aes(label = lab), position = position_stack(vjust = 0.5), size = 2.5) +
  geom_text(data = tot_df, aes(x = Location, y = 1.04, label = tot_lab),
            inherit.aes = FALSE, size = 2.8, colour = "grey30") +
  facet_wrap(~ Dataset, nrow = 1) +
  scale_fill_manual(values = model_cols, name = "Model") +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, 1.07), clip = "off") +
  labs(x = NULL, y = "Proportion of CpG\u2013gene pairs") +
  guides(fill = guide_legend(nrow = 1)) +
  base_theme +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.45, "cm"),
        axis.text.x = element_text(angle = 40, hjust = 1))

# ---- 2B: % of each model per location (dodged), easier to compare locations ----
## Same numbers as panel A, rearranged: x = model, one bar per location.
## In a stacked bar only the bottom segment shares a baseline, so changes
## in middle segments are hard to judge; here every bar starts at 0, which
## makes location differences within each model easy to compare
## (e.g. is M1.2 more common in the gene body than at TSS200?).
## Bars of the same colour across models add up to 100% within a cancer.
p2B <- ggplot(prop_df, aes(x = Model, y = 100 * prop, fill = Location)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8,
           colour = "white", linewidth = 0.15) +
  facet_wrap(~ Dataset, ncol = 2) +
  scale_fill_manual(values = loc_cols, name = "Location") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "% of pairs at location") +
  guides(fill = guide_legend(nrow = 1)) +
  base_theme +
  theme(legend.position = "bottom")

# ---- combine & save ----
## A on top, B below (B a bit taller because it has two rows of panels).
## The caption explains the numbers above the bars in panel A.
fig2 <- p2A / p2B +
  plot_layout(heights = c(1, 1.3)) +
  plot_annotation(
    title = "Causal-model distribution by genomic location across cancers",
    caption = "Numbers above bars in A: total CpG\u2013gene pairs per location",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.caption = element_text(size = 9, colour = "grey30", hjust = 0)))

fig2
ggsave("Fig2_model_distribution_by_location_6regions.pdf", fig2,
       width = 14, height = 11, device = cairo_pdf)

# ---- Figure 2 statistics ----
## (a) The % behind both panels: one row per cancer x model, one column
##     per location. Use it to quote numbers in the text
##     (e.g. "M1.2 rose from x% at TSS200 to y% in the gene body").
print(dcast(prop_df[, .(Dataset, Model, Location, pct = round(100 * prop, 1))],
            Dataset + Model ~ Location, value.var = "pct"))

# Dataset  Model TSS1500 TSS200 5'UTR 1st exon Gene body 3'UTR
# <fctr> <fctr>   <num>  <num> <num>    <num>     <num> <num>
# 1:     BLCA   M0.1    31.2   28.3  26.7     26.3      25.0  27.3
# 2:     BLCA   M0.2     6.6    5.0   4.5      5.6       5.4   5.9
# 3:     BLCA   M1.1    10.8    9.9  14.2     10.8      15.4  14.2
# 4:     BLCA   M1.2     2.4    2.1   3.0      2.6       4.4   4.4
# 5:     BLCA   M2.1     2.1    2.7   4.3      3.0       5.7   4.8
# 6:     BLCA   M2.2     0.6    0.5   0.9      0.6       2.1   1.8
# 7:     BLCA     M3    19.1   25.0  19.5     20.7      13.7  15.2
# 8:     BLCA     M4     5.0    7.1   5.8      5.2       5.1   4.8
# 9:     BLCA  Other    22.2   19.4  21.0     25.2      23.1  21.6
# 10:     LIHC   M0.1    35.0   32.4  30.1     29.1      28.6  31.8
# 11:     LIHC   M0.2     5.6    4.6   4.8      5.2       5.4   5.9
# 12:     LIHC   M1.1     8.1    7.6  10.7      7.8      12.5  10.3
# 13:     LIHC   M1.2     1.7    1.7   2.2      2.1       3.8   2.7
# 14:     LIHC   M2.1     1.0    1.3   2.0      1.5       2.7   2.0
# 15:     LIHC   M2.2     0.4    0.3   0.6      0.5       1.2   1.0
# 16:     LIHC     M3    15.3   22.2  16.8     17.9      12.7  13.6
# 17:     LIHC     M4     2.8    3.6   3.2      2.8       2.7   2.2
# 18:     LIHC  Other    30.2   26.3  29.7     33.0      30.4  30.5
# 19: BRCA ER−   M0.1    42.9   45.9  42.2     40.0      40.1  41.8
# 20: BRCA ER−   M0.2     5.5    5.0   5.1      5.1       5.9   6.3
# 21: BRCA ER−   M1.1     3.3    3.5   4.4      3.1       4.1   3.4
# 22: BRCA ER−   M1.2     0.5    0.6   0.6      0.6       0.8   0.7
# 23: BRCA ER−   M2.1     0.4    0.7   0.9      0.9       0.7   0.5
# 24: BRCA ER−   M2.2     0.0    0.1   0.1      0.1       0.4   0.1
# 25: BRCA ER−     M3     8.1   10.0   8.2      6.9       6.8   8.0
# 26: BRCA ER−     M4     0.3    0.3   0.3      0.2       0.3   0.3
# 27: BRCA ER−  Other    39.2   33.9  38.1     43.2      40.8  38.8
# 28: BRCA ER+   M0.1    29.0   27.2  24.5     25.1      23.6  25.4
# 29: BRCA ER+   M0.2     7.8    6.0   6.5      7.7       7.1   8.5
# 30: BRCA ER+   M1.1    10.8    9.9  12.8     10.8      14.4  11.8
# 31: BRCA ER+   M1.2     2.5    2.4   3.0      3.2       3.5   3.6
# 32: BRCA ER+   M2.1     2.4    2.9   3.9      3.1       4.2   3.1
# 33: BRCA ER+   M2.2     0.5    0.6   0.8      0.7       1.5   1.1
# 34: BRCA ER+     M3    18.1   22.3  18.5     17.8      16.7  19.2
# 35: BRCA ER+     M4     7.3    9.6   8.4      6.8       7.4   7.4
# 36: BRCA ER+  Other    21.6   19.0  21.5     24.7      21.6  19.8
# Dataset  Model TSS1500 TSS200 5'UTR 1st exon Gene body 3'UTR
# <fctr> <fctr>   <num>  <num> <num>    <num>     <num> <num>
  
## (b) Does the model mix differ by location? One test per cancer.
##     - table(): location x model contingency table; droplevels() removes
##       empty levels so the test has no all-zero rows/columns
##     - chi-square test of independence; suppressWarnings() hides the
##       "approximation may be incorrect" warning from rare cells
##     - Cramer's V = sqrt(chi2 / (n * (min(rows, cols) - 1))), from 0 (no
##       association) to 1; with this many trios p is always tiny, so V is
##       the number to report (~0.1 small, ~0.3 medium, ~0.5 large)
##     As in Fig 1, split trios are counted twice, so n is slightly inflated.
print(dat_m[, {
  tab <- table(droplevels(Location), droplevels(Model))
  ct  <- suppressWarnings(chisq.test(tab))
  .(chisq_p   = ct$p.value,
    cramers_v = round(sqrt(ct$statistic / (sum(tab) * (min(dim(tab)) - 1))), 3))
}, by = Dataset])

# Dataset chisq_p cramers_v
# <fctr>   <num>     <num>
#   1: BRCA ER+       0     0.050
# 2: BRCA ER−       0     0.038
# 3:     BLCA       0     0.078
# 4:     LIHC       0     0.067




##############################################################
# Figure 3. E-M correlation by inferred causal direction
# Figure 4. E-M regression coefficient (coef12) by inferred causal direction
#           (both: BLCA, LIHC, BRCA ER-, BRCA ER+; all directed trios)
#
# Input: the Location table made by trio_probe_gene_annot_all_datasets.R
#        (trio_probe_gene_annot_loc_split_allCancerALL.rds)
#   - one row per trio x Location (6 Illumina regions)
#   - multi-region trios with >1 accession were removed
#   - 5'UTR;1stExon and 1stExon;3'UTR trios (1 accession) were split,
#     so they appear in BOTH locations
#
# Two kinds of tables are used below:
#   * location rows (one row per trio x Location): for the plots, which are
#     shown BY location, so a split trio contributes to each of its regions
#   * trio rows (one row per trio, made with unique(..., by = trio key)):
#     for the trio counts, so split trios are not counted twice
#
# Models (T1 = expression E, T2 = methylation M, V1 = CNA):
#   M1.1  CNA -> E -> M        chain,    direction E -> M  (E-mediated)
#   M1.2  CNA -> M -> E        chain,    direction M -> E  (M-mediated)
#   M2.1  CNA -> E <- M        collider, direction M -> E
#   M2.2  CNA -> M <- E        collider, direction E -> M
#   M4 (fully connected) is not used: its E-M direction is not identifiable.
#
# cor_exp_meth = marginal E-M correlation (Figure 3)
# coef12       = CNA- and confounder-adjusted E-M t-statistic (Figure 4)
#
# Distributions are drawn as binned densities, not kernel densities.
# Directed trios need a significant E-M edge, so values near zero
# (|t| < ~2 for coef12) are absent by construction. A kernel density
# would smooth across that empty region; binned lines show it as zero.
#
# Grey area = overlap of the two directions (in each bin, the smaller of
# the two densities). The % in each panel is the overlap coefficient:
# the total grey area, 100% = identical distributions, 0% = no overlap.
##############################################################

## Packages:
##   data.table : grouping, joins
##   ggplot2    : plots
library(data.table)
library(ggplot2)

# ======================= 0. Input and output =======================
combined_dir <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_combined"
loc_df <- readRDS(file.path(combined_dir, "trio_probe_gene_annot_loc_split_allCancerALL.rds"))
setDT(loc_df)                      # tibble -> data.table, in place

fig_dir <- file.path(combined_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)
setwd(fig_dir)

## Columns needed for Figures 3-4.
## Stop now with a clear message if any is missing or named differently.
need_cols <- c("Dataset", "trio_row", "Location",
               "Inferred.Model.BH_fdr", "cor_exp_meth", "coef12")
miss <- setdiff(need_cols, names(loc_df))
if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

## Unique trio key. A split trio has the same Dataset + trio_row on two rows.
trio_key <- c("Dataset", "trio_row")

# ======================= 1. Shared settings =======================
cancer_levels <- c("BLCA", "LIHC", "BRCA ER\u2212", "BRCA ER+")

## Six locations in gene order (upstream -> downstream), as in Figs 1-2.
## loc_levels = values in the data, loc_labels = what is shown on the axis.
loc_levels <- c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR")
loc_labels <- c("TSS1500", "TSS200", "5'UTR", "1st exon", "Gene body", "3'UTR")

## Direction colours and the model -> direction lookup table.
## "\u2192" is the arrow character.
dir_cols <- c("M \u2192 E" = "#D55E00", "E \u2192 M" = "#0072B2")
dir_tab  <- data.table(Model     = c("M1.1", "M1.2", "M2.1", "M2.2"),
                       Direction = c("E \u2192 M", "M \u2192 E", "M \u2192 E", "E \u2192 M"))

overlap_col <- "grey90"

## Convert raw dataset names to display names (BLCA / LIHC unchanged)
relabel_dataset <- function(x) fcase(x == "negER_BRCA", "BRCA ER\u2212",
                                     x == "posER_BRCA", "BRCA ER+",
                                     default = x)

## Shared plot theme (same as Figs 1-2)
base_theme <- theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# ======================= 2. Shared data: all directed trios =======================
## All M1.1 / M1.2 / M2.1 / M2.2 trios at the six locations (location rows)
dir_df <- loc_df[Inferred.Model.BH_fdr %in% dir_tab$Model & Location %in% loc_levels]

## Keep the original Dataset (part of trio_key) and add a Cancer column
## for plotting: negER_BRCA and posER_BRCA stay separate
dir_df[, Cancer := relabel_dataset(Dataset)]
dir_df <- dir_df[Cancer %in% cancer_levels]
dir_df[, Cancer   := factor(Cancer, levels = cancer_levels)]
dir_df[, Location := factor(Location, levels = loc_levels, labels = loc_labels)]

## Direction from the lookup table (data.table update-join)
dir_df[dir_tab, on = .(Inferred.Model.BH_fdr = Model), Direction := i.Direction]
dir_df[, Direction := factor(Direction, levels = names(dir_cols))]

# ======================= 3. Binning, overlap and plot functions =======================
## Binned density of one E-M measure, per cancer x location x direction.
## Bins are [k*binwidth, (k+1)*binwidth), so 0 is always a bin edge and no
## bin straddles the empty region around zero. Empty bins are kept as 0,
## so the lines drop to zero in the gap.
## density = n / (N * binwidth): each direction's area sums to 1 per panel.
bin_dir_dist <- function(dt, xvar, binwidth) {
  b   <- dt[, .(Cancer, Location, Direction, idx = floor(get(xvar) / binwidth))]
  cnt <- b[, .N, by = .(Cancer, Location, Direction, idx)]
  
  ## Full grid: every cell x direction x bin (one empty bin added at each end)
  fac  <- function(f) factor(levels(f), levels = levels(f))
  grid <- CJ(Cancer    = fac(dt$Cancer),
             Location  = fac(dt$Location),
             Direction = fac(dt$Direction),
             idx       = seq(min(b$idx) - 1, max(b$idx) + 1))
  
  res <- cnt[grid, on = .(Cancer, Location, Direction, idx)]
  res[is.na(N), N := 0L]
  res[, dens := if (sum(N) > 0) N / (sum(N) * binwidth) else 0,
      by = .(Cancer, Location, Direction)]
  res[, x := (idx + 0.5) * binwidth]            # bin midpoint
  res[]
}

## Overlap per bin (the smaller of the two densities) and per panel
## (overlap coefficient = total grey area).
overlap_dir_dist <- function(res, binwidth) {
  ov_bin  <- res[, .(ov = min(dens)), by = .(Cancer, Location, idx, x)]
  ov_cell <- ov_bin[, .(OVL = sum(ov) * binwidth), by = .(Cancer, Location)]
  list(bin = ov_bin, cell = ov_cell)
}

## Plot: grey overlap area, one coloured line per direction, overlap % label.
##   xlim   : optional zoom c(lo, hi); coord_cartesian keeps all data in
##            the bins and only changes the visible range
##   breaks : x-axis breaks (waiver() = ggplot default)
plot_dir_dist <- function(res, ov, xlab, title, xlim = NULL, breaks = waiver()) {
  ggplot() +
    geom_ribbon(data = ov$bin, aes(x = x, ymin = 0, ymax = ov),
                fill = overlap_col, colour = NA) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line(data = res, aes(x = x, y = dens, colour = Direction),
              linewidth = 0.7) +
    geom_text(data = ov$cell,
              aes(x = Inf, y = Inf, label = sprintf("overlap %.0f%%", 100 * OVL)),
              hjust = 1.1, vjust = 1.5, size = 3, colour = "grey25") +
    facet_grid(Cancer ~ Location) +
    scale_colour_manual(values = dir_cols) +
    scale_x_continuous(breaks = breaks) +
    coord_cartesian(xlim = xlim) +
    labs(title    = title,
         subtitle = "Trios with a directed E\u2013M edge (M1.1, M1.2, M2.1, M2.2)",
         x = xlab, y = "Density", colour = "Inferred direction",
         caption  = paste("Grey: overlap between the two directions;",
                          "% = overlap coefficient (100% = identical distributions)")) +
    base_theme +
    theme(legend.position = "top",
          plot.title    = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10, colour = "grey30"),
          plot.caption  = element_text(size = 8.5, colour = "grey40", hjust = 0))
}


##############################################################
# Figure 3. E-M correlation by inferred causal direction
##############################################################
## One line per direction in each cancer x location cell. If the two
## directions overlap, the marginal E-M correlation does not tell us
## whether M -> E or E -> M.

fig3_df <- dir_df[!is.na(cor_exp_meth)]

cat("\nFig 3 trios per cancer x direction:\n")
print(table(unique(fig3_df, by = trio_key)[, .(Cancer, Direction)]))

bw3   <- 0.05
res3  <- bin_dir_dist(fig3_df, "cor_exp_meth", bw3)
ov3   <- overlap_dir_dist(res3, bw3)

cat("\nFig 3 overlap coefficient (%):\n")
print(dcast(ov3$cell[, .(Cancer, Location, OVL = round(100 * OVL))],
            Cancer ~ Location, value.var = "OVL"))

fig3 <- plot_dir_dist(res3, ov3,
                      xlab   = "E\u2013M correlation",
                      title  = "E\u2013M correlation by inferred causal direction",
                      xlim   = c(-1, 1),
                      breaks = c(-1, -0.5, 0, 0.5, 1))
fig3
ggsave("Fig3.1_EM_correlation_vs_direction_4cancers_6regions.pdf", fig3,
       width = 14, height = 10, device = cairo_pdf)


##############################################################
# Figure 4. E-M regression coefficient (coef12) by inferred causal direction
##############################################################
## Same layout as Figure 3, but with coef12, the CNA- and confounder-
## adjusted E-M t-statistic, instead of the marginal correlation.
## coef12 is not bounded to [-1, 1] and can have long tails, so the x-axis
## is zoomed to the central 98% of values (1st-99th percentile).

fig4_df <- dir_df[!is.na(coef12)]

cat("\nFig 4 trios per cancer x direction:\n")
print(table(unique(fig4_df, by = trio_key)[, .(Cancer, Direction)]))

cat("\nFig 4 coef12 summary:\n")
print(summary(fig4_df$coef12))

bw4   <- 0.5
res4  <- bin_dir_dist(fig4_df, "coef12", bw4)
ov4   <- overlap_dir_dist(res4, bw4)

cat("\nFig 4 overlap coefficient (%):\n")
print(dcast(ov4$cell[, .(Cancer, Location, OVL = round(100 * OVL))],
            Cancer ~ Location, value.var = "OVL"))

coef_lim <- quantile(fig4_df$coef12, c(0.01, 0.99))

fig4 <- plot_dir_dist(res4, ov4,
                      xlab  = "Adjusted E\u2013M effect (coef12, t)",
                      title = "Adjusted E\u2013M effect by inferred causal direction",
                      xlim  = coef_lim)
fig4
ggsave("Fig4.1_EM_coef12_vs_direction_4cancers_6regions.pdf", fig4,
       width = 14, height = 10, device = cairo_pdf)


