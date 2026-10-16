all_annot <- readRDS("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_combined/ALL_trio_results_with_stats_location.rds")
df=all_annot[,c(1:14,58:60)]
head(df)

library(data.table)
library(ggplot2)
library(ggpubr)

##########################################

library(data.table)
library(ggplot2)
library(ggpubr)   # for stat_compare_means

df <- all_annot[, c(1:14, 58:60)]

# Check what's in Group first (there may be NA, "", 1stExon, IGR, etc.)
table(df$Group, df$Dataset, useNA = "ifany")

loc_levels <- c("TSS", "5'/3'UTR", "Body")

plot_df <- df[!is.na(cor_exp_meth) & Group %in% loc_levels]
plot_df[, Group := factor(Group, levels = loc_levels,
                          labels = c("TSS", "3'/5' UTR", "Gene body"))]
plot_df[, Dataset := factor(Dataset)]   # reorder levels here if you want a specific order

## Main Figure:
if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
library(ggpubr)
loc_cols <- c("TSS" = "#E64B35", "3'/5' UTR" = "#4DBBD5", "Gene body" = "#00A087")

my_comparisons <- list(c("TSS", "3'/5' UTR"),
                       c("TSS", "Gene body"),
                       c("3'/5' UTR", "Gene body"))


################


########################
plot_df[, Dataset := factor(Dataset,
                            levels = c("BLCA", "LIHC", "BRCA ER\u2212", "BRCA ER+"))]

cancer_cols <- c("BLCA"        = "#66C2A5",   # green
                 "LIHC"        = "#FC8D62",   # orange
                 "BRCA ER\u2212" = "#8DA0CB", # blue
                 "BRCA ER+"    = "#E78AC3")   # pink

p2 <- ggplot(plot_df, aes(x = Group, y = cor_exp_meth, fill = Dataset)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8),
               width = 0.7, linewidth = 0.4) +
  scale_fill_manual(values = cancer_cols) +
  scale_y_continuous(breaks = round(seq(-0.6, 0.6, 0.2), 1)) +
  coord_cartesian(ylim = c(-0.6, 0.65)) +
  labs(x = NULL, y = "Correlation (Expression vs. Methylation)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
p2
# ggsave("cor_EM_by_location.pdf", p1, width = 10, height = 4.5)


###############################
#############################
# ---- data ----
grid_df <- plot_df[!is.na(Inferred.Model.BH_fdr) & Inferred.Model.BH_fdr != "" &
                     !is.na(cor_exp_meth)]

# models: natural order, "Other" last
model_levels <- sort(unique(grid_df$Inferred.Model.BH_fdr))
model_levels <- c(setdiff(model_levels, "Other"), intersect("Other", model_levels))
grid_df[, Model := factor(Inferred.Model.BH_fdr, levels = model_levels)]   # <- the missing step

# locations: match by pattern so any earlier label works
grid_df[, Group := fcase(
  grepl("UTR", Group),       "5'/3' UTR",
  grepl("TSS", Group),       "TSS",
  grepl("Body|body", Group), "Gene body"
)]
grid_df <- grid_df[!is.na(Group)]
grid_df[, Group := factor(Group, levels = c("5'/3' UTR", "TSS", "Gene body"))]

# cancers (use "BRCA ER-" here if you switched to a plain hyphen)
grid_df[, Dataset := factor(as.character(Dataset),
                            levels = c("BLCA", "LIHC", "BRCA ER\u2212", "BRCA ER+"))]

# counts
n_grid <- grid_df[, .N, by = .(Dataset, Model, Group)]
n_grid[, lab := ifelse(N >= 1000, paste0(round(N / 1000, 1), "k"), as.character(N))]

loc_cols <- c("5'/3' UTR" = "#A6CEE3", "TSS" = "#E78AC3", "Gene body" = "#A6D854")

# quick checks: no NA rows should appear
table(grid_df$Group, useNA = "ifany")
table(grid_df$Dataset, useNA = "ifany")

p_grid <- ggplot(grid_df, aes(x = Group, y = cor_exp_meth, fill = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.35) +
  geom_text(data = n_grid, aes(x = Group, y = -0.95, label = lab),
            inherit.aes = FALSE, size = 2, colour = "grey30") +
  facet_grid(Dataset ~ Model) +
  scale_fill_manual(values = loc_cols) +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5)) +
  coord_cartesian(ylim = c(-1, 0.9)) +
  labs(x = "Location", y = "Correlation (Expression vs. Methylation)",
       title = "E–M correlation by genomic location, cancer type and inferred model") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(0.3, "lines"))

p_grid
ggsave("cor_EM_location_by_cancer_model.pdf", p_grid, width = 15, height = 9)


#######################################
library(data.table)
library(ggplot2)

# ---- data ----
comb_df <- plot_df[!is.na(Inferred.Model.BH_fdr) & Inferred.Model.BH_fdr != ""]

model_levels <- sort(unique(comb_df$Inferred.Model.BH_fdr))   # M0.1 ... M4, Other
model_levels <- c(setdiff(model_levels, "Other"), "Other")   # keep "Other" last
comb_df[, Model := factor(Inferred.Model.BH_fdr, levels = model_levels)]

comb_df[, Dataset := factor(Dataset,
                            levels = c("BLCA", "LIHC", "BRCA ER\u2212", "BRCA ER+"))]

comb_df[, Group := fcase(
  grepl("TSS",  Group),                 "TSS",
  grepl("UTR",  Group),                 "5'/3' UTR",
  grepl("Body|body", Group),            "Gene body"
)]
comb_df <- comb_df[!is.na(Group)]
comb_df[, Group := factor(Group, levels = c("TSS", "5'/3' UTR", "Gene body"))]

n_comb <- comb_df[, .N, by = .(Dataset, Group, Model)]   # recompute after fixing

# Set2 colours as in your attached figure, white for "Other"
model_cols <- setNames(
  c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
    "#FFD92F", "#E5C494", "#B3B3B3", "white")[seq_along(model_levels)],
  model_levels)

# ---- plot ----
dodge <- position_dodge(width = 0.85)

p_comb <- ggplot(comb_df, aes(x = Group, y = cor_exp_meth, fill = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_boxplot(position = dodge, width = 0.8, outlier.shape = NA,
               linewidth = 0.3) +
  geom_text(data = n_comb,
            aes(x = Group, y = -1.05, label = N, group = Model),
            position = dodge, angle = 90, hjust = 0, size = 2.2,
            inherit.aes = FALSE) +
  facet_wrap(~ Dataset, ncol = 2) +
  scale_fill_manual(values = model_cols) +
  scale_y_continuous(breaks = seq(-1, 1, 0.5)) +
  coord_cartesian(ylim = c(-1.05, 1)) +
  labs(x = NULL, y = "Correlation (expression vs. methylation)", fill = "Model") +
  guides(fill = guide_legend(nrow = 1)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

p_comb
ggsave("cor_EM_location_model_combined.pdf", p_comb, width = 13, height = 8)

#############################
#### Figure 2
############################

library(data.table)
library(ggplot2)

# ---- proportions per cancer x location ----
prop_df <- comb_df[, .N, by = .(Dataset, Group, Model)]
prop_df[, total := sum(N), by = .(Dataset, Group)]
prop_df[, prop := N / total]

tot_df <- unique(prop_df[, .(Dataset, Group, total)])
tot_df[, lab := paste0(round(total / 1000), "k")]

# ---- plot ----
p_dist <- ggplot(prop_df, aes(x = Group, y = prop, fill = Model)) +
  geom_col(width = 0.75, colour = "grey30", linewidth = 0.2) +
  geom_text(data = tot_df, aes(x = Group, y = 1.04, label = lab),
            inherit.aes = FALSE, size = 3, colour = "grey30") +
  facet_wrap(~ Dataset, nrow = 1) +
  scale_fill_manual(values = model_cols) +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, 1.07), clip = "off") +
  labs(x = NULL, y = "Proportion of CpG\u2013gene pairs", fill = "Model",
       title = "Causal-model distribution by genomic location",
       caption = "Numbers above bars: total CpG\u2013gene pairs per location") +
  guides(fill = guide_legend(nrow = 1)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.caption = element_text(size = 9, colour = "grey30", hjust = 0))
p_dist <- p_dist +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.key.size = unit(0.45, "cm"),
        legend.text = element_text(size = 9),
        plot.margin = margin(5, 10, 5, 5))

model_cols <- c(M0.1 = "#A6CEE3", M0.2 = "#1F78B4", M1.1 = "#B2DF8A", M1.2 = "#33A02C",
                M2.1 = "#FB9A99", M2.2 = "#E31A1C", M3 = "#FDBF6F", M4 = "#FF7F00",
                Other = "grey75")

p_dist <- p_dist +
  scale_fill_manual(values = model_cols, name = "Model") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))   # also fixes the cut-off "Other"

p_dist
ggsave("model_distribution_by_location.pdf", p_dist,
       width = 11, height = 5.5, device = cairo_pdf)

############################################################################
