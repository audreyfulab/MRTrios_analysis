setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/raw_Data_Methyl")
cna         <- fread('data_CNA.txt', data.table = F)
gene.exp    <- fread('data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt', data.table = F)
TCGA.meth   <- fread("split.names.TCGA.meth.logit.txt", data.table = F)
trios       <- fread("trio.final.protein.coding.txt")
clinical    <- fread("names.pos.patient2.txt", header = FALSE)


com.ind <- intersect(
  intersect(colnames(gene.exp)[3:ncol(gene.exp)],
            colnames(TCGA.meth)[5:ncol(TCGA.meth)]),
  colnames(cna)[3:ncol(cna)])

unique_id <- intersect(com.ind, clinical$V1)

clinical_filter <- clinical %>%
  arrange(match(V1, unique_id)) %>%
  filter(V1 %in% unique_id)

age  <- clinical_filter$V2
race <- clinical_filter$V3

# ── Treat empty/blank/placeholder values as NA ─────────────
# fread() already turns truly empty numeric cells into NA, but
# values stored as blank strings ("", " ") or literal "NA"/"NULL"
# text survive as character data. Normalize all of these to NA
# before they're used as confounders.
age  <- trimws(as.character(age))
race <- trimws(as.character(race))
age[age   %in% c("", "NA", "NULL", "N/A")]  <- NA
race[race %in% c("", "NA", "NULL", "N/A")]  <- NA
age <- suppressWarnings(as.numeric(age))

cna_filter <- cna %>%
  mutate(cna.row = rownames(cna)) %>%
  filter(cna.row %in% unique(trios$cna.row)) %>%
  select(cna.row, all_of(unique_id))

exp_filter <- gene.exp %>%
  mutate(gene.row = rownames(gene.exp)) %>%
  filter(gene.row %in% unique(trios$gene.row)) %>%
  select(gene.row, all_of(unique_id))

meth_filter <- TCGA.meth %>%
  mutate(meth.row = rownames(TCGA.meth)) %>%
  filter(meth.row %in% unique(trios$meth.row)) %>%
  select(meth.row, all_of(unique_id))

# cna_filter <- fread("pos_cna.txt")
# exp_filter <- fread("pos_exp.txt")
# meth_filter <- fread("pos_meth.txt")

# ======================================================================
###For m0.1(C -> E):###
#Fliter the M0.1 model datasets
# mod0_1 <- pos_model %>% filter(Inferred.Model=="M0.1",Inferred.Model.BH_fdr=="M0.1",Inferred.Model.qval=="M0.1")
mod0_1 <- pos_model %>% filter(Inferred.Model=="M0.1",Inferred.Model.BH_fdr=="M0.1")

head(mod0_1,3)
#Running the trio_row number#270001
# trio (270001) isn't M0.1 model when check the Inferred.Model.qval
trios[270001,]
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
#make sure these values are not in "lists"
#Find the correlation value of 2 varibales
# dat_m0_1 <- data.frame(CNA=unlist(positive_CNA[264,]),Meth=unlist(positive_meth[462699,]),GE=unlist(positive_GE[1197,]))
dat_m0_1 <- data.frame(CNA=unlist(cna_filter %>% filter(cna.row==264) %>% select(-cna.row)),
                       Meth=unlist(meth_filter %>% filter(meth.row==462699) %>% select(-meth.row)),
                       GE=unlist(exp_filter %>% filter(gene.row==1197) %>% select(-gene.row)))

dim(dat_m0_1)
cor(dat_m0_1$GE,dat_m0_1$Meth, use = "complete.obs")
cor(dat_m0_1$CNA,dat_m0_1$GE, use = "complete.obs")
cor(dat_m0_1$CNA,dat_m0_1$Meth, use = "complete.obs")

##ggplot2##

#plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8)
library(ggplot2)
library(patchwork)
# Scatterplot GE & Meth
p01_1=ggplot(dat_m0_1,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M0.1 Gene Expression v. Methylation',subtitle = "(r=0.118)", x= 'Gene Expression', y = 'Methylation')
#Boxplot CNA & GE
p01_2 <- ggplot(dat_m0_1, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8))+labs(title = 'M0.1 CNA v. GE', subtitle = "(r=0.327)",x = 'Copy Number Alternation', y = 'Gene Expression')
# Boxplot CNA & Meth
p01_3 <- ggplot(dat_m0_1, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M0.1 CNA v. Methylation', subtitle = "(r=0.120)",x = 'Copy Number Alternation', y = 'Methylation')
(p01_1)/(p01_2|p01_3)
ggsave("Pos_M0.1.pdf")


##  read the saved filter cna,exp,meth data
setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026")
cna_filter <- fread("pos_cna.txt")
exp_filter <- fread("pos_exp.txt")
meth_filter <- fread("pos_meth.txt")

# ── Paths (defined once) ──────────────────────────────────────────────────────
DATA_DIR <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA"
RAW_DIR  <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/raw_Data_Methyl"
FIG_DIR  <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Figure"

# ── Load model results (loaded once) ─────────────────────────────────────────
model_pos     <- fread(file.path(DATA_DIR, "Output_posER_BRCA/posER_BRCA_trio_Model_results_ALL_with_BH_fdr_qval_byLZ.txt"))
model_neg     <- fread(file.path(DATA_DIR, "Output_negER_BRCA/negER_BRCA_trio_Model_results_ALL_with_BH_fdr_qval_byLZ.txt"))
model_pos_loc <- readRDS(file.path(DATA_DIR, "Output_posER_BRCA/LZ_trio_results_pos_ALL_with_BH_fdr_qval_with_location.rds"))
model_neg_loc <- readRDS(file.path(DATA_DIR, "Output_negER_BRCA/LZ_trio_results_neg_ALL_with_BH_fdr_qval_with_location.rds"))


## M1.1
### trio: 135939 : 4 infer model are M1.1
model_pos %>% filter(trio_row==135939)

trio_df <- data.frame(CNA=unlist(cna_filter %>% filter(cna.row==model_pos %>% filter(trio_row==135939) %>% pull(cna.row)) %>% select(-cna.row)),
                       Meth=unlist(meth_filter %>% filter(meth.row==model_pos %>% filter(trio_row==135939) %>% pull(meth.row)) %>% select(-meth.row)),
                       GE=unlist(exp_filter %>% filter(gene.row==model_pos %>% filter(trio_row==135939) %>% pull(gene.row)) %>% select(-gene.row)))

dim(trio_df)
cor(trio_df$GE,trio_df$Meth, use = "complete.obs")
cor(trio_df$CNA,trio_df$GE, use = "complete.obs")
cor(trio_df$CNA,trio_df$Meth, use = "complete.obs")

# cor(trio_df$GE,trio_df$Meth, use = "complete.obs")  -0.6762893
# cor(trio_df$CNA,trio_df$GE, use = "complete.obs")    0.08559773
# cor(trio_df$CNA,trio_df$Meth, use = "complete.obs")  0.02036697

# Scatterplot GE & Meth
p01_1=ggplot(trio_df,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M0.1 Gene Expression v. Methylation',subtitle = "(r=-0.6762893)", x= 'Gene Expression', y = 'Methylation')
#Boxplot CNA & GE
p01_2 <- ggplot(trio_df, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8))+labs(title = 'M0.1 CNA v. GE', subtitle = "(r=0.08559773)",x = 'Copy Number Alternation', y = 'Gene Expression')
# Boxplot CNA & Meth
p01_3 <- ggplot(trio_df, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M0.1 CNA v. Methylation', subtitle = "(r=0.02036697)",x = 'Copy Number Alternation', y = 'Methylation')
(p01_1)/(p01_2|p01_3)
ggsave("Pos_M1.1.pdf")

# ====================================== Noconfounder models ======================================================

library(data.table)


## Local path:
input_dir <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/NoConfounder_Output_posER_BRCA_part"
output_dir <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_posER_BRCA"

# ── Find all shard files (exclude merged ALL file) ─────────
shard_files <- sort(list.files(input_dir,
                               pattern = "NoConfounder_trio_posER_BRCA_results_part.*\\.txt",
                               full.names = TRUE))

cat(sprintf("Found %d shard files. Merging...\n", length(shard_files)))
##  599 shard files : missing part57

# ── Merge all shards ───────────────────────────────────────
all_results <- rbindlist(lapply(shard_files, fread), fill = TRUE)

cat(sprintf("Total rows: %d | Total cols: %d\n",
            nrow(all_results), ncol(all_results)))

## Total rows: 292043 | Total cols: 25

# ── Save merged file ───────────────────────────────────────
out_file <- file.path(output_dir, "NoConfounder_trio_posER_BRCA_results_ALL_missPart57.txt")
fwrite(all_results, file = out_file, sep = "\t")

cat(sprintf("Saved: %s\n", out_file))


## 
Noconfounder_pos2 <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_posER_BRCA/NoConfounder_trio_posER_BRCA_results_ALL_missPart57.txt")

Noconfounder_pos2_filter <- Noconfounder_pos2 %>% mutate(InferModel_Noconfounder2=Inferred.Model) %>% select(trio_row,cor_CNA_Exp,cor_CNA_Meth,cor_Exp_Meth,InferModel_Noconfounder2) 

model_po2_loc_full <- Noconfounder_pos2_filter %>% left_join(model_pos_loc,by="trio_row")

m1.1_interset2 <- model_po2_loc_full %>% filter(Inferred.Model=="M1.1",Inferred.Model.BH_fdr=="M1.1",Inferred.Model.qval=="M1.1",InferModel_Noconfounder2=="M1.1")

m1.1_interset2_filter <- m1.1_interset2 %>% select(trio_row,cna.row,gene.row,meth.row,Total.PC.Count,Confounders,Inferred.Model,Inferred.Model.BH_fdr,Inferred.Model.qval,InferModel_Noconfounder2,Name,Gene.name,UCSC_RefGene_Group,Group,cor_CNA_Exp,cor_CNA_Meth,cor_Exp_Meth) %>% distinct()

dff <- m1.1_interset2_filter %>% select(-c(UCSC_RefGene_Group,Group)) %>% distinct()

dff2=dff %>% filter(abs(cor_CNA_Exp) < abs(cor_Exp_Meth)) %>%
  arrange(desc(abs(cor_Exp_Meth) - abs(cor_CNA_Exp)))
## 4430 rows

# ================================================ Noconfounder models  com id with clinical ======================================================

library(data.table)


## Local path:
input_dir <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/NoConfounder_Clinical_Output_posER_BRCA_part"
output_dir <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_posER_BRCA"

# ── Find all shard files (exclude merged ALL file) ─────────
shard_files <- sort(list.files(input_dir,
                               pattern = "NoConfounder_Clinical_trio_posER_BRCA_results_part.*\\.txt",
                               full.names = TRUE))

cat(sprintf("Found %d shard files. Merging...\n", length(shard_files)))
##  599 shard files : missing part57

# ── Merge all shards ───────────────────────────────────────
all_results <- rbindlist(lapply(shard_files, fread), fill = TRUE)

cat(sprintf("Total rows: %d | Total cols: %d\n",
            nrow(all_results), ncol(all_results)))

## Total rows: 292043 | Total cols: 25

# ── Save merged file ───────────────────────────────────────
out_file <- file.path(output_dir, "NoConfounder_Clinical_trio_posER_BRCA_results_ALL_missPart57.txt")
fwrite(all_results, file = out_file, sep = "\t")

cat(sprintf("Saved: %s\n", out_file))


## 
Noconfounder_pos <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_posER_BRCA/NoConfounder_Clinical_trio_posER_BRCA_results_ALL_missPart57.txt")

Noconfounder_pos_filter <- Noconfounder_pos %>% mutate(InferModel_Noconfounder=Inferred.Model) %>% select(trio_row,cor_CNA_Exp,cor_CNA_Meth,cor_Exp_Meth,InferModel_Noconfounder) 

model_pos_loc_full <- model_pos_loc %>% left_join(Noconfounder_pos_filter,by="trio_row")

m1.1_interset <- model_pos_loc_full %>% filter(Inferred.Model=="M1.1",Inferred.Model.BH_fdr=="M1.1",Inferred.Model.qval=="M1.1",InferModel_Noconfounder=="M1.1")

m1.1_interset_filter <- m1.1_interset %>% select(trio_row,cna.row,gene.row,meth.row,Total.PC.Count,Confounders,Inferred.Model,Inferred.Model.BH_fdr,Inferred.Model.qval,InferModel_Noconfounder,Name,Gene.name,UCSC_RefGene_Group,Group,cor_CNA_Exp,cor_CNA_Meth,cor_Exp_Meth) %>% distinct()
library(tidyverse)
df <- m1.1_interset_filter %>% select(-c(UCSC_RefGene_Group,Group)) %>% distinct()

df_interest=df %>% filter(abs(cor_CNA_Exp) < abs(cor_CNA_Meth)) %>%
  arrange(desc(abs(cor_CNA_Meth) - abs(cor_CNA_Exp)))
## trio: 213934,213938,213944 only 3 trios have cor_C_M: 0.2 > cor_C_E :0.1
## around 10~11k trios result m1.1 result affected and not affected by confounders.
# ==================================================================================================================================================================

trioRow=213934

trio_df <- data.frame(CNA=unlist(cna_filter %>% filter(cna.row==model_pos %>% filter(trio_row==trioRow) %>% pull(cna.row)) %>% select(-cna.row)),
                      Meth=unlist(meth_filter %>% filter(meth.row==model_pos %>% filter(trio_row==trioRow) %>% pull(meth.row)) %>% select(-meth.row)),
                      GE=unlist(exp_filter %>% filter(gene.row==model_pos %>% filter(trio_row==trioRow) %>% pull(gene.row)) %>% select(-gene.row)))

dim(trio_df)
cor(trio_df$GE,trio_df$Meth, use = "complete.obs")
cor(trio_df$CNA,trio_df$GE, use = "complete.obs")
cor(trio_df$CNA,trio_df$Meth, use = "complete.obs")

# Scatterplot GE & Meth
p01_1=ggplot(trio_df,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 Gene Expression v. Methylation',subtitle = "(r= -0.5890477)", x= 'Gene Expression', y = 'Methylation')
#Boxplot CNA & GE
p01_2 <- ggplot(trio_df, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. GE', subtitle = "(r= 0.1220043)",x = 'Copy Number Alternation', y = 'Gene Expression')
# Boxplot CNA & Meth
p01_3 <- ggplot(trio_df, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. Methylation', subtitle = "(r= -0.2474198)",x = 'Copy Number Alternation', y = 'Methylation')
(p01_1)/(p01_2|p01_3)

# ==================================================================================================================================================================

# ==================================================================================================================================================================

trioRow=213938

trio_df <- data.frame(CNA=unlist(cna_filter %>% filter(cna.row==model_pos %>% filter(trio_row==trioRow) %>% pull(cna.row)) %>% select(-cna.row)),
                      Meth=unlist(meth_filter %>% filter(meth.row==model_pos %>% filter(trio_row==trioRow) %>% pull(meth.row)) %>% select(-meth.row)),
                      GE=unlist(exp_filter %>% filter(gene.row==model_pos %>% filter(trio_row==trioRow) %>% pull(gene.row)) %>% select(-gene.row)))

dim(trio_df)
cor(trio_df$GE,trio_df$Meth, use = "complete.obs")
cor(trio_df$CNA,trio_df$GE, use = "complete.obs")
cor(trio_df$CNA,trio_df$Meth, use = "complete.obs")

# Scatterplot GE & Meth
p01_1=ggplot(trio_df,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 Gene Expression v. Methylation',subtitle = "(r= -0.5194308)", x= 'Gene Expression', y = 'Methylation')
#Boxplot CNA & GE
p01_2 <- ggplot(trio_df, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. GE', subtitle = "(r= 0.1220043)",x = 'Copy Number Alternation', y = 'Gene Expression')
# Boxplot CNA & Meth
p01_3 <- ggplot(trio_df, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. Methylation', subtitle = "(r= --0.2077591)",x = 'Copy Number Alternation', y = 'Methylation')
(p01_1)/(p01_2|p01_3)

# ==================================================================================================================================================================


# ==================================================================================================================================================================

trioRow=213944

trio_df <- data.frame(CNA=unlist(cna_filter %>% filter(cna.row==model_pos %>% filter(trio_row==trioRow) %>% pull(cna.row)) %>% select(-cna.row)),
                      Meth=unlist(meth_filter %>% filter(meth.row==model_pos %>% filter(trio_row==trioRow) %>% pull(meth.row)) %>% select(-meth.row)),
                      GE=unlist(exp_filter %>% filter(gene.row==model_pos %>% filter(trio_row==trioRow) %>% pull(gene.row)) %>% select(-gene.row)))

dim(trio_df)
cor(trio_df$GE,trio_df$Meth, use = "complete.obs")
cor(trio_df$CNA,trio_df$GE, use = "complete.obs")
cor(trio_df$CNA,trio_df$Meth, use = "complete.obs")

# Scatterplot GE & Meth
p01_1=ggplot(trio_df,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 Gene Expression v. Methylation',subtitle = "(r= -0.4814745)", x= 'Gene Expression', y = 'Methylation')
#Boxplot CNA & GE
p01_2 <- ggplot(trio_df, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. GE', subtitle = "(r= 0.1220043)",x = 'Copy Number Alternation', y = 'Gene Expression')
# Boxplot CNA & Meth
p01_3 <- ggplot(trio_df, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. Methylation', subtitle = "(r= -0.2038931)",x = 'Copy Number Alternation', y = 'Methylation')
(p01_1)/(p01_2|p01_3)
## saved: 5*7 Portrait trio*_cor_scatter_plot.pdf
# ==================================================================================================================================================================

m1.2_interset <- model_pos_loc_full %>% filter(Inferred.Model=="M1.2",Inferred.Model.BH_fdr=="M1.2",Inferred.Model.qval=="M1.2",InferModel_Noconfounder=="M1.2")

m1.2_interset_filter <- m1.2_interset %>% select(trio_row,cna.row,gene.row,meth.row,Total.PC.Count,Confounders,Inferred.Model,Inferred.Model.BH_fdr,Inferred.Model.qval,InferModel_Noconfounder,Name,Gene.name,UCSC_RefGene_Group,Group,cor_CNA_Exp,cor_CNA_Meth,cor_Exp_Meth) %>% distinct()

df_m1.2 <- m1.2_interset_filter %>% select(-c(UCSC_RefGene_Group,Group)) %>% distinct()

df_interest_m1.2=df_m1.2 %>% filter(abs(cor_CNA_Exp) > abs(cor_CNA_Meth)) %>%
  arrange(desc(abs(cor_CNA_Exp) - abs(cor_CNA_Meth)))
### df_interest_m1.2 : 0 row
# ==================================================================================================================================================================
# ==================================================================================================================================================================
## Below are inferred to M1.1 without confounder, inferred to M1.2 with confounder effect
## 196 trios 
m1.2_interset2 <- model_pos_loc_full %>% filter(Inferred.Model=="M1.2",Inferred.Model.BH_fdr=="M1.2",Inferred.Model.qval=="M1.2",InferModel_Noconfounder=="M1.1")

m1.2_interset_filter2 <- m1.2_interset2 %>% select(trio_row,cna.row,gene.row,meth.row,Total.PC.Count,Confounders,Inferred.Model,Inferred.Model.BH_fdr,Inferred.Model.qval,InferModel_Noconfounder,Name,Gene.name,UCSC_RefGene_Group,Group,cor_CNA_Exp,cor_CNA_Meth,cor_Exp_Meth) %>% distinct()

df2_m1.2 <- m1.2_interset_filter2 %>% select(-c(UCSC_RefGene_Group,Group)) %>% distinct()

df2_interest_m1.2=df2_m1.2 %>% filter(abs(cor_CNA_Exp) > abs(cor_CNA_Meth)) %>%
  arrange(desc(abs(cor_CNA_Exp) - abs(cor_CNA_Meth)))

# ==================================================================================================================================================================

trioRow=89761
# cor(C_E)

trio_df <- data.frame(CNA=unlist(cna_filter %>% filter(cna.row==model_pos %>% filter(trio_row==trioRow) %>% pull(cna.row)) %>% select(-cna.row)),
                      Meth=unlist(meth_filter %>% filter(meth.row==model_pos %>% filter(trio_row==trioRow) %>% pull(meth.row)) %>% select(-meth.row)),
                      GE=unlist(exp_filter %>% filter(gene.row==model_pos %>% filter(trio_row==trioRow) %>% pull(gene.row)) %>% select(-gene.row)))

dim(trio_df)
cor(trio_df$GE,trio_df$Meth, use = "complete.obs") #0.1494763
cor(trio_df$CNA,trio_df$GE, use = "complete.obs") #0.2807822
cor(trio_df$CNA,trio_df$Meth, use = "complete.obs") #-0.04046334

# Scatterplot GE & Meth
p01_1=ggplot(trio_df,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.2 Gene Expression v. Methylation',subtitle = "(r= 0.149)", x= 'Gene Expression', y = 'Methylation')
#Boxplot CNA & GE
p01_2 <- ggplot(trio_df, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8))+labs(title = 'M1.2 CNA v. GE', subtitle = "(r= 0.281)",x = 'Copy Number Alternation', y = 'Gene Expression')
# Boxplot CNA & Meth
p01_3 <- ggplot(trio_df, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.2 CNA v. Methylation', subtitle = "(r= -0.040)",x = 'Copy Number Alternation', y = 'Methylation')
(p01_1)/(p01_2|p01_3)
## saved: 5*7 Portrait M1.2_trio89761_NoconfounderM1.1_cor_plot.pdf

# ==================================================================================================================================================================
## Below are inferred to M1.1 without confounder, inferred to M1.2 with confounder effect
## 164  trios 
m1.1_interset2 <- model_pos_loc_full %>% filter(Inferred.Model=="M1.1",Inferred.Model.BH_fdr=="M1.1",Inferred.Model.qval=="M1.1",InferModel_Noconfounder=="M1.2")

m1.1_interset_filter2 <- m1.1_interset2 %>% select(trio_row,cna.row,gene.row,meth.row,Total.PC.Count,Confounders,Inferred.Model,Inferred.Model.BH_fdr,Inferred.Model.qval,InferModel_Noconfounder,Name,Gene.name,UCSC_RefGene_Group,Group,cor_CNA_Exp,cor_CNA_Meth,cor_Exp_Meth) %>% distinct()

df2_m1.1 <- m1.1_interset_filter2 %>% select(-c(UCSC_RefGene_Group,Group)) %>% distinct()

df2_interest_m1.1=df2_m1.1 %>% filter(abs(cor_CNA_Exp) < abs(cor_CNA_Meth)) %>%
  arrange(desc(abs(cor_CNA_Meth) - abs(cor_CNA_Exp)))

# ==================================================================================================================================================================

trioRow=208243
# 

trio_df <- data.frame(CNA=unlist(cna_filter %>% filter(cna.row==model_pos %>% filter(trio_row==trioRow) %>% pull(cna.row)) %>% select(-cna.row)),
                      Meth=unlist(meth_filter %>% filter(meth.row==model_pos %>% filter(trio_row==trioRow) %>% pull(meth.row)) %>% select(-meth.row)),
                      GE=unlist(exp_filter %>% filter(gene.row==model_pos %>% filter(trio_row==trioRow) %>% pull(gene.row)) %>% select(-gene.row)))

dim(trio_df)
cor(trio_df$GE,trio_df$Meth, use = "complete.obs") #-0.2417344
cor(trio_df$CNA,trio_df$GE, use = "complete.obs") #0.00283551
cor(trio_df$CNA,trio_df$Meth, use = "complete.obs") #0.1824511

# Scatterplot GE & Meth
p01_1=ggplot(trio_df,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 Gene Expression v. Methylation',subtitle = "(r= -0.2417)", x= 'Gene Expression', y = 'Methylation')
#Boxplot CNA & GE
p01_2 <- ggplot(trio_df, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. GE', subtitle = "(r= 0.0028)",x = 'Copy Number Alternation', y = 'Gene Expression')
# Boxplot CNA & Meth
p01_3 <- ggplot(trio_df, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. Methylation', subtitle = "(r= 0.1825)",x = 'Copy Number Alternation', y = 'Methylation')
(p01_1)/(p01_2|p01_3)
## saved: 5*7 Portrait M1.1_trio208243_NoconfounderM1.2_cor_plot.pdf



# ==================================================================================================================================================================

# MRGN inferred raw and adjusted models are all M1.1
m1.1_interset <- model_pos_loc_full %>% filter(Inferred.Model=="M1.1",Inferred.Model.BH_fdr=="M1.1",Inferred.Model.qval=="M1.1")

m1.1_interset_filter <- m1.1_interset %>% select(trio_row,cna.row,gene.row,meth.row,Total.PC.Count,Confounders,Inferred.Model,Inferred.Model.BH_fdr,Inferred.Model.qval,InferModel_Noconfounder,Name,Gene.name,UCSC_RefGene_Group,Group,cor_CNA_Exp,cor_CNA_Meth,cor_Exp_Meth) %>% distinct()
df_m1.1 <- m1.1_interset_filter %>% select(-c(UCSC_RefGene_Group,Group)) %>% distinct()

table(df_m1.1$InferModel_Noconfounder)

# M0.1  M0.2  M1.1  M1.2  M2.1  M2.2    M3    M4 Other 
# 4540    48 11497   164   670    15  1070   680  3152 
# length(unique(df_m1.1$trio_row))
# 21864

# ==================================================================================================================================================================

# MRGN inferred raw and adjusted models are all M1.2
m1.2_interset <- model_pos_loc_full %>% filter(Inferred.Model=="M1.2",Inferred.Model.BH_fdr=="M1.2",Inferred.Model.qval=="M1.2")

m1.2_interset_filter <- m1.2_interset %>% select(trio_row,cna.row,gene.row,meth.row,Total.PC.Count,Confounders,Inferred.Model,Inferred.Model.BH_fdr,Inferred.Model.qval,InferModel_Noconfounder,Name,Gene.name,UCSC_RefGene_Group,Group,cor_CNA_Exp,cor_CNA_Meth,cor_Exp_Meth) %>% distinct()
df_m1.2 <- m1.2_interset_filter %>% select(-c(UCSC_RefGene_Group,Group)) %>% distinct()

table(df_m1.2$InferModel_Noconfounder)
# M0.1  M0.2  M1.1  M1.2  M2.1  M2.2    M3    M4 Other 
# 41   384   196  1399    49    90   116   264  1432

# length(unique(df_m1.2$trio_row))
# 3989

# ==================================================================================================================================================================
# trio_76: posER_brca_M0.1_M0.2_byMRGN_vs_M3_M4_by_baycn

trioRow=76
# cor(C_E)

## trio:76, MRGN_InferModel_raw = "M0.1",MRGN_InferModel.BH = "M0.1", baycn_InferModel="M3"
trio_df <- data.frame(CNA=unlist(cna_filter %>% filter(cna.row==model_pos %>% filter(trio_row==trioRow) %>% pull(cna.row)) %>% select(-cna.row)),
                      Meth=unlist(meth_filter %>% filter(meth.row==model_pos %>% filter(trio_row==trioRow) %>% pull(meth.row)) %>% select(-meth.row)),
                      GE=unlist(exp_filter %>% filter(gene.row==model_pos %>% filter(trio_row==trioRow) %>% pull(gene.row)) %>% select(-gene.row)))

dim(trio_df)
cor(trio_df$GE,trio_df$Meth, use = "complete.obs") #0.1928472
cor(trio_df$CNA,trio_df$GE, use = "complete.obs") #0.7435656
cor(trio_df$CNA,trio_df$Meth, use = "complete.obs") #0.2120911

# Scatterplot GE & Meth
p01_1=ggplot(trio_df,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M0.1 Gene Expression v. Methylation',subtitle = "(r= 0.1928)", x= 'Gene Expression', y = 'Methylation')
#Boxplot CNA & GE
p01_2 <- ggplot(trio_df, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8))+labs(title = 'M0.1 CNA v. GE', subtitle = "(r= 0.7436)",x = 'Copy Number Alternation', y = 'Gene Expression')
# Boxplot CNA & Meth
p01_3 <- ggplot(trio_df, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M0.1 CNA v. Methylation', subtitle = "(r= 0.2121)",x = 'Copy Number Alternation', y = 'Methylation')
(p01_1)/(p01_2|p01_3)

## Saved: /Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Figure/trio_76_posER_brca_M0.1_MRGN_baycn_M3.pdf

# ==================================================================================================================================================================

