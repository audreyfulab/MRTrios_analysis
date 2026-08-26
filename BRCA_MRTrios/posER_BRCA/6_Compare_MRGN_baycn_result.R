library(data.table)
library(tidyverse)

# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## models for BLCA and LIHC:
model_BLCA <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BLCA/Output_BLCA/BLCA_trio_Model_result_baycns_ALL_with_BH_fdr_qval_byLZ.txt")
model_LIHC <- fread("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_LIHC/Output_LIHC/LIHC_trio_Model_result_baycns_ALL_with_BH_fdr_qval_byLZ.txt")

table(model_BLCA$Inferred.Model.BH_fdr)/sum(table(model_BLCA$Inferred.Model.BH_fdr))
table(model_LIHC$Inferred.Model.BH_fdr)/sum(table(model_LIHC$Inferred.Model.BH_fdr))
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

## models for posER BRCA :

setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/raw_Data_Methyl")
cna         <- fread('data_CNA.txt', data.table = F)
gene.exp    <- fread('data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt', data.table = F)
TCGA.meth   <- fread("split.names.TCGA.meth.logit.txt", data.table = F)
trios       <- fread("trio.final.protein.coding.txt")
clinical    <- fread("names.pos.patient2.txt", header = FALSE)

# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## process cna, gene.exp, meth data:
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

## Save:
setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026")
fwrite(cna_filter,"pos_cna.txt")
fwrite(exp_filter,"pos_exp.txt")
fwrite(meth_filter,"pos_meth.txt")
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

##  read the saved filter cna,exp,meth data
# setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/output_Data_2026")
cna_filter <- fread("pos_cna.txt")
exp_filter <- fread("pos_exp.txt")
meth_filter <- fread("pos_meth.txt")

# ── Paths (defined once) ──────────────────────────────────────────────────────
DATA_DIR <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA"
#RAW_DIR  <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_analysis/raw_Data_Methyl"
#FIG_DIR  <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Figure"

# ── Load model result_MRGN (loaded once) ─────────────────────────────────────────
model_pos     <- fread(file.path(DATA_DIR, "Output_posER_BRCA/posER_BRCA_trio_Model_result_baycns_ALL_with_BH_fdr_qval_byLZ.txt"))
#model_neg     <- fread(file.path(DATA_DIR, "Output_negER_BRCA/negER_BRCA_trio_Model_result_baycns_ALL_with_BH_fdr_qval_byLZ.txt"))
#model_pos_loc <- readRDS(file.path(DATA_DIR, "Output_posER_BRCA/LZ_trio_result_baycns_pos_ALL_with_BH_fdr_qval_with_location.rds"))
#model_neg_loc <- readRDS(file.path(DATA_DIR, "Output_negER_BRCA/LZ_trio_result_baycns_neg_ALL_with_BH_fdr_qval_with_location.rds"))

# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## M1.1 exampe:  Figure :M1.1, trio: 135939
model_pos %>% filter(trio_row==135939)

trio_df <- data.frame(CNA=unlist(cna_filter %>% filter(cna.row==model_pos %>% filter(trio_row==135939) %>% pull(cna.row)) %>% select(-cna.row)),
                       Meth=unlist(meth_filter %>% filter(meth.row==model_pos %>% filter(trio_row==135939) %>% pull(meth.row)) %>% select(-meth.row)),
                       GE=unlist(exp_filter %>% filter(gene.row==model_pos %>% filter(trio_row==135939) %>% pull(gene.row)) %>% select(-gene.row)))

dim(trio_df) # 564,3
cor(trio_df$GE,trio_df$Meth) # -0.68
cor(trio_df$CNA,trio_df$GE)  # 0.08
cor(trio_df$CNA,trio_df$Meth)  # 0.02

# Scatterplot GE & Meth
p01_1=ggplot(trio_df,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 Gene Expression v. Methylation',subtitle = "(r=-0.6762893)", x= 'Gene Expression', y = 'Methylation')
#Boxplot CNA & GE
p01_2 <- ggplot(trio_df, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. GE', subtitle = "(r=0.08559773)",x = 'Copy Number Alternation', y = 'Gene Expression')
# Boxplot CNA & Meth
p01_3 <- ggplot(trio_df, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. Methylation', subtitle = "(r=0.02036697)",x = 'Copy Number Alternation', y = 'Methylation')
(p01_1)/(p01_2|p01_3)
ggsave("Pos_M1.1_trio135939.pdf")
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────




## ============ Compare baycn with MRGN ============================
## Load baycn result with 4312 sample trios:
result_baycn <- fread("/Users/lianzuo/Downloads/trios_sample_baycn.csv")
result_MRGN=head(model_pos,4213)

## models distribution:
sum(table(model_pos$Inferred.Model.BH_fdr)) # 292525
table(model_pos$Inferred.Model.BH_fdr)/292525

## confusion matrix:
a3 <- table(result_baycn$Inferred.Model, result_MRGN$Inferred.Model.BH_fdr)

         M0.1 M0.2 M1.1 M1.2 M2.1 M2.2  M3  M4 Other
  Error   16   60   11   19    6    1   3   0   509
  M0.1   808    4   74    0    9    0  99   5    67
  M0.2     1  111    0   16    0    1   2   0    97
  M1.1    74    2  163    3    9    2  15   9    19
  M1.2     0    8    0   19    0    2   0   0    26
  M2.1    70    3  156    0   75    4  28  18    27
  M2.2     0    9    3   14    0   16   4   0    18
  M3     241   27   35    5    9    0 529  85    38
  M4      25    7   98   14   42   49  78 188    28

a3_df <- as.data.frame.matrix(a3)
a3_df <- cbind(baycn_model = rownames(a3_df), a3_df)
fwrite(a3_df, "/Users/lianzuo/Desktop/posER_brca_baycn_MRGN_confusion_matrix.csv")

# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────


## ========= M1.1 by MRGN but M2.1 by baycn ===================
## Combine the two dataframe:
colnames(result_MRGN)[24]="Inferred.Model_raw"
df_result_combine=result_MRGN %>% left_join(result_baycn,by=c("trio_row","meth.row","cna.row","gene.row","Total.PC.Count"))
colnames(df_result_combine)[51]="baycn_Inferred.Model"

## filter to M1.1 by MRGN but M2.1 by baycn
a1=df_result_combine %>% filter(Inferred.Model.BH_fdr=="M1.1",baycn_Inferred.Model=="M2.1")
fwrite(a1,"trio_156samples_MRGN_M1.1_baycn_M2.1.csv")

## ========= Figure with two trio correlation plot =============
## check several examples:
### ======== trio: 15
Noconfounder_pos %>% filter(trio_row==15)
model_pos %>% filter(trio_row==15)

trioRow==15
trio_df <- data.frame(CNA=unlist(cna_filter %>% filter(cna.row==model_pos %>% filter(trio_row==trioRow) %>% pull(cna.row)) %>% select(-cna.row)),
                      Meth=unlist(meth_filter %>% filter(meth.row==model_pos %>% filter(trio_row==trioRow) %>% pull(meth.row)) %>% select(-meth.row)),
                      GE=unlist(exp_filter %>% filter(gene.row==model_pos %>% filter(trio_row==trioRow) %>% pull(gene.row)) %>% select(-gene.row)))

dim(trio_df)
cor(trio_df$GE,trio_df$Meth)
cor(trio_df$CNA,trio_df$GE)
cor(trio_df$CNA,trio_df$Meth)

# Scatterplot GE & Meth
p01_1=ggplot(trio_df,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 Gene Expression v. Methylation',subtitle = "(r=)", x= 'Gene Expression', y = 'Methylation')
#Boxplot CNA & GE
p01_2 <- ggplot(trio_df, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. GE', subtitle = "(r=)",x = 'Copy Number Alternation', y = 'Gene Expression')
# Boxplot CNA & Meth
p01_3 <- ggplot(trio_df, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. Methylation', subtitle = "(r=)",x = 'Copy Number Alternation', y = 'Methylation')
(p01_1)/(p01_2|p01_3)

### ======== trio: 99
model_pos %>% filter(trio_row==99)

trioRow=99
trio_df <- data.frame(CNA=unlist(cna_filter %>% filter(cna.row==model_pos %>% filter(trio_row==trioRow) %>% pull(cna.row)) %>% select(-cna.row)),
                      Meth=unlist(meth_filter %>% filter(meth.row==model_pos %>% filter(trio_row==trioRow) %>% pull(meth.row)) %>% select(-meth.row)),
                      GE=unlist(exp_filter %>% filter(gene.row==model_pos %>% filter(trio_row==trioRow) %>% pull(gene.row)) %>% select(-gene.row)))

dim(trio_df)
cor(trio_df$GE,trio_df$Meth)
cor(trio_df$CNA,trio_df$GE)
cor(trio_df$CNA,trio_df$Meth)

# Scatterplot GE & Meth
p01_1=ggplot(trio_df,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 Gene Expression v. Methylation',subtitle = "(r=)", x= 'Gene Expression', y = 'Methylation')
#Boxplot CNA & GE
p01_2 <- ggplot(trio_df, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. GE', subtitle = "(r=)",x = 'Copy Number Alternation', y = 'Gene Expression')
# Boxplot CNA & Meth
p01_3 <- ggplot(trio_df, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. Methylation', subtitle = "(r=)",x = 'Copy Number Alternation', y = 'Methylation')
(p01_1)/(p01_2|p01_3)

# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

## ============== M0.1 or M0.2 by MRGN but M3 or M4 by baycn ==================

# a2=posER_brca_M0.1_M0.2_byMRGN_vs_M3_M4_by_baycn.csv is the 300 subset of combined MRGN and baycn 4213 samples result : 
# MRGN classified the model (M0.1 or M0.2), while baycn independently classified as the model (M3 or M4)

a2=df_result_combine %>%
  filter(Inferred.Model.BH_fdr %in% c("M0.1", "M0.2")) %>% filter(
    baycn_Inferred.Model %in% c("M3", "M4"))

fwrite(a2, "/Users/lianzuo/Desktop/posER_brca_M0.1_M0.2_byMRGN_vs_M3_M4_by_baycn.csv")
a2=fread("/Users/lianzuo/Desktop/posER_brca_M0.1_M0.2_byMRGN_vs_M3_M4_by_baycn.csv")
table(a2$baycn_Inferred.Model)

# M3  M4 
# 268  32 

## ========= trios rerun but in baycn for 100,000 iterations ================
baycn <- fread("/Users/lianzuo/Downloads/trio.reruns.lian.csv")
table(baycn$Inferred.Model)

# Error  M0.1  M0.2  M1.1  M1.2  M2.1    M3    M4 
#    1    35     7     7     3     1   219    27 

### Cofusion matrix
table(a2$baycn_Inferred.Model,baycn$Inferred.Model)
    
#      Error M0.1  M0.2  M1.1  M1.2  M2.1   M3    M4
#  M3     1   35    7     4     0      0    215    6
#  M4     0    0    0     3     3      1     4    21

# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────







