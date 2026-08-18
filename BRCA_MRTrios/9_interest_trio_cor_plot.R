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
cor(dat_m0_1$GE,dat_m0_1$Meth)
cor(dat_m0_1$CNA,dat_m0_1$GE)
cor(dat_m0_1$CNA,dat_m0_1$Meth)

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
