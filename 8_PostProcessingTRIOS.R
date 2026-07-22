### week1 (10/31/24-11/6/2024) ###
# set work directory
setwd('/Users/lianzuo/LianZuo/WayneSate_Research/data_methylation')
# install packages and load the library
library(data.table)
# library(splitstackshape)
library(tidyverse)
#Loading in causal model interference of positive ER+(code reference: methylation_data)
df1 <-read.delim('model.trio.MRGN.all.posER.reclassify2.txt')

#Loading in list of TRIOS consisting of gene expression (GE), methylation and   (code reference: trios)
trios <- fread('trio.final.protein.coding.txt',data.table = F)

###For all trios unique genes
# uni_trios <- trios %>% mutate(comb=paste(Gene.name,gene.row,sep = "_")) %>% distinct(comb,Gene.name,gene.row)
# uni_GEname <- trios %>% group_by(Gene.name) %>% summarise(count=n()) %>% filter(count>1)
### KRTAP5-1 (count:2) ,PRR20E(count:4),  QSOX1 (count:2) ###
# uni_generow <- unige %>% group_by(gene.row) %>% summarise(count=n()) %>% filter(count>1)
### 5 gene.row: 1573,6055,9350,12731,14901,17732, with count of 2.###

#Finding overall unique genes within ER+ Data
positive_genes <- trios[df1$trio_row,]
# results <-unique(positive_genes$Gene.name)
# length(results) # = 16331
# length(unique(positive_genes$gene.row)) #16328? (off by 3)

## ER+ move gene.row= c(6482,14474,1573,6055,14901,17732))


# add a new character column which combine the Gene.name and gene.row,remove duplicate rows based on the new column;
positive_genes_unique <- positive_genes %>% mutate(comb=paste(Gene.name,gene.row,sep = "_")) %>% distinct(comb,Gene.name,gene.row)
dim(positive_genes_unique)

#extract the count of gene.name >1 
extra_Gene.name <- positive_genes_unique %>% group_by(Gene.name) %>% summarise(count=n()) %>% filter(count>1)
extra_Gene.name

#  x <- positive_genes[,c(1,4)]
#  x2 <- x%>%mutate(comb=paste(Gene.name, gene.row, sep="_"))%>%distinct(comb, .keep_all=T)
# dim(x2)
# dim(x)
# head(x2)

#extract the count of gene.row >1 
extra_gene.row <- positive_genes_unique %>% group_by(gene.row) %>% summarise(count=n()) %>% filter(count>1)
extra_gene.row

# summ <- x2%>%group_by(Gene.name)%>%summarise(ny=n(), .groups="drop")
# dim(summ)
# summ%>%filter(ny>1)
#  summ <- x2%>%group_by(Gene.name)%>%summarise(ny=n(), .groups="drop")%>%ungroup()
# summ%>%filter(ny>1)%>%as.data.frame()
# length(unique(x$comb))
# length(unique(x2$gene.row))
# summ2 <- x2%>%group_by(gene.row)%>%summarise(ny=n(), .groups = "drop")%>%ungroup()
# summ2%>%filter(ny>1)

#find the unique gene
extra_Gename <- positive_genes_unique %>% filter(Gene.name=="QSOX1")
extra_Gename

extra_Gerow_1 <- positive_genes_unique %>% filter(gene.row==1573)
extra_Gerow_2 <- positive_genes_unique %>% filter(gene.row==6055)
extra_Gerow_3 <- positive_genes_unique %>% filter(gene.row==14901)
extra_Gerow_4 <- positive_genes_unique %>% filter(gene.row==17732)
extra_Gerow_1
extra_Gerow_2
extra_Gerow_3
extra_Gerow_4


index_pos <- which(positive_genes$gene.row %in% c(6482,14474,1573,6055,14901,17732))
length(index_pos)
uni_pos_gene <- positive_genes[-index_pos,]


### week1:11/5/2024 ###

##################### Pos ER+ #############################################
# set work directory
setwd('/Users/lianzuo/LianZuo/WayneSate_Research/data_methylation')
#Loading in raw data for TRIOS (use 'fread' to load in large & heavily filtered datasets)
CNA_rawdata <- fread('data_CNA.txt',data.table = F) #CNA
GE_rawdata <- fread('data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt',data.table = F) #gene_expression
meth_rawdata <- fread('split.names.TCGA.meth.logit.txt',data.table = F) #methylation

#Loading in all intermediate files
##ER+ intermediate files (PCA files will have genes as rows instead of columns)
ERpos_patients <- read.delim('names.pos.patient2.txt', sep = ',',header = F)#ER+ patient list
pos_GE_PCA <- read.delim('PCA.gene.exp.posER.txt',row.names = 1) #gene expression ER+ PCA scores
pos_meth_PCA <- read.delim('PCA.meth.posER.txt',row.names = 1) #methylation ER+ PCA scores

#GOAL: Generate scatter plots from causal inferences (code reference: "methylation_data") by aligning patient's TRIOS
#PURPOSE: finding each patient in each trios raw dataset by matching first column of PCA scores in ER+/- to column names of raw data which contain patients 

#Avoiding errors of "list" by turning into data frame
#GE_rawdata_tmp <- as.data.frame(GE_rawdata) #made with Dr. Fu
#ME_rawdata_tmp <- as.data.frame(meth_rawdata)
#CNA_rawdata_tmp <- as.data.frame(CNA_rawdata)

#Extracting the row names of ER+ methylation PCA intermediate file 
pos_meth_PCA_row <- rownames(pos_meth_PCA)

#Matching columns of GE_raw data to extracted row names of ER+ methylation PCA intermediate file
final_matching <- match(pos_meth_PCA_row, colnames(GE_rawdata))
#final_matching
GE_ERpos <- GE_rawdata[,final_matching]
# GE_ERpos_1 <- GE_rawdata_tmp[,pos_GE_PCA_row]

#Extracting the row names of ER+ GE PCA intermediate file
pos_GE_PCA_row <- rownames(pos_GE_PCA)

#Matching columns of GE_raw data to extracted row names of ER+ GE PCA intermediate file 
matching <- match(pos_GE_PCA_row , colnames(GE_rawdata))
#matching
ER_pos <- GE_rawdata[,matching]

##OBSERVATION: using the GE and methylation PCA intermediate files will give you the same result
#Code below proves this point:
table(pos_meth_PCA_row == pos_GE_PCA_row)  #all elements are TRUE
##CONCLUSION: do not have to match all raw datas to both GE and methylation PCA files



######
#Match the row names of meth_rawdata with the row names of ER+ methylation PCA 
more_matching <- match(pos_meth_PCA_row, colnames(meth_rawdata))
more_ER_pos <- meth_rawdata[, more_matching] 


#Match the column names of CNA_rawdata with the row names of pos_meth_PCA 
moree_matchingg <- match(pos_meth_PCA_row, colnames(CNA_rawdata))

##Why it did not work?:
#Check "moree_matchingg" vector with code: moree_matchingg[1:10]
#Conclusion: matching function is not a fan of "NA" and dimensions of GE_rawdata & CNA rawdata are different

##Question: How many NAs/Which NAs do we remove from CNA_rawdata?
#Code: sum(is.na(moree_matchingg))  = 8

#Revising pos_meth_PCA_row for CNA
revised_patients <- pos_meth_PCA_row[ -which (is.na (moree_matchingg))] #basically keeps "first_column" - all the NAs

#Matching extracted row names of ER- methylation PCA intermediate file to columns of CNA_rawdata
moree_matchingg <- match(revised_patients, colnames(CNA_rawdata))
positive_CNA <- CNA_rawdata[,moree_matchingg]


#Now the number of patients is different between each raw data set (572 variables v. 564 variables). We must eliminate the same eight patients for all matching scenarios

#NEW matching with GE data using "revised_patients"
matching <-match(revised_patients, colnames(GE_rawdata))
positive_GE<- GE_rawdata[,matching]

#NEW matching with meth raw data using "revised_patients"
more_matching <- match(revised_patients, colnames(meth_rawdata))
positive_meth <- meth_rawdata[, more_matching] 

###For m0.1(C -> E):###
#Fliter the M0.1 model datasets
mod0_1 <- df1 %>% filter(Inferred.Model=="M0.1",Inferred.Model2=="M0.1",Inferred.Model3=="M0.1")
head(mod0_1,3)
#Running the trio_row number#270001
trios[270001,]
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
#make sure these values are not in "lists"
#Find the correlation value of 2 varibales
dat_m0_1 <- data.frame(CNA=unlist(positive_CNA[264,]),Meth=unlist(positive_meth[462699,]),GE=unlist(positive_GE[1197,]))
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


###For m0.2(C -> M):###
#Fliter the M0.1 model datasets
mod0_2 <- df1 %>% filter(Inferred.Model=="M0.2",Inferred.Model2=="M0.2",Inferred.Model3=="M0.2")
head(mod0_2,3)
#Running the trio_row number#270025
trios[270025,]
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
#make sure these values are not in "lists"
#Find the correlation value of 2 varibales
dat_m0_2 <- data.frame(CNA=unlist(positive_CNA[20630,]),Meth=unlist(positive_meth[252319,]),GE=unlist(positive_GE[5411,]))
dim(dat_m0_2)
cor(dat_m0_2$GE,dat_m0_2$Meth)
cor(dat_m0_2$CNA,dat_m0_2$GE)
cor(dat_m0_2$CNA,dat_m0_2$Meth)


# Scatterplot GE & Meth
p02_1 <- ggplot(dat_m0_2,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M0.2 Gene Expression v. Methylation',subtitle = "(r = - 0.002)", x= 'Gene Expression', y = 'Methylation')
#Boxplot CNA & GE
p02_2 <- ggplot(dat_m0_2, aes(x=factor(CNA), y=GE))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M0.2 CNA v. GE',subtitle = "(r = - 0.116)", x = 'Copy Number Alternation', y = 'Gene Expression')
# Boxplot CNA & Meth
p02_3 <- ggplot(dat_m0_2, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M0.2 CNA v. Methylation',subtitle = "(r = - 0.126)", x = 'Copy Number Alternation', y= 'Methylation')
(p02_1)/(p02_2|p02_3)
ggsave("Pos_M0.2.pdf")

#For m1.1(C -> E -> M):
#Fliter the M1.1 model datasets
mod1_1 <- df1 %>% filter(Inferred.Model=="M1.1",Inferred.Model2=="M1.1",Inferred.Model3=="M1.1")
head(mod1_1,3)
#Running the trio_row number#270016
trios[270016,]
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
#make sure these values are not in "lists"
#Find the correlation value of 2 varibales
dat_m1_1 <- data.frame(CNA=unlist(positive_CNA[15092,]),Meth=unlist(positive_meth[442663,]),GE=unlist(positive_GE[9100,]))
dim(dat_m1_1)
cor(dat_m1_1$GE,dat_m1_1$Meth)
cor(dat_m1_1$CNA,dat_m1_1$GE)
cor(dat_m1_1$CNA,dat_m1_1$Meth)

# Scatter plot GE & Meth
p11_1 <- ggplot(dat_m1_1,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 Gene Expression v. Methylation',subtitle = "(r = NA)", x= 'Gene Expression', y = 'Methylation')

#Boxplot CNA & GE
p11_2 <- ggplot(dat_m1_1, aes(x=factor(CNA), y=GE))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. GE',subtitle = "(r = 0.217)", x = 'Copy Number Alternation', y = 'Gene Expression')

# Boxplot CNA & Meth
p11_3 <- ggplot(dat_m1_1, aes(x=factor(CNA), y=Meth ))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. Methylation',subtitle = "(r = NA)", x = 'Copy Number Alternation', y = 'Methylation')
(p11_1)/(p11_2|p11_3)
ggsave("Pos_M1.1.pdf")

###For m1.2(C -> M -> E):###
#Fliter the M1.2 model datasets
mod1_2 <- df1 %>% filter(Inferred.Model=="M1.2",Inferred.Model2=="M1.2",Inferred.Model3=="M1.2")
head(mod1_2,3)
#Running the trio_row number#270023
trios[270023,]
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
#make sure these values are not in "lists"
#Find the correlation value of 2 varibales
dat_m1_2 <- data.frame(CNA=unlist(positive_CNA[20630,]),Meth=unlist(positive_meth[222723,]),GE=unlist(positive_GE[5411,]))
dim(dat_m1_2)

cor(dat_m1_2$GE,dat_m1_2$Meth)
cor(dat_m1_2$CNA,dat_m1_2$GE)
cor(dat_m1_2$CNA,dat_m1_2$Meth)

# Scatter plot GE & Meth
p12_1 <- ggplot(dat_m1_2,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.2 Gene Expression v. Methylation',subtitle = "(r = 0.317)", x= 'Gene Expression', y = 'Methylation')

#Boxplot CNA & GE
p12_2 <- ggplot(dat_m1_2, aes(x=factor(CNA), y=GE))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.2 CNA v. GE',subtitle = "(r = -0.116)", x = 'Copy Number Alternation', y = 'Gene Expression')

# Boxplot CNA & Meth
p12_3 <- ggplot(dat_m1_2, aes(x=factor(CNA), y=Meth ))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.2 CNA v. Methylation',subtitle = "(r = 0.018)", x = 'Copy Number Alternation', y = 'Methylation')
(p12_1)/(p12_2|p12_3)
ggsave("Pos_M1.2.pdf")




###For m2.1(C -> E <- M):###
#Fliter the M0.1 model datasets
mod2_1 <- df1 %>% filter(Inferred.Model=="M2.1",Inferred.Model2=="M2.1",Inferred.Model3=="M2.1")
head(mod2_1,3)
#Running the trio_row number#270001
trios[270510,]
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
#make sure these values are not in "lists"
#Find the correlation value of 2 varibales
dat_m2_1 <- data.frame(CNA=unlist(positive_CNA[13439,]),Meth=unlist(positive_meth[481639,]),GE=unlist(positive_GE[3194,]))
dim(dat_m2_1)

cor(dat_m2_1$GE,dat_m2_1$Meth)
cor(dat_m2_1$CNA,dat_m2_1$GE)
cor(dat_m2_1$CNA,dat_m2_1$Meth)

# Scatter plot GE & Meth
p21_1 <- ggplot(dat_m2_1,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M2.1 Gene Expression v. Methylation',subtitle = "(r = -0.010)", x= 'Gene Expression', y = 'Methylation')

#Boxplot CNA & GE
p21_2 <- ggplot(dat_m2_1, aes(x=factor(CNA), y=GE))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M2.1 CNA v. GE',subtitle = "(r = 0.379)", x = 'Copy Number Alternation', y = 'Gene Expression')

# Boxplot CNA & Meth
p21_3 <- ggplot(dat_m2_1, aes(x=factor(CNA), y=Meth ))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M2.1 CNA v. Methylation',subtitle = "(r = -0.042)", x = 'Copy Number Alternation', y = 'Methylation')
(p21_1)/(p21_2|p21_3)
ggsave("Pos_M2.1.pdf")


###11/7/24###
###For m2.2(C -> M <- E):###
#Fliter the M0.1 model datasets
mod2_2 <- df1 %>% filter(Inferred.Model=="M2.2",Inferred.Model2=="M2.2",Inferred.Model3=="M2.2")
head(mod2_2,3)
#Running the trio_row number#271127
trios[271127,]
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
#make sure these values are not in "lists"
#Find the correlation value of 2 varibales
dat_m2_2 <- data.frame(CNA=unlist(positive_CNA[16315,]),Meth=unlist(positive_meth[422513,]),GE=unlist(positive_GE[18832,]))
dim(dat_m2_2)

cor(dat_m2_2$GE,dat_m2_2$Meth)
cor(dat_m2_2$CNA,dat_m2_2$GE)
cor(dat_m2_2$CNA,dat_m2_2$Meth)

# Scatter plot GE & Meth
p22_1 <- ggplot(dat_m2_2,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M2.2 Gene Expression v. Methylation',subtitle = "(r = 0.043)", x= 'Gene Expression', y = 'Methylation')

#Boxplot CNA & GE
p22_2 <- ggplot(dat_m2_2, aes(x=factor(CNA), y=GE))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M2.2 CNA v. GE',subtitle = "(r = 0.078)", x = 'Copy Number Alternation', y = 'Gene Expression')

# Boxplot CNA & Meth
p22_3 <- ggplot(dat_m2_2, aes(x=factor(CNA), y=Meth ))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M2.2 CNA v. Methylation',subtitle = "(r = -0.221)", x = 'Copy Number Alternation', y = 'Methylation')
(p22_1)/(p22_2|p22_3)
ggsave("Pos_M2.2.pdf")

##For m3 model(E ← C → M):###
#Filter the m3 model datasets
mod3 <- df1 %>% filter(Inferred.Model=='M3',Inferred.Model2=='M3',Inferred.Model3 == 'M3')
dim(mod3)
head(mod3)
mod3[100:105,]
#Running the trio_row number#270700
trios[270700,]
dat_m3 <-  data.frame(CNA=unlist(positive_CNA[16443,]),Meth=unlist(positive_meth[146234,]),GE=unlist(positive_GE[17015,]))
dim(dat_m3)

cor(dat_m3$GE,dat_m3$Meth)
cor(dat_m3$CNA,dat_m3$GE)
cor(dat_m3$CNA,dat_m3$Meth)

# Scatter plot GE & Meth
p3_1 <- ggplot(dat_m3,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M3 Gene Expression v. Methylation',subtitle = "(r = -0.235)", x= 'Gene Expression', y = 'Methylation')

#Boxplot CNA & GE
p3_2 <- ggplot(dat_m3, aes(x=factor(CNA), y=GE))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M3 CNA v. GE',subtitle = "(r = 0.338)", x = 'Copy Number Alternation', y = 'Gene Expression')

# Boxplot CNA & Meth
p3_3 <- ggplot(dat_m3, aes(x=factor(CNA), y=Meth ))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M3 CNA v. Methylation',subtitle = "(r = -0.209)", x = 'Copy Number Alternation', y = 'Methylation')
(p3_1)/(p3_2|p3_3)
ggsave("Pos_M3.pdf")



### M4 model()###
#Filter the m4 model datasets
mod4 <-  df1 %>% filter(Inferred.Model == 'M4',Inferred.Model2 == 'M4',Inferred.Model3 == 'M4')
dim(mod4)
mod4[50:53,]
trios[270550,]
dat_m4 <-  data.frame(CNA=unlist(positive_CNA[23279,]),GE=unlist(positive_GE[18505,]),Meth=unlist(positive_meth[367555,]))
dim(dat_m4)
head(dat_m4)

cor(dat_m4$GE,dat_m4$Meth)
cor(dat_m4$CNA,dat_m4$GE)
cor(dat_m4$CNA,dat_m4$Meth)

# Scatter plot GE & Meth
p4_1 <- ggplot(dat_m4,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M4 Gene Expression v. Methylation',subtitle = "(r = NA)", x= 'Gene Expression', y = 'Methylation')

#Boxplot CNA & GE
p4_2 <- ggplot(dat_m4, aes(x=factor(CNA), y=GE))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M4 CNA v. GE',subtitle = "(r = 0.369)", x = 'Copy Number Alternation', y = 'Gene Expression')

# Boxplot CNA & Meth
p4_3 <- ggplot(dat_m4, aes(x=factor(CNA), y=Meth ))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M4 CNA v. Methylation',subtitle = "(r = NA)", x = 'Copy Number Alternation', y = 'Methylation')
(p4_1)/(p4_2|p4_3)
ggsave("Pos_M4.pdf")


### Other model###
modOther <-  df1 %>% filter(Inferred.Model == 'Other',Inferred.Model2 == 'Other',Inferred.Model3 == 'Other')
dim(modOther)
head(modOther)
trios[270027,]
dat_Other <-  data.frame(CNA=unlist(positive_CNA[20630,]),GE=unlist(positive_GE[5411,]),Meth=unlist(positive_meth[275322,]))
dim(dat_Other)
head(dat_Other)

#cor(GE,Meth) =0.267855
#cor(GE,CNA) = -0.1155507
#cor(Meth,CNA) = 0.002762125


cor(dat_Other$GE,dat_Other$Meth)
cor(dat_Other$CNA,dat_Other$GE)
cor(dat_Other$CNA,dat_Other$Meth)

# Scatter plot GE & Meth
p_other_1 <- ggplot(dat_Other,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'Other Gene Expression v. Methylation',subtitle = "(r = 0.268)", x= 'Gene Expression', y = 'Methylation')

#Boxplot CNA & GE
p_other_2 <- ggplot(dat_Other, aes(x=factor(CNA), y=GE))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'Other CNA v. GE',subtitle = "(r = - 0.116)", x = 'Copy Number Alternation', y = 'Gene Expression')

# Boxplot CNA & Meth
p_other_3 <- ggplot(dat_Other, aes(x=factor(CNA), y=Meth ))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'Other CNA v. Methylation',subtitle = "(r = 0.003)", x = 'Copy Number Alternation', y = 'Methylation')
(p_other_1)/(p_other_2|p_other_3)
ggsave("Pos_M_Other.pdf")


##################### Neg ER- #############################################
# install packages and load the library
library(data.table)
library(tidyverse)
df2 <-read.delim('model.trio.MRGN.all.negER.reclassify2.txt')

#Loading in list of TRIOS consisting of gene expression (GE), methylation and   (code reference: trios)
trios <- data.frame(fread('trio.final.protein.coding.txt'))
#Finding overall unique genes within ER- Data
neg_genes <- trios[df2$trio_row,]
# add a new character column which combine the Gene.name and gene.row,remove duplicate rows based on the new column;
neg_genes_unique <- neg_genes %>% mutate(comb=paste(Gene.name,gene.row,sep = "_")) %>% distinct(comb,Gene.name,gene.row)
dim(neg_genes_unique)

#extract the count of gene.name >1 
neg_extra_Gene.name <- neg_genes_unique %>% group_by(Gene.name) %>% summarise(count=n()) %>% filter(count>1)
neg_extra_Gene.name

#extract the count of gene.row >1 
neg_extra_gene.row <- neg_genes_unique %>% group_by(gene.row) %>% summarise(count=n()) %>% filter(count>1)
neg_extra_gene.row

#find the unique gene
neg_extra_Gename <- neg_genes_unique %>% filter(Gene.name=="QSOX1")
neg_extra_Gename

neg_extra_Gerow_1 <- neg_genes_unique %>% filter(gene.row==1573)
neg_extra_Gerow_2 <- neg_genes_unique %>% filter(gene.row==6055)
neg_extra_Gerow_3 <- neg_genes_unique %>% filter(gene.row==14901)
neg_extra_Gerow_4 <- neg_genes_unique %>% filter(gene.row==17732)
neg_extra_Gerow_1
neg_extra_Gerow_2
neg_extra_Gerow_3
neg_extra_Gerow_4
### week1:11/12/2024 ###


#Loading in ER- intermediate files
ERneg_patients <- read.delim('names.neg.patient2.txt', sep = ',') #ER- patient list
neg_GE_PCA <- read.delim('PCA.gene.exp.negER.txt', row.names = 1) #gene expression ER- PCA scores
neg_meth_PCA <-read.delim('PCA.meth.negER.txt', row.names = 1) #methylation ER- PCA scores


#GOAL: Generate scatter plots from causal inferences (code reference: "methylation_data") by aligning patient's TRIOS
#PURPOSE: finding each patient in each trios raw dataset by matching first column of PCA scores in ER- to column names of raw data which contain patients 

# #Avoiding errors of "list" by turning into data frame
# GE_rawdata_tmp <- as.data.frame(GE_rawdata) #made with Dr. Fu
# ME_rawdata_tmp <- as.data.frame(meth_rawdata)
# CNA_rawdata_tmp <- as.data.frame(CNA_rawdata)

#Extracting the row names of ER+ methylation PCA intermediate file 
neg_meth_PCA_row <- rownames(neg_meth_PCA)

#Matching columns of GE_raw data to extracted row names of ER+ methylation PCA intermediate file
neg_final_matching <- match(neg_meth_PCA_row, colnames(GE_rawdata))
#final_matching
GE_ERneg <- GE_rawdata[,neg_final_matching]
# GE_ERpos_1 <- GE_rawdata_tmp[,pos_GE_PCA_row]

#Extracting the row names of ER+ GE PCA intermediate file
neg_GE_PCA_row <- rownames(neg_GE_PCA)

#Matching columns of GE_raw data to extracted row names of ER+ GE PCA intermediate file 
neg_matching <- match(neg_GE_PCA_row , colnames(GE_rawdata))
#matching
ER_neg <- GE_rawdata[,matching]

##OBSERVATION: using the GE and methylation PCA intermediate files will give you the same result
#Code below proves this point:
table(neg_meth_PCA_row == neg_GE_PCA_row)  #all elements are TRUE
##CONCLUSION: do not have to match all raw datas to both GE and methylation PCA files


######
#Match the row names of meth_rawdata with the row names of ER+ methylation PCA 
neg_more_matching <- match(neg_meth_PCA_row, colnames(meth_rawdata))
more_ER_neg <- meth_rawdata[, neg_more_matching] 


#Match the column names of CNA_rawdata with the row names of pos_meth_PCA 
neg_moree_matchingg <- match(neg_meth_PCA_row, colnames(CNA_rawdata))

##Why it did not work?:
#Check "moree_matchingg" vector with code: moree_matchingg[1:10]
#Conclusion: matching function is not a fan of "NA" and dimensions of GE_rawdata & CNA rawdata are different

##Question: How many NAs/Which NAs do we remove from CNA_rawdata?
#Code: sum(is.na(moree_matchingg))  = 8

#Revising pos_meth_PCA_row for CNA
neg_revised_patients <- neg_meth_PCA_row[ -which (is.na (neg_moree_matchingg))] #basically keeps "first_column" - all the NAs

#Matching extracted row names of ER- methylation PCA intermediate file to columns of CNA_rawdata
neg_moree_matchingg2 <- match(neg_revised_patients, colnames(CNA_rawdata))
neg_CNA <- CNA_rawdata[,neg_moree_matchingg2]


#Now the number of patients is different between each raw data set (572 variables v. 564 variables). We must eliminate the same eight patients for all matching scenarios

#NEW matching with GE data using "revised_patients"
neg_matching2 <-match(neg_revised_patients, colnames(GE_rawdata))
neg_GE<- GE_rawdata[,neg_matching2]

#NEW matching with meth raw data using "revised_patients"
neg_more_matching2 <- match(neg_revised_patients, colnames(meth_rawdata))
neg_meth <- meth_rawdata[, neg_more_matching2] 

###For m0.1(C -> E):###
#Fliter the M0.1 model datasets
neg_mod0_1 <- df2 %>% filter(Inferred.Model=="M0.1",Inferred.Model2=="M0.1",Inferred.Model3=="M0.1")
head(neg_mod0_1,3)
#Running the trio_row number#270001
trios[270001,]
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
#make sure these values are not in "lists"
#Find the correlation value of 2 varibales
neg_dat_m0_1 <- data.frame(CNA=unlist(neg_CNA[264,]),Meth=unlist(neg_meth[462699,]),GE=unlist(neg_GE[1197,]))
dim(neg_dat_m0_1)
cor(neg_dat_m0_1$GE,neg_dat_m0_1$Meth)
cor(neg_dat_m0_1$CNA,neg_dat_m0_1$GE)
cor(neg_dat_m0_1$CNA,neg_dat_m0_1$Meth)

##ggplot2##

# Scatterplot GE & Meth
neg_p01_1=ggplot(neg_dat_m0_1,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = '(ER-) M0.1 GE v. Meth',subtitle = "(r=0.130)", x= 'Gene Expression', y = 'Methylation')
#Boxplot CNA & GE
neg_p01_2 <- ggplot(neg_dat_m0_1, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M0.1 CNA v. GE', subtitle = "(r=0.405)",x = 'Copy Number Alternation', y = 'Gene Expression')+theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))

# Boxplot CNA & Meth
neg_p01_3 <- ggplot(neg_dat_m0_1, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M0.1 CNA v. Meth', subtitle = "(r=0.144)",x = 'CNA', y = 'Methylation')+theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
(neg_p01_1)/(neg_p01_2|neg_p01_3)
ggsave("Neg_M0.1.pdf")



###For m0.2(C -> M):###
#Fliter the M0.1 model datasets
neg_mod0_2 <- df2 %>% filter(Inferred.Model=="M0.2",Inferred.Model2=="M0.2",Inferred.Model3=="M0.2")
head(neg_mod0_2,3)
#Running the trio_row number#270025
trios[270025,]
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
#make sure these values are not in "lists"
#Find the correlation value of 2 varibales
neg_dat_m0_2 <- data.frame(CNA=unlist(neg_CNA[20630,]),Meth=unlist(neg_meth[252319,]),GE=unlist(neg_GE[5411,]))
dim(neg_dat_m0_2)
cor(neg_dat_m0_2$GE,neg_dat_m0_2$Meth)
cor(neg_dat_m0_2$CNA,neg_dat_m0_2$GE)
cor(neg_dat_m0_2$CNA,neg_dat_m0_2$Meth)

# Scatterplot GE & Meth
neg_p02_1=ggplot(neg_dat_m0_2,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = '(ER-) M0.2 GE v. Meth',subtitle = "(r = -0.098)", x= 'Gene Expression', y = 'Methylation')

#Boxplot CNA & GE
neg_p02_2 <- ggplot(neg_dat_m0_2, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M0.2 CNA v. GE', subtitle = "(r = 0.044)",x = 'Copy Number Alternation', y = 'Gene Expression')+theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))

# Boxplot CNA & Meth
neg_p02_3 <- ggplot(neg_dat_m0_2, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M0.2 CNA v. Meth', subtitle = "(r = -0.288)",x = 'CNA', y = 'Methylation')+theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))

(neg_p02_1)/(neg_p02_2|neg_p02_3)
ggsave("Neg_M0.2.pdf")



#For m1.1(C -> E -> M):
#Fliter the M1.1 model datasets
neg_mod1_1 <- df2 %>% filter(Inferred.Model=="M1.1",Inferred.Model2=="M1.1",Inferred.Model3=="M1.1")
head(neg_mod1_1,3)
#Running the trio_row number#270284
trios[270284,]
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
#make sure these values are not in "lists"
#Find the correlation value of 2 varibales
neg_dat_m1_1 <- data.frame(CNA=unlist(neg_CNA[7607,]),Meth=unlist(neg_meth[157650,]),GE=unlist(neg_GE[12220,]))
dim(neg_dat_m1_1)
cor(neg_dat_m1_1$GE,neg_dat_m1_1$Meth)
cor(neg_dat_m1_1$CNA,neg_dat_m1_1$GE)
cor(neg_dat_m1_1$CNA,neg_dat_m1_1$Meth)

# Scatterplot GE & Meth
neg_p11_1=ggplot(neg_dat_m1_1,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = '(ER-) M1.1 GE v. Meth',subtitle = "(r = -0.14)", x= 'Gene Expression', y = 'Methylation')

#Boxplot CNA & GE
neg_p11_2 <- ggplot(neg_dat_m1_1, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. GE', subtitle = "(r = 0.54)",x = 'Copy Number Alternation', y = 'Gene Expression')+theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))

# Boxplot CNA & Meth
neg_p11_3 <- ggplot(neg_dat_m1_1, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.1 CNA v. Meth', subtitle = "(r = -0.25)",x = 'CNA', y = 'Methylation')+theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))

(neg_p11_1)/(neg_p11_2|neg_p11_3)
ggsave("Neg_M1.1.pdf")



###For m1.2(C -> M -> E):###
#Fliter the M1.2 model datasets
neg_mod1_2 <- df2 %>% filter(Inferred.Model=="M1.2",Inferred.Model2=="M1.2",Inferred.Model3=="M1.2")
head(neg_mod1_2,3)
#Running the trio_row number#270381
trios[270381,]
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
#make sure these values are not in "lists"
#Find the correlation value of 2 varibales
neg_dat_m1_2 <- data.frame(CNA=unlist(neg_CNA[7605,]),Meth=unlist(neg_meth[342302,]),GE=unlist(neg_GE[2942,]))
dim(neg_dat_m1_2)

cor(neg_dat_m1_2$GE,neg_dat_m1_2$Meth)
cor(neg_dat_m1_2$CNA,neg_dat_m1_2$GE)
cor(neg_dat_m1_2$CNA,neg_dat_m1_2$Meth)

# Scatter plot GE & Meth
neg_p12_1 <- ggplot(neg_dat_m1_2,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = '(ER-) M1.2 GE v. Methylation',subtitle = "(r = -0.35)", x= 'Gene Expression', y = 'Methylation')

#Boxplot CNA & GE
neg_p12_2 <- ggplot(neg_dat_m1_2, aes(x=factor(CNA), y=GE))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.2 CNA v. GE',subtitle = "(r = 0.37)", x = 'Copy Number Alternation', y = 'Gene Expression')+theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))

# Boxplot CNA & Meth
neg_p12_3 <- ggplot(neg_dat_m1_2, aes(x=factor(CNA), y=Meth ))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M1.2 CNA v. Meth',subtitle = "(r = -0.44)", x = 'CNA', y = 'Methylation')+theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))

(neg_p12_1)/(neg_p12_2|neg_p12_3)
ggsave("Neg_M1.2.pdf")




###For m2.1(C -> E <- M):###
#Fliter the M0.1 model datasets
neg_mod2_1 <- df2 %>% filter(Inferred.Model=="M2.1",Inferred.Model2=="M2.1",Inferred.Model3=="M2.1")
head(neg_mod2_1,3)
#Running the trio_row number#270118
trios[270118,]
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
#make sure these values are not in "lists"
#Find the correlation value of 2 varibales
neg_dat_m2_1 <- data.frame(CNA=unlist(neg_CNA[20855,]),Meth=unlist(neg_meth[59464,]),GE=unlist(neg_GE[20245,]))
dim(neg_dat_m2_1)

cor(neg_dat_m2_1$GE,neg_dat_m2_1$Meth)
cor(neg_dat_m2_1$CNA,neg_dat_m2_1$GE)
cor(neg_dat_m2_1$CNA,neg_dat_m2_1$Meth)

# Scatterplot GE & Meth
neg_p21_1=ggplot(neg_dat_m2_1,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = '(ER-) M2.1 GE v. Meth',subtitle = "(r = 0.253)", x= 'Gene Expression', y = 'Methylation')

#Boxplot CNA & GE
neg_p21_2 <- ggplot(neg_dat_m2_1, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M2.1 CNA v. GE', subtitle = "(r = 0.519)",x = 'Copy Number Alternation', y = 'Gene Expression')+theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))

# Boxplot CNA & Meth
neg_p21_3 <- ggplot(neg_dat_m2_1, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M2.1 CNA v. Meth', subtitle = "(r = -0.126)",x = 'CNA', y = 'Methylation')+theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))

(neg_p21_1)/(neg_p21_2|neg_p21_3)
ggsave("Neg_M2.1.pdf")



###11/7/24###
###For m2.2(C -> M <- E):###
#Fliter the M0.1 model datasets
neg_mod2_2 <- df2 %>% filter(Inferred.Model=="M2.2",Inferred.Model2=="M2.2",Inferred.Model3=="M2.2")
head(neg_mod2_2,3)
#Running the trio_row number#272239
trios[272239,]
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
#make sure these values are not in "lists"
#Find the correlation value of 2 varibales
neg_dat_m2_2 <- data.frame(CNA=unlist(neg_CNA[240,]),Meth=unlist(neg_meth[183166,]),GE=unlist(neg_GE[4742,]))
dim(neg_dat_m2_2)

cor(neg_dat_m2_2$GE,neg_dat_m2_2$Meth)
cor(neg_dat_m2_2$CNA,neg_dat_m2_2$GE)
cor(neg_dat_m2_2$CNA,neg_dat_m2_2$Meth)

# Scatterplot GE & Meth
neg_p22_1=ggplot(neg_dat_m2_2,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = '(ER-) M2.2 GE v. Meth',subtitle = "(r = 0.183)", x= 'Gene Expression', y = 'Methylation')

#Boxplot CNA & GE
neg_p22_2 <- ggplot(neg_dat_m2_2, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M2.2 CNA v. GE', subtitle = "(r = 0.113)",x = 'Copy Number Alternation', y = 'Gene Expression')+theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))

# Boxplot CNA & Meth
neg_p22_3 <- ggplot(neg_dat_m2_2, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M2.2 CNA v. Meth', subtitle = "(r = 0.500)",x = 'CNA', y = 'Methylation')+theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))

(neg_p22_1)/(neg_p22_2|neg_p22_3)

ggsave("Neg_M2.2.pdf")


##For m3 model(E ← C → M):###
#Filter the m3 model datasets
neg_mod3 <- df2 %>% filter(Inferred.Model=='M3',Inferred.Model2=='M3',Inferred.Model3 == 'M3')
dim(neg_mod3)
head(neg_mod3)
neg_mod3[100:105,]
#Running the trio_row number#270860
trios[270860,]
neg_dat_m3 <-  data.frame(CNA=unlist(neg_CNA[19832,]),Meth=unlist(neg_meth[453714,]),GE=unlist(neg_GE[11321,]))
dim(neg_dat_m3)

cor(neg_dat_m3$GE,neg_dat_m3$Meth)
cor(neg_dat_m3$CNA,neg_dat_m3$GE)
cor(neg_dat_m3$CNA,neg_dat_m3$Meth)

# Scatterplot GE & Meth
neg_p3_1=ggplot(neg_dat_m3,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = '(ER-) M3 GE v. Meth',subtitle = "(r = NA)", x= 'Gene Expression', y = 'Methylation')

#Boxplot CNA & GE
neg_p3_2 <- ggplot(neg_dat_m3, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M3 CNA v. GE', subtitle = "(r = 0.502)",x = 'Copy Number Alternation', y = 'Gene Expression')+theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))

# Boxplot CNA & Meth
neg_p3_3 <- ggplot(neg_dat_m3, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M3 CNA v. Meth', subtitle = "(r = NA)",x = 'CNA', y = 'Methylation')+theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))

(neg_p3_1)/(neg_p3_2|neg_p3_3)
ggsave("Neg_M3.pdf")


### M4 model()###
#Filter the m4 model datasets
neg_mod4 <-  df2 %>% filter(Inferred.Model == 'M4',Inferred.Model2 == 'M4',Inferred.Model3 == 'M4')
dim(neg_mod4)
neg_mod4[50:53,]
trios[270440,]
neg_dat_m4 <-  data.frame(CNA=unlist(neg_CNA[21367,]),GE=unlist(neg_GE[19537,]),Meth=unlist(neg_meth[146150,]))
dim(neg_dat_m4)
head(neg_dat_m4)

cor(neg_dat_m4$GE,neg_dat_m4$Meth)
cor(neg_dat_m4$CNA,neg_dat_m4$GE)
cor(neg_dat_m4$CNA,neg_dat_m4$Meth)

# Scatterplot GE & Meth
neg_p4_1=ggplot(neg_dat_m4,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = '(ER-) M4 GE v. Meth',subtitle = "(r = 0.274)", x= 'Gene Expression', y = 'Methylation')

#Boxplot CNA & GE
neg_p4_2 <- ggplot(neg_dat_m4, aes(x=factor(CNA), y=GE, color=factor(CNA)))+geom_boxplot()+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M4 CNA v. GE', subtitle = "(r = 0.470)",x = 'Copy Number Alternation', y = 'Gene Expression')+theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))

# Boxplot CNA & Meth
neg_p4_3 <- ggplot(neg_dat_m4, aes(x=factor(CNA), y=Meth))+geom_boxplot(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'M4 CNA v. Meth', subtitle = "(r = 0.313)",x = 'CNA', y = 'Methylation')+theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))

(neg_p4_1)/(neg_p4_2|neg_p4_3)
ggsave("Neg_M4.pdf")

### Other model###
neg_modOther <-  df2 %>% filter(Inferred.Model == 'Other',Inferred.Model2 == 'Other',Inferred.Model3 == 'Other')
dim(neg_modOther)
neg_modOther[14:18,]
trios[270027,]
neg_dat_Other <-  data.frame(CNA=unlist(neg_CNA[20630,]),GE=unlist(neg_GE[5411,]),Meth=unlist(neg_meth[275322,]))
dim(neg_dat_Other)
head(neg_dat_Other)

cor(neg_dat_Other$GE,neg_dat_Other$Meth)
cor(neg_dat_Other$CNA,neg_dat_Other$GE)
cor(neg_dat_Other$CNA,neg_dat_Other$Meth)

# Scatter plot GE & Meth
neg_p_other_1 <- ggplot(neg_dat_Other,aes(GE,Meth))+geom_point(aes(color=factor(CNA)))+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = '(ER-) Other Gene Expression v. Methylation',subtitle = "(r = 0.294)", x= 'Gene Expression', y = 'Methylation')

#Boxplot CNA & GE
neg_p_other_2 <- ggplot(neg_dat_Other, aes(x=factor(CNA), y=GE))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'Other CNA v. GE',subtitle = "(r = 0.044)", x = 'Copy Number Alternation', y = 'Gene Expression')

# Boxplot CNA & Meth
neg_p_other_3 <- ggplot(neg_dat_Other, aes(x=factor(CNA), y=Meth ))+geom_boxplot(aes(color=factor(CNA)),outlier.shape=NA)+theme_bw()+theme(legend.position="none",plot.title = element_text(size = 10),axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))+labs(title = 'Other CNA v. Methylation',subtitle = "(r = 0.084)", x = 'Copy Number Alternation', y = 'Methylation')
(neg_p_other_1)/(neg_p_other_2|neg_p_other_3)

ggsave("Neg_M_other.pdf")

