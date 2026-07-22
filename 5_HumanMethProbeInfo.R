#load the packages
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(data.table)

#import the data
df1 <- read.delim("/mnt/ceph/kark6289/PCandTrioAnalysis/output.10.14/model.trio.MRGN.all.posER.reclassify2.txt", header = TRUE)

df2 <- read.delim("/mnt/ceph/kark6289/PCandTrioAnalysis/output.10.14/model.trio.MRGN.all.negER.reclassify2.txt", header = TRUE)

df1 <- read.delim("/mnt/ceph/kark6289/PCandTrioAnalysis/downsampling/pos_ER/output/model.trio.MRGN.all.posER.reclassify.txt", header = TRUE)

df.loc <- read.csv("/mnt/ceph/kark6289/test_trio/trios/trio.final.protein.coding.location.csv", sep = "\t")

# this probe info file is identical to GPL13534_HumanMethylation450_15017482_v.1.1 2.csv
humanmeth <- read.csv("/mnt/ceph/kark6289/humanmethylation450_15017482_v1-2.csv", skip = 7, header = TRUE)

TCGA.meth <- as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/split.names.TCGA.meth.logit.txt"))
dim(TCGA.meth)

#read in the original datasets
gene.exp <- as.data.frame(fread("/mnt/ceph/fu_lab/TCGA/cbioportal/brca_tcga/data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt"))
dim(gene.exp)

cna <- as.data.frame(fread("/mnt/ceph/fu_lab/TCGA/cbioportal/brca_tcga/data_CNA.txt"))
dim(cna)

#reading in the Trios data
trios <- data.frame(fread("/mnt/ceph/kark6289/test_trio/trios/trio.final.protein.coding.txt"))

#read in the bio mart data
biomart <- read.delim("/mnt/ceph/kark6289/PCandTrioAnalysis/output.10.14/ensembl37_genes_p13_biomart.txt", header = TRUE)

#reading in the neg and pos ER individuals data
clinical.neg <- fread("/mnt/ceph/kark6289/TCGA_analysis/names.neg.patient2.txt", header = FALSE)
dim(clinical.neg)

neg.ind <- clinical.neg[,1]

clinical.pos <- fread("/mnt/ceph/kark6289/TCGA_analysis/names.pos.patient2.txt", header = FALSE)
dim(clinical.pos)

pos.ind <- clinical.pos[,1]


M0.1_pos = extractHumanMethProbeInfo(df1, "M0.1", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, pos.ind)
M0.2_pos = extractHumanMethProbeInfo(df1, "M0.2", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, pos.ind)

M0.1_neg = extractHumanMethProbeInfo(df2, "M0.1", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, neg.ind)
M0.2_neg = extractHumanMethProbeInfo(df2, "M0.2", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, neg.ind)

M1.1_pos = extractHumanMethProbeInfo(df1, "M1.1", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, pos.ind)
M1.2_pos = extractHumanMethProbeInfo(df1, "M1.2", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, pos.ind)

M1.1_neg = extractHumanMethProbeInfo(df2, "M1.1", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, neg.ind)
M1.2_neg = extractHumanMethProbeInfo(df2, "M1.2", TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, 5, 3, 3, neg.ind)
