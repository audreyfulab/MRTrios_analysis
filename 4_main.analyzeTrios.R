library(data.table)
library(na.tools)
library(MRGN, lib = "/mnt/ceph/kark6289/Rlibs")

#read in the original datasets
gene.exp <- as.data.frame(fread("/mnt/ceph/fu_lab/TCGA/cbioportal/brca_tcga/data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt"))
dim(gene.exp)

cna <- as.data.frame(fread("/mnt/ceph/fu_lab/TCGA/cbioportal/brca_tcga/data_CNA.txt"))
dim(cna)

TCGA.meth <- as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/split.names.TCGA.meth.logit.txt"))
dim(TCGA.meth)

#Read in the PC score matrix
pc.meth.pos <- read.table("/mnt/ceph/kark6289/PCandTrioAnalysis/PCA.meth.posER.txt", row.names = 1)
pc.meth.neg <- read.table("/mnt/ceph/kark6289/PCandTrioAnalysis/PCA.meth.negER.txt", row.names = 1)

pc.gene.pos <- read.table("/mnt/ceph/kark6289/PCandTrioAnalysis/PCA.gene.exp.posER.txt", row.names = 1)
pc.gene.neg <- read.table("/mnt/ceph/kark6289/PCandTrioAnalysis/PCA.gene.exp.negER.txt", row.names = 1)


#read in the neg and pos ER individuals data
clinical.neg <- fread("/mnt/ceph/kark6289/TCGA_analysis/names.neg.patient2.txt", header = FALSE)
dim(clinical.neg)

clinical.pos <- fread("/mnt/ceph/kark6289/TCGA_analysis/names.pos.patient2.txt", header = FALSE)
dim(clinical.pos)

#reading in the Trios data
trios <- data.frame(fread("/mnt/ceph/kark6289/test_trio/trios/trio.final.protein.coding.txt"))

#read in the indices table
meth.table.pos <- fread("/mnt/ceph/kark6289/PCandTrioAnalysis/meth.posER.table.txt", drop = 1)
gene.table.pos <- fread("/mnt/ceph/kark6289/PCandTrioAnalysis/gene.exp.posER.table.txt", drop = 1)

meth.table.neg <- fread("/mnt/ceph/kark6289/PCandTrioAnalysis/meth.negER.table.txt", drop = 1)
gene.table.neg <- fread("/mnt/ceph/kark6289/PCandTrioAnalysis/gene.exp.negER.table.txt", drop = 1)

#read in the sig pcs data
meth.sig.asso.pcs.pos <- readRDS("/mnt/ceph/kark6289/PCandTrioAnalysis/meth.posER.sig.asso.pcs.RData")
gene.sig.asso.pcs.pos <- readRDS("/mnt/ceph/kark6289/PCandTrioAnalysis/gene.exp.posER.sig.asso.pcs.RData")
  
meth.sig.asso.pcs.neg <- readRDS("/mnt/ceph/kark6289/PCandTrioAnalysis/meth.negER.sig.asso.pcs.RData")
gene.sig.asso.pcs.neg <- readRDS("/mnt/ceph/kark6289/PCandTrioAnalysis/gene.exp.negER.sig.asso.pcs.RData")

analyzeTrios(TCGA.meth, gene.exp, cna, trios, pc.meth.pos, pc.gene.pos, meth.sig.asso.pcs.pos, gene.sig.asso.pcs.pos, clinical.pos, meth.table.pos, gene.table.pos, path = "/mnt/ceph/kark6289/PCandTrioAnalysis/output1.posER/model.trio.MRGN.txt")
analyzeTrios(TCGA.meth, gene.exp, cna, trios, pc.meth.neg, pc.gene.neg, meth.sig.asso.pcs.neg, gene.sig.asso.pcs.neg, clinical.neg, meth.table.neg, gene.table.neg, path = "/mnt/ceph/kark6289/PCandTrioAnalysis/output1.negER/model.trio.MRGN.txt")
