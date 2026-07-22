library(data.table)
library(na.tools)
library(MRGN, lib = "/mnt/ceph/kark6289/Rlibs")

#read in the datasets
TCGA.meth <- as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/split.names.TCGA.meth.logit.txt"))
dim(TCGA.meth)

#read in the datasets
gene.exp <- as.data.frame(fread("/mnt/ceph/fu_lab/TCGA/cbioportal/brca_tcga/data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt"))
dim(gene.exp)


#reading in the neg and pos ER individuals data
clinical.neg <- fread("/mnt/ceph/kark6289/TCGA_analysis/names.neg.patient2.txt", header = FALSE)
dim(clinical.neg)

neg.ind <- clinical.neg[,1]

clinical.pos <- fread("/mnt/ceph/kark6289/TCGA_analysis/names.pos.patient2.txt", header = FALSE)
dim(clinical.pos)

pos.ind <- clinical.pos[,1]

#finding common individuals among the 3 datasets
com.ind = intersect(colnames(gene.exp)[3:ncol(gene.exp)], colnames(TCGA.meth)[5:ncol(TCGA.meth)]) #length 787

#Positive ER
meth.pos.final <- findPCs(TCGA.meth, startCol = 5, GeneNameCol = 2, pos.ind, com.ind, "pos")
gene.pos.final <- findPCs(gene.exp, startCol = 3, GeneNameCol = 1, pos.ind, com.ind, "pos")

#Negative ER
meth.neg.final <- findPCs(TCGA.meth, startCol = 5, GeneNameCol = 2, neg.ind, com.ind, "neg")
gene.neg.final <- findPCs(gene.exp, startCol = 3, GeneNameCol = 1, neg.ind, com.ind, "neg")


#Meth pos ER
write.table(meth.pos.final[[1]], file = paste("/mnt/ceph/kark6289/PCandTrioAnalysis/pcdata2/PCA.meth.posER.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = NA, append = TRUE, quote=FALSE)

saveRDS(meth.pos.final[[2]], file="/mnt/ceph/kark6289/PCandTrioAnalysis/pcdata2/meth.posER.sig.asso.pcs.RData")

write.table(meth.pos.final[[3]], file = paste("/mnt/ceph/kark6289/PCandTrioAnalysis/pcdata2/meth.posER.table.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = NA, append = TRUE, quote=FALSE)

write.table(meth.pos.final[[4]], file = paste("/mnt/ceph/kark6289/PCandTrioAnalysis/pcdata2/step3.meth.posER.data.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = TRUE, append = FALSE, quote=FALSE)

#Gene exp pos ER
write.table(gene.pos.final[[1]], file = paste("/mnt/ceph/kark6289/PCandTrioAnalysis/pcdata2/PCA.gene.exp.posER.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = NA, append = TRUE, quote=FALSE)

saveRDS(gene.pos.final[[2]], file="/mnt/ceph/kark6289/PCandTrioAnalysis/pcdata2/gene.exp.posER.sig.asso.pcs.RData")

write.table(gene.pos.final[[3]], file = paste("/mnt/ceph/kark6289/PCandTrioAnalysis/pcdata2/gene.exp.posER.table.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = NA, append = TRUE, quote=FALSE)

write.table(gene.pos.final[[4]], file = paste("/mnt/ceph/kark6289/PCandTrioAnalysis/pcdata2/step3.gene.exp.posER.data.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = TRUE, append = FALSE, quote=FALSE)

#Meth neg ER
write.table(meth.neg.final[[1]], file = paste("/mnt/ceph/kark6289/PCandTrioAnalysis/pcdata2/PCA.meth.negER.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = NA, append = TRUE, quote=FALSE)

saveRDS(meth.neg.final[[2]], file="/mnt/ceph/kark6289/PCandTrioAnalysis/pcdata2/meth.negER.sig.asso.pcs.RData")

write.table(meth.neg.final[[3]], file = paste("/mnt/ceph/kark6289/PCandTrioAnalysis/pcdata2/meth.negER.table.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = NA, append = TRUE, quote=FALSE)

write.table(meth.neg.final[[4]], file = paste("/mnt/ceph/kark6289/PCandTrioAnalysis/pcdata2/step3.meth.negER.data.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = TRUE, append = FALSE, quote=FALSE)

#Gene exp neg ER
write.table(gene.neg.final[[1]], file = paste("/mnt/ceph/kark6289/PCandTrioAnalysis/pcdata2/PCA.gene.exp.negER.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = NA, append = TRUE, quote=FALSE)

saveRDS(gene.neg.final[[2]], file="/mnt/ceph/kark6289/PCandTrioAnalysis/pcdata2/gene.exp.negER.sig.asso.pcs.RData")

write.table(gene.neg.final[[3]], file = paste("/mnt/ceph/kark6289/PCandTrioAnalysis/pcdata2/gene.exp.negER.table.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = NA, append = TRUE, quote=FALSE)

write.table(gene.neg.final[[4]], file = paste("/mnt/ceph/kark6289/PCandTrioAnalysis/pcdata2/step3.gene.exp.negER.data.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = TRUE, append = FALSE, quote=FALSE)
