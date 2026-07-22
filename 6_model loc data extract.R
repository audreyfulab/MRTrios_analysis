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

mypath = "/mnt/ceph/kark6289/PCandTrioAnalysis/output13.posER"
setwd(mypath)

# Create list of text files
txt_files_ls = list.files(path=mypath, pattern="*.txt") 
# Read the files in, assuming comma separator
txt_files_df <- lapply(txt_files_ls, function(x) {read.table(file = x, header = F, sep ="")})
# Combine them
combined_df1 <- do.call("rbind", lapply(txt_files_df, as.data.frame))

table(combined_df1$V2)

#310412

mypath = "/mnt/ceph/kark6289/PCandTrioAnalysis/output13.negER"
setwd(mypath)

# Create list of text files
txt_files_ls = list.files(path=mypath, pattern="*.txt") 
# Read the files in, assuming comma separator
txt_files_df <- lapply(txt_files_ls, function(x) {read.table(file = x, header = F, sep ="")})
# Combine them
combined_df2 <- do.call("rbind", lapply(txt_files_df, as.data.frame))

table(combined_df2$V2)

df.loc <- read.csv("/mnt/ceph/kark6289/test_trio/trios/trio.final.protein.coding.location.csv", sep = "\t")

rows1 <- combined_df1$V1[which(combined_df1$V2 == "M1.1")]
rows2 <- combined_df1$V1[which(combined_df1$V2 == "M1.2")]


#create a summary table for the location
table(unlist(strsplit(as.character(df.loc[rows1,2]), ';')))

table(unlist(strsplit(as.character(df.loc[rows2,2]), ';')))


#function to extract row numbers of trios that have a specific MRGN model and location
extract_trios_row <- function(model, location, combined_df, df.loc){
  
  #find the rows in trios that have a certain model type
  rows <- combined_df$V1[which(combined_df$V2 == as.character(model))]
  
  #summary of the locations of that model type
  print(table(unlist(strsplit(as.character(df.loc[rows,2]), ';'))))
  
  # g1 <- rep(seq_along(split.loc), sapply(split.loc, length))
  # g1[1:5]
  # 
  # df.loc.row <- g1[which(unlist(split.loc) == as.character(location))]
  # df.loc.row
  
  #find which of the trios with that model type are located in a specific location
  row.with.loc <- grep(as.character(location), df.loc[rows,2])
  
  #find the row number in trios data
  row.in.trios <- df.loc[rows[row.with.loc],1]
  df.loc[row.in.trios[1:10],]
  
  #the trio data with a specific model result and location
  trios[row.in.trios[1:5],]
  
  #return the row number
  return(row.in.trios)
}

#Pos
rows_in_trios_pos <- extract_trios_row("M2", "TSS1500", combined_df1, df.loc)

#Neg
rows_in_trios_neg <- extract_trios_row("M1.1", "TSS1500", combined_df2, df.loc)
