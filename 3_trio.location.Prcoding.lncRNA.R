#load the package
library(data.table)

#read in the datasets
# probe info data downloaded from GEO:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534
df <- read.csv("/mnt/ceph/kark6289/GPL13534_HumanMethylation450_15017482_v.1.2.csv", skip = 7, header = TRUE)
dim(df)

trios <- data.frame(fread("/mnt/ceph/kark6289/test_trio/trios/trio.final.protein.coding.txt"))

TCGA.meth <- as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/split.names.TCGA.meth.logit.txt"))
dim(TCGA.meth)


#create a text file for location
write.table(t(c("trios.row", "location")), file = "/mnt/ceph/kark6289/test_trio/trios/trio.final.protein.coding.location.csv", sep = "\t", row.names = FALSE,
            col.names = FALSE, append = FALSE,quote=FALSE)

#loop to go over each gene
for(i in 1:nrow(trios)){
  
  print(i)
  
  #find if the gene is in the Human Methylation data
  row <- which(TCGA.meth[trios[i,2],1] == df$Name)
  row
  
  #check if the length is greater than 0
  if(length(row)>0){
  
  #use the row number to find the location
  loc <- df$UCSC_RefGene_Group[row]
  loc
  
  #save the gene row in trio and location in HM data to a new variable
  final <- c(i, loc)
  final
  
  #write final to the location file
  write.table(t(final), file = "/mnt/ceph/kark6289/test_trio/trios/trio.final.protein.coding.location.csv", sep = "\t", row.names = FALSE,
              col.names = FALSE, append = TRUE, quote=FALSE)
  
  }
}

#read in the created dataset for location
df.loc <- read.csv("/mnt/ceph/kark6289/test_trio/trios/trio.final.protein.coding.location.csv", header = TRUE, sep = "\t")
dim(df.loc)

df.loc[1:5,]

#create a summary table for the location
table(unlist(strsplit(as.character(df.loc$location), ';')))

