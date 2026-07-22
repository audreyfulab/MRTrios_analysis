

############ changing "." between ID names to "-" ################################

#set the directory
setwd("/mnt/ceph/megheib/Cancer_Genomics/TCGA_data_results")
getwd()

#load the library
library(data.table)
library(splitstackshape)

load("TCGA.meth.RData")
dim(TCGA.meth)
meth = TCGA.meth

for(i in 5:ncol(meth)){

  print(i)

  #split the name in the TCGA data
  split.ID <- unlist(strsplit(as.character(colnames(meth)[i]), "[.]"))
  split.ID

  #now merge the first 4 parts
  new.rowname <- paste(split.ID[1],"-",split.ID[2],"-",split.ID[3],"-01", sep = "")
  new.rowname

  #replace the old column name to the new column name
  colnames(meth)[i] <- new.rowname


}

#save it to a new txt file
write.table(meth, file = "/mnt/ceph/kark6289/TCGA_analysis/split.names.TCGA.meth.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE, append = FALSE, quote=FALSE)


##################### logit transformation of methylation data ############################

TCGA.meth <- as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/split.names.TCGA.meth.txt"))
dim(TCGA.meth)

data <- TCGA.meth[rowSums(is.na(TCGA.meth[,-(1:4)])) != ncol(TCGA.meth[,-(1:4)]), ]
dim(data)

# split data into numbers and probe info
data.info <- data[,1:4]

data.nona <- log(data[,-(1:4)]/(1-data[,-(1:4)]))

final <- cbind(data.info, data.nona)

#save it to a new txt file
write.table(final, file = "/mnt/ceph/kark6289/TCGA_analysis/split.names.TCGA.meth.logit.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE, append = FALSE, quote=FALSE)

############################ split the patients into ER+ and ER- patients ##################################

setwd("/Users/12083/OneDrive - University of Idaho/DOCUMENTS/AUDREY FU LAB/TCGA")
# getwd()

library(data.table)
TCGA.meth <- as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/split.names.TCGA.meth.logit.txt"))
dim(TCGA.meth)

data <- TCGA.meth[rowSums(is.na(TCGA.meth[,5:ncol(TCGA.meth)])) != ncol(TCGA.meth[,5:ncol(TCGA.meth)]), ]
dim(data)
data[1:5,1:5]

dim(data)

# split data into numbers and probe info
data.only <- t(data[,-(1:4)])
data.info <- data[,1:4]

clinical.data <- fread("brca_tcga_clinical_data.txt")
clinical.data[1:5,1:5]
dim(clinical.data)

library(tidyr)
tmp <- clinical.data %>% drop_na(`ER Status By IHC`)
tmp[1:5,1:5]
dim(tmp)

# rows <- grep("TCGA-A7-A26I",rownames(data.only))
# rows[1:5]

for(i in 1:nrow(data.only)){
  
#split the name in the TCGA data  
split.ID <- unlist(strsplit(as.character(rownames(data.only)[i]), '-'))
split.ID

#now use this ID to search in the clinical data
new.rowname <- paste(split.ID[1],"-",split.ID[2],"-",split.ID[3], sep = "")
new.rowname

#find the patient ID in the clinical data
row.num <- which(as.character(new.rowname) == as.character(tmp$`Patient ID`))
row.num

if(length(row.num)>0){
  
if(tmp$`ER Status By IHC`[row.num] == "Positive"){
  
  ER.pos <- c(new.rowname,data.only[i,],tmp$`ER Status By IHC`[row.num])
  
  write.table(t(ER.pos), file = "/mnt/ceph/kark6289/TCGA_analysis/ERpos.csv", sep = ",", row.names = FALSE,
              col.names = FALSE, append = TRUE,quote=FALSE)
}else if(tmp$`ER Status By IHC`[row.num] == "Negative"){
  
  ER.neg <- c(new.rowname,data.only[i,],tmp$`ER Status By IHC`[row.num])
  write.table(t(ER.neg), file = "/mnt/ceph/kark6289/TCGA_analysis/ERneg.csv", sep = ",", row.names = FALSE,
              col.names = FALSE, append = TRUE,quote=FALSE)

}
}
}


########## histograms for few methylation probes (transformed data) ###########

data <- as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/split.names.TCGA.meth.txt"))
dim(data)

pdf(file = "/mnt/ceph/kark6289/TCGA_analysis/newhist/unscaled/TCGA.pdf", onefile=TRUE)

for(i in 1:100) {
  
  print(i)
  
  data.nona <- log(data[i,5:899]/(1-data[i,5:899]))
  
  hist(as.numeric(data.nona), xlim = c(-2,3.5), main = paste("Histogram of (TCGA) ",data$Row.names[i], sep = ""), xlab = "Methylation")
  
}

dev.off()

