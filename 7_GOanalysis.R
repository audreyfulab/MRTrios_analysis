##################################################################
#                                                                #
#                                                                #
#                           MAIN CODE                            #
#                                                                #
#                                                                #
##################################################################

#Note: you can just use the function below (other codes before the function are just trying to see how it works)
library(gprofiler2)

#set directory
setwd("/Users/kark6289/OneDrive - University of Idaho/DOCUMENTS/AUDREY FU LAB/Probe info/addGeneMeans")

#read in the dataset containing genes of interest
M1.1.pos = read.csv("M1.1_pos.csv")

M1.2.pos = read.csv("M1.2_pos.csv")

M1.1.neg = read.csv("M1.1_neg.csv")

M1.2.neg = read.csv("M1.2_neg.csv")

#Extract the genes of interest
up.genes = unique(M1.2.pos$Biomart_GeneName)
#length(up.genes)
#head(up.genes)

####### Trying different methods #########
gostres1 = gost(query = up.genes, organism = "hsapiens", sources = "GO:BP")
head(gostres1$result)
dim(gostres1$result)

gostplot(gostres1)

gostres2 = gost(query = up.genes, organism = "hsapiens", correction_method = "fdr", sources = "GO:BP")
head(gostres2$result)
dim(gostres2$result)

gostplot(gostres2)

gostres3 = gost(query = up.genes, organism = "hsapiens", correction_method = "bonferroni", sources = "GO:BP")
head(gostres3$result)

######### function to perform the analysis
GOanalysis <- function(data, source){
  
  #Extract the genes of interest
  up.genes = unique(data$Biomart_GeneName)
  
  #use the gost function to get the analysis result
  gostres1 = gost(query = up.genes, organism = "hsapiens", sources = source)
  
  #head(gostres1$result)
  #dim(gostres1$result)
  
  #changing the result to a dataframe so it can be saved as csv file
  df <- apply(gostres1$result,2,as.character)
  
  return(df)
  
}

M1.1.pos_GO = GOanalysis(M1.1.pos, "GO:BP")
M1.2.pos_GO = GOanalysis(M1.2.pos, "GO:BP")
M1.1.neg_GO = GOanalysis(M1.1.neg, "GO:BP")
M1.2.neg_GO = GOanalysis(M1.2.neg, "GO:BP")

write.csv(M1.1.pos_GO, "/Users/kark6289/OneDrive - University of Idaho/DOCUMENTS/AUDREY FU LAB/GO enrichment analysis/GOtable_M1.1_ERpos_BP_gSCS.csv", row.names=FALSE)
write.csv(M1.2.pos_GO, "/Users/kark6289/OneDrive - University of Idaho/DOCUMENTS/AUDREY FU LAB/GO enrichment analysis/GOtable_M1.2_ERpos_BP_gSCS.csv", row.names=FALSE)
write.csv(M1.1.neg_GO, "/Users/kark6289/OneDrive - University of Idaho/DOCUMENTS/AUDREY FU LAB/GO enrichment analysis/GOtable_M1.1_ERneg_BP_gSCS.csv", row.names=FALSE)
write.csv(M1.2.neg_GO, "/Users/kark6289/OneDrive - University of Idaho/DOCUMENTS/AUDREY FU LAB/GO enrichment analysis/GOtable_M1.2_ERneg_BP_gSCS.csv", row.names=FALSE)


##################################################################
#                                                                #
#                                                                #
#                      Extra (Other codes)                       #
#                                                                #
#                                                                #
##################################################################

library(topGO)
library(ALL)
library(org.Hs.eg.db)
library(biomaRt)

data(ALL)
data(geneList)

affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)

sum(topDiffGenes(geneList))

#sample_list = geneList[1:100]

sampleGOdata <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = geneList, geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.db, affyLib = affyLib)

sampleGOdata

resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultFisher

resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

allRes


########################################################################

setwd("/Users/kark6289/OneDrive - University of Idaho/DOCUMENTS/AUDREY FU LAB/Probe info/addGeneMeans")

M1.1.pos = read.csv("M1.1_pos.csv")
dim(M1.1.pos)
M1.1.pos[1,]

# Read in genes of interest
candidate_list = unique(M1.1.pos$Biomart_GeneName)
length(candidate_list)
head(candidate_list)

# create GO db for genes to be used using biomaRt - please note that this takes a while
db= useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl', host="www.ensembl.org")
go_ids= getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), filters='external_gene_name', values=candidate_list, mart=db)


listAttributes(db)

# build the gene 2 GO annotation list (needed to create topGO object)
gene_2_GO=unstack(go_ids[,c(1,2)])

# remove any candidate genes without GO annotation
keep = candidate_list %in% go_ids[,2]
keep =which(keep==TRUE)
candidate_list=candidate_list[keep]

###????????????????????

# make named factor showing which genes are of interest
geneList=factor(as.integer(candidate_list %in% candidate_list))
names(geneList)= candidate_list


GOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO)


################################

library(clusterProfiler)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

# Entrez gene ID
head(gene)

ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)

ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

########################################################

library(data.table)

#read the datasets containing all the genes and their gene type
gene.info.data <- data.frame(fread("/mnt/ceph/jarredk/Reg_Net/mart_export_merged_lncRNA_fixed.txt"))

#extract the rows of protein coding genes and lnRNAs
rows = which(gene.info.data$Gene.type == "protein_coding" | gene.info.data$Gene.type == "lncRNA" )

#extract the gene names and gene types from those rows
genes = gene.info.data$Gene.name[rows]
genetype = gene.info.data$Gene.type[rows]

#combine the gene names and types
data = cbind(genes, genetype)

#save to file
write.table(data, file = "/mnt/ceph/kark6289/GOanalysis/geneList.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE, append = TRUE, quote=FALSE)

############################################################

library("data.table")
library("topGO")
library("org.Hs.eg.db")
library("ggplot2")
library("knitr")

#set directory
setwd("/Users/kark6289/OneDrive - University of Idaho/DOCUMENTS/AUDREY FU LAB/Probe info/addGeneMeans")

#read in the dataset containing genes of interest
M1.1.pos = read.csv("M1.1_pos.csv")
#dim(M1.1.pos)
#M1.1.pos[1,]

#Extract the genes of interest
up.genes = unique(M1.1.pos$Biomart_GeneName)
#length(up.genes)
#head(up.genes)

#directory for the data with the universal list
setwd("/Users/kark6289/OneDrive - University of Idaho/DOCUMENTS/AUDREY FU LAB/GO enrichment analysis")

#read in the dataset containing the universal gene list
geneList = fread("geneList.txt")

#extract the universal gene list
all.genes = geneList$genes
#head(all.genes)


#Method = "classic" and Ontology = "BP" (Biological Process)
  
#specify the ontology (BP, MF, or CC)
ontology <- "BP"

#specify the algorithm (classic or weight01)
algorithm <- "classic"

#specify the test statistic
statistic <- "fisher"

#identify the genes in universal gene list that are in the genes of interest list
#so if a gene in universal gene list contains a gene that is in the genes of interest list
#we assign values 1
#if not 0
upList <- factor(as.integer(all.genes %in% up.genes))
names(upList) <- all.genes

head(upList, 30)

table(upList)

#create topG0 data with the GO IDs for the specified ontology
upGOdata <- new("topGOdata", ontology = ontology, allGenes = upList, geneSel = function(x)(x == 1), 
                nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "SYMBOL")

#extract the GO term table using the specified algorithm and test statistic
upRes <- runTest(upGOdata, algorithm = algorithm, statistic = statistic)
upRes

#sort the table by the pvalues
up.tab <- GenTable(upGOdata, Pval = upRes, topNodes = 6873)

#generates a fancy table
kable(up.tab, caption = "Significance of GO terms according to classic method and BP ontology")

# write.csv(up.tab, "/Users/kark6289/OneDrive - University of Idaho/DOCUMENTS/AUDREY FU LAB/GO enrichment analysis/GOtable_classic_BP.csv", row.names=TRUE)

#Method = "weight01" and Ontology = "BP" (Biological Process)
  
#specify the algorithm
algorithm <- "weight01"

#extract the GO term table using the specified algorithm and test statistic
upRes <- runTest(upGOdata, algorithm = algorithm, statistic = statistic)
upRes

#sort the table by the pvalues
up.tab <- GenTable(upGOdata, Pval = upRes, topNodes = 6873)

#generates a fancy table
kable(up.tab, caption = "Significance of GO terms according to weight01 method and BP ontology")

# write.csv(up.tab, "/Users/kark6289/OneDrive - University of Idaho/DOCUMENTS/AUDREY FU LAB/GO enrichment analysis/GOtable_weight01_BP.csv", row.names=TRUE)


#Method = "classic" and Ontology = "MF" (Molecular Function)

#specify the ontology
ontology <- "MF"

#specify the algorithm
algorithm <- "classic"

#create topG0 data with the GO IDs for the specified ontology
upGOdata <- new("topGOdata", ontology = ontology, allGenes = upList, geneSel = function(x)(x == 1), 
                nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "SYMBOL")

#extract the GO term table using the specified algorithm and test statistic
upRes <- runTest(upGOdata, algorithm = algorithm, statistic = statistic)
upRes

#sort the table by the pvalues
up.tab <- GenTable(upGOdata, Pval = upRes, topNodes = 1312)

#generates a fancy table
kable(up.tab, caption = "Significance of GO terms according to classic method and MF ontology")

# write.csv(up.tab, "/Users/kark6289/OneDrive - University of Idaho/DOCUMENTS/AUDREY FU LAB/GO enrichment analysis/GOtable_classic_MF.csv", row.names=TRUE)

#Method = "weight01" and Ontology = "MF" (Molecular Function)

#specify the algorithm
algorithm <- "weight01"

#extract the GO term table using the specified algorithm and test statistic
upRes <- runTest(upGOdata, algorithm = algorithm, statistic = statistic)
upRes

#sort the table by the pvalues
up.tab <- GenTable(upGOdata, Pval = upRes, topNodes = 1312)

#generates a fancy table
kable(up.tab, caption = "Significance of GO terms according to weight01 method and MF ontology")

# write.csv(up.tab, "/Users/kark6289/OneDrive - University of Idaho/DOCUMENTS/AUDREY FU LAB/GO enrichment analysis/GOtable_weight01_MF.csv", row.names=TRUE)

#Method = "classic" and Ontology = "CC" (Cellular Component)

#specify the ontology
ontology <- "CC"

#specify the algorithm
algorithm <- "classic"

#create topG0 data with the GO IDs for the specified ontology
upGOdata <- new("topGOdata", ontology = ontology, allGenes = upList, geneSel = function(x)(x == 1), 
                nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "SYMBOL")

#extract the GO term table using the specified algorithm and test statistic
upRes <- runTest(upGOdata, algorithm = algorithm, statistic = statistic)
upRes

#sort the table by the pvalues
up.tab <- GenTable(upGOdata, Pval = upRes, topNodes = 883)

#generates a fancy table
kable(up.tab, caption = "Significance of GO terms according to classic method and CC ontology")

# write.csv(up.tab, "/Users/kark6289/OneDrive - University of Idaho/DOCUMENTS/AUDREY FU LAB/GO enrichment analysis/GOtable_classic_CC.csv", row.names=TRUE)

#Method = "weight01" and Ontology = "CC" (Cellular Component)

#specify the algorithm
algorithm <- "weight01"

#extract the GO term table using the specified algorithm and test statistic
upRes <- runTest(upGOdata, algorithm = algorithm, statistic = statistic)
upRes

#sort the table by the pvalues
up.tab <- GenTable(upGOdata, Pval = upRes, topNodes = 883)

#generates a fancy table
kable(up.tab, caption = "Significance of GO terms according to weight01 method and CC ontology")

# write.csv(up.tab, "/Users/kark6289/OneDrive - University of Idaho/DOCUMENTS/AUDREY FU LAB/GO enrichment analysis/GOtable_weight01_CC.csv", row.names=TRUE)

#############################################################

