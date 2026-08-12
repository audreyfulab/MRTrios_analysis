###################################################
##       Load datasets                           ##
###################################################
library(data.table)
library(tidyverse)

# Read the datasets from saved ones
# posER_BRCA
setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_posER_BRCA/")

M0.1_pos <- fread("M0.1_pos_extract.txt")
M0.2_pos <- fread("M0.2_pos_extract.txt")
M1.1_pos <- fread("M1.1_pos_extract.txt")
M1.2_pos <- fread("M1.2_pos_extract.txt")

# negER_BRCA
setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_negER_BRCA/")

M0.1_neg <- fread("M0.1_neg_extract.txt")
M0.2_neg <- fread("M0.2_neg_extract.txt")
M1.1_neg <- fread("M1.1_neg_extract.txt")
M1.2_neg <- fread("M1.2_neg_extract.txt")
