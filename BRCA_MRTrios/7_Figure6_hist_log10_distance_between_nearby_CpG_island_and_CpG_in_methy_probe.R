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

##########################################################################################################################
######  Distribution of the log10 (distance) between the nearby CpG island and the CpG in the methylation probe ##########
##########################################################################################################################

setwd("/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Figure")
library(ggplot2)
library(patchwork)
library(cowplot)

#hist_M0.1_pos_log_distance=ggplot(M0.1_pos %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() , aes(x = log10(abs(diff_cpG_mapinfo)))) +
#  geom_histogram(aes(y = ..density..), bins = 40, fill = "#f781bf", color = "black", alpha = 0.7) +
#  labs(title = "M0.1",x = "log10 (distance)",y = "Density") + theme_bw()+
#  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"),
#        axis.title.x = element_text(size = 8),  
#        axis.title.y = element_text(size = 8))+xlim(-2,4)+ylim(0,1)

hist_M0.1_pos_log_distance=ggplot(M0.1_pos %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() , aes(x = log10(abs(diff_cpG_mapinfo)))) +
  geom_histogram(aes(y = ..density..), bins = 34, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+xlim(-2,4)+ylim(0,1)

hist_M0.2_pos_log_distance=ggplot(M0.2_pos%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 34, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

hist_M1.1_pos_log_distance=ggplot(M1.1_pos%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 34, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

hist_M1.2_pos_log_distance=ggplot(M1.2_pos%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 34, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

hist_M0.1_neg_log_distance=ggplot(M0.1_neg%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 34, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

hist_M0.2_neg_log_distance=ggplot(M0.2_neg%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 34, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

hist_M1.1_neg_log_distance=ggplot(M1.1_neg%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 34, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "log10 (distance)", y = "Density") +theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+ xlim(-2,4)+ylim(0,1)

hist_M1.2_neg_log_distance=ggplot(M1.2_neg%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 34, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

log_distance_layout_1=(hist_M0.1_pos_log_distance|hist_M1.1_pos_log_distance)/(hist_M0.2_pos_log_distance|hist_M1.2_pos_log_distance)+
  plot_annotation(title = "ER+")& theme(plot.title = element_text(hjust = 0.5,size = 12))
log_distance_layout_1

#ggsave("hist_log_distance_1.pdf") (5*7inches) Landscape

log_distance_layout_2=(hist_M0.1_neg_log_distance|hist_M1.1_neg_log_distance)/(hist_M0.2_neg_log_distance|hist_M1.2_neg_log_distance)+
  plot_annotation(title = "ER-")& theme(plot.title = element_text(hjust = 0.5,size = 12))
log_distance_layout_2

#ggsave("hist_log_distance_2.pdf")(5*7inches) Landscape

plot_grid(log_distance_layout_1,log_distance_layout_2,labels = "AUTO",label_size = 15, ncol = 1)
#ggsave(6*8inches) potrait
#plot_grid(log_distance_layout_1,log_distance_layout_2,labels = c("A (ER+)","B (ER-)"), ncol = 1)

## Figure6_hist_log10_distance_between_nearby_CpG_island_and_CpG_in_methy_probe.pdf



