
# Read the datasets from saved ones
M1.2_neg <- read.delim("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/M1.2_neg_extract.txt",sep = "\t")
M1.2_pos <- read.delim("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/M1.2_pos_extract.txt",sep = "\t")
M1.1_neg <- read.delim("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/M1.1_neg_extract.txt",sep = "\t")
M1.1_pos <- read.delim("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/M1.1_pos_extract.txt",sep = "\t")

M0.2_neg <- read.delim("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/M0.2_neg_extract.txt",sep = "\t")
M0.2_pos <- read.delim("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/M0.2_pos_extract.txt",sep = "\t")
M0.1_neg <- read.delim("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/M0.1_neg_extract.txt",sep = "\t")
M0.1_pos <- read.delim("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/M0.1_pos_extract.txt",sep = "\t")


##########################################################################################################################
######  Distribution of the log10 (distance) between the nearby CpG island and the CpG in the methylation probe ##########
##########################################################################################################################

setwd("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Plot_Result_Methyl_LZ")
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
  geom_histogram(aes(y = ..density..), bins = 40, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+xlim(-2,4)+ylim(0,1)

hist_M0.2_pos_log_distance=ggplot(M0.2_pos%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 40, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

hist_M1.1_pos_log_distance=ggplot(M1.1_pos%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 40, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

hist_M1.2_pos_log_distance=ggplot(M1.2_pos%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 40, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

hist_M0.1_neg_log_distance=ggplot(M0.1_neg%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

hist_M0.2_neg_log_distance=ggplot(M0.2_neg%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

hist_M1.1_neg_log_distance=ggplot(M1.1_neg%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "log10 (distance)", y = "Density") +theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+ xlim(-2,4)+ylim(0,1)

hist_M1.2_neg_log_distance=ggplot(M1.2_neg%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

log_distance_layout_1=(hist_M0.1_pos_log_distance|hist_M1.1_pos_log_distance)/(hist_M0.2_pos_log_distance|hist_M1.2_pos_log_distance)+
  plot_annotation(title = "ER+")& theme(plot.title = element_text(hjust = 0.5,size = 12))
log_distance_layout_1

#ggsave("hist_log_distance_1.pdf")

log_distance_layout_2=(hist_M0.1_neg_log_distance|hist_M1.1_neg_log_distance)/(hist_M0.2_neg_log_distance|hist_M1.2_neg_log_distance)+
  plot_annotation(title = "ER-")& theme(plot.title = element_text(hjust = 0.5,size = 12))
log_distance_layout_2

#ggsave("hist_log_distance_2.pdf")

plot_grid(log_distance_layout_1,log_distance_layout_2,labels = "AUTO",label_size = 15, ncol = 1)
#ggsave(6*8inches) potrait
#plot_grid(log_distance_layout_1,log_distance_layout_2,labels = c("A (ER+)","B (ER-)"), ncol = 1)
##########


######################################################################################################################
### Distribution of the log 10 (distance) between the CpG in the methylation probe and the start position of a gene.
######################################################################################################################

hist_M0.1_pos_log_distance_start=ggplot(M0.1_pos %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs(), aes(x = log10(abs(diff_mapinfo_geneStart)))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.1", x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+
  ylim(0,0.9)+xlim(-2.5,7.5)

hist_M0.2_pos_log_distance_start=ggplot(M0.2_pos %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs(), aes(x = log10(abs(diff_mapinfo_geneStart)))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+
  ylim(0,0.9)+xlim(-2.5,7.5)

hist_M1.1_pos_log_distance_start=ggplot(M1.1_pos %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs(), aes(x = log10(abs(diff_mapinfo_geneStart)))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.1", x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+
  ylim(0,0.9)+xlim(-2.5,7.5)

hist_M1.2_pos_log_distance_start=ggplot(M1.2_pos %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs(), aes(x = log10(abs(diff_mapinfo_geneStart)))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "log10 (distance)",y = "Density")  + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+
  ylim(0,0.9)+xlim(-2.5,7.5)


hist_M0.1_neg_log_distance_start=ggplot(M0.1_neg %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs(), aes(x = log10(abs(diff_mapinfo_geneStart)))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "log10 (distance)", y = "Density")  + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+
  ylim(0,0.9)+xlim(-2.5,7.5)

hist_M0.2_neg_log_distance_start=ggplot(M0.2_neg %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs(), aes(x = log10(abs(diff_mapinfo_geneStart)))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "log10 (distance)",y = "Density")  + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+
  ylim(0,0.9)+xlim(-2.5,7.5)

hist_M1.1_neg_log_distance_start=ggplot(M1.1_neg %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs(), aes(x = log10(abs(diff_mapinfo_geneStart)))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.1", x = "log10 (distance)",y = "Density")  + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+
  ylim(0,0.9)+xlim(-2.5,7.5)

hist_M1.2_neg_log_distance_start=ggplot(M1.2_neg %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs(), aes(x = log10(abs(diff_mapinfo_geneStart)))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.2", x = "log10 (distance)",y = "Density")  + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+
  ylim(0,0.9)+xlim(-2.5,7.5)

log_distance_start_layout_1=(hist_M0.1_pos_log_distance_start|hist_M1.1_pos_log_distance_start)/(hist_M0.2_pos_log_distance_start|hist_M1.2_pos_log_distance_start)+
  plot_annotation(title = "ER+")& theme(plot.title = element_text(hjust = 0.5,size = 12))
log_distance_start_layout_1

#ggsave("hist_log_distance_start_1.pdf")

log_distance_start_layout_2=(hist_M0.1_neg_log_distance_start|hist_M1.1_neg_log_distance_start)/(hist_M0.2_neg_log_distance_start|hist_M1.2_neg_log_distance_start)+
  plot_annotation(title = "ER-")& theme(plot.title = element_text(hjust = 0.5,size = 12))
log_distance_start_layout_2

#ggsave("hist_log_distance_start_2.pdf")

plot_grid(log_distance_start_layout_1,log_distance_start_layout_2,labels = "AUTO", label_size = 15,ncol = 1)
#ggsave(6*8inches)
############


######################################################################################################################
### Distribution of the log 10 (length) between the start position and the end position of a gene. 
######################################################################################################################

hist_M0.1_pos_logLength=ggplot(M0.1_pos%>% select(gene_length) %>% distinct(), aes(x = log10(gene_length))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "log10 (length)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+ylim(0,0.8)+xlim(1,7)

hist_M0.2_pos_logLength=ggplot(M0.2_pos %>% select(gene_length) %>% distinct(), aes(x = log10(gene_length))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "log10 (length)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+ylim(0,0.8)+xlim(1,7)

hist_M1.1_pos_logLength=ggplot(M1.1_pos %>% select(gene_length) %>% distinct(), aes(x = log10(gene_length))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "log10 (length)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+ylim(0,0.8)+xlim(1,7)

hist_M1.2_pos_logLength=ggplot(M1.2_pos %>% select(gene_length) %>% distinct(), aes(x = log10(gene_length))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "log10 (length)",y = "Density")+ theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+ylim(0,0.8)+xlim(1,7)

hist_M0.1_neg_logLength=ggplot(M0.1_neg %>% select(gene_length) %>% distinct(), aes(x = log10(gene_length))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.1", x = "log10 (length)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+ylim(0,0.8)+xlim(1,7)

hist_M0.2_neg_logLength=ggplot(M0.2_neg %>% select(gene_length) %>% distinct(), aes(x = log10(gene_length))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "log10 (length)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+ylim(0,0.8)+xlim(1,7)

hist_M1.1_neg_logLength=ggplot(M1.1_neg %>% select(gene_length) %>% distinct(), aes(x = log10(gene_length))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "log10 (length)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+ylim(0,0.8)+xlim(1,7)

hist_M1.2_neg_logLength=ggplot(M1.2_neg %>% select(gene_length) %>% distinct(), aes(x = log10(gene_length))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "log10 (length)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+ylim(0,0.8)+xlim(1,7)

logLength_layout_1=(hist_M0.1_pos_logLength|hist_M1.1_pos_logLength)/(hist_M0.2_pos_logLength|hist_M1.2_pos_logLength)+
  plot_annotation(title = "ER+")& theme(plot.title = element_text(hjust = 0.5,size = 12))
logLength_layout_1

#ggsave("hist_logLength_1.pdf")

logLength_layout_2=(hist_M0.1_neg_logLength|hist_M1.1_neg_logLength)/(hist_M0.2_neg_logLength|hist_M1.2_neg_logLength)+
  plot_annotation(title = "ER-")& theme(plot.title = element_text(hjust = 0.5,size = 12))
logLength_layout_2

#ggsave("hist_logLength_2.pdf")

plot_grid(logLength_layout_1,logLength_layout_2,labels = "AUTO", label_size = 15,ncol = 1)
#6*8 potrait
###########


######################################################################################################################
# Distribution of the relation to CpG island :five sections: Island, N shore, S shore, N shelf, S shelf, and No Info 
######################################################################################################################

## pie_M0.1_pos
pie_dat_M0.1_pos <- M0.1_pos %>% count(Relation_to_UCSC_CpG_Island) %>% mutate(percentage = n / sum(n) * 100)
pie_dat_M0.1_pos[1,1]="No Info"

## pie_M0.2_pos
pie_dat_M0.2_pos <- M0.2_pos %>% count(Relation_to_UCSC_CpG_Island) %>% mutate(percentage = n / sum(n) * 100)
pie_dat_M0.2_pos[1,1]="No Info"

## pie_M1.1_pos
pie_dat_M1.1_pos <- M1.1_pos %>% count(Relation_to_UCSC_CpG_Island) %>% mutate(percentage = n / sum(n) * 100)
pie_dat_M1.1_pos[1,1]="No Info"

## pie_M1.2_pos
pie_dat_M1.2_pos <- M1.2_pos %>% count(Relation_to_UCSC_CpG_Island) %>% mutate(percentage = n / sum(n) * 100)
pie_dat_M1.2_pos[1,1]="No Info"

pie_dat_M0.1_neg <- M0.1_neg %>% count(Relation_to_UCSC_CpG_Island) %>% mutate(percentage = n / sum(n) * 100)
pie_dat_M0.1_neg[1,1]="No Info"

pie_dat_M0.2_neg <- M0.2_neg %>% count(Relation_to_UCSC_CpG_Island) %>% mutate(percentage = n / sum(n) * 100)
pie_dat_M0.2_neg[1,1]="No Info"

pie_dat_M1.1_neg <- M1.1_neg %>% count(Relation_to_UCSC_CpG_Island) %>% mutate(percentage = n / sum(n) * 100)
pie_dat_M1.1_neg[1,1]="No Info"

pie_dat_M1.2_neg <- M1.2_neg %>% count(Relation_to_UCSC_CpG_Island) %>% mutate(percentage = n / sum(n) * 100)
pie_dat_M1.2_neg[1,1]="No Info"

#########funtion########
create_pie <- function(data_pie, title) {pie(data_pie$n, labels = paste0(data_pie$Relation_to_UCSC_CpG_Island," (",round(data_pie$percentage, 1), "%)"), col = c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f'), main = title,cex = 1,cex.main=1.5 )}

par(mfrow = c(2, 2), mar = c(2, 2, 2, 2), oma = c(2, 2, 2, 2))
#create_pie(pie_dat_M0.1_pos,"M0.1 pos")

pie_M0.1_pos <- create_pie(pie_dat_M0.1_pos,"M0.1")
pie_M1.1_pos <- create_pie(pie_dat_M1.1_pos,"M1.1")
pie_M0.2_pos <- create_pie(pie_dat_M0.2_pos,"M0.2")
pie_M1.2_pos <- create_pie(pie_dat_M1.2_pos,"M1.2")
mtext("ER+", outer = TRUE, cex = 1.4, font = 2)
mtext("A", side = 3, adj = 0, outer = TRUE, line = -1.5, cex = 2, font = 1)

#"pie_pos.pdf" 5*7 landscape
#pie_pos=(pie_M0.1_pos|pie_M1.1_pos)/(pie_M0.2_pos|pie_M1.2_pos)+
#  plot_annotation(title = "ER+")& theme(plot.title = element_text(hjust = 0.5,size = 12))
#pie_pos
#ggsave("pie_1_pos.pdf")

par(mfrow = c(2, 2), mar = c(2, 2, 2, 2), oma = c(2, 2, 2, 2))
create_pie(pie_dat_M0.1_neg,"M0.1")
create_pie(pie_dat_M1.1_neg,"M1.1")
create_pie(pie_dat_M0.2_neg,"M0.2")
create_pie(pie_dat_M1.2_neg,"M1.2")
mtext("ER-", outer = TRUE, cex = 1.4, font = 2)
mtext("B", side = 3, adj = 0, outer = TRUE, line = -1.5, cex = 2, font = 1)
#"pie_neg.pdf" 5*7 landscape
#ggsave("pie_2_neg.pdf")
############



######################################################################################################################
# Distribution of methylation mean values for the baseline and mediation models
######################################################################################################################

hist_M0.1_pos_MethyMean=ggplot(M0.1_pos, aes(x = Methyl_mean)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "Mean methylation",y = "Density") + theme_bw()+xlim(-5.0,5.0)+ylim(0,0.65)

hist_M0.2_pos_MethyMean=ggplot(M0.2_pos, aes(x = Methyl_mean)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "Mean methylation",y = "Density") + theme_bw()+xlim(-5.0,5.0)+ylim(0,0.65)

hist_M1.1_pos_MethyMean=ggplot(M1.1_pos, aes(x = Methyl_mean)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "Mean methylation",y = "Density") + theme_bw()+xlim(-5.0,5.0)+ylim(0,0.65)

hist_M1.2_pos_MethyMean=ggplot(M1.2_pos, aes(x = Methyl_mean)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "Mean methylation",y = "Density") + theme_bw()+xlim(-5.0,5.0)+ylim(0,0.65)

hist_M0.1_neg_MethyMean=ggplot(M0.1_neg, aes(x = Methyl_mean)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "Mean methylation",y = "Density") + theme_bw()+xlim(-5.0,5.0)+ylim(0,0.65)

hist_M0.2_neg_MethyMean=ggplot(M0.2_neg, aes(x = Methyl_mean)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "Mean methylation",y = "Density") + theme_bw()+xlim(-5.0,5.0)+ylim(0,0.65)

hist_M1.1_neg_MethyMean=ggplot(M1.1_neg, aes(x = Methyl_mean)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "Mean methylation",y = "Density") + theme_bw()+xlim(-5.0,5.0)+ylim(0,0.65)

hist_M1.2_neg_MethyMean=ggplot(M1.2_neg, aes(x = Methyl_mean)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "Mean methylation", y = "Density") + theme_bw()+xlim(-5.0,5.0)+ylim(0,0.65)

MethyMean_layout_1=(hist_M0.1_pos_MethyMean|hist_M1.1_pos_MethyMean)/(hist_M0.2_pos_MethyMean|hist_M1.2_pos_MethyMean)+
  plot_annotation(title = "ER+")& theme(plot.title = element_text(hjust = 0.5,size = 12))
MethyMean_layout_1

#ggsave("hist_MethylMean_1.pdf")

MethyMean_layout_2=(hist_M0.1_neg_MethyMean|hist_M1.1_neg_MethyMean)/(hist_M0.2_neg_MethyMean|hist_M1.2_neg_MethyMean)+
  plot_annotation(title = "ER-")& theme(plot.title = element_text(hjust = 0.5,size = 12))
MethyMean_layout_2

#ggsave("hist_MethylMean_2.pdf")

plot_grid(MethyMean_layout_1,MethyMean_layout_2,labels = "AUTO",  label_size = 15,ncol = 1)
#6*8 potrait ("hist_MethylMean.pdf")
###########




######################################################################################################################
# Distribution of the GC content
######################################################################################################################
#The figure shows the distribution of the GC content which is the percentage of C and G among the four nucleotide base pairs (A, T, C, and G) in the baseline and mediation models across the cancer subtypes.


p_M0.1_pos=ggplot(M0.1_pos %>% select(Gene...GC.content), aes(x = Gene...GC.content)) +
  geom_histogram(aes(y = ..density..), binwidth = 2.5, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "GC content",y = "Density") + theme_bw()+xlim(28,75)+ylim(0,0.1)

p_M0.2_pos=ggplot(M0.2_pos %>% select(Gene...GC.content) , aes(x = Gene...GC.content)) +
  geom_histogram(aes(y = ..density..), binwidth = 2.5, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "GC content",y = "Density") + theme_bw()+xlim(28,75)+ylim(0,0.1)

p_M1.1_pos=ggplot(M1.1_pos %>% select(Gene...GC.content) , aes(x = Gene...GC.content)) +
  geom_histogram(aes(y = ..density..), binwidth = 2.5, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "GC content",y = "Density") + theme_bw()+xlim(28,75)+ylim(0,0.1)

p_M1.2_pos=ggplot(M1.2_pos %>% select(Gene...GC.content) , aes(x = Gene...GC.content)) +
  geom_histogram(aes(y = ..density..), binwidth = 2.5, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "GC content",y = "Density") + theme_bw()+xlim(28,75)+ylim(0,0.1)

p_M0.1_neg=ggplot(M0.1_neg %>% select(Gene...GC.content) , aes(x = Gene...GC.content)) +
  geom_histogram(aes(y = ..density..), binwidth = 2.5, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "GC content",y = "Density") + theme_bw()+xlim(28,75)+ylim(0,0.1)

p_M0.2_neg=ggplot(M0.2_neg %>% select(Gene...GC.content) , aes(x = Gene...GC.content)) +
  geom_histogram(aes(y = ..density..), binwidth = 2.5, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "GC content",y = "Density") + theme_bw()+xlim(28,75)+ylim(0,0.25)

p_M1.1_neg=ggplot(M1.1_neg %>% select(Gene...GC.content) , aes(x = Gene...GC.content)) +
  geom_histogram(aes(y = ..density..), binwidth = 2.5, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "GC content", y = "Density") + theme_bw()+xlim(28,75)+ylim(0,0.1)

p_M1.2_neg=ggplot(M1.2_neg %>% select(Gene...GC.content) , aes(x = Gene...GC.content)) +
  geom_histogram(aes(y = ..density..), binwidth = 2.5, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "GC content",y = "Density") + theme_bw()+xlim(28,75)+ylim(0,0.25)

layout_1=(p_M0.1_pos|p_M1.1_pos)/(p_M0.2_pos|p_M1.2_pos)+
  plot_annotation(title = "ER+")& theme(plot.title = element_text(hjust = 0.5,size = 12))
layout_1
#ggsave("Figure_GC_content_1.pdf")

layout_2=(p_M0.1_neg|p_M1.1_neg)/(p_M0.2_neg|p_M1.2_neg)+
  plot_annotation(title = "ER-")& theme(plot.title = element_text(hjust = 0.5,size = 12))
layout_2
#ggsave("Figure_GC_content_2.pdf")

plot_grid(layout_1,layout_2,labels = "AUTO", label_size = 15,ncol = 1)
#6*8 potrait

#########



################### calculate average distance ############
# Mean of the log10 (distance) between the nearby CpG island and the CpG in the methylation probe.

M0.1_pos %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_cpG_mapinfo > -2) %>% pull() %>% summary()
M0.2_pos %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_cpG_mapinfo > -2) %>% pull() %>% summary()
M1.1_pos %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_cpG_mapinfo > -2) %>% pull() %>% summary()
M1.2_pos %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_cpG_mapinfo > -2) %>% pull() %>% summary()

M0.1_neg %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_cpG_mapinfo > -2) %>% pull() %>% summary()
M0.2_neg %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_cpG_mapinfo > -2) %>% pull() %>% summary()
M1.1_neg %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_cpG_mapinfo > -2) %>% pull() %>% summary()
M1.2_neg %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_cpG_mapinfo > -2) %>% pull() %>% summary()

# Mean of the log 10 (distance) between the CpG in the methylation probe and the start position of a gene.
M0.1_pos %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_mapinfo_geneStart > -2) %>% pull() %>% summary()
M0.2_pos %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_mapinfo_geneStart > -2) %>% pull() %>% summary()
M1.1_pos %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_mapinfo_geneStart > -2) %>% pull() %>% summary()
M1.2_pos %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_mapinfo_geneStart > -2) %>% pull() %>% summary()

M0.1_neg %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_mapinfo_geneStart > -2) %>% pull() %>% summary()
M0.2_neg %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_mapinfo_geneStart > -2) %>% pull() %>% summary()
M1.1_neg %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_mapinfo_geneStart > -2) %>% pull() %>% summary()
M1.2_neg %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_mapinfo_geneStart > -2) %>% pull() %>% summary()

# Mean of the log 10 (length) between the start position and the end position of a gene.
M0.1_pos %>% select(gene_length) %>% drop_na() %>% abs() %>% log10() %>% filter(gene_length > -2) %>% pull() %>% summary()
M0.2_pos %>% select(gene_length) %>% drop_na() %>% abs() %>% log10() %>% filter(gene_length > -2) %>% pull() %>% summary()
M1.1_pos %>% select(gene_length) %>% drop_na() %>% abs() %>% log10() %>% filter(gene_length > -2) %>% pull() %>% summary()
M1.2_pos %>% select(gene_length) %>% drop_na() %>% abs() %>% log10() %>% filter(gene_length > -2) %>% pull() %>% summary()

M0.1_neg %>% select(gene_length) %>% drop_na() %>% abs() %>% log10() %>% filter(gene_length > -2) %>% pull() %>% summary()
M0.2_neg %>% select(gene_length) %>% drop_na() %>% abs() %>% log10() %>% filter(gene_length > -2) %>% pull() %>% summary()
M1.1_neg %>% select(gene_length) %>% drop_na() %>% abs() %>% log10() %>% filter(gene_length > -2) %>% pull() %>% summary()
M1.2_neg %>% select(gene_length) %>% drop_na() %>% abs() %>% log10() %>% filter(gene_length > -2) %>% pull() %>% summary()

# Average value of methylation mean values in the baseline and mediation models 
#M0.1_pos %>% select(Methyl_mean) %>% summary()

summary(M0.1_pos$Methyl_mean)
summary(M0.2_pos$Methyl_mean)
summary(M1.1_pos$Methyl_mean)
summary(M1.2_pos$Methyl_mean)

## GC content
summary(M0.1_pos$Gene...GC.content)
summary(M0.2_pos$Gene...GC.content)
summary(M1.1_pos$Gene...GC.content)
summary(M1.2_pos$Gene...GC.content)

summary(M0.1_neg$Gene...GC.content)
summary(M0.2_neg$Gene...GC.content)
summary(M1.1_neg$Gene...GC.content)
summary(M1.2_neg$Gene...GC.content)

##################################
#
# Location_AllmodelType Figure #
#
##################################


######## 1. Load data files and packages  #########
#set the directory
library(data.table)
setwd("/Users/lianzuo/LZ/WayneSU/raw_Data_Methyl")
#read in the datasets
HumanMethInfo <- read.csv("GPL13534_HumanMethylation450_15017482_v.1.1 2.csv",skip = 7, header = TRUE)
dim(HumanMethInfo)
trios <- fread("trio.final.protein.coding.txt", data.table=F)
TCGA_meth <- fread("split.names.TCGA.meth.logit.txt", data.table=F)

#TC0 <- TCGA_meth%>%dplyr::select(Name=Row.names)
#df2 <- df%>%dplyr::select(Name, UCSC_RefGene_Group) ## 486428 rows and unique-485,592
#df2 <- df2%>%dplyr::filter(Name%in%TC0$Name) ##485,577 rows and unique-485,577
#TC0 <- TC0%>%left_join(df2, by="Name") ## Row.name and Loc 

####### 2. clean data #########

# find the meth rows in trios
trios_methrow <- trios$meth.row
#sort the meth data with the same rows in trios
meth_trio <- TCGA_meth[trios_methrow,]
# define the Rownames in sorted meth data as Name, using for merge later.
Name <- meth_trio$Row.names
# add the Name to tiros data
trios_add_Name=cbind(trios,Name)
head(trios_add_Name)
# add the trio rows as a column 
trios_add_Name$trios_row=c(1:nrow(trios_add_Name))

#select the two coloums in HumanMethInfo data,it is easy to calculate
info_1 <- HumanMethInfo %>% select(Name,UCSC_RefGene_Group)
# filter the matched Rownames in df_1 data
info_2 <- info_1 %>% filter(Name %in% trios_add_Name$Name)
#merge two dataframes together by the same Name
merge_df <- merge(trios_add_Name,info_2,by="Name",all=F)
####merge_df <- left_join(trios_add_Name,info_2,by="Name")
#sort the data in trios_row order(from 1 to 295492)
merge_trio_loc=merge_df %>% arrange(trios_row)
head(merge_trio_loc)
#df data using the v1.1 version,the problem is: trios_add_Name : count(Name)=295492,but has duplicates,there are 29462 Names not match in df_2 Name;
#df_2 :count(Name)=246070, no duplicates,the final location csv:nrow()=266030 .not equal to nrow(trios)
# 29462 Rownames in meth data are not in df(humanmethylation) data,that means they don't has UCSC_RefGene_Group data.
#trios[(221:230,1200:1219),]
# df data using the v1.2 version, the rownames in meth data are all in df(Humanmethylation) data.
dim(merge_df)# [295492,2]

table(is.na(match(trios_add_Name$Name,info_2$Name)))
#before:TRUE :29462(not match names) #after:FALSE :295492(all matched)

##create a summary table for the location
table(unlist(strsplit(as.character(merge_trio_loc$UCSC_RefGene_Group), ';')))

#write final to the location file
write.table(merge_trio_loc, file = "/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/final_protein_coding_location_LZ_03_19_25.csv", sep = "\t", row.names = FALSE,col.names = T, append = F, quote=FALSE)

# load the data


######## 3/16/2025 ###########
loc <- fread("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/final_protein_coding_location_LZ_03_19_25.csv",data.table = F)
pos_modeltype<-read.delim('/Users/lianzuo/LZ/WayneSU/raw_Data_Methyl/model.trio.MRGN.all.posER.reclassify2.txt')
#combine the pos modeltype with the location infos by trio_row
merge_pos_df=left_join(pos_modeltype,loc,by=c("trio_row"="trios_row"))
pos_df=merge_pos_df %>% separate_rows(UCSC_RefGene_Group, sep = ";")

#M1.1=dff_1 %>% filter(Inferred.Model3=="M1.1")
#table(M1.1$UCSC_RefGene_Group)

## Combine the 6 locations into 3
pos_df_2=pos_df %>% mutate(Group= case_when(
  UCSC_RefGene_Group %in% c("TSS200","TSS1500") ~ "TSS",
  UCSC_RefGene_Group %in% c("Body","1stExon") ~ "Body",
  UCSC_RefGene_Group %in% c("3'UTR","5'UTR") ~ "5'/3'UTR",
  TRUE ~NA_character_))
## barplot
#ggplot(dff_2)+geom_bar(aes(x=Inferred.Model3,fill=Group,color=Group),position = position_dodge(0.8),width = 0.7)
#+theme_bw()+theme(legend.position = "right")
#+labs(title = "ER+",x="Inferrence Models",y="Count (thousands)")
#+scale_y_continuous(limits = c(0, 180000), labels=seq(0,15,5),breaks = seq(0, 150000, 50000))


pos_df_counts <- pos_df_2 %>%
  group_by(Inferred.Model3, Group) %>%
  summarise(Count = n(), .groups = 'drop')
table(dff_2$Group)
df_counts$amount=rep(c(121627,362129,191201),9)
df_counts$percent=round((df_counts$Count/df_counts$amount)*100,digits = 1)

# Create the bar plot (Count)
ggplot(pos_df_counts, aes(x = Inferred.Model3, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "(ER+) Count of Locations by Model Type",x = "Model Type",y = "Count (thousand)",fill = "Group") +
  theme_bw()+scale_y_continuous(limits = c(0, 180000), labels=seq(0,15,5),breaks = seq(0, 150000, 50000))+
  theme(legend.position = "top",plot.title = element_text(size=14,hjust = 0.5))+
  scale_fill_manual(name = "Location",  values = c("TSS" = "#377eb8", "Body" = "#4daf4a", "5'/3'UTR" = "#e7298a"))  # Optional: Customize colors)

#bar plot _ log10(Count) 
ggplot(pos_df_counts, aes(x = Inferred.Model3, y = log10(Count), fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "(ER+) Count of Locations by Model Type",x = "Model Type",y = "log10(Count)",fill = "Group") +
  theme_bw()+scale_y_continuous(limits = c(0, 6), labels=seq(0,6,1),breaks = seq(0, 6, 1))+
  theme(legend.position = "top",plot.title = element_text(size=14,hjust = 0.5))+
  scale_fill_manual(name = "Location",  values = c("TSS" = "#377eb8", "Body" = "#4daf4a", "5'/3'UTR" = "#e7298a"))  # Optional: Customize colors)

# Create the bar plot (Percentage)
ggplot(pos_df_counts, aes(x = Inferred.Model3, y = percent, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "( ER+ ) Count of Locations by Model Type",x = "Model Type",y = "Percentage",fill = "Group") +theme_bw()+
  theme(legend.position = "top",plot.title = element_text(size=14,hjust = 0.5))+
  scale_fill_manual(name = "Location",  # Change the legend title here
                    values = c("TSS" = "#377eb8", "Body" = "#4daf4a", "5'/3'UTR" = "#e7298a")  # Optional: Customize colors
  )
##add labels to a dodged barplot
#+geom_text(aes(label=Count,group=Group),position = position_dodge(0.8),vjust=-0.3,size=3.5)


##### neg ######

neg_df0 <-read.delim('/Users/lianzuo/LZ/WayneSU/raw_Data_Methyl/model.trio.MRGN.all.negER.reclassify2.txt')
neg_df=left_join(neg_df0,loc,by=c("trio_row"="trios_row"))
neg_dff_2=neg_df %>% separate_rows(UCSC_RefGene_Group, sep = ";")

## Combine the 6 locations into 3
neg_dff_2=neg_dff_2 %>% mutate(Group= case_when(
  UCSC_RefGene_Group %in% c("TSS200","TSS1500") ~ "TSS",
  UCSC_RefGene_Group %in% c("Body","1stExon") ~ "Body",
  UCSC_RefGene_Group %in% c("3'UTR","5'UTR") ~ "5'/3'UTR",
  TRUE ~NA_character_))

neg_df_counts <- neg_dff_2 %>%
  group_by(Inferred.Model3, Group) %>%
  summarise(Count = n(), .groups = 'drop')

table(neg_dff_2$Group)
neg_df_counts$amount=rep(c(121627,362079,191129),9)
neg_df_counts$percent=round((neg_df_counts$Count/neg_df_counts$amount)*100,digits = 1)

# ER- Create the bar plot
# Create the bar plot
ggplot(neg_df_counts, aes(x = Inferred.Model3, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "ER- :Count of Locations by Model Type",x = "Model Type",y = "Count (thousand)",fill = "Group") +
  theme_bw()+scale_y_continuous(limits = c(0, 150000), labels=seq(0,15,5),breaks = seq(0, 150000, 50000))+
  theme(legend.position = "top",plot.title = element_text(size=14,hjust = 0.5))+
  scale_fill_manual(name = "Location",  # Change the legend title here
                    values = c("TSS" = "#377eb8", "Body" = "#4daf4a", "5'/3'UTR" = "#e7298a")  # Optional: Customize colors
  )

# Create the bar plot(Pecentage)
ggplot(neg_df_counts, aes(x = Inferred.Model3, y = percent, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "( ER- ) Count of Locations by Model Type",x = "Model Type",y = "Percentage",fill = "Group") +theme_bw()+
  theme(legend.position = "top",plot.title = element_text(size=14,hjust = 0.5))+
  scale_fill_manual(name = "Location",  # Change the legend title here
                    values = c("TSS" = "#377eb8", "Body" = "#4daf4a", "5'/3'UTR" = "#e7298a")  # Optional: Customize colors
  )

