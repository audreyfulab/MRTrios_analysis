## Figure10_log10_length_between_startposition_and_endposition_of_a_gene.pdf


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

