

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

## Figure8_hist_log10_distance_between_CpG_in_methy_probe_and_startposition_of_a_gene.pdf



