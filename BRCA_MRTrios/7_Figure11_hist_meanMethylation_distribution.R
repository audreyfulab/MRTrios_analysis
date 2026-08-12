## Figure11_hist_meanMethylation_distribution.pdf

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


