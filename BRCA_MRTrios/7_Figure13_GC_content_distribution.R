## Figure13_GC_content_distribution.pdf

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
