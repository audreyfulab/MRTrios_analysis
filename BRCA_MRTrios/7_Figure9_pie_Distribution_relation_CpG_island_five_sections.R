## Figure9_pie_Distribution_relation_CpG_island_five_sections

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
