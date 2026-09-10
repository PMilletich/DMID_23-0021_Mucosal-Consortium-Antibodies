# ---- Script Metadata ---- 
# Title: Supplemental_Figure_2.R
# Author: Trish Milletich, PhD
# Data: 2026-07-31
# --------------------------------

library(ggplot2); library(ggpubr)
library(tidyr); library(dplyr)

All_data = read.csv("Saliva_Serum_Nasal.csv")
All_data = subset(All_data, All_data$Ag == "Spike")
All_data = subset(All_data, All_data$Days == "Day 1")
All_data$Sample.Type = gsub("PooledNostril", "NLF", All_data$Sample.Type)
All_data = subset(All_data, All_data$Sample.Type %in% c("Saliva", "NLF", 'Serum'))

####################################################### #
# ---- IgG/IgA for all Participants ----
####################################################### #
summary_ID = All_data %>% 
  group_by(Sample.Type, SN) %>%
  dplyr::summarise(.groups = "drop_last", 
                   Total_IgG.IgA = 
                     Total.ug.mL[Ab == "IgG"]/Total.ug.mL[Ab == "IgA"],
                   Unadjusted_IgG.IgA = 
                     Specific.AU.mL[Ab == "IgG"]/Specific.AU.mL[Ab == "IgA"], 
                   Normalized_IgG.IgA = 
                     Specific.Total.AU.ug[Ab == "IgG"]/Specific.Total.AU.ug[Ab == "IgA"], ) 

summary_ID_long = pivot_longer(summary_ID, 
                          cols = c('Total_IgG.IgA', 'Unadjusted_IgG.IgA', 'Normalized_IgG.IgA'),
                          names_to = "Measurement", 
                          values_to = "Titer")

summary_ID_long = data.frame(summary_ID_long)
head(summary_ID_long)

summary_ID_long$ST_Measurement = paste(summary_ID_long$Sample.Type, summary_ID_long$Measurement)
summary_ID_long = subset(summary_ID_long, ! summary_ID_long$ST_Measurement %in% 
                           c('Serum Total_IgG.IgA', 'Serum Normalized_IgG.IgA'))
summary_ID_long$Measurement = gsub("_IgG.IgA", "", summary_ID_long$Measurement)
  
summary_ID_long$Measurement = factor(summary_ID_long$Measurement, 
                                     levels = c("Total", "Unadjusted", "Normalized"))

####################################################### #
# Average across Samples and Measurements  ----
####################################################### #
label_data = summary_ID %>% 
  group_by(Sample.Type) %>%
  dplyr::summarise(.groups = "drop_last", 
                   Total_IgG.IgA = 
                     paste(round(mean((Total_IgG.IgA), na.rm = T), 2),
                           round(sd((Total_IgG.IgA), na.rm = T),2), sep = "±"),
                   Unadjusted_IgG.IgA = 
                     paste(round(mean((Unadjusted_IgG.IgA), na.rm = T), 2),
                           round(sd((Unadjusted_IgG.IgA), na.rm = T),2), sep = "±"), 
                   Normalized_IgG.IgA = 
                     paste(round(mean((Normalized_IgG.IgA), na.rm = T), 2),
                           round(sd((Normalized_IgG.IgA), na.rm = T),2), sep = "±") )

label_data = pivot_longer(label_data, 
                          cols = c('Total_IgG.IgA', 'Unadjusted_IgG.IgA', 'Normalized_IgG.IgA'),
                          values_to = "Label", 
                          names_to = "Measurement")
label_data = data.frame(label_data)
label_data$Measurement = gsub("_IgG.IgA", "", label_data$Measurement)
label_data$ST_Measurement = paste(label_data$Sample.Type, label_data$Measurement)
label_data = subset(label_data, ! label_data$ST_Measurement %in% 
                           c('Serum Total', 'Serum Normalized'))
label_data$Measurement = factor(label_data$Measurement, 
                                     levels = c("Total", "Unadjusted", "Normalized"))


jpeg("Supplemental_2.jpeg", res = 600, height = 3000, width = 6000)
ggplot(summary_ID_long, aes(x = Sample.Type, y = (Titer), fill = Sample.Type)) +
  geom_hline(yintercept = 1, linetype = "dashed") + 
  facet_wrap(~Measurement, ncol = 3, scales = "free_x") + 
  geom_boxplot(color = "black") + 
  geom_label(data = label_data, aes(x = Sample.Type, y = 1000, label = Label), 
             fill = "white", vjust = 1) + 
  theme_bw() + 
  ylab("IgG/IgA") + 
  theme(axis.title.x = element_blank(), 
        strip.background = element_rect(fill = "white"), 
        legend.position = "bottom") + 
  scale_fill_manual(breaks = c("NLF", "Saliva", "Serum"), 
                    values = c("#88CCEE", "#CC6677", "grey30"),  
                    name = "Specimen") + 
  scale_y_continuous(
    transform = pseudo_log_trans(base = 10, sigma = 1), 
    breaks = c( 0, 1, 10, 100, 1000))
dev.off()
