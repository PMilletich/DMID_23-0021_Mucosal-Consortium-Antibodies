# ---- Script Metadata ---- 
# Title: Figure_6.R
# Author: Trish Milletich, PhD
# Data: 2026-08-03
# --------------------------------

library(ggplot2);library(ggpubr)
library(pheatmap)

High_Binding = c("CV06FH003", "CV06BW010", "CV06FH011", "CV06FH005", "CV06HC003", "CV06UM009")
Low_Binding = c( "CV06HC001", "CV06HC004", "CV06HC014", "CV06UM014")

##########
# Binding Antibodies 
Binding_data = read.csv("Saliva_Serum_Nasal.csv")

##########
All_Neut = read.csv("Neutralizing_Saliva_Nasal.csv")
All_Neut$Trend = ifelse(All_Neut$SN %in% High_Binding, "High", "Low")
All_Neut$PseudoSpike_ID50_2 = ifelse(All_Neut$PseudoSpike_ID50 == "<10", 0,
                                      ifelse(All_Neut$PseudoSpike_ID50 == "nt", NA,  
                                             All_Neut$PseudoSpike_ID50))
All_Neut$PseudoSpike_ID50_2 = as.numeric(All_Neut$PseudoSpike_ID50_2)

All_Neut$Non.Specific = ifelse(All_Neut$Non.Specific == "Y", "Yes", "No")

All_Neut$Day_NSB = paste(All_Neut$Days, All_Neut$Non.Specific)
All_Neut$ST_NSB = paste(All_Neut$Sample.Type, All_Neut$Non.Specific)

################################### # 
# ---- Scatter Plot ----
Raw_Count = ggplot() + 
  geom_jitter(data = All_Neut, aes(x = Trend, y = PseudoSpike_ID50_2, 
                                   shape = Days, color = Sample.Type, 
                                   fill = ST_NSB), 
              height = 0, width = 0.3, size = 3,
              alpha = 0.7) + 
  theme_bw() + 
  facet_wrap(~Sample.Type, scales = "free_x") + 
  scale_color_manual(breaks = c("Saliva", "NLF"), 
                     values = c("#CC6677", "#88CCEE"),  
                    name = "Specimen") + 
  
  scale_fill_manual(breaks = c("Saliva No",  "Saliva Yes", "NLF No" ), 
                     values = c("#CC6677", "white", "#88CCEE"),  
                     guide = "none") + 
  scale_shape_manual(breaks = c("Day 1", "Day 3", "Day 15"), 
                     values = c(21, 22,  24), 
                     name = "Day") + 
  ylab("Pseudospike ID50") + 
  theme(axis.title.x = element_blank(), 
        strip.background = element_rect(fill = "white")); Raw_Count

#############################################################
## Day comparisons 
Clean_data_IC50 = All_Neut[All_Neut$Non.Specific == "No",
                           c("Days", "SN", "Sample.Type", "PseudoSpike_ID50")]
Clean_data_IC50 = subset(Clean_data_IC50, ! Clean_data_IC50$PseudoSpike_ID50 %in% c("nt", "<10"))
Clean_data_IC50$PseudoSpike_ID50 = as.numeric(Clean_data_IC50$PseudoSpike_ID50)

Clean_data_IC50 = reshape(Clean_data_IC50, idvar = c("SN", "Sample.Type"), timevar = "Days", direction = "wide")
colnames(Clean_data_IC50) = gsub("PseudoSpike_ID50.", "", colnames(Clean_data_IC50))
colnames(Clean_data_IC50) = gsub(" ", "", colnames(Clean_data_IC50))

Clean_data_IC50$Day.1.3 = Clean_data_IC50$Day1/Clean_data_IC50$Day3
Clean_data_IC50$Day.1.15 = Clean_data_IC50$Day1/Clean_data_IC50$Day15

Clean_data_IC50 = Clean_data_IC50[,c("SN", "Sample.Type", "Day.1.3", "Day.1.15")]
Clean_data_IC50_long = reshape(Clean_data_IC50, direction = "long", 
                               varying = c("Day.1.3", "Day.1.15"),
                               idvar = c("SN", "Sample.Type"),
                               timevar = "day", 
                               v.names = "FC",
                               times = c("Day1/Day3", "Day1/Day15"))

Clean_data_IC50_long = subset(Clean_data_IC50_long, is.na(Clean_data_IC50_long$FC) == F)

Clean_data_IC50_long$day = factor(Clean_data_IC50_long$day, 
                                  levels = c("Day1/Day3", "Day1/Day15"))
t.test_labels = data.frame(Clean_data_IC50_long %>%
  group_by(day, Sample.Type) %>%
  dplyr::summarise(P = t.test(log2(FC), mu = 0)$p.value,
                   .groups = "drop_last"))
t.test_labels$P2 = paste("p=", round(t.test_labels$P,2), sep = "")
#Remove Saliva since there's only two participants 
t.test_labels = subset(t.test_labels, t.test_labels$Sample.Type == "NLF")

Time = ggplot(Clean_data_IC50_long, aes(x = day, y = log2(FC), fill = Sample.Type)) + 
  geom_hline(yintercept = 0, color = "grey45") + 
  geom_boxplot(position = position_dodge(width = 0.5), 
               width = 0.3, color = "black") + 
  geom_point(position = position_dodge(width = 0.5)) + 
  theme_bw() + 
  geom_label(data = t.test_labels, 
             aes(x = day, y = 1, label = P2),
             show.legend = F, vjust = 1, fill = "white")  + 
  facet_wrap(~Sample.Type) + 
  scale_fill_manual(breaks = c("NLF", "Saliva"), 
                     values = c("#88CCEE", "#CC6677"),  
                     name = "Specimen") + 
  theme(axis.title.x = element_blank(), 
        strip.background = element_rect(fill = "white")); Time

#############################
#Correlate Nasal Nasal_Neut with other Nasal Antibodies 
Binding_Total = subset(Binding_data, Binding_data$Sample.Type %in% c("PooledNostril", "Saliva"))
Binding_Total = subset(Binding_Total, Binding_Total$Ag == "Spike")
Binding_Total = subset(Binding_Total, Binding_Total$SN %in% All_Neut$SN)

Binding_Total$Ag_Ab = paste(Binding_Total$Ab, Binding_Total$Sample.Type, sep = ":")
Binding_Total= Binding_Total[,c("Ag_Ab", 'SN', "Days", "Specific.Total.AU.ug", "Specific.AU.mL")]

Clean_data_IC50 = All_Neut[All_Neut$Non.Specific == "No",
                           c("Days", "SN", "Sample.Type", "PseudoSpike_ID50")]
Clean_data_IC50 = subset(Clean_data_IC50, ! Clean_data_IC50$PseudoSpike_ID50 %in% c("nt", "<10"))
Clean_data_IC50$PseudoSpike_ID50 = as.numeric(Clean_data_IC50$PseudoSpike_ID50)
Clean_data_IC50 = reshape(Clean_data_IC50, idvar = c("SN", "Days"), timevar = "Sample.Type", direction = "wide")
Neut_Binding = merge(Clean_data_IC50, Binding_Total)
head(Neut_Binding)

R_data = matrix(ncol = 2, nrow = 4)
rownames(R_data) = c("Unadjusted IgA", "Unadjusted IgG", 
                           "Normalized IgA", "Normalized IgG")
colnames(R_data) = c("NLF", "Saliva")
P_data = R_data

current_Group = unique(Neut_Binding$Ag_Ab)[4]
for (current_Group in unique(Neut_Binding$Ag_Ab)) {
  current_subset = subset(Neut_Binding, Neut_Binding$Ag_Ab == current_Group)
  current_ab = ifelse(grepl("IgG", current_Group), "IgG", "IgA")

  if (grepl("Nostril", current_Group)) {
    Unadjusted_Spearman = cor.test(log2(current_subset$PseudoSpike_ID50.NLF), 
                                        log2(current_subset$Specific.AU.mL), 
                                   method = "pearson")
    
    Normalized_Spearman = cor.test(log2(current_subset$PseudoSpike_ID50.NLF), 
                                        log2(current_subset$Specific.Total.AU.ug), 
                                   method = "pearson")
    
    R_data[paste("Unadjusted", current_ab), "NLF"] = round(Unadjusted_Spearman$estimate,3)
    P_data[paste("Unadjusted", current_ab), "NLF"] = 
      paste("r=",round(Unadjusted_Spearman$estimate,3), 
            "\np=",   round(Unadjusted_Spearman$p.value,3), sep = "")
    
    R_data[paste("Normalized", current_ab), "NLF"] = round(Normalized_Spearman$estimate,3)
    P_data[paste("Normalized", current_ab), "NLF"] =  
      paste("r=",round(Normalized_Spearman$estimate,3), 
            "\np=",   round(Normalized_Spearman$p.value,3), sep = "")
    
  } else {
    Unadjusted_Spearman = cor.test(log2(current_subset$PseudoSpike_ID50.Saliva), 
                                        log2(current_subset$Specific.AU.mL), 
                                   method = "pearson")
    Normalized_Spearman = cor.test(log2(current_subset$PseudoSpike_ID50.Saliva), 
                                        log2(current_subset$Specific.Total.AU.ug), 
                                   method = "pearson")
    
    R_data[paste("Unadjusted", current_ab), "Saliva"] = round(Unadjusted_Spearman$estimate,3)
    P_data[paste("Unadjusted", current_ab), "Saliva"] = paste("r=",round(Unadjusted_Spearman$estimate,3), 
                                                              "\np=",  
                                                              round(Unadjusted_Spearman$p.value,3),
                                                              sep = "")
    
    R_data[paste("Normalized", current_ab), "Saliva"] = round(Normalized_Spearman$estimate,3)
    P_data[paste("Normalized", current_ab), "Saliva"] =  paste("r=",round(Normalized_Spearman$estimate,3), 
                                                               "\np=", 
                                                               round(Normalized_Spearman$p.value,3),
                                                               sep = "")
  }
}

# Only used to see what measures are significant to neutralizing
pheatmap(R_data, cluster_rows = F, cluster_cols = F, display_numbers = P_data, 
         angle_col = 0, border_color = "black", gaps_col = 1,
         gaps_row = 2, number_color = "black")

############################# #
# Correlation scatter of significant measurements 
############################# #

Spec_Nasal = ggplot(Neut_Binding[Neut_Binding$Ag_Ab == "IgA:PooledNostril",],
       aes(x = log2(PseudoSpike_ID50.NLF), y = log2(Specific.AU.mL))) +
  geom_smooth(method = "lm", formula = "y ~ x", color = "#88CCEE") +
  stat_cor( p.accuracy = 0.001, size = 4.5, 
            aes(label = gsub("R", "r", 
                             paste(..r.label.., ..p.label.., sep = "~`,`~"))),
            method = "pearson") + 
  geom_point(color = "#88CCEE", aes(shape = Days),
             size = 3) +
  theme_bw() +
  scale_shape_manual(breaks = c("Day 1", "Day 3", "Day 15"), 
                     values = c(16, 15,  17), 
                     name = "NLF") + 
  xlab("PseudoSpike ID50 [log2]") +
  ylab("Unadjusted IgA [log2]"); Spec_Nasal

Spec_IgG_Nasal = ggplot(Neut_Binding[Neut_Binding$Ag_Ab == "IgG:PooledNostril",],
                    aes(x = log2(PseudoSpike_ID50.NLF), y = log2(Specific.AU.mL))) +
  geom_smooth(method = "lm", formula = "y ~ x", color = "#88CCEE") +
  stat_cor( p.accuracy = 0.001, size = 4.5, 
            aes(label = gsub("R", "r", 
                             paste(..r.label.., ..p.label.., sep = "~`,`~"))),
            method = "pearson") + 
  geom_point(color = "#88CCEE", aes(shape = Days),
             size = 3) +
  theme_bw() +
  scale_shape_manual(breaks = c("Day 1", "Day 3", "Day 15"), 
                     values = c(16, 15,  17), 
                     name = "NLF") + 
  xlab("PseudoSpike ID50 [log2]") +
  ylab("Unadjusted IgG [log2]"); Spec_IgG_Nasal

jpeg("Figure6.jpeg", res = 600, height = 3000, width = 6500)
ggarrange(ggarrange(Raw_Count, plot.new(), Time, 
                    nrow = 3, heights = c(1, 0.1, 1), common.legend = T), 
          plot.new(),
          ggarrange(Spec_Nasal, Spec_IgG_Nasal, 
                    nrow = 2, common.legend = T), 
          ncol = 3, widths = c(1, 0.1, 1))
dev.off()
