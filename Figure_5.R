# ---- Script Metadata ---- 
# Title: Figure_5.R
# Author: Trish Milletich, PhD
# Data: 2026-07-31
# --------------------------------

library(ggplot2); library(ggpubr)
library(dplyr); library(tidyr)

#Load Data 
sIgA = read.csv("sIgA_NLF_Saliva.csv")

#Create Ratios to compare Days 
sIgA_2 = sIgA %>% 
  group_by(Subject.ID, Specimen.Type) %>%
  summarize(    
    Normalized_sIgA = { 
      day3  <- Spike.sIgA.Normalized[Days == "Day 3"]
      day15 <- Spike.sIgA.Normalized[Days == "Day 15"]
      
      if (length(day3) == 1 && length(day15) == 1) {
        day3 / day15
      } else { NA_real_  }
    },
    
    Unadjusted_sIgA = { 
      day3  <- Spike.sIgA.Unadjusted[Days == "Day 3"]
      day15 <- Spike.sIgA.Unadjusted[Days == "Day 15"]
      
      if (length(day3) == 1 && length(day15) == 1) {
        day3 / day15
      } else { NA_real_  }
    },
    
    Total_IgA = { 
      day3  <- Total.IgA[Days == "Day 3"]
      day15 <- Total.IgA[Days == "Day 15"]
      
      if (length(day3) == 1 && length(day15) == 1) {
        day3 / day15
      } else { NA_real_     }
    },
    
    Unadjusted_IgG = { 
      day3  <- Spike.IgG[Days == "Day 3"]
      day15 <- Spike.IgG[Days == "Day 15"]
      
      if (length(day3) == 1 && length(day15) == 1) {
        day3 / day15
      } else { NA_real_     }
    },
    
    .groups = "drop"
  )


sIgA_2 = subset(sIgA_2, is.na(sIgA_2$Unadjusted_sIgA) == F)

sIgA_2 = sIgA_2%>%
  pivot_longer(
    cols = c("Normalized_sIgA", "Unadjusted_sIgA", "Total_IgA", "Unadjusted_IgG"),
    names_to = "Sample",  values_to = "FC"   ) 

#Clean Data
sIgA_2$Specimen.Type = gsub("PooledNostril", "NLF", sIgA_2$Specimen.Type)
sIgA_2$Sample = gsub("_", "\n", sIgA_2$Sample)
sIgA_2$Sample_2 = paste(sIgA_2$Sample, sIgA_2$Specimen.Type)
sIgA_2$Sample_2 = gsub("Unadjusted\n", "", sIgA_2$Sample_2)
sIgA_2$Sample_2 = gsub("Normalized\n", "", sIgA_2$Sample_2)
sIgA_2$Sample_2 = gsub("Total\n", "", sIgA_2$Sample_2)

sIgA_2$Sample = factor(sIgA_2$Sample, levels = 
                         c("Total\nIgA", "Unadjusted\nsIgA", "Normalized\nsIgA", 'Unadjusted\nIgG'))

######################## #
# ---- T Test ---- 
TTest_pvalues = sIgA_2 %>%
  group_by(Sample, Specimen.Type) %>%
  summarise(P = round(t.test(log2(FC), mu = 0)$p.value,2),
            .groups = "drop_last")

TTest_pvalues$P2 = paste("P=", TTest_pvalues$P, sep = "")
TTest_pvalues$Specimen.Type = gsub("PooledNostril", "NLF", TTest_pvalues$Specimen.Type)
sIgA_2$Sample_2 = gsub("sIgA", "IgA", sIgA_2$Sample_2)

sIgA_2$Sample = factor(sIgA_2$Sample, levels = c("Normalized\nsIgA", "Unadjusted\nsIgA", "Unadjusted\nIgG", "Total\nIgA") )
TTest_pvalues$Sample = factor(TTest_pvalues$Sample, levels = c("Normalized\nsIgA", "Unadjusted\nsIgA", "Unadjusted\nIgG", "Total\nIgA") )

FC_graph = ggplot(sIgA_2, aes(x = Sample, y = log2(FC), fill = Sample_2)) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_boxplot(color = "black", outlier.shape = NA) + 
  geom_jitter(width = 0.1, height = 0, alpha = 0.5) + 
  facet_wrap(~Specimen.Type, nrow = 2, scales = "free") + 
  theme_bw() + 
  scale_fill_manual(breaks = c("sIgA NLF", "IgA NLF", "IgG NLF",
                               "sIgA Saliva", "IgA Saliva", "IgG Saliva"), 
                    values = c("#CC6677", "#CC6677", "#7B0E45", 
                               "#88CCEE", "#88CCEE", "#065882"),
                    name = "")+ 
  geom_label(data = TTest_pvalues,
             aes(x = Sample, y = 3, label = P2), 
             fill = "white", vjust = 1) + 
  theme(strip.background = element_rect(fill = "white"), 
        legend.position = "top", 
        axis.title.x = element_blank()) + 
  ylab("Day3/Day15 [log2]")

###########################################################
# Load and Clean Data 
All_data = read.csv("Saliva_Serum_Nasal.csv")
All_data = subset(All_data, All_data$Ag == "Spike")
All_data = subset(All_data, All_data$Days != "Day 1" &
                    All_data$Sample.Type != "Serum")

All_data$Ag = NULL; All_data$sample.ID = NULL; All_data$Visit = NULL
colnames(All_data) = c("Ab", "Subject.ID", "Unadjusted", "Total", "Normalized", "Specimen.Type", "Days" )

All_data_wide = reshape(All_data, 
                        idvar = c("Subject.ID", "Days", "Specimen.Type"), 
                        timevar = "Ab", 
                        direction = "wide")

All_data_wide$Total.IgG = NULL; All_data_wide$Normalized.IgG = NULL
All_data_wide$Specimen.Type = gsub("PooledNostril", "NLF", All_data_wide$Specimen.Type)

#Clean sIgA to allow merging
colnames(sIgA) = c("Subject.ID", "Days", "Specimen.Type","Ab_ST_day",
                   "Normalized.sIgA", "Unadjusted.sIgA","Total.IgA.ELISA", "Unadjusted.IgG.ELISA" )
sIgA$Ab_ST_day = NULL

Merge_data= merge(All_data_wide, sIgA, all = T)
Merge_data = data.frame(Merge_data)

comparison_list = list(
  #Adjusted sIgA vs IgA 
  c("Normalized.IgA", "Normalized.sIgA"),
  #Unadjusted sIgA vs IgA 
  c("Unadjusted.IgA", "Unadjusted.sIgA"),
  #Total Ab
  c("Total.IgA", "Total.IgA.ELISA"),
  #IgG 
  c( "Unadjusted.IgG", "Unadjusted.IgG.ELISA")
)

current_sample = "Saliva";current_compare = comparison_list[1]
Total_plot_list = c()

for (current_compare in comparison_list) {
  #Format data based on current comparison
  Merge_data$Xvalue = Merge_data[,current_compare[1]]
  Merge_data$Yvalue = Merge_data[,current_compare[2]]
  
  if (current_compare[1] == "Unadjusted.IgG") {
    color_list = c( "#7B0E45",  "#065882")
    Ylabel = "Unadjusted IgG ECLIA"
    Xlabel = current_compare[2]
    shape = 17
  } else if (current_compare[1] == "Total.IgA") {
    color_list = c("#CC6677","#88CCEE")
    Ylabel = "Total IgA ECLIA"
    Xlabel = current_compare[2]
    shape = 16
  } else  {
    color_list = c("#CC6677","#88CCEE")
    Ylabel = paste(current_compare[1], "ECLIA")
    Xlabel = paste(current_compare[2], "ELISA")
    shape = 16
  }
  
  plot_data = subset(Merge_data, is.na(Merge_data$Xvalue) == F & 
                       is.na(Merge_data$Yvalue) == F)

  plot_data <- plot_data %>%
    group_by(Subject.ID, Specimen.Type) %>%
    summarise(
      Xvalue = if(n() > 1) mean(Xvalue, na.rm = TRUE) else Xvalue,
      Yvalue = if(n() > 1) mean(Yvalue, na.rm = TRUE) else Yvalue,
      .groups = "drop"
    )

  
  current_plot = ggplot(plot_data,
                        aes(x = log2(Xvalue),
                            y = log2(Yvalue),
                            color = Specimen.Type) ) +
    geom_point(data = plot_data,
               aes(x = log2(Xvalue),
                   y = log2(Yvalue)),
               shape = shape) +
    geom_smooth(method = "lm", formula = "y~x", color = "black") +
    stat_cor( aes(label = gsub("R", "r",
                               paste(..r.label.., ..p.label.., sep = "~`,`~"))),
              method = "pearson", p.accuracy = 0.001, size = 5, color = "black",
              vjust = 0.5) +
    theme_bw() +
    facet_wrap(~Specimen.Type, scales = "free", nrow = 2) +
    ylab(paste(gsub("\\.", " ", Ylabel), "[log2]")) +
    xlab(paste(gsub("\\.", " ", Xlabel), "[log2]")) +
    scale_color_manual(breaks = c("NLF", "Saliva"),
                       values = color_list) +
    theme(strip.background = element_rect(fill = "white")); current_plot
  
  Total_plot_list = c(Total_plot_list, current_plot)

}

jpeg("Figure5.jpeg", res = 600, height = 4000, width = 8500)
ggarrange(FC_graph, plot.new(), Total_plot_list[[1]], Total_plot_list[[2]], 
          Total_plot_list[[4]], Total_plot_list[[3]],
          ncol = 6, widths = c(1.5, 0.1, 1, 1, 1, 1), 
          common.legend = T,
          align = "h")
dev.off()
