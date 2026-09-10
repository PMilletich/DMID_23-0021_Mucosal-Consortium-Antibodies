# ---- Script Metadata ---- 
# Title: Figure_2.R
# Author: Trish Milletich, PhD
# Data: 2026-07-31
# - ------------------------------ -

library(ggplot2); library(ggpubr)
library(tidyr); library(dplyr)

#Upload data and subset to nasal baseline samples for the Spike Antigen
All_data = read.csv("Saliva_Serum_Nasal.csv")
All_data = subset(All_data, All_data$Days == "Day 1" & 
                    All_data$Sample.Type %in% c("Left", "Right", "PooledNostril") &
                    All_data$Ag == "Spike")
All_data$Sample_Ab = paste(All_data$Sample.Type, All_data$Ab, sep = ":")

#Create list for plots 
Boxplot_list = c()
Left.Right_list = c() 
Pooled_list = c()

#For each measurement type....
for (current_measurement in c("Total", "Unadjusted", "Normalized")) {
  if (current_measurement == "Total") { 
    current_subset = All_data[,c("SN", "Ab", "Sample.Type",  "Sample_Ab", "Total.ug.mL")]
    titer = "(ug/mL)"
  } else if (current_measurement == "Unadjusted") { 
    current_subset = All_data[,c("SN", "Ab", "Sample.Type",  "Sample_Ab", "Specific.AU.mL")]
    titer = "(AU/mL)"
  } else if (current_measurement == "Normalized") { 
    current_subset = All_data[,c("SN", "Ab", "Sample.Type", "Sample_Ab", "Specific.Total.AU.ug")]
    titer = "(AU/ug)"
  } else {
    print("uh oh")
  }
  # Rename columns and clean data 
  colnames(current_subset)= c("SN",  "Ab", "Sample.Type", "Sample_Ab", "Titer")
  current_subset$Sample.Type = gsub("Nostril", "", current_subset$Sample.Type)
  current_subset$Sample.Type = factor(current_subset$Sample.Type, 
                                levels = c("Left", "Right", "Pooled"))
  current_subset$Ab = paste(current_measurement, current_subset$Ab)
  
  # ################# #
  # ---- Boxplots of means ----
  # ################# #
  All3 = ggplot(current_subset, aes(x = Sample.Type, y = log2(Titer))) + 
    geom_boxplot(color = "black") + 
    theme_bw() + 
    facet_wrap(~Ab) + 
    scale_fill_manual(breaks = c("Left", "Right", "Pooled"), 
                       values = c("#CC79A7", "#009E73", "#88CCEE")) + 
    stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))), 
                       vjust = 1, size = 5) + 
    theme(strip.background = element_rect(fill = "white"), 
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 12)); All3

  Boxplot_list = c(Boxplot_list, All3)
  
  # ################# #
  #Pooled_Correlation ----
  # ################# #
  for (current_group in c("IgA", "IgG")) {
    current_subset_long_Ab = subset(current_subset, current_subset$Ab == 
                                      paste(current_measurement, current_group))
    current_subset_long_Ab = current_subset_long_Ab[,c("SN", "Sample.Type", "Titer")]
    current_subset_long_Ab = reshape(direction = "wide", data = current_subset_long_Ab, 
                                       idvar = c("SN"), timevar = "Sample.Type")
    colnames(current_subset_long_Ab) = gsub("Titer.", "", colnames(current_subset_long_Ab))
    current_subset_long_Ab = pivot_longer(current_subset_long_Ab, 
                                          cols= c("Right", "Left"), 
                                          names_to = "Side", 
                                          values_to = "Individual.Nostril")
    
    current_subset_long_Ab = data.frame(current_subset_long_Ab)
    current_subset_long_Ab$Ab = current_group
    
    if (current_measurement == "Normalized") {
      Vjust_1 = 1; Vjust_2 = 2.2
    } else{
      Vjust_1 = 0.5; Vjust_2 = 1.7
    }
  
    Pooled = ggplot(current_subset_long_Ab, 
           aes(x = log2(Pooled), y = log2(Individual.Nostril), color = Side, shape = Ab)) + 
      geom_point(size = 3) + 
      geom_smooth(method = "lm", formula = "y~x") + 
      stat_cor( data = current_subset_long_Ab[current_subset_long_Ab$Side == "Right",],
                aes(x = log2(Pooled), y = log2(Individual.Nostril), color = Side,
                    label = gsub("R", "r", paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~"))),
                p.accuracy = 0.001, size = 4.5,
                vjust = Vjust_1, method = "pearson", show.legend = F) +
      stat_cor( data = current_subset_long_Ab[current_subset_long_Ab$Side == "Left",], 
                aes(x = log2(Pooled), y = log2(Individual.Nostril), color = Side,
                    label = gsub("R", "r", paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~"))),
                p.accuracy = 0.001, size = 4.5,
                vjust = Vjust_2,
                method = "pearson", show.legend = F) + 
      theme_bw() + 
      scale_shape_manual(breaks = c("IgA", "IgG"),
                         values = c(16, 17), guide = "none") + 
      scale_color_manual(breaks = c("Left", "Right"), 
                         values = c("#CC79A7", "#009E73"), 
                         name = "Nostril Side") + 
      xlab(paste("Pooled", titer))+ 
      ylab(paste("Individual Nostril", titer))
    
    Pooled_list = c(Pooled_list, Pooled)
    }

  # ################# #
  #Left Right IgA ----
  # ################# #
  current_subset = current_subset[,c("SN", "Sample_Ab", "Titer")]
  current_subset = reshape(direction = "wide", current_subset, idvar = c("SN"), timevar = "Sample_Ab")
  
  LR_IgA = ggplot(current_subset, aes ( x = log2(`Titer.Left:IgA`), y = log2(`Titer.Right:IgA`))) + 
    geom_point(shape = 16, size = 3, alpha = 0.7) + 
    geom_smooth(method = "lm", formula = "y~x") +
    ylab(paste("Right Nostril IgA", titer, "[log2]"))+
    xlab(paste("Left Nostril IgA", titer, "[log2]")) + 
    stat_cor( aes(label = gsub("R", "r", paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~"))),
              p.accuracy = 0.001, size = 5,
              method = "pearson") + 
    theme_bw()
  
  # ################# #
  #Left Right IgG ----
  # ################# #  
  LR_IgG = ggplot(current_subset, aes ( x = log2(`Titer.Left:IgG`), y = log2(`Titer.Right:IgG`))) + 
    geom_point(shape = 17, size = 3, alpha = 0.7) + 
    stat_cor( aes(label = gsub("R", "r", paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~"))), 
              p.accuracy = 0.001, size = 5, method = "pearson") + 
    geom_smooth(method = "lm", formula = "y~x") + 
    ylab(paste("Right Nostril IgG", titer, "[log2]"))+
    xlab(paste("Left Nostril IgG", titer, "[log2]")) + 
    theme_bw()

  Left.Right_list = c(Left.Right_list, LR_IgA, LR_IgG) 
  
}

jpeg("Figure2.jpeg", res = 600, height = 5000, width = 9000)
ggarrange(ggarrange(Boxplot_list[[1]], plot.new(), 
                    Boxplot_list[[2]], plot.new(), 
                    Boxplot_list[[3]],
                    ncol = 5, widths = c(1, 0.1, 1, 0.1, 1), common.legend = T),  
          plot.new(), 
          ggarrange(Left.Right_list[[1]], Left.Right_list[[2]], plot.new(),
                    Left.Right_list[[3]], Left.Right_list[[4]], plot.new(),
                    Left.Right_list[[5]], Left.Right_list[[6]],
                    ncol = 8, widths = c(1,1,0.1,1,1,0.1,1,1)),  
          plot.new(), 
          ggarrange(Pooled_list[[1]], Pooled_list[[2]], plot.new(),
                    Pooled_list[[3]], Pooled_list[[4]], plot.new(),
                    Pooled_list[[5]], Pooled_list[[6]],
                    ncol = 8, widths = c(1,1,0.1,1,1,0.1,1,1), 
                    common.legend = T, legend = "bottom"),
          nrow = 5, heights = c(1, 0.2, 1, 0.2, 1))
dev.off()
