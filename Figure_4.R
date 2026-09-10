# ---- Script Metadata ---- 
# Title: Figure_4.R
# Author: Trish Milletich, PhD
# Data: 2026-07-31
# --------------------------------

library(ggplot2); library(ggpubr)

#Load and Clean Data 
All_data = read.csv("Saliva_Serum_Nasal.csv")
All_data = subset(All_data, All_data$Ag == "Spike")
All_data = subset(All_data, All_data$Days == "Day 1")
All_data = subset(All_data, ! All_data$Sample.Type %in% c("Left", "Right"))
All_data$Sample.Type = gsub("PooledNostril", "NLF", All_data$Sample.Type)
All_data = All_data[order(All_data$Ab),]
All_data$Sample_Ab = paste(All_data$Sample.Type,  All_data$Ab)

comparison_list = list(c("IgA", "Normalized"), c("IgG", "Normalized"), 
                       c("IgA", "Unadjusted"), c("IgG", "Unadjusted"))

for (current_compare in comparison_list) {
  data_subset = subset(All_data, All_data$Ab == current_compare[1])
  
  if (current_compare[2] == "Normalized") {
    data_subset$Titer = data_subset$Specific.Total.AU.ug
    value_list = c(16, 17) #values for geom_point shape
  } else if (current_compare[2] == "Unadjusted") {
    data_subset$Titer = data_subset$Specific.AU.mL
    value_list = c(1, 2) #values for geom_point shape
  } else {
    break
  }
  
  data_subset = data_subset[,c("SN", "Titer", "Sample.Type")]
  
  #Reshape to wide to allow comparisons 
  data_subset = reshape(data_subset, direction = "wide", 
                        idvar = "SN", timevar = "Sample.Type")
  colnames(data_subset) = gsub("Titer.", "", colnames(data_subset))
  
  data_subset$Ab = current_compare[1]

  Saliva.NLF = ggplot(data_subset, aes(x = log2(Saliva), y = log2(NLF), shape = Ab)) + 
    geom_point(size = 2) +
    geom_smooth( method = "lm", formula = "y~x") +
    stat_cor(aes(label = gsub("R", "r", paste(..r.label.., ..p.label.., sep = "~`,`~"))),
              p.accuracy = 0.001, size = 4.5,
             method = "pearson", show.legend = F) + 
    theme_bw() + 
    scale_shape_manual(breaks = c("IgA", "IgG"),
                       values = value_list, guide = "none") + 
    xlab(paste(current_compare[2], "Saliva", current_compare[1], "[log2]")) +
    ylab(paste(current_compare[2], "NLF", current_compare[1], "[log2]")) 
    
  Serum.NLF = ggplot(data_subset, aes(x = log2(NLF), y = log2(Serum), shape = Ab)) + 
    geom_point(size = 2) +
    geom_smooth( method = "lm", formula = "y~x") +
    stat_cor(aes(label = gsub("R", "r", paste(..r.label.., ..p.label.., sep = "~`,`~"))),
             p.accuracy = 0.001, size = 4.5,
             method = "pearson", 
             vjust = -0.75,
             show.legend = F) + 
    theme_bw() + 
    scale_shape_manual(breaks = c("IgA", "IgG"),
                       values = value_list, guide = "none") + 
    xlab(paste(current_compare[2], "NLF", current_compare[1], "[log2]"))+ 
    ylab(paste("Serum", current_compare[1], "[log2]")); Serum.NLF
  
  Serum.Saliva = ggplot(data_subset, aes(x = log2(Saliva), y = log2(Serum), shape = Ab)) + 
    geom_point(size = 2) +
    geom_smooth( method = "lm", formula = "y~x") +
    stat_cor(aes(label = gsub("R", "r", paste(..r.label.., ..p.label.., sep = "~`,`~"))),
             p.accuracy = 0.001, size = 4.5,
             method = "pearson", show.legend = F) + 
    theme_bw() + 
    scale_shape_manual(breaks = c("IgA", "IgG"),
                       values = value_list, guide = "none") + 
    xlab(paste(current_compare[2], "Saliva", current_compare[1], "[log2]"))+ 
    ylab(paste("Serum", current_compare[1], "[log2]"))
  
  all_plots = ggarrange(Serum.NLF, Serum.Saliva, Saliva.NLF, ncol = 3)
  assign(paste(current_compare[1], current_compare[2], sep = "_"), all_plots)
}

jpeg("/Figure4.jpeg", res = 600, height = 3500, width = 9500)
ggarrange(IgG_Normalized, plot.new(), IgA_Normalized, 
          plot.new(), plot.new(), plot.new(),
          IgG_Unadjusted, plot.new(), IgA_Unadjusted, 
          nrow = 3, ncol = 3, 
          widths = c(1, 0.03, 1), 
          heights = c(1, 0.2, 1))
dev.off()


