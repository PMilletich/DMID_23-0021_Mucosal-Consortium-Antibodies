#---- Script Metadata #----
# Title: Antigen_Heatmap.R
# Author: Trish Milletich, PhD
# Date: 2026-04-17
#--------------------------

library(ggplot2); library(ggpubr)
library(pheatmap); library(ggplotify)
library(tidyr)

#Data Load in and Subset
All_data = read.csv("Saliva_Serum_Nasal.csv")
All_data = subset(All_data, All_data$Days == "Day 1")
All_data = subset(All_data, All_data$Sample.Type %in% c("Serum", "Saliva", "PooledNostril"))
All_data$Ab_ST = paste(All_data$Sample.Type, All_data$Ab)

#Antigen List 
Antigen_list = c("Alpha", "Beta", "Delta", "Gamma", "Omicron.BA.1", "Omicron.BA.5", "RBD", "Spike")

#Default Value 
current_sample = "Saliva IgG"
for (current_sample in unique(All_data$Ab_ST)) {
  current_subset = All_data[All_data$Ab_ST == current_sample,]
  if (grepl("Serum", current_sample) == T) {
    current_subset$Titer = current_subset$Specific.AU.mL
    current_title = paste("Unadjusted", current_sample)
  } else {
    current_subset$Titer = current_subset$Specific.Total.AU.ug
    current_title = paste("Normalized", current_sample)
    current_title = gsub("PooledNostril", "NLF", current_title)
  }
  current_subset = current_subset[,c("SN", "Ag", "Ab", "Titer")]

  #Create Starting Matrices
  R_Matrix = matrix(ncol = 8, nrow = 8)
  colnames(R_Matrix) = Antigen_list
  rownames(R_Matrix) = Antigen_list
  P_Matrix = R_Matrix
  
  for (i in 1:8) {
    for (j in i:8) {
      ij_subset = current_subset[current_subset$Ag %in% c(Antigen_list[i], Antigen_list[j]), ]
      ij_subset = data.frame(ij_subset[,c("SN", "Ag", "Titer")])
      ij_subset = reshape(direction = "wide", data = ij_subset, 
                          idvar = c("SN"), timevar = "Ag") 
      
      if (ncol(ij_subset) == 2) {
        R = NA; P = NA
      } else {
        Pearson = cor.test(log2(ij_subset[,2]), log2(ij_subset[,3]), method = "pearson")
        R = Pearson$estimate
        P = Pearson$p.value
      }
      R_Matrix[Antigen_list[i],Antigen_list[j]] = R
      P_Matrix[Antigen_list[i],Antigen_list[j]] = P
    } #j
  } #i

  # P Matrix shows all are significant 
  R_Matrix_2 = round(R_Matrix, 2)
  R_Matrix_2[is.na(R_Matrix_2)] = ""

  colnames(R_Matrix) = ifelse(colnames(R_Matrix) %in% c("Spike", "RBD"),
                              paste("Wuhan-", colnames(R_Matrix), sep = ""), colnames(R_Matrix))
  rownames(R_Matrix) = ifelse(rownames(R_Matrix) %in% c("Spike", "RBD"),
                              paste("Wuhan-", rownames(R_Matrix), sep = ""), rownames(R_Matrix))
  
  colnames(R_Matrix_2) = ifelse(colnames(R_Matrix_2) %in% c("Spike", "RBD"),
                                paste("Wuhan-", colnames(R_Matrix_2), sep = ""), colnames(R_Matrix_2))
  rownames(R_Matrix_2) = ifelse(rownames(R_Matrix_2) %in% c("Spike", "RBD"),
                                paste("Wuhan-", rownames(R_Matrix_2), sep = ""), rownames(R_Matrix_2))

  current_heatmap = as.ggplot(pheatmap(R_Matrix, cluster_rows = F, cluster_cols = F, 
                                       display_numbers = R_Matrix_2, number_color = "black", fontsize_number = 10,
                                       legend = F,
                                       border_color = "black", na_col = "white",fontsize = 10,
                                       breaks = c(seq(0, 1, length.out = 100)), angle_col = 315))
  
  current_heatmap = current_heatmap + ggtitle(current_title)
  assign(paste(gsub(" ", "_", current_sample), "plot", sep = "_"), current_heatmap)
  print(paste(gsub(" ", "_", current_sample), "plot", sep = "_"))
  dev.off()
}


jpeg("Supplemental_Figure_1.jpeg", res = 600, height = 4000, width = 8000)
ggarrange(
  Serum_IgG_plot, Saliva_IgG_plot, PooledNostril_IgG_plot, 
  Serum_IgA_plot, Saliva_IgA_plot, PooledNostril_IgA_plot,
  ncol = 3, nrow = 2)
dev.off()
