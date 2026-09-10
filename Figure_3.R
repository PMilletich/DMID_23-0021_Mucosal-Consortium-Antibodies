# ---- Script Metadata ---- 
# Title: Figure_3.R
# Author: Trish Milletich, PhD
# Data: 2026-07-31
# ------------------------------- -  

library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

# Upload and clean data and subset to Spike antigen and mucosal samples  
All_data = read.csv("Saliva_Serum_Nasal.csv")
Spike = subset(All_data, All_data$Ag == "Spike")
Spike = subset(Spike, Spike$Sample.Type %in% c("Saliva", "PooledNostril"))
Spike$Sample.Type = ifelse(Spike$Sample.Type == "PooledNostril", "NLF", Spike$Sample.Type)

############################################################ #
# Append NAs to missing data to all samples to have the same count
ID_Count = data.frame(table(Spike$SN))

Total_data = Spike[Spike$SN == ID_Count[ID_Count$Freq == 12,"Var1"][1],]
Total_data[,c("Specific.AU.mL", "Total.ug.mL", "Specific.Total.AU.ug")] = NA
Total_data$ST_Ab_Day = paste(Total_data$Ab, Total_data$Sample.Type, Total_data$Days)

Missing_data = subset(ID_Count, ID_Count$Freq < 12)

for (current_ID in Missing_data$Var1) {
  current_subset = subset(Spike, Spike$SN == current_ID)
  current_subset$ST_Ab_Day = paste(current_subset$Ab, current_subset$Sample.Type, current_subset$Days)
  Total_data_subset = subset(Total_data, ! Total_data$ST_Ab_Day %in% current_subset$ST_Ab_Day)
  Total_data_subset$SN = current_ID
  Total_data_subset$ST_Ab_Day = NULL
  Spike = rbind(Spike, Total_data_subset)
}

Spike$Ab_Day = paste(Spike$Ab, Spike$Days, sep = "_")

#Clean up dataframes no longer needed
rm(Total_data_subset, Total_data, ID_Count, Missing_data, current_subset)

################################################################ #
# ---- Data Transformation: Ratios and %CV ---- 
################################################################ #
CV_data = data.frame()

#Default Values 
current_ST = "NLF"; current_method = "Normalized"
for (current_ST in c("Saliva", "NLF")) {
  for (current_method in c("Total", "Unadjusted", "Normalized")) {
    #Select data based on Sample type and measurement type 
    Spike_Subset = subset(Spike, Spike$Sample.Type == current_ST)
    
    if (current_method == "Total") {
      Spike_Subset_2 = Spike_Subset[,c("Ab", "Ab_Day", "SN","Sample.Type", "Total.ug.mL")]
      Spike_Subset_2$Measurement = Spike_Subset_2$Total.ug.mL 
      Spike_Subset_2$Total.ug.mL = NULL
    } else if (current_method == "Unadjusted") {
      Spike_Subset_2 = Spike_Subset[,c("Ab", "Ab_Day", "SN","Sample.Type", "Specific.AU.mL")]
      Spike_Subset_2$Measurement = Spike_Subset_2$Specific.AU.mL 
      Spike_Subset_2$Specific.AU.mL = NULL
    } else {
      Spike_Subset_2 = Spike_Subset[,c("Ab", "Ab_Day", "SN","Sample.Type", "Specific.Total.AU.ug")]
      Spike_Subset_2$Measurement = Spike_Subset_2$Specific.Total.AU.ug 
      Spike_Subset_2$Specific.Total.AU.ug = NULL
    }
    
    ############################ #
    # ---- %CV by Participant ----
    ############################ #
    #Create dataframe of CV by participant and Ab type 
    
    CV_df = Spike_Subset_2 %>%
      group_by(SN, Ab) %>% 
      summarise(CV = (sd(Measurement, na.rm = TRUE) / mean(Measurement, na.rm = TRUE)) * 100,
                count = length(Measurement),
                .groups = "drop_last")
    CV_df$ST = current_ST 
    CV_df$Measurement = current_method
    CV_data = rbind(CV_data, CV_df)
    #clean up unused dataframe
    remove(CV_df)
    
    ############################# #
    # ---- Reshape Data and Format ratios ----
    ############################# #
    #Rehape to wide to allow ratios
    Spike_Wide = reshape(Spike_Subset_2[,c("Ab_Day", "SN","Sample.Type", "Measurement")],
                         idvar = c("SN", "Sample.Type"), timevar = "Ab_Day", direction = "wide")
    colnames(Spike_Wide) = gsub("Measurement.", "", colnames(Spike_Wide))
    head(Spike_Wide)
    
    Spike_Wide$IgA_D1.3 = Spike_Wide$`IgA_Day 1`/Spike_Wide$`IgA_Day 3`
    Spike_Wide$IgA_D1.15 = Spike_Wide$`IgA_Day 1`/Spike_Wide$`IgA_Day 15`

    Spike_Wide$IgG_D1.3 = Spike_Wide$`IgG_Day 1`/Spike_Wide$`IgG_Day 3`
    Spike_Wide$IgG_D1.15 = Spike_Wide$`IgG_Day 1`/Spike_Wide$`IgG_Day 15`

    Spike_Wide = Spike_Wide[,c("SN", "IgG_D1.3", "IgG_D1.15",
                               "IgA_D1.3", "IgA_D1.15")]
    
    #Reshape to long for graphing
    Spike_Long <- Spike_Wide %>% 
      pivot_longer(cols = `IgG_D1.3`:`IgA_D1.15`, 
                   names_to = "Ab",
                   values_to = "value" )
    Spike_Long$Ab = gsub("_", "\n", Spike_Long$Ab)
    Spike_Long$Ab = gsub("D1.", "Day1/Day", Spike_Long$Ab)

    Spike_Long$Main = ifelse(grepl("IgA", Spike_Long$Ab), "IgA", "IgG")
    
    Spike_Long$Ab = factor(Spike_Long$Ab,
                           levels = c("IgA\nDay1/Day3", "IgG\nDay1/Day3", 
                                      "IgA\nDay1/Day15", "IgG\nDay1/Day15" ))
    
    Spike_Long$ST_Main = paste(current_ST, Spike_Long$Main)
    Spike_Long$ST_Main = factor(Spike_Long$ST_Main, 
                                levels = c("NLF IgA", "NLF IgG",
                                           "Saliva IgA", "Saliva IgG"))
    
    Spike_Long$AB_ST = paste(current_ST, Spike_Long$Main, sep = ": ")
    Spike_Long = Spike_Long[is.na(Spike_Long$value) == F,]
    
    #################################### #
    # T-test of Day differences
    t_test_df = data.frame() 
    
    for (current_pair in unique(Spike_Long$Ab)) {
      p = t.test(log2(Spike_Long[Spike_Long$Ab == current_pair, "value"]), mu = 0)
      p_value = p$p.value
      t_test_df = rbind(t_test_df, 
                        data.frame( "Ab" = current_pair, 
                                    "P" = round(p_value, 3)))
    }
    t_test_df$Main = ifelse(grepl("IgA", t_test_df$Ab), "IgA", "IgG")
    t_test_df$AB_ST = paste(current_ST, t_test_df$Main, sep = ": ")
    
    #Plot Ratios with labelled T-test p-values
    t_plot = ggplot(Spike_Long, 
                    aes(x = Ab, y = log2(value), fill = ST_Main, shape = Main)) + 
      geom_hline(yintercept = 0, linetype = "dashed") + 
      facet_grid(~AB_ST, scales = "free_x") + 
      geom_boxplot(outlier.shape = "", color = "black", show.legend = F) +
      geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
      ylab(paste(current_method, "log2(FC)")) + 
      theme_bw() + 
      geom_label(data = t_test_df,
                 aes(x = Ab, y = Inf, label = paste("P=", round(P,2), sep = "")),
                 fill = "white", vjust = 2) +
      scale_fill_manual(values = c("NLF IgA"="#CC6677", "NLF IgG" = "#7B0E45",
                                   "Saliva IgA" = "#88CCEE", "Saliva IgG" = "#065882"),
                        limits = c("NLF IgA", "NLF IgG", "Saliva IgA", "Saliva IgG"),
                        drop = F) +
      scale_shape_manual(breaks = c("IgA", "IgG"), 
                         values = c(16, 17), 
                         guide = "none") + 
      theme(legend.position = "none",
            strip.background = element_rect(fill = "white"), 
            legend.margin=margin(0,0,0,0),
            axis.title.x = element_blank(), 
            strip.text = element_text(size = 13,
                                      margin = margin(1,1,1,1, "mm")),
            legend.box.margin=margin(0,0,0,0)) + 
      guides(fill = guide_legend(override.aes = list(shape = 22, size = 5)));t_plot
    
    assign(paste(current_ST, current_method, "t.Plot", sep = "_"), t_plot)
    print(paste(current_ST, current_method, "t.Plot", sep = "_"))
  }
}

###################################### #
# --- CV Plot ----
###################################### #
#Clean CV data 
CV_data$ST_ab = paste(CV_data$ST, CV_data$Ab)
CV_data$ST_ab = gsub(" Ig", "\nIg", CV_data$ST_ab)
CV_data$ST_ab = factor(CV_data$ST_ab, 
                       levels = c("Saliva\nIgA", "NLF\nIgA", "Saliva\nIgG", "NLF\nIgG"))
CV_data$Measurement = factor(CV_data$Measurement, 
                             levels = c("Total", "Unadjusted", "Normalized"))

#Create dataframe of means for labels 
CV_data_Mean = CV_data %>% 
  group_by(Measurement, ST_ab) %>% 
  summarise(Mean = mean(CV, na.rm = T))

CV_plot = ggplot(CV_data, aes(x = ST_ab, y = CV, fill = ST_ab)) + 
  facet_wrap(~Measurement, scales = "free_x") +
  geom_hline(yintercept = 30, color = "grey") + 
  geom_boxplot(color = "black") + theme_bw() + 
  scale_fill_manual(values = c("NLF\nIgA"="#CC6677", "NLF\nIgG" = "#7B0E45",
                               "Saliva\nIgA" = "#88CCEE", "Saliva\nIgG" = "#065882")) + 
  geom_label(data = CV_data_Mean, aes(x = ST_ab, y = 120, label = round(Mean,2)), fill = "white") + 
  ylab("%CV across time") + 
  theme(axis.title.x = element_blank(), 
        legend.position = "top", 
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 13,
                                  margin = margin(1,1,1,1, "mm")),
        legend.margin=margin(0,0,0,0)); CV_plot


################################################## #
# ---- Mean Line Graphs ----
################################################## #
Spike$Sample.Day = paste(Spike$Sample.Type, Spike$Days, sep = ": ")
#Clean Data 
Spike$Sample.Day = gsub("Pooled", "Pooled ", Spike$Sample.Day)
Spike$Sample.Day = gsub("Left", "Left Nostril", Spike$Sample.Day)
Spike$Sample.Day = gsub("Right", "Right Nostril", Spike$Sample.Day)
Spike$Sample.Day = factor(Spike$Sample.Day, 
                          levels = c("Saliva: Day 1","Saliva: Day 3","Saliva: Day 15",
                                     "Left Nostril: Day 1", "Right Nostril: Day 1",
                                     "Pooled Nostril: Day 1", "Pooled Nostril: Day 3" , "Pooled Nostril: Day 15",
                                     "Serum: Day 1"))
Spike$Ab_SD = paste(Spike$Sample.Day, Spike$Ab, sep = "\n")
Spike$ST_Ab = paste(Spike$Sample.Type, Spike$Ab, sep = ": ")
Spike = Spike[,c("SN", "Sample.Type", 'Days', "Ab", 
                 "Specific.AU.mL", "Total.ug.mL", "Specific.Total.AU.ug")]
Spike$Days = factor(Spike$Days, levels = c("Day 1", 'Day 3', "Day 15"))
Spike$Group = paste(Spike$Days, Spike$Ab, sep = "_")

# Format to Long to allow T-test comparison
Spike_Long = Spike %>% 
  pivot_longer(
    cols = c("Specific.AU.mL","Total.ug.mL", "Specific.Total.AU.ug"), 
    names_to = "Measurement", 
    values_to = "Values"
  )
Spike_Long = data.frame(Spike_Long)
Spike_Long$Sample.Type_Ab= paste(Spike_Long$Sample.Type, Spike_Long$Ab, sep = ": ")
Spike_Long$Sample.Type_Ab = gsub("Nost", " Nost", Spike_Long$Sample.Type_Ab)
Spike_Long$Sample.Type_Ab_Day = paste(Spike_Long$Sample.Type, Spike_Long$Ab, Spike_Long$Days, sep = ": ")
Spike_Long$log2_Values = log2(Spike_Long$Values)

current_measure = "Normalized"
for (current_measure in c("Total", "Unadjusted", "Normalized")) {
  if (current_measure == "Total") {
    ylab_current = "Total [log2]"
    Measure_2 = "Total.ug.mL"
    ymax = 200
  } else  if (current_measure == "Unadjusted") {
    ylab_current = "Unadjusted [log2]"
    Measure_2 = "Specific.AU.mL"
    ymax = 9000
  } else  if (current_measure == "Normalized") {
    ylab_current = "Normalized [log2]"
    Measure_2 = "Specific.Total.AU.ug"
    ymax = 120
  } else {
    print ("uh oh")
  }
  
  line_graph =   ggplot(Spike_Long[Spike_Long$Measurement == Measure_2, ], 
                        aes(x = Days,
                            y = log2_Values,
                            color = Sample.Type_Ab,
                            shape = Sample.Type_Ab,
                            linetype = Sample.Type,
                            group = Sample.Type_Ab)) +
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar", 
                 linetype = "solid", 
                 width = 0.1, position = position_dodge(width = 0.07)) + 
    stat_summary(fun = mean, geom = "line", 
                 linewidth = 0.8, 
                 position = position_dodge(width = 0.07)) +
    stat_summary(fun = mean, geom = "point",
                 size = 4, 
                 position = position_dodge(width = 0.07)) +
    theme_bw() + 
    scale_color_manual(values = c("NLF: IgA"="#CC6677", "NLF: IgG" = "#7B0E45",
                                  "Saliva: IgA" = "#88CCEE", "Saliva: IgG" = "#065882"),
                       name = "Mucosal Ab")+ 
    scale_shape_manual(values = c("NLF: IgA"=16, 
                                  "NLF: IgG" = 17,
                                  "Saliva: IgA" = 16, 
                                  "Saliva: IgG" = 17),
                       name = "Mucosal Ab") + 
    scale_linetype_manual(values = c("NLF"="solid", 
                                     "Saliva" = "dashed"),
                          name = "Sample Type") + 
    ylab(ylab_current) + 
    theme(strip.background = element_rect(fill = "white"), 
          legend.position = "bottom", 
          axis.title.x = element_blank()); line_graph
  
  assign(paste(current_measure, "_line", sep = ""), line_graph)
  print(paste(current_measure, "_line", sep = ""))
}


jpeg("Figure3.jpeg", res = 750, height = 6500, width = 8000)
ggarrange(ggarrange(Total_line, Unadjusted_line, Normalized_line, 
          plot.new(), plot.new(), plot.new(), 
          Saliva_Total_t.Plot, Saliva_Unadjusted_t.Plot, Saliva_Normalized_t.Plot,
          plot.new(), plot.new(), plot.new(), 
          NLF_Total_t.Plot, NLF_Unadjusted_t.Plot, NLF_Normalized_t.Plot,
          ncol = 3, nrow = 5, 
          heights = c(1,0.1,1.5,0.1,1.5),
          common.legend = T, legend = "top"), 
          plot.new(), 
          CV_plot, nrow = 3, heights = c(3,0.1,1), legend = "none")
dev.off()

