library(ggplot2)
library(reshape2)
setwd("C:/Users/au634287/OneDrive - Aarhus universitet/PhD/Rstudio/Data/")
miRBase <- read.table(file="high_confidence_miRBase21-master.tsv", header = TRUE, fill = TRUE, sep =  " ")
miRBase <- miRBase[which(miRBase$high_conf==T),]

setwd("C:/Users/au634287/OneDrive - Aarhus universitet/PhD/Rstudio/Output/DICER_length_ratio_FC/fold_change_cumulative/")
MIRNA <- read.csv(file = "MIRNA_expr.csv", sep = ",")
unique(MIRNA$SAMPLE)

setwd("C:/Users/au634287/OneDrive - Aarhus universitet/PhD/Rstudio/Output/plots/FIGURES/FC/")
arms <- unique(miRBase[, c(1, 16)]) 

for (sample_type in unique(MIRNA$SAMPLE)) {
  
  # Filter data for the current sample
  WT_MUT <- MIRNA[MIRNA$SAMPLE == sample_type, ]
  
  # Add strand information
  WT_MUT <- merge(WT_MUT, arms, by = "MIRNA")
  
  # Filter out entries where CPM_mut is 0
  WT_MUT <- WT_MUT[which(WT_MUT$CPM_mut > 0), ]
  
  # Calculate log2 Fold Change
  WT_MUT$log2FC <- log2(WT_MUT$CPM_mut / WT_MUT$CPM_wt)
  
  # Separate data by strand type
  WT_MUT_5P <- WT_MUT[which(WT_MUT$STRAND == "5P"), ]
  WT_MUT_3P <- WT_MUT[which(WT_MUT$STRAND == "3P"), ]
  

  # Calculate p-value
  p_val_less <- wilcox.test(WT_MUT_5P$log2FC, WT_MUT_3P$log2FC, alternative = "less")$p.value
  
  # Optionally, format the p-value in scientific notation for display purposes (but keep it numeric)
  p_val_less_formatted <- format(p_val_less, scientific = TRUE)
  

  
  # Create cumulative fraction data
  breaks <- seq(-12, 12, by = 0.01)
  set_5P <- as.numeric(WT_MUT_5P$log2FC)
  set_5P_cumulative <- c(0, cumsum(table(cut(set_5P, breaks, right = FALSE)) / nrow(WT_MUT_5P)))
  set_3P <- as.numeric(WT_MUT_3P$log2FC)
  set_3P_cumulative <- c(0, cumsum(table(cut(set_3P, breaks, right = FALSE)) / nrow(WT_MUT_3P)))
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    breaks = breaks,
    `5P` = set_5P_cumulative[1:length(breaks)],
    `3P` = set_3P_cumulative[1:length(breaks)]
  )
  
  # Convert data to long format for ggplot
  plot_data_long <- melt(plot_data, id.vars = "breaks", variable.name = "Strand", value.name = "CumulativeFraction")
  
  # Rename the levels of the Strand variable to ensure they are correct
  plot_data_long$Strand <- gsub("^X", "", plot_data_long$Strand)
  
  # Plot using ggplot2 with a white background and vertical line at 0
  p <- ggplot(plot_data_long, aes(x = breaks, y = CumulativeFraction, color = Strand)) +
    geom_line(size = 1.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +  # Add vertical line at 0
    scale_color_manual(values = c("black", "red")) +
    labs(title = sample_type, 
         x = "Log2FC CPM", 
         y = "Cumulative fraction",
         color = "Strand") +
       # Use theme_minimal to start with a simple base
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # Set background to white
      plot.background = element_rect(fill = "white", color = NA),   # Set plot background to white
      panel.grid.major = element_blank(),  # Remove major gridlines
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      axis.line = element_line(color = "black"),# Keep axis lines for clarity
      plot.title = element_text(hjust = 0.5)
      )
  
  print(p)

  # Save the plot
  ggsave(filename = paste0("plot_", sample_type, ".png"), plot = p, width = 8, height = 6)
}

setwd("C:/Users/au634287/OneDrive - Aarhus universitet/PhD/Rstudio/Output/plots/Cumulative_FC/")

file_Name <- paste0("Weighted_LEN_",MUT_select,".tsv")
write.table(rbind(breaks, set_5P_cumulative, set_3P_cumulative), file_Name, sep="\t", append = FALSE)
rm(set_3P, set_3P_cumulative, set_5P,set_5P_cumulative, MUT, WT_MUT, WT_MUT_3P, WT_MUT_5P, p_val_less, breaks, )

write.table(x  = WT_MUT_3P, file = "WT_MUT_3P.tsv", row.names = F, col.names = T, sep = ",")
write.table(x  = WT_MUT_5P, file = "WT_MUT_5P.tsv", row.names = F, col.names = T, sep = ",")
write.table(x  = combinations, file = "combinations.csv", row.names = F, col.names = T, sep = ",")


