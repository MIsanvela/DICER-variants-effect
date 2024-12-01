library(ggplot2)
library(dplyr)

setwd("C:/Users/au634287/OneDrive - Aarhus universitet/PhD/Rstudio/Data")
miRBase <- read.table(file="high_confidence_miRBase21-master.tsv", header = TRUE, fill = TRUE, sep =  " ")
#To delete noise and take only sequences that have been checked as a microRNAs there is a colum of high fidelity. There has to be true 
miRBase <- miRBase[which(miRBase$high_conf==T),]
arms <- unique(miRBase[, c(1, 16)]) 
#setwd("C:/Users/au634287/OneDrive - Aarhus universitet/PhD/Rstudio/Data")
#reads <- read.csv(file = "READS_to_WT_47MUT.tsv", sep = "\t")

#WTs <- reads[which(reads$SAMPLE=="WT1"|reads$SAMPLE=="WT2"|reads$SAMPLE=="WT3"), ]
#WTs <- aggregate(WTs$CPM, by=list(WTs$SAMPLE,WTs$MIRNA),sum)
#WTs$count <- 1
#colnames(WTs)
#colnames(WTs)<- c("WT","MIRNA", "CPM", "count")
#WT_miR <- aggregate(WTs$count, by=list(WTs$MIRNA),sum)
#WT_miR <- WT_miR[which(WT_miR$x>=2),]
#WTs <- WTs[WTs$MIRNA %in% WT_miR$Group.1,]
#WTs <- aggregate(WTs$CPM, by=list(WTs$MIRNA),mean)
#colnames(WTs)
#colnames(WTs) <- c("MIRNA","CPM")

#templated <- unique(reads[,c(2,3)])

#samples <- unique(reads$SAMPLE)
#miRNAs <- unique(WTs$Group.1)

#combinations <- expand.grid(samples,miRNAs)
#combinations$CPM <- NA
#combinations$wLEN <- NA

#i <- 1
#for (i in 1:nrow(combinations)) {
  test_set <-reads[which(reads$SAMPLE==combinations$Var1[i] & reads$MIRNA==combinations$Var2[i]), ]
  combinations$CPM[i] <- sum(test_set$CPM)
  test_set$ratio2 <- test_set$RATIO/100
  test_set$wLEN <- test_set$ratio2*test_set$LEN_READ
  wLEN_test <- sum(test_set$wLEN)
  combinations$wLEN[i] <- wLEN_test
  rm(wLEN_test,test_set)
  print(i/nrow(combinations)*100)
#}

setwd("C:/Users/au634287/OneDrive - Aarhus universitet/PhD/Rstudio/Output/")
combinations <- read.csv(file = "combinations.csv", sep = ",")
colnames(combinations)
colnames(combinations) <- c("counts","SAMPLE","MIRNA","CPM","wLEN")
WTs <- combinations[which(combinations$SAMPLE=="WT1"|combinations$SAMPLE=="WT2"|combinations$SAMPLE=="WT3"), ]
WTs <- WTs[which(WTs$CPM>0),]
WT_wL <- aggregate(WTs$wLEN, by=list(WTs$MIRNA),mean)
colnames(WT_wL) <- c("MIRNA","wLEN_WT")
WT_CPM <- aggregate(WTs$CPM, by=list(WTs$MIRNA),mean)
colnames(WT_CPM) <- c("MIRNA","CPM_WT")

WT_wL_SD <- aggregate(WTs$wLEN, by=list(WTs$MIRNA),sd)
colnames(WT_wL_SD) <- c("MIRNA","wL_SD")

wLENs <- merge(WT_CPM,WT_wL,by="MIRNA")
wLENs <- merge(wLENs,WT_wL_SD,by="MIRNA")
rm(WT_wL,WT_CPM,WT_wL_SD)

# Define the mutations in a single vector


# Predefine breaks for the cumulative plot
breaks <- seq(-10, 10, by = 0.01)
setwd("C:/Users/au634287/OneDrive - Aarhus universitet/PhD/Rstudio/Output/plots/FIGURES/LENGTH/")
# Loop through each mutation
for (MUT_select in MUT_selects) {
  
  # Extract the data for the current mutation
  MUT <- combinations[which(combinations$SAMPLE == MUT_select), ]
  
  # Rename columns
  colnames(MUT) <- c("counts", "MUT", "MIRNA", "CPM_MUT", "wL_MUT")
  
  # Merge with WT data
  WT_MUT <- merge(wLENs, MUT, by = "MIRNA")
  
  # Calculate delta
  WT_MUT$Delta <- WT_MUT$wL_MUT - WT_MUT$wLEN_WT
  
  # Merge with arms information
  WT_MUT <- merge(WT_MUT, arms, by = "MIRNA")
  
  # Filter for CPM_MUT > 0
  WT_MUT <- WT_MUT[which(WT_MUT$CPM_MUT > 0), ]
  
  ### Split into 5P and 3P
  WT_MUT_5P <- WT_MUT[which(WT_MUT$STRAND == "5P"), ]
  WT_MUT_3P <- WT_MUT[which(WT_MUT$STRAND == "3P"), ]
  
  ### Statistical Test
  p_val_less <- wilcox.test(WT_MUT_5P$Delta, WT_MUT_3P$Delta, alternative = "greater")$p.value
  
  ### Prepare data for ggplot
  set_5P <- as.numeric(WT_MUT_5P$Delta)
  set_5P_cumulative <- c(0, cumsum(table(cut(set_5P, breaks, right = FALSE)) / nrow(WT_MUT_5P)))
  
  set_3P <- as.numeric(WT_MUT_3P$Delta)
  set_3P_cumulative <- c(0, cumsum(table(cut(set_3P, breaks, right = FALSE)) / nrow(WT_MUT_3P)))
  
  # Create a long-format data frame for ggplot
  plot_data <- data.frame(
    breaks = rep(breaks, 2),
    CumulativeFraction = c(set_5P_cumulative, set_3P_cumulative),
    Strand = rep(c("5P", "3P"), each = length(breaks))
  )
  
  ### ggplot
  p <- ggplot(plot_data, aes(x = breaks, y = CumulativeFraction, color = Strand)) +
    geom_line(size = 1.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +  # Add vertical line at 0
    scale_color_manual(values = c("black", "red")) +
    labs(title = MUT_select, 
         x = "Delta weighted length (MUT-WT)", 
         y = "Cumulative fraction",
         color = "Strand") +
    theme_minimal(base_size = 15) +  # Use theme_minimal to start with a simple base
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # Set background to white
      plot.background = element_rect(fill = "white", color = NA),   # Set plot background to white
      panel.grid.major = element_blank(),  # Remove major gridlines
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      axis.line = element_line(color = "black"),  # Keep axis lines for clarity
      plot.title = element_text(hjust = 0.5)
      )
  
  print(p)  # Print the plot for each mutation
  
  ### Optionally save the plot to a file (e.g., PNG)
   ggsave(filename = paste0("plot_", MUT_select, ".png"), plot = p, width = 8, height = 6)
}



