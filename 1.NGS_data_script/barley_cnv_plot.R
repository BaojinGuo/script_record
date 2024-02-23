#!/usr/bin/Rscript
##Rscript barley_cnv_plot.R --cnv_data $file1 --chromosome_lengths $file2

# Load necessary libraries
library(ggplot2)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check the number of arguments
if (length(args) != 2) {
  cat("Usage: ./plot_cnv_script.R <CNV data file> <chromosome lengths file>\n")
  quit(status = 1)
}

# Read chromosome lengths data
chromosome_lengths <- read.table(args[2], header = FALSE, col.names = c('Chromosome', 'Length'))

# Read CNV data
cnv_data <- read.table(args[1], header = FALSE, col.names = c('Chromosome', 'Start', 'End', 'Type'))

# Specify the desired order of CNV types
desired_order <- c("SGR", "PHR", "Deletion_both", "Duplication_both", "BASSdeletion", "BASSduplication", "Flindersdeletion", "Flindersduplication")

# Convert Type to factor
cnv_data$Type <- as.factor(cnv_data$Type)

# Define color palette
color_palette <- c("#22BBEE","#111148", "#007046","#AA0000","#60D6A9","#FF8B73","#9FEE00","#FFD9E8")

# Merge chromosome lengths and CNV data
merged_data <- merge(cnv_data, chromosome_lengths, by='Chromosome')
merged_data$Type <- factor(merged_data$Type, levels = desired_order)
# Plot the chromosome graph with Type information and custom colors
p <- ggplot(merged_data, aes(x = Start, y = Chromosome, fill = Type)) +
  geom_tile(height=0.8,alpha=1) +
  scale_fill_manual(values = color_palette) +
  labs(x = "Position", y = "Chromosome", title = "Distribution of CNV Data on Chromosomes", fill = "CNV Type") +
  theme_minimal()+
  scale_x_continuous(breaks = c(0, 100000000, 200000000, 300000000, 400000000, 500000000, 600000000, 700000000, 800000000), labels = c("0Mbp", "", "200Mbp", "", "400Mbp", "", "600Mbp", "", "800Mbp"))+
  theme(text = element_text(family = "Times New Roman", face = "bold", size = 20))+
  theme(axis.text = element_text(size = 16,face = "bold"))+
  theme(legend.text = element_text(size = 16))
# Save the plot as a PDF file
ggsave("output_plot.tiff", plot = p, width = 20, height = 6)

