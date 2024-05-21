##This script is used to generate separate plots for each chromosome and each type of structural variation (SV) in the data.
##By running this script, you can generate a set of plots that provide insights into the distribution and characteristics of structural variations across different chromosomes and types within a barley dataset.

#here there are four methods to settle down the order

library(ggplot2)
custom_colors <- c("DEL" = "#66c2a5", "DUP" = "#fc8d62", "INS" = "#8da0cb", "INV" = "#e78ac3")

data<-read.table("barley.all.SV.plot.txt",header = T)
#Accession Chromosome Position Type
#Baudin         1H   180730  INS 
#Baudin         1H   182233  DEL
#Baudin         1H   185665  INS
#Baudin         1H   188906  INS
#Baudin         1H   241466  DEL
#Baudin         1H   241963  DEL
data$Type <- factor(data$Type, levels = c("INS", "DEL", "DUP", "INV"))
for (chromosome in unique(data$Chromosome)) {
    # Create a directory to save the plots
    dir.create(paste0("plots_", chromosome), showWarnings = FALSE)
    
    # Subset data for the current chromosome
    chromosome_data <- subset(data, Chromosome == chromosome)
    
    # Loop through each type of variation
    for (type in levels(data$Type)) {
        # Subset data for the current type of variation
        type_data <- subset(chromosome_data, Type == type)
        
        # Create a plot for the current type of variation
        p <- ggplot(type_data, aes(x = Position, y = as.factor(Accession), color = Type)) +
            geom_point(shape = '|', size = 5, alpha = 0.7) +
            labs(title = paste("Chromosomal Variants -", chromosome, "-", type),
                 x = "Position", y = "Accession") +
            theme_minimal() +
            theme(
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 16, hjust = 0.5),
                text = element_text(family = "Times New Roman", face = "bold", size = 20),
                legend.text = element_text(size = 16),
                panel.background = element_rect(fill = "white") # Set panel background to white
            ) +
            scale_color_manual(values = custom_colors) +
            scale_x_continuous(breaks = c(0, 100000000, 200000000, 300000000, 400000000, 500000000, 600000000, 700000000, 800000000), labels = c("0Mbp", "", "200Mbp", "", "400Mbp", "", "600Mbp", "", "800Mbp"))
        
        # Save the plot in TIFF format with white background
        ggsave(paste0("plots_", chromosome, "/plot_", type, ".tiff"), plot = p, width = 10, height = 6, units = "in", bg = "white")
    }
}

##################################################################################################################################################
##the following is based on the order of genetic tree 

# Load required library
library(ggplot2)

# Define the order of accessions
accession_order <- c("Schooner", "Yambla", "Buloke", "Clipper", "Beecher", "Eth69", "Sahara3771", "Prior", "Fathom", "Compass", "Leabrook", "Yeti", "Laperouse", "Stirling", "Vlamingh", "Baudin", "LaTrobe", "Maximus", "Buff", "Mundah", "Betzes", "Harrington", "RGT_Planet", "Gairdner", "Halcyon", "ILAN15", "Galleon", "Flagship")

# Read the data
data <- read.table("barley.all.SV.plot.txt", header = TRUE)

# Order Accession column according to defined order
data$Accession <- factor(data$Accession, levels = accession_order)

# Set the order of types
data$Type <- factor(data$Type, levels = c("INS", "DEL", "DUP", "INV"))

# Loop through each chromosome
for (chromosome in unique(data$Chromosome)) {
    # Create a directory to save the plots
    dir.create(paste0("plots_", chromosome), showWarnings = FALSE)
    
    # Subset data for the current chromosome
    chromosome_data <- subset(data, Chromosome == chromosome)
    
    # Loop through each type of variation
    for (type in levels(data$Type)) {
        # Subset data for the current type of variation
        type_data <- subset(chromosome_data, Type == type)
        
        # Create a plot for the current type of variation
        p <- ggplot(type_data, aes(x = Position, y = Accession, color = Type)) +
            geom_point(shape = '|', size = 5, alpha = 0.7) +
            labs(title = paste("Chromosomal Variants -", chromosome, "-", type),
                 x = "Position", y = "Accession") +
            theme_minimal() +
            theme(
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 16, hjust = 0.5),
                text = element_text(family = "Times New Roman", face = "bold", size = 20),
                legend.text = element_text(size = 16),
                panel.background = element_rect(fill = "white") # Set panel background to white
            ) +
            scale_color_manual(values = custom_colors) +
            scale_x_continuous(breaks = c(0, 100000000, 200000000, 300000000, 400000000, 500000000, 600000000, 700000000, 800000000), labels = c("0Mbp", "", "200Mbp", "", "400Mbp", "", "600Mbp", "", "800Mbp"))
        
        # Save the plot in TIFF format with white background
        ggsave(paste0("plots_", chromosome, "/plot_", type, ".tiff"), plot = p, width = 10, height = 6, units = "in", bg = "white")
    }
}

##################################################################################################################################################
##split chromosome cluster results
# Load required libraries
library(ggplot2)
library(dendextend)  # For hierarchical clustering and dendrogram manipulation

# Read the data
data <- read.table("barley.all.SV.plot.txt", header = TRUE)

# Set the order of accessions
accession_order <- c("Schooner", "Yambla", "Buloke", "Clipper", "Beecher", "Eth69", "Sahara3771", "Prior", "Fathom", "Compass", "Leabrook", "Yeti", "Laperouse", "Stirling", "Vlamingh", "Baudin", "LaTrobe", "Maximus", "Buff", "Mundah", "Betzes", "Harrington", "RGT_Planet", "Gairdner", "Halcyon", "ILAN15", "Galleon", "Flagship")

# Loop through each chromosome
for (chromosome in unique(data$Chromosome)) {
    # Subset data for the current chromosome
    chromosome_data <- subset(data, Chromosome == chromosome)
    
    # Perform hierarchical clustering for the current chromosome
    accession_matrix <- table(chromosome_data$Accession)
    dend <- as.dendrogram(hclust(dist(accession_matrix)))
    ordered_accessions <- labels(dend)
    
    # Reorder Accession column based on clustering
    chromosome_data$Accession <- factor(chromosome_data$Accession, levels = ordered_accessions)
    
    # Set the order of types
    chromosome_data$Type <- factor(chromosome_data$Type, levels = c("INS", "DEL", "DUP", "INV"))
    
    # Create a directory to save the plots
    dir.create(paste0("plots_", chromosome), showWarnings = FALSE)
    
    # Loop through each type of variation
    for (type in levels(chromosome_data$Type)) {
        # Subset data for the current type of variation
        type_data <- subset(chromosome_data, Type == type)
        
        # Create a plot for the current type of variation
        p <- ggplot(type_data, aes(x = Position, y = Accession, color = Type)) +
            geom_point(shape = '|', size = 5, alpha = 0.7) +
            labs(title = paste("Chromosomal Variants -", chromosome, "-", type),
                 x = "Position", y = "Accession") +
            theme_minimal() +
            theme(
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 16, hjust = 0.5),
                text = element_text(family = "Times New Roman", face = "bold", size = 20),
                legend.text = element_text(size = 16),
                panel.background = element_rect(fill = "white") # Set panel background to white
            ) +
            scale_color_manual(values = custom_colors) +
            scale_x_continuous(breaks = c(0, 100000000, 200000000, 300000000, 400000000, 500000000, 600000000, 700000000, 800000000), labels = c("0Mbp", "", "200Mbp", "", "400Mbp", "", "600Mbp", "", "800Mbp"))
        
        # Save the plot in TIFF format with white background
        ggsave(paste0("plots_", chromosome, "/plot_", type, ".tiff"), plot = p, width = 10, height = 6, units = "in", bg = "white")
    }
}

##################################################################################################################################################
# Load required libraries
library(ggplot2)
library(dendextend)  # For hierarchical clustering and dendrogram manipulation

# Read the data
data <- read.table("barley.all.SV.plot.txt", header = TRUE)

# Set the order of accessions
accession_order <- c("Schooner", "Yambla", "Buloke", "Clipper", "Beecher", "Eth69", "Sahara3771", "Prior", "Fathom", "Compass", "Leabrook", "Yeti", "Laperouse", "Stirling", "Vlamingh", "Baudin", "LaTrobe", "Maximus", "Buff", "Mundah", "Betzes", "Harrington", "RGT_Planet", "Gairdner", "Halcyon", "ILAN15", "Galleon", "Flagship")

# Perform hierarchical clustering on the entire dataset
accession_matrix <- table(data$Accession)
dend <- as.dendrogram(hclust(dist(accession_matrix)))
ordered_accessions <- labels(dend)

# Reorder Accession column based on clustering across all chromosomes
data$Accession <- factor(data$Accession, levels = ordered_accessions)

# Loop through each chromosome
for (chromosome in unique(data$Chromosome)) {
    # Subset data for the current chromosome
    chromosome_data <- subset(data, Chromosome == chromosome)
    
    # Set the order of types
    chromosome_data$Type <- factor(chromosome_data$Type, levels = c("INS", "DEL", "DUP", "INV"))
    
    # Create a directory to save the plots
    dir.create(paste0("plots_", chromosome), showWarnings = FALSE)
    
    # Loop through each type of variation
    for (type in levels(chromosome_data$Type)) {
        # Subset data for the current type of variation
        type_data <- subset(chromosome_data, Type == type)
        
        # Create a plot for the current type of variation
        p <- ggplot(type_data, aes(x = Position, y = Accession, color = Type)) +
            geom_point(shape = '|', size = 5, alpha = 0.7) +
            labs(title = paste("Chromosomal Variants -", chromosome, "-", type),
                 x = "Position", y = "Accession") +
            theme_minimal() +
            theme(
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 16, hjust = 0.5),
                text = element_text(family = "Times New Roman", face = "bold", size = 20),
                legend.text = element_text(size = 16),
                panel.background = element_rect(fill = "white") # Set panel background to white
            ) +
            scale_color_manual(values = custom_colors) +
            scale_x_continuous(breaks = c(0, 100000000, 200000000, 300000000, 400000000, 500000000, 600000000, 700000000, 800000000), labels = c("0Mbp", "", "200Mbp", "", "400Mbp", "", "600Mbp", "", "800Mbp"))
        
        # Save the plot in TIFF format with white background
        ggsave(paste0("plots_", chromosome, "/plot_", type, ".tiff"), plot = p, width = 10, height = 6, units = "in", bg = "white")
    }
}




