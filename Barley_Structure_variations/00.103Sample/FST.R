#########################lineplot
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(scales)

file_name <- "Cultiar.AUS_OTH.SV.Fst.windowed.weir.fst"

fst_data <- read_tsv(file_name, col_names = TRUE)

# Filter data to include only chromosomes 1H to 7H
fst_data_filtered <- fst_data %>%
    filter(CHROM %in% c("1H", "2H", "3H", "4H", "5H", "6H", "7H"))

p <- ggplot(fst_data_filtered, aes(x = BIN_START)) +
    #geom_line(aes(y = WEIGHTED_FST, color = "Weighted FST"), size = 1, alpha = 0.7) +
    geom_line(aes(y = MEAN_FST, color = "Mean FST"), size = 0.5, alpha = 0.7) +
    facet_wrap(~ CHROM, scales = "free_x", ncol = 1) +  # Arrange in one column
    #scale_color_manual(values = c("Weighted FST" = "blue", "Mean FST" = "red")) +
    scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = "Mb")) +  # Convert to Mb
    theme_minimal() +
    labs(title = "",
         x = "",
         y = "FST",
         ) +
    theme(legend.position = "none", 
          strip.text = element_text(size = 12),
          axis.text.x = element_text( hjust = 1)) + # Increase spacing between facets
    ylim(-0.1,0.7)

ggsave("FST_AUS_OTH_lineplot.pdf", plot = p, width = 12, height = 18, dpi = 300, bg = "white")

###############################
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)

# Define file name
file_name <- "Cultiar.AUS_OTH.SV.Fst.windowed.weir.fst"  

# Read in the FST data
fst_data <- read_tsv(file_name, col_names = TRUE)

fst_data_filtered <- fst_data %>%
    select(MEAN_FST) %>%
    mutate(Chromosome = "All Chromosomes")# Add a column to indicate the group

p <- ggplot(fst_data_filtered, aes(x = Chromosome, y = MEAN_FST)) +
    geom_violin(fill = "lightblue", color = "black") +
    geom_boxplot(width = 0.1, color = "black") +
    labs(title = "",
         x = "AUS vs OTH",
         y = "FST") +
    theme_minimal() +
    theme(element_text(size=16),axis.text.x = element_blank(),axis.text.y = element_text(size = 12))

ggsave("FST_AUS_OTH_violinplot.pdf", plot = p, width = 4, height = 8, dpi = 300, bg = "white")

####################################
data<-read.table("Cultiar.AUS_OTH.SV.Fst.windowed.weir.fst",header = T)
chromosome_colors <- c("1H" = "#1f77b4", "2H" = "#ff7f0e", "3H" = "#2ca02c",
                       "4H" = "#d62728", "5H" = "#9467bd", "6H" = "#8c564b",
                       "7H" = "#e377c2")
data <- data %>%
    filter(CHROM %in% c("1H", "2H", "3H", "4H", "5H", "6H", "7H"))
data <- data %>%
    group_by(CHROM) %>%
    mutate(threshold = quantile(MEAN_FST, 0.99))
p<-ggplot(data, aes(x = BIN_START, y = MEAN_FST)) +
    geom_point(aes(color = CHROM, shape = MEAN_FST > threshold, size = MEAN_FST > threshold), alpha = 0.8) +
    scale_color_manual(values = chromosome_colors) +
    scale_shape_manual(values = c(16, 17)) +  # Circle and triangle for regular and highlighted points
    scale_size_manual(values = c(1, 2)) +  # Enlarge top 1% points
    labs(x = NULL, y = "Fst") +  # Remove the plot title
    theme_minimal() +
    theme(
        axis.text.x = element_blank(),  # Remove x-axis tick labels
        axis.title.x = element_blank(), # Remove x-axis title
        panel.grid.major.x = element_blank(),  # Remove vertical grid lines
        panel.grid.minor.x = element_blank(),  # Remove vertical grid lines
        legend.position = "none",  # Remove legend
        strip.text = element_text(size = 12, face = "bold"),  # Customize facet labels
        strip.placement = "outside",  # Place facet labels outside of the plot
        strip.background = element_rect(fill = "white", color = "black"),  # Background color of facet labels
        panel.spacing = unit(0.1, "lines")  # Adjust spacing between panels
    ) +
    facet_wrap(~ CHROM, scales = "free_x", nrow = 1)+theme(element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size = 16))
ggsave("FST_AUS_OTH_manhattan.pdf", plot = p, width = 4, height = 8, dpi = 300, bg = "white")
