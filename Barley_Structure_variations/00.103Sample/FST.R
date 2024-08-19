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

