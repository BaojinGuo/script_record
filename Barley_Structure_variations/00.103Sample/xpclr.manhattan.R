# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read the data
data <- read.table("xpclr.manhattan.input", header = TRUE, sep = " ")



# Define colors for each chromosome
chromosome_colors <- c("1H" = "#1f77b4", "2H" = "#ff7f0e", "3H" = "#2ca02c",
                       "4H" = "#d62728", "5H" = "#9467bd", "6H" = "#8c564b",
                       "7H" = "#e377c2")

# Find the top 1% threshold for each chromosome
data <- data %>%
    group_by(Chromosome) %>%
    mutate(threshold = quantile(Xpclr, 0.99))

# Plot the data
p<-ggplot(data, aes(x = Position, y = Xpclr)) +
    geom_point(aes(color = Chromosome, shape = Xpclr > threshold, size = Xpclr > threshold), alpha = 0.8) +
    scale_color_manual(values = chromosome_colors) +
    scale_shape_manual(values = c(16, 17)) +  # Circle and triangle for regular and highlighted points
    scale_size_manual(values = c(1, 2)) +  # Enlarge top 1% points
    labs(x = NULL, y = "Xpclr") +  # Remove the plot title
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
    facet_wrap(~ Chromosome, scales = "free_x", nrow = 1)+theme(element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size = 16))
ggsave("Xpclr_manhattan.pdf", plot = p, width = 6, height = 6, dpi = 300, bg = "white")




