## Draw a phenotypic correlation map
setwd("d:/0Murdoch/00.Barley/01.figure-make/")
library(PerformanceAnalytics)

# Read the dataset from CSV file
data <- read.csv("phenotype-cor.csv", header = TRUE)

# Use chart.Correlation to create a phenotypic correlation map
chart.Correlation(data)

## Draw QTL LOD picture
library(ggplot2)

# Create a scatter plot and line plot for QTL LOD
ggplot(data, aes(x = Position, y = LOD, color = TraitName)) +
  geom_point() +  # Add points
  geom_line() +   # Add lines
  theme_minimal() +  # Apply a minimal theme
  geom_hline(yintercept = 3, linetype = "dashed", color = "red") +  # Add a dashed line at y = 3 in red
  theme(
    text = element_text(family = "Times New Roman", face = "bold", size = 20),  # Modify font style
  ) +
  geom_text(x = 49, y = 3.5, label = "LOD = 3", hjust = -0.1, vjust = 0.5, color = "red", size = 6) +  # Add text annotation
  labs(x = "chr6H (cM)")  # Add x-axis label

