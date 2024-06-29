## Draw a phenotypic correlation map
setwd("d:/0Murdoch/00.Barley/01.figure-make/")
library(PerformanceAnalytics)

# Read the dataset from CSV file
data <- read.csv("phenotype-cor.csv", header = TRUE)
# head(data)

#BDT21.TKW.obs.mean	MEI21.TKW.obs.mean	MER21.TKWAs.obs.mean	WHS21.TKW.obs.mean
#53.8	46.5	47.2	44.6
#50.7	42.7	40	42.9
#49.9	44.1	42.5	42.9
#53	44	47.1	46.5
#53.2	43	47.15	46.5
#55.8	46.8	43.2	46.6
#NA	44.4	45.7	43.9
#57	45.05	46	47.3

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

