


# Load required libraries
library(ggplot2)
library(plotly)
library(reshape2)


data<-read.csv("Book1.csv")

shapes <- c("Asia" = 24, "Europe" = 22, "Africa" = 21, "North_America" = 23, "Australia" = 25)
colors <- c("wild" = "#1B9E77", "landrace" = "#D95F02", "cultivar" = "#7570B3")

# Create 2D PCA plots
p1 <- ggplot(data, aes(x = PC1, y = PC2, color = Status, shape = Region)) +
    geom_point(size = 3,stroke = 1.5) +
    scale_shape_manual(values = shapes) +
    scale_color_manual(values = colors) +
    labs(title = "2D PCA Plot: PC1 vs PC2", x = "PC1", y = "PC2") + theme_minimal()+
    theme(text = element_text(family = "Times New Roman", face = "bold", size = 20))

p2 <- ggplot(data, aes(x = PC1, y = PC3, color = Status, shape = Region)) +
    geom_point(size = 3,stroke = 1.5) +
    scale_shape_manual(values = shapes) +
    scale_color_manual(values = colors) +
    labs(title = "2D PCA Plot: PC1 vs PC3", x = "PC1", y = "PC3") + theme_minimal()+
    theme(text = element_text(family = "Times New Roman", face = "bold", size = 20))
   

p3 <- ggplot(data, aes(x = PC2, y = PC3, color = Status, shape = Region)) +
    geom_point(size = 3,stroke = 1.5) +
    scale_shape_manual(values = shapes) +
    scale_color_manual(values = colors) +
    labs(title = "2D PCA Plot: PC2 vs PC3", x = "PC2", y = "PC3") + theme_minimal()+
    theme(text = element_text(family = "Times New Roman", face = "bold", size = 20))

# Create 3D PCA plot
p4 <- plot_ly(data, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Status, symbol = ~Region, symbols = c('circle', 'square', 'diamond')) %>%
    add_markers() %>%
    layout(title = "3D PCA Plot",
           scene = list(xaxis = list(title = 'PC1'),
                        yaxis = list(title = 'PC2'),
                        zaxis = list(title = 'PC3')))

ggsave("PC1_2.tiff", plot = p1, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("PC1_3.tiff", plot = p2, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("PC2_3.tiff", plot = p3, width = 10, height = 6, dpi = 300, bg = "white")

library(scatterplot3d)

data$Status <- as.factor(data$Status)
data$Region <- as.factor(data$Region)

# Define the shapes and colors
shapes <- c("Asia" = 24, "Europe" = 22, "Africa" = 21, "North_America" = 23, "Australia" = 25)
colors <- c("wild" = "#1B9E77", "landrace" = "#D95F02", "cultivar" = "#7570B3")

# Map shapes and colors to the data
data$shape <- shapes[as.character(data$Region)]
data$color <- colors[as.character(data$Status)]

# Set text parameters for the plot
par(family = "Times", font = 1, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)

# Plot 3D PCA using scatterplot3d
s3d <- scatterplot3d(data$PC1, data$PC2, data$PC3, pch = as.numeric(data$shape), color = data$color,
                     xlab = "PC1", ylab = "PC2", zlab = "PC3", main = "3D PCA Plot",
                     angle = 55, scale.y = 0.7)

# Reset par to default to avoid affecting other plots
par(family = "Times New Roman", font = 20, cex.axis = 1, cex.lab = 1, cex.main = 1)

# Add legend
legend("topleft", legend = levels(data$Status), col = colors[levels(data$Status)], pch = 16, title = "Status")
legend("topright", legend = names(shapes), pch = shapes, title = "Region")
