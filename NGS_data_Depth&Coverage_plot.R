##The input file is from bedtools coverage
##example:
##chr1H	0	100000	28392	88991	100000	0.8899100
##chr1H	100000	200000	15141	73078	100000	0.7307800
##chr1H	200000	300000	11218	69491	100000	0.6949100
##chr1H	300000	400000	11049	50241	100000	0.5024100
##chr1H	400000	500000	26269	64803	100000	0.6480300

# Load ggplot2 library
library(ggplot2)

# Read data from files
data1 <- read.table("BASS_morexv1.coverage", header = FALSE)
data2 <- read.table("Flinders_morexv1.coverage", header = FALSE)

# Rename columns for better clarity
colnames(data1) <- c("Chromosome", "Start", "End", "Reads_number", "Length", "Window_size", "Coverage")

# Calculate Depth and Locus for data1
data1$Depth <- data1$Reads_number * 100 / data1$Length
data1$Locus <- (data1$Start + data1$End) / 2

# Rename columns for better clarity
colnames(data2) <- c("Chromosome", "Start", "End", "Reads_number", "Length", "Window_size", "Coverage")

# Calculate Depth and Locus for data2
data2$Depth <- data2$Reads_number * 100 / data2$Length
data2$Locus <- (data2$Start + data2$End) / 2




# Create a scatter plot for data1 and data2 with Depth values
ggplot() +
  geom_point(data = data1, mapping = aes(Locus, Depth), color = "blue") +  # Scatter plot for data1
  geom_point(data = data2, mapping = aes(Locus, -1 * Depth), color = "orange") +  # Scatter plot for data2

  facet_grid(Chromosome ~ .) +  # Facet the plot based on Chromosome variable

  ylim(-60, 60) +  # Set the y-axis limits

  theme_light() +  # Use a light theme

  labs(x = "Chromosomes", y = "Depth", title = "BASS vs Flinders Depth") +  # Add labels to axes and the title

  scale_x_continuous(
    breaks = c(0, 100000000, 200000000, 300000000, 400000000, 500000000, 600000000, 700000000, 800000000),
    labels = c("0Mbp", "", "200Mbp", "", "400Mbp", "", "600Mbp", "", "800Mbp")
  ) +  # Set breaks and labels for the x-axis

  theme(
    text = element_text(family = "Times New Roman", face = "bold", size = 20),  # Set text properties
    axis.text.x = element_text(angle = 30, hjust = 1),  # Rotate x-axis labels for better readability
    axis.text.y = element_text(margin = margin(r = -40)),  # Adjust the margin for y-axis labels
    strip.text = element_text(color = "black"),  # Set the color of chromosome labels to black   
    plot.title = element_text(hjust = 0.5)  # Center the plot title
  ) + # Set various theme elements
 # Add red dashed lines at y = -50 and y = 50, solid line at y = 0
  annotate("segment", x = 0, xend = 800000000, y = -30, yend = -30, color = "red", linetype = "dashed") +
  annotate("segment", x = 0, xend = 800000000, y = 30, yend = 30, color = "red", linetype = "dashed")+
  annotate("segment",x = 0, xend = 800000000, y = 0, yend = 0, color = "red", linetype = "solid")




# Create a scatter plot for data1 and data2 with Coverage values multiplied by 100
ggplot() +
  geom_point(data = data1, mapping = aes(Locus, 100 * Coverage), color = "blue") +
  geom_point(data = data2, mapping = aes(Locus, -100 * Coverage), color = "orange") +

  # Facet the plot based on Chromosome variable
  facet_grid(Chromosome ~ .) +

  # Set y-axis limits
  ylim(-100, 100) +

  # Use a light theme
  theme_light() +

  # Add labels to axes and the title
  labs(x = "Chromosomes", y = "Coverage(%)", title = "BASS vs Flinders Coverage") +

  # Set breaks and labels for the x-axis
  scale_x_continuous(
    breaks = c(0, 100000000, 200000000, 300000000, 400000000, 500000000, 600000000, 700000000, 800000000),
    labels = c("0Mbp", "", "200Mbp", "", "400Mbp", "", "600Mbp", "", "800Mbp")
  ) +

  # Set various theme elements
  theme(
    text = element_text(family = "Times New Roman", face = "bold", size = 20),
    axis.text.x = element_text(angle = 30, hjust = 1),  # Rotate x-axis labels for better readability
    axis.text.y = element_text(margin = margin(r = -40)),  # Adjust the margin for y-axis labels
    strip.text = element_text(color = "black"),  # Set the color of chromosome labels to black   
    plot.title = element_text(hjust = 0.5)  # Center the plot title
  ) +

  # Add red dashed lines at y = -50 and y = 50, solid line at y = 0
  annotate("segment", x = 0, xend = 800000000, y = -50, yend = -50, color = "red", linetype = "dashed") +
  annotate("segment", x = 0, xend = 800000000, y = 50, yend = 50, color = "red", linetype = "dashed") +
  annotate("segment",x = 0, xend = 800000000, y = 0, yend = 0, color = "red", linetype = "solid")
