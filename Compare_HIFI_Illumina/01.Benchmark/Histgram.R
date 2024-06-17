# Load necessary libraries
library(ggplot2)
library(tidyr)
library(RColorBrewer)

# Read the data
data<-read.delim("Illumina.benchmark.length",col.names = c("Hockett","Igri","OUN333","RGT Planet","Vlamingh"))
# Convert the data to long format for ggplot2
data_long <- gather(data, key = "Sample", value = "Length", Hockett, Igri, OUN333, RGT.Planet, Vlamingh)
# Plot the data
p<-ggplot(data_long, aes(x = Length, fill = Sample, color = Sample)) +
    geom_histogram(alpha = 1,bins = 100) +  # Use density plot with increased transparency
    labs(title = "", x = "INDEL Length (bp)", y = "") +
    theme_minimal() +
    theme(legend.position = "top") +
    scale_fill_brewer(palette = "Set3") +
    scale_color_brewer(palette = "Set3") +
    scale_y_continuous(labels = scales::comma) +  # Avoid scientific notation
    facet_wrap(~ Sample, ncol = 1, scales = "free_y") + # Arrange panels in one column with five rows
    theme(legend.title = element_blank(),
        text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5))+  
    xlim(-50,50)
ggsave("Illumina_indel_length_distribution.tiff", plot = p, width = 8, height = 6, dpi = 300, bg = "white")
