# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)

# Define file names
file_names <- c("Cultivar_AUS.SV.1M.Tajima.D",
                "Cultivar_OTH.SV.1M.Tajima.D",
                "Cultivar.SV.1M.Tajima.D",
                "Landrace.SV.1M.Tajima.D",
                "Wild.SV.1M.Tajima.D")

# Function to read and prepare data
read_tajima_data <- function(file) {
    data <- read_tsv(file, col_names = TRUE) %>%
        mutate(Sample = gsub(".SV.1M.Tajima.D", "", file))  # Add sample name
    return(data)
}

# Read and combine data from all files
tajima_data <- bind_rows(lapply(file_names, read_tajima_data))

# Filter the data to include only chromosomes 1H to 7H
tajima_data_filtered <- tajima_data %>%
    filter(CHROM %in% c("1H", "2H", "3H", "4H", "5H", "6H", "7H")
tajima_data_filtered$Sample<-factor(tajima_data_filtered$Sample,levels = c("Cultivar_AUS","Cultivar_OTH","Cultivar","Landrace","Wild"))
p <- ggplot(tajima_data_filtered, aes(x = BIN_START / 1e6, y = TajimaD, shape = Sample, fill = Sample)) +
    geom_point(size = 2, stroke = 0.5, color = "black", alpha = 0.5) +  # Add black borders and transparency
    scale_fill_brewer(palette = "Set3") +  # Use Set3 palette for fill colors
    scale_shape_manual(values = c(21, 22, 23, 24, 25)) +  # Assign different shapes to each sample
    facet_wrap(~ CHROM, ncol = 1) +  # Separate plots for each chromosome (1H to 7H)
    labs(title = "",
         x = "",
         y = "Tajima's D") +
    theme_minimal() +
    theme(axis.text.x = element_text( hjust = 1),
          legend.position = "top",  # Position legend at the top
          strip.text = element_text(size = 10)) +
    scale_x_continuous(expand = c(0.01, 0), labels = scales::number_format(scale = 1, suffix = " Mb"))+theme(legend.position = "bottom")
ggsave("TajimaD_point_plot.pdf", plot = p, width = 12, height = 14, dpi = 300, bg = "white")









           
