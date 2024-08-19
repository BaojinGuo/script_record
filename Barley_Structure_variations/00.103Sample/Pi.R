##############violin+boxplot
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# Define file names and corresponding sample names
file_names <- c("Cultivar_AUS.SV.1M500k.windowed.pi",
                "Cultivar_OTH.SV.1M500k.windowed.pi",
                "Cultivar.SV.1M500k.windowed.pi",
                "Landrace.SV.1M500k.windowed.pi",
                "Wild.SV.1M500k.windowed.pi")

sample_names <- c("Cultivar_AUS", "Cultivar_OTH", "Cultivar", "Landrace", "Wild")

# Initialize an empty data frame to store all data
all_data <- data.frame()

# Read in data and combine into one data frame
for (i in seq_along(file_names)) {
  data <- read_tsv(file_names[i], col_names = TRUE)
  data <- data %>% 
    mutate(Sample = sample_names[i])
  all_data <- bind_rows(all_data, data)
}
all_data$Sample<-factor(all_data$Sample,levels = c("Cultivar_AUS","Cultivar_OTH","Cultivar","Landrace","Wild"))
# Plot violin-boxplot
p<-ggplot(all_data, aes(x = Sample, y = PI, fill = Sample)) +
    geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.8) +
    scale_fill_brewer(palette = "Set3") +
    theme_minimal() +
    labs(title = "",
         x = "",
         y = expression(pi)) +
    theme(legend.position = "none")+theme(text = element_text(size = 16))
ggsave("PI_violin_boxplot.pdf", plot = p, width = 8, height = 6, dpi = 300, bg = "white")

#######################lineplot
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(scales)  # For scale formatting

# Define file names and corresponding sample names
file_names <- c("Cultivar_AUS.SV.1M500k.windowed.pi",
                "Cultivar_OTH.SV.1M500k.windowed.pi",
                "Cultivar.SV.1M500k.windowed.pi",
                "Landrace.SV.1M500k.windowed.pi",
                "Wild.SV.1M500k.windowed.pi")

sample_names <- c("Cultivar_AUS", "Cultivar_OTH", "Cultivar", "Landrace", "Wild")

# Initialize an empty data frame to store all data
all_data <- data.frame()

# Read in data and combine into one data frame
for (i in seq_along(file_names)) {
    data <- read_tsv(file_names[i], col_names = TRUE)
    data <- data %>% 
        mutate(Sample = sample_names[i])
    all_data <- bind_rows(all_data, data)
}
all_data$Sample<-factor(all_data$Sample,levels = c("Cultivar_AUS","Cultivar_OTH","Cultivar","Landrace","Wild"))
# Filter the data to include only chromosomes 1H to 7H
all_data_filtered <- all_data %>%
    filter(CHROM %in% c("1H", "2H", "3H", "4H", "5H", "6H", "7H"))

# Identify high Pi values for Cultivar_AUS
highlight_points <- all_data_filtered %>%
    filter(Sample == "Cultivar_AUS" & PI == max(PI, na.rm = TRUE))

# Plot line plot for chromosomes 1H to 7H
p <- ggplot(all_data_filtered, aes(x = BIN_START, y = PI, color = Sample, group = Sample)) +
    geom_line(size = 0.5, alpha = 0.7) +  # Increased line transparency
    geom_point(data = highlight_points, aes(x = BIN_START, y = PI), color = "red", size = 3, shape = 21, fill = "red") +
    scale_color_brewer(palette = "Set3") +
    facet_wrap(~ CHROM, scales = "free_x", ncol = 1) +  # Arrange in one column
    scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = "Mb")) +  # Convert to Mb and avoid scientific notation
    theme_minimal() +
    labs(title = "",
         x = "",
         y = expression(pi)) +
    theme(legend.title = element_blank(), legend.position = "bottom",  
          legend.justification = c(1, 1),  
          strip.text = element_text(size = 12),
          axis.text.x = element_text( hjust = 1))

# Save the plot to a PDF file
ggsave("PI_lineplot.pdf", plot = p, width = 12, height = 18, dpi = 300, bg = "white")


######################################################################################################################################################
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(scales)  # For scale formatting

# Define file names and corresponding sample names
file_names <- c("Cultivar_AUS.SV.1M500k.windowed.pi",
                "Cultivar_OTH.SV.1M500k.windowed.pi",
                "Cultivar.SV.1M500k.windowed.pi"
          )

sample_names <- c("Cultivar_AUS", "Cultivar_OTH", "Cultivar")

# Initialize an empty data frame to store all data
all_data <- data.frame()

# Read in data and combine into one data frame
for (i in seq_along(file_names)) {
    data <- read_tsv(file_names[i], col_names = TRUE)
    data <- data %>% 
        mutate(Sample = sample_names[i])
    all_data <- bind_rows(all_data, data)
}
all_data$Sample<-factor(all_data$Sample,levels = c("Cultivar_AUS","Cultivar_OTH","Cultivar"))
# Filter the data to include only chromosomes 1H to 7H
all_data_filtered <- all_data %>%
    filter(CHROM %in% c("1H", "2H", "3H", "4H", "5H", "6H", "7H"))

# Identify high Pi values for Cultivar_AUS
highlight_points <- all_data_filtered %>%
    filter(Sample == "Cultivar_AUS" & PI == max(PI, na.rm = TRUE))

# Plot line plot for chromosomes 1H to 7H
p <- ggplot(all_data_filtered, aes(x = BIN_START, y = PI, color = Sample, group = Sample)) +
    geom_line(size = 0.5, alpha = 0.7) +  # Increased line transparency
    #geom_point(data = highlight_points, aes(x = BIN_START, y = PI), color = "red", size = 3, shape = 21, fill = "red") +
    scale_color_brewer(palette = "Set3") +
    facet_wrap(~ CHROM, scales = "free_x", ncol = 1) +  # Arrange in one column
    scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = "Mb")) +  # Convert to Mb and avoid scientific notation
    theme_minimal() +
    labs(title = "",
         x = "",
         y = expression(pi)) +
    theme(legend.title = element_blank(), legend.position = "bottom",  
          legend.justification = c(1, 1),  
          strip.text = element_text(size = 12),
          axis.text.x = element_text( hjust = 1))

# Save the plot to a PDF file
ggsave("PI_lineplot_cultivar.pdf", plot = p, width = 12, height = 18, dpi = 300, bg = "white")
