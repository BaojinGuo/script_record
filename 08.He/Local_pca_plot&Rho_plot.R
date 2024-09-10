
library(ggplot2)

data<-read.csv("mds_coords.csv",header = T)

library(dplyr)

chromosomes <- unique(data$Chromosome)


for (chr in chromosomes) {
    chr_data <- filter(data, Chromosome == chr)
    
    # 图 1: Position vs MDS1
    p1 <- ggplot(chr_data, aes(x = Position, y = MDS1)) +
        geom_point(size=0.5) +
        labs(title = paste(chr, "- Position vs MDS1"),
             x = "Position (Mb)", y = "MDS1") +
        theme_minimal()+ylim(-0.85,0.35)
    
    # 图 2: Position vs MDS2
    p2 <- ggplot(chr_data, aes(x = Position, y = MDS2)) +
        geom_point(size=0.5) +
        labs(title = paste(chr, "- Position vs MDS2"),
             x = "Position (Mb)", y = "MDS2") +
        theme_minimal()+ylim(-0.65,0.45)
    
    # 图 3: MDS1 vs MDS2
    p3 <- ggplot(chr_data, aes(x = MDS1, y = MDS2)) +
        geom_point(size=0.5) +
        labs(title = paste(chr, "- MDS1 vs MDS2"),
             x = "MDS1", y = "MDS2") +
        theme_minimal()+xlim(-0.85,0.35)+ylim(-0.65,0.45)
    
    # 保存图表为 PDF 格式
    ggsave(filename = paste0(chr, "_Position_vs_MDS1.pdf"), plot = p1, width = 6, height = 4)
    ggsave(filename = paste0(chr, "_Position_vs_MDS2.pdf"), plot = p2, width = 6, height = 4)
    ggsave(filename = paste0(chr, "_MDS1_vs_MDS2.pdf"), plot = p3, width = 6, height = 4)
}

###################################
library(gridExtra)

p1_list <- list()
p2_list <- list()

for (chr in chromosomes) {
    chr_data <- filter(data, Chromosome == chr)
    
    # 图 1: Position vs MDS1
    p1 <- ggplot(chr_data, aes(x = Position, y = MDS1)) +
        geom_point(size=0.5) +
        labs(title = "",
             x = paste(chr, " (Mb)"), y = "") +
        theme_minimal()+ylim(-1,0.4)+xlim(0,650)
    p1_list[[chr]] <- p1
    # 图 2: Position vs MDS2
    p2 <- ggplot(chr_data, aes(x = Position, y = MDS2)) +
        geom_point(size=0.5) +
        labs(title = "",
             x = paste(chr," (Mb)"), y = "") +
        theme_minimal()+ylim(-0.8,0.6)+xlim(0,650)
    p2_list[[chr]] <- p2
}

p1_combined <- marrangeGrob(grobs = p1_list, ncol = 1, nrow = 7, top = "Position vs MDS1")
p2_combined <- marrangeGrob(grobs = p2_list, ncol = 1, nrow = 7, top = "Position vs MDS2")
ggsave("Combined_Position_vs_MDS1.pdf", p1_combined, width = 10, height = 12, dpi=300)
ggsave("Combined_Position_vs_MDS2.pdf", p2_combined, width = 10, height = 12, dpi=300)

###################################################################################################
##FASTEPRR Rho value plot
# Load necessary libraries
library(ggplot2)
library(gridExtra)

# Define the folder path where the files are stored
folder_path <- "path/to/your/files/"


# File names and corresponding chromosome labels
files <- list("chrchr1H.txt", "chrchr2H.txt", "chrchr3H.txt", "chrchr4H.txt", 
              "chrchr5H.txt", "chrchr6H.txt", "chrchr7H.txt")
chr_names <- c("1H", "2H", "3H", "4H", "5H", "6H", "7H")

plot_rho <- function(file, chr_name) {
    # Construct the full file path
    file_path <- paste0(folder_path, file)
    
    # Read the file
    data <- read.table(file_path, header = TRUE)
    
    # Create position by taking the midpoint of the Start and End, and convert to Mb
    data$Position <- (data$Start + data$End) / 2 / 1000000  # Convert to Mb
    
    # Create the plot
    p <- ggplot(data, aes(x = Position, y = Rho)) +
        geom_point(size=0.5) +
        theme_minimal() +
        xlab(paste(chr_name, "(Mb)")) +  # Update the x-axis label
        ylab("Rho") +xlim(0,680)+ylim(0,650)
        theme(plot.title = element_blank())  # Correctly remove the plot title
    
    return(p)
}



plots <- lapply(1:7, function(i) plot_rho(files[i], chr_names[i]))

# Arrange the plots in one row with 7 columns
grid.arrange(grobs = plots, ncol = 1)

ggsave("rho_plots_mb.pdf", grid.arrange(grobs = plots, ncol = 1), width = 10, height = 12,dpi=300)

