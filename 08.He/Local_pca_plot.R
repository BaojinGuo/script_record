
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
        labs(title = paste(chr),
             x = "Position (Mb)", y = "MDS1") +
        theme_minimal()+ylim(-0.85,0.35)
    p1_list[[chr]] <- p1
    # 图 2: Position vs MDS2
    p2 <- ggplot(chr_data, aes(x = Position, y = MDS2)) +
        geom_point(size=0.5) +
        labs(title = paste(chr),
             x = "Position (Mb)", y = "MDS2") +
        theme_minimal()+ylim(-0.65,0.45)
    p2_list[[chr]] <- p2
}

p1_combined <- marrangeGrob(grobs = p1_list, ncol = 1, nrow = 7, top = "Position vs MDS1")
p2_combined <- marrangeGrob(grobs = p2_list, ncol = 1, nrow = 7, top = "Position vs MDS2")
ggsave("Combined_Position_vs_MDS1.pdf", p1_combined, width = 10, height = 12, units = "in")
ggsave("Combined_Position_vs_MDS2.pdf", p2_combined, width = 10, height = 12, units = "in")

