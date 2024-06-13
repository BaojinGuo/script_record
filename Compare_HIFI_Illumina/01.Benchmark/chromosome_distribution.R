library(ggplot2)
library(dplyr)
library(RColorBrewer)

H <- read.delim("Hockett.summary", header = FALSE)
R <- read.delim("RGT.summary", header = FALSE)
O <- read.delim("OUN333.summary", header = FALSE)
V <- read.delim("Vlamingh.summary", header = FALSE)
I <- read.delim("Igri.summary", header = FALSE)
# 添加文件标识
H$file <- "Hockett"
R$file <- "RGT"
O$file <- "OUN333"
V$file <- "Vlamingh"
I$file <- "Igri" 
data <- rbind(H, R, O, V, I)
colnames(data)<-c("Chromosome","POS","Len","TYPE","1","2","3","Accession")
data$Chromosome <- as.factor(data$Chromosome)
data$POS <- as.numeric(data$POS)

# 只选择1H到7H染色体的数据
data <- data %>%
  filter(Chromosome %in% c("1H", "2H", "3H", "4H", "5H", "6H", "7H"))
# Filter for only deletions (DEL)


#data <- data %>%
#  filter(TYPE == "DEL")


data <- data %>%
    mutate(Chromosome_Accession = paste0(Chromosome, "-", Accession))

color_palette <- brewer.pal(n = 5, name = "Set3")

# 绘图
p<-ggplot(data, aes(x = POS, y = Chromosome_Accession, color = Accession)) +
    geom_point(shape = '|', size = 5, alpha = 0.7) +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    theme(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        text = element_text(family = "Times New Roman", face = "bold", size = 20),
        axis.text = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16)
    ) +
    scale_color_manual(values = color_palette) +
    scale_x_continuous(
        breaks = c(0, 100000000, 200000000, 300000000, 400000000, 500000000, 600000000, 700000000, 800000000),
        labels = c("0Mbp", "", "200Mbp", "", "400Mbp", "", "600Mbp", "", "800Mbp")
    ) +
    theme(
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5, angle = 0),
        panel.spacing.y = unit(1, "lines")
    )+ theme(legend.position = "none")
ggsave("SVdistribution_plot.tiff", plot = p, width = 18, height = 8,  dpi = 300,bg = "white")


