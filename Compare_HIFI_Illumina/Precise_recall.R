# Load necessary libraries
library(ggplot2)
library(dplyr)
library(scales)

data<-read.csv("Book1.csv")
#Variant_type	Samples	Depth	Reads_type	Precision	Recall	F-score
#SNP	Hockett	5X	Short reads	84.71%	84.39%	84.55%
#SNP	Hockett	10X	Short reads	83.79%	92.89%	88.10%
#SNP	Hockett	15X	Short reads	82.89%	93.89%	88.05%
#SNP	Igri	5X	Short reads	90.45%	89.39%	89.92%
#SNP	Igri	10X	Short reads	88.10%	97.55%	92.59%
#SNP	Igri	15X	Short reads	86.54%	98.69%	92.22%
#SNP	OUN333	5X	Short reads	89.22%	86.87%	88.03%
#SNP	OUN333	10X	Short reads	86.44%	97.41%	91.60%
#SNP	OUN333	15X	Short reads	84.75%	98.65%	91.18%
#SNP	RGT	5X	Short reads	90.44%	87.03%	88.70%
#SNP	RGT	10X	Short reads	87.84%	97.49%	92.41%
#SNP	RGT	15X	Short reads	86.23%	98.65%	92.03%
#SNP	Vlamingh	5X	Short reads	78.78%	84.74%	81.65%
#SNP	Vlamingh	10X	Short reads	73.73%	95.33%	83.15%
#SNP	Vlamingh	15X	Short reads	70.59%	97.48%	81.88%
#INDEL	Hockett	5X	Short reads	81.16%	73.54%	77.16%
#INDEL	Hockett	10X	Short reads	80.88%	85.88%	83.30%
#INDEL	Hockett	15X	Short reads	80.23%	88.43%	84.13%
#INDEL	Igri	5X	Short reads	88.39%	77.73%	82.72%
#INDEL	Igri	10X	Short reads	86.77%	89.40%	88.07%
#INDEL	Igri	15X	Short reads	85.47%	92.31%	88.76%
#INDEL	OUN333	5X	Short reads	87.09%	75.11%	80.66%
#INDEL	OUN333	10X	Short reads	85.19%	89.41%	87.25%
#INDEL	OUN333	15X	Short reads	83.77%	92.47%	87.91%
#INDEL	RGT	5X	Short reads	88.41%	75.59%	81.50%
#INDEL	RGT	10X	Short reads	86.63%	89.59%	88.09%


snp_data <- data %>% filter(Variant_type == "SNP")
snp_data$Sample_Depth <- paste(snp_data$Samples, snp_data$Depth, sep = "_")
snp_data$Group <- paste(snp_data$Samples, snp_data$Reads_type, sep = "_")
snp_data$Sample_Depth <- factor(snp_data$Sample_Depth, levels = c("Hockett_5X", "Hockett_10X", "Hockett_15X","Igri_5X", "Igri_10X", "Igri_15X","OUN333_5X", "OUN333_10X", "OUN333_15X","RGT_5X", "RGT_10X", "RGT_15X","Vlamingh_5X", "Vlamingh_10X", "Vlamingh_15X"))

shape_palette <- c(
  "Short reads" = 16, # Circle
  "Long reads" = 17   # Triangle
)



my_palette <- c(
    Hockett1 = "#1F77B4", Hockett2 = "#69A4CD", Hockett3 = "#B4D1E6",
    Igri1 = "#FF7F0E", Igri2 = "#FFA95E", Igri3 = "#FFD4AE",
    OUN3331 = "#2CA02C", OUN3332 = "#72BF72", OUN3333 = "#B8DFB8",
    RGT1 = "#D62728", RGT2 = "#E36F6F", RGT3 = "#F1B7B7",
    Vlamingh1 = "#9467BD", Vlamingh2 = "#B799D3", Vlamingh3 = "#DBCCE9"
)

sample_depth_levels <- unique(snp_data$Sample_Depth)
colors <- setNames(my_palette[1:length(sample_depth_levels)],sample_depth_levels)

snp_data$Recall <- as.numeric(gsub("%", "", snp_data$Recall)) / 100
snp_data$Precision <- as.numeric(gsub("%", "", snp_data$Precision)) / 100


snp<-ggplot(snp_data, aes(x = Recall, y = Precision, color = Sample_Depth, shape = Reads_type)) +
    geom_point(size = 3) +
    geom_line(aes(group = Group), linetype = "solid") +
    scale_color_manual(values = colors, labels = function(x) gsub("_", " Depth: ", x)) +
    scale_shape_manual(values = shape_palette) +
    theme_minimal() +
    labs(
        title = "SNP Precision vs Recall",
        x = "Recall",
        y = "Precision",
        color = "Sample + Depth",
        shape = "Reads Type"
    ) +
    theme(
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
    ) +  xlim(0.5, 1) + ylim(0, 1)

ggsave("SNPrecall.tiff", plot = snp, width = 10, height = 6, dpi = 300, bg = "white")



snp_data$F.score <- as.numeric(gsub("%", "", snp_data$F.score)) / 100
snp_data$Depth<-factor(snp_data$Depth,levels = c("5X","10X","15X"))
snpf<-ggplot(snp_data,aes(x=Depth,y=F.score,fill=Reads_type))+geom_boxplot()+theme_minimal() +scale_fill_manual(values = brewer.pal(n = 3, name = "Set3")) +theme(text = element_text(size = 16),legend.position = "none")
ggsave("SNPfscore.tiff", plot = snpf, width = 8, height = 6, dpi = 300, bg = "white")





                       
                       
