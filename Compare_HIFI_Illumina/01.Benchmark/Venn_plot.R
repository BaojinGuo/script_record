#####

library(ggVennDiagram)
library(ggplot2)
data<-read.delim("venn.input",header = T)
##input file
#Hockett.cuteSV Hockett.sniffles  Hockett.SVIM   Igri.cuteSV Igri.sniffles     Igri.SVIM OUN333.cuteSV OUN333.sniffles  OUN333.SVIM   RGT.cuteSV RGT.sniffles     RGT.SVIM Vlamingh.cuteSV Vlamingh.sniffles Vlamingh.SVIM
#1H_6150_DEL     1H_16986_DEL  1H_16986_DEL  1H_16985_DEL  1H_16986_DEL  1H_16985_DEL  1H_55827_DEL    1H_55828_DEL 1H_55828_DEL 1H_24767_INS 1H_24767_INS 1H_32488_DEL     1H_6150_DEL       1H_6150_DEL   1H_6150_DEL
#1H_16985_DEL     1H_64815_INS  1H_64815_INS  1H_64815_INS  1H_64815_INS  1H_64815_INS  1H_64028_DEL    1H_64028_DEL 1H_64028_DEL 1H_50844_INS 1H_50844_INS 1H_50844_INS   1H_180730_INS     1H_180730_INS 1H_180730_INS
#1H_64815_INS     1H_71483_DEL  1H_71483_DEL  1H_71482_DEL  1H_71483_DEL  1H_71483_DEL  1H_87645_DEL    1H_87518_DEL 1H_87520_DEL 1H_56036_DEL 1H_81805_INS 1H_55898_DEL   1H_182278_DEL     1H_182273_DEL 1H_182279_DEL
#1H_71483_DEL    1H_107760_INS  1H_88423_DEL  1H_88402_DEL 1H_161434_INS  1H_88403_DEL  1H_91183_INS    1H_88199_DEL 1H_88258_DEL 1H_81805_INS 1H_83690_DEL 1H_81805_INS   1H_183486_DEL     1H_183487_DEL 1H_183487_DEL
#1H_88423_DEL    1H_161434_INS  1H_96621_DEL 1H_107760_INS 1H_163845_DEL 1H_163845_DEL 1H_141169_INS    1H_88782_DEL 1H_88782_DEL 1H_83689_DEL 1H_83926_DEL 1H_83690_DEL   1H_185665_INS     1H_185665_INS 1H_185665_INS
#1H_96691_DEL    1H_163845_DEL 1H_163845_DEL 1H_161415_INS 1H_165695_DEL 1H_165695_DEL 1H_142042_DEL    1H_91183_INS 1H_89448_DEL 1H_84060_INS 1H_84060_INS 1H_83926_DEL   1H_188906_INS     1H_188906_INS 1H_188906_INS


x<-list(cuteSV=data$Hockett.cuteSV,sniffles=data$Hockett.sniffles,SVIM=data$Hockett.SVIM)

p<-ggVennDiagram(x,label_alpha = 0,set_size = 8,label_size = 8)+
  scale_fill_distiller(palette = "Reds",direction = 1)+ 
  ggtitle("Hockett") + 
  theme(plot.title = element_text(hjust = 0.5, size = 28,face = "bold"))

ggsave("Hockett_venn_diagram.tiff", plot = p, width = 10, height = 8, dpi = 600,bg="white")


####the following is the loop for all sample
# Script Purpose: This script reads data from a file containing outputs from SV caller detection tools, groups the data into predefined sets, generates Venn diagrams for each set, and saves the generated plots to files.

# Import necessary libraries
library(ggVennDiagram)
library(ggplot2)
library(dplyr)

# Read the data file
data <- read.delim("venn.input", header = TRUE)

# Define comparisons for different groups, with each group consisting of three columns
groups <- list(
  Hockett = c("Hockett.cuteSV", "Hockett.sniffles", "Hockett.SVIM"),
  Igri = c("Igri.cuteSV", "Igri.sniffles", "Igri.SVIM"),
  OUN333 = c("OUN333.cuteSV", "OUN333.sniffles", "OUN333.SVIM"),
  RGT = c("RGT.cuteSV", "RGT.sniffles", "RGT.SVIM"),
  Vlamingh = c("Vlamingh.cuteSV", "Vlamingh.sniffles", "Vlamingh.SVIM")
)

# Loop through each group to generate and save Venn diagrams
for (group_name in names(groups)) {
  # Get the column names for the current group
  group <- groups[[group_name]]
  
  # Create data list for the Venn diagram
  x <- list(
    cuteSV = data[[group[1]]],
    sniffles = data[[group[2]]],
    SVIM = data[[group[3]]]
  )
  
  # Generate the Venn diagram
  p <- ggVennDiagram(x, label_alpha = 0, set_size = 8, label_size = 8) +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    ggtitle(group_name) +
    theme(plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
          panel.background = element_rect(fill = "white", color = "white"),
          plot.background = element_rect(fill = "white", color = "white"))
  
  # Save the plot
  ggsave(paste0(group_name, "_venn_diagram.tiff"), plot = p, width = 10, height = 8, dpi = 600, bg = "white")
}
