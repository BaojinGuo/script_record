
library(ggplot2)
data<-read.csv("multi.sniffles.sv.csv",header = T)
head(data)
#Chromosome    POS TYPE    length
#        1H 182233  DEL 516505932
#        1H 185665  INS 516505932
#        1H 188906  INS 516505932
#        1H 241466  DEL 516505932
#        1H 241963  DEL 516505932
str(data)
#'data.frame':	65917 obs. of  4 variables:
# $ Chromosome: Factor w/ 7 levels "1H","2H","3H",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ POS       : num  182233 185665 188906 241466 241963 ...
# $ TYPE      : Factor w/ 4 levels "INS","DEL","DUP",..: 2 1 1 2 2 1 2 1 2 2 ...
# $ length    : int  516505932 516505932 516505932 516505932 516505932 516505932 516505932 516505932 516505932 516505932 ...

data$Chromosome<-as.factor(data$Chromosome)
data$POS<-as.numeric(data$POS)
data$TYPE<-as.factor(data$TYPE)
desired_order <-c("INS","DEL","DUP","INV")
data$TYPE<-factor(data$TYPE,levels=desired_order)
custom_colors <- c("DEL" = "#66c2a5", "DUP" = "#fc8d62", "INS" = "#8da0cb", "INV" = "#e78ac3")
ggplot(data, aes(x = POS, y = Chromosome, color = TYPE)) +
    geom_point(shape = '|', size = 5, alpha = 0.7) +
    labs(title = "Chromosomal Variants", x = "Position", y = "Chromosome") +
    theme_minimal() +
    theme(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5)
    ) +
    scale_color_manual(values = custom_colors)+ scale_x_continuous(breaks = c(0, 100000000, 200000000, 300000000, 400000000, 500000000, 600000000, 700000000, 800000000), labels = c("0Mbp", "", "200Mbp", "", "400Mbp", "", "600Mbp", "", "800Mbp"))+
    theme(text = element_text(family = "Times New Roman", face = "bold", size = 20))+
    theme(axis.text = element_text(size = 16,face = "bold"))+
    theme(legend.text = element_text(size = 16))


