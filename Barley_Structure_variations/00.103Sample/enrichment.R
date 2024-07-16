library(ggplot2)
library(paletteer)

data<-read.csv("AUS_special_go.csv",header = T)
data %>% filter(pvalue<0.05) -> dat01
p<-ggplot(dat01[c(1:12),],aes(y=reorder(Description,GeneRatio),x=GeneRatio))+geom_point(aes(size=Count,color=pvalue))+scale_color_gradient(low="red",high="blue")+labs(y="")+theme_light()+theme(
    text = element_text(family = "Times New Roman", face = "bold", size = 20),  # Set text properties
    strip.text = element_text(color = "black"),  # Set the color of chromosome labels to black   
    plot.title = element_text(hjust = 0.5)  # Center the plot title
)
ggsave("GO.tiff", plot = p, width = 10, height = 6, dpi = 300, bg = "white")

data<-read.csv("AUS_special_kegg.csv",header = T)
data$GeneRatio1<-data$Count/138
data %>% filter(pvalue<0.05) -> dat01
p<-ggplot(dat01,aes(y=reorder(Description,GeneRatio1),x=GeneRatio1))+geom_point(aes(size=Count,color=pvalue))+scale_color_gradient(low="red",high="blue")+labs(y="",x="GeneRatio")+theme_light()+theme(
    text = element_text(family = "Times New Roman", face = "bold", size = 20),  # Set text properties
    strip.text = element_text(color = "black"),  # Set the color of chromosome labels to black   
    plot.title = element_text(hjust = 0.5)  # Center the plot title
)

ggsave("KEGG.tiff", plot = p, width = 10, height = 6, dpi = 300, bg = "white")
