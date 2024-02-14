##draw a phenotypic correlation map
setwd("d:/0Murdoch/00.Barley/01.figure-make/")
library(PerformanceAnalytics)
data<-read.csv("phenotype-cor.csv",header = T)
chart.Correlation(data)

##draw qtl LOD picture
library(ggplot2)
ggplot(data,aes(x=Position,y = LOD,color=TraitName))+geom_point()+geom_line()+theme_minimal()+geom_hline(yintercept = 3, linetype = "dashed", color = "red")+theme(text = element_text(family = "Times New Roman",face = "bold",size = 20))+geom_text(x= 49,y = 3.5, label = "LOD = 3", hjust = -0.1, vjust = 0.5, color = "red", size = 6)+labs(x="chr6H (cM)")


