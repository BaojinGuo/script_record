
library(ggplot2)
data<-read.csv("anova.csv",header = T)

head(data)
  Environment        Genotype TKW
1         BDT 17MUB002-SSD002 199
2         BDT 17MUB002-SSD003 164
3         BDT 17MUB002-SSD004 157
4         BDT 17MUB002-SSD005 191
5         BDT 17MUB002-SSD006 193
6         BDT 17MUB002-SSD007 219

data$Environment<-as.factor(data$Environment)
data$TKW<-as.numeric(data$TKW)

ggplot(data , aes(x = Environment, y = TKW, fill = Environment)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", color = "black") +
    labs(x = "Environment", y = "Thousand Kernel Weight (g)") +
    theme_minimal()
    +theme(text = element_text(family = "Times New Roman", size = 16, face = "bold"))
