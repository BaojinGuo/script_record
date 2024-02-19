
library(ggplot2)
data<-read.csv("anova.csv",header = T)

head(data)
#  Environment  Genotype  TKW
#  BDT  17MUB002-SSD002  199
#  BDT  17MUB002-SSD003  164
#  BDT  17MUB002-SSD004  157
#  BDT  17MUB002-SSD005  191
#  BDT  17MUB002-SSD006  193
#  BDT  17MUB002-SSD007  219

# Convert the "Environment" column to a factor variable
data$Environment <- as.factor(data$Environment)

# Convert the "TKW" column to numeric
data$TKW <- as.numeric(data$TKW)

# Use ggplot2 to create a violin plot and boxplot
ggplot(data, aes(x = Environment, y = TKW, fill = Environment)) +
  geom_violin(trim = FALSE) +  # Add a violin plot without trimming the tails
  geom_boxplot(width = 0.1, fill = "white", color = "black") +  # Add a boxplot with specified width and styling
  
  # Add labels to the axes
  labs(x = "Environment", y = "Thousand Kernel Weight (g)") +
  
  # Apply a minimal theme to the plot
  theme_minimal() +
  
  # Modify font style using Times New Roman, size 16, and bold
  theme(text = element_text(family = "Times New Roman", size = 16, face = "bold"))
