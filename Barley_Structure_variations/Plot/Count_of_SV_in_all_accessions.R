library(ggplot2)
data<-read.table("barley.all.SV.plot.txt",header = T)
head(data)
#  Accession Chromosome Position Type
#  Baudin         1H   180730  INS
#  Baudin         1H   182233  DEL
#  Baudin         1H   185665  INS
# Set the order of types
data$Type <- factor(data$Type, levels = c("INS", "DEL", "DUP", "INV"))

# Plot
ggplot(data, aes(x = Accession, fill = Accession)) +
    geom_bar(position = "dodge") + theme_minimal() +
    geom_text(stat = "count", aes(label = after_stat(count)), position = position_dodge(width = 0.9), hjust = 1, vjust = 0, angle = 90) +  # Add count labels on top
    labs(x = "Accession", y = "", title = "Count of Structure Vriations") +
    facet_grid(Type ~ ., scales = "free", switch = "y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ # Rotate x-axis labels for better readability
    theme(text = element_text(family = "Times New Roman", face = "bold", size = 20))+
    theme(axis.text = element_text(size = 16,face = "bold"))+
    theme(legend.text = element_text(size = 16))+
    theme(plot.title = element_text(hjust = 0.5))+
    guides(fill = guide_legend(ncol = 1))  # Arrange legend items in one column
