

library(ggplot2)
library(treemapify)


data<-read.csv("S2_cons_go_all_brdLJ_treemap.csv",header = T)
ggplot(data, aes(area = Count, fill = Category, subgroup = Category, label = Description)) +
    geom_treemap(colour ="black", size = 0.5) +
    geom_treemap_subgroup_border(colour = "black", size = 4) +
    geom_treemap_text(colour = "black", place = "centre", reflow = TRUE) +  # 只显示每个块的 Description
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") 
