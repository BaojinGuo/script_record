####for one K value
library(ggplot2)
library(reshape2)



tbl=read.table("filter.snp.3.Q") 

out.dist=dist(tbl,method="euclidean")
out.hclust=hclust(out.dist,method = "ward.D")
order_s<-out.hclust$order
df_sample <- read.table("filter.snp.nosex", header = F, sep = "\t")
df_sample  <- df_sample [,1]

df_sample <- df_sample[order_s]
i <- 3
col_panel <- rainbow(i)
file <- paste("filter.snp.", i, ".Q", sep = "")
df <- read.table(file, header = F, sep = " ")
df <- df[order_s,]
df_new <- transform.data.frame(samples=df_sample, df)
aql <- melt(df_new, id.vars = "samples")
aql$samples <- factor(x=aql$samples, levels = df_sample)
y_lab <- paste("K=", i, sep = "")
p4 <- ggplot(aql) + geom_bar(aes(x=samples, weight=value, fill=variable), position = "stack", width = 1) +
    scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) + 
    scale_fill_manual(values = col_panel ) + ####可根据需要调整颜色
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, size = 10),
          axis.title.x = element_blank(),
          axis.ticks = element_blank(), axis.title.y = element_text(size = 13),
          panel.grid = element_blank()) + ylab(y_lab)

###################################################################################################
##a loop to plot K=2 to K=10
library(ggplot2)
library(reshape2)

# 读取样本文件
df_sample <- read.table("filter.snp.nosex", header = F, sep = "\t")
df_sample <- df_sample[,1]

for (i in 2:10) {
    # 读取对应K值的文件
    file <- paste("filter.snp.", i, ".Q", sep = "")
    tbl <- read.table(file, header = F, sep = " ")
    
    # 计算样本顺序
    out.dist <- dist(tbl, method = "euclidean")
    out.hclust <- hclust(out.dist, method = "ward.D")
    order_s <- out.hclust$order
    df_sample_ordered <- df_sample[order_s]
    
    # 重新排序数据
    df <- tbl[order_s,]
    df_new <- transform.data.frame(samples = df_sample_ordered, df)
    aql <- melt(df_new, id.vars = "samples")
    aql$samples <- factor(x = aql$samples, levels = df_sample_ordered)
    
    # 生成颜色面板
    col_panel <- rainbow(i)
    
    # 生成图表标签
    y_lab <- paste("K=", i, sep = "")
    
    # 绘制图表
    p <- ggplot(aql) + 
        geom_bar(aes(x = samples, weight = value, fill = variable), position = "stack", width = 1) +
        scale_x_discrete(expand = c(0, 0)) + 
        scale_y_continuous(expand = c(0, 0)) + 
        scale_fill_manual(values = col_panel) + 
        theme(
            legend.position = "none",
            panel.background = element_blank(),
            axis.text.x = element_text(angle = 90, size = 10, family = "Times New Roman", face = "bold", hjust = 1,vjust = 0.4),
            axis.title.x = element_blank(),
            axis.ticks = element_blank(),
            axis.title.y = element_text(size = 13, family = "Times New Roman", face = "bold"),
            panel.grid = element_blank()
        ) + 
        ylab(y_lab)+
        theme(element_text(face="bold"))
    
    # 保存图像
    ggsave(paste("structure_plot_K", i, ".png", sep = ""), plot = p, width = 12, height = 4, dpi = 300)
}


