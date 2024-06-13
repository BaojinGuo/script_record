


# 设置文件名和颜色
file_names <- c("Hockett.summary", "RGT.summary", "OUN333.summary", "Vlamingh.summary", "Igri.summary")
colors <- brewer.pal(n = 5, name = "Set3")

# 读取每个文件并处理第四列数据
data_list <- lapply(file_names, function(file) {
    data <- read.delim(file, header = FALSE)
    colnames(data) <- c("V1", "V2", "V3", "Type", "V5", "V6", "V7")
    type_counts <- table(data$Type)
    type_counts <- data.frame(Type = names(type_counts), Freq = as.numeric(type_counts), file = file)  # 转换为数据框并添加文件名列
    type_counts <- type_counts[type_counts$Type %in% c("DEL", "INS", "INV", "DUP"), ]  # 筛选出指定的类型
    return(type_counts)
})

# 合并数据
combined_data <- do.call(rbind, data_list)

# 设置文件列的因子变量，确保顺序和文件名一致
combined_data$file <- factor(combined_data$file, levels = file_names)
combined_data$file<-gsub("\\.summary", "", combined_data$file)
# 创建左右两个Y轴的范围
y1_max <- max(combined_data$Freq[combined_data$Type %in% c("DEL", "INS")])
y2_max <- max(combined_data$Freq[combined_data$Type %in% c("INV", "DUP")])

# 绘制第一个图表（DEL和INS）
p1 <- ggplot(combined_data, aes(x = Type, y = Freq, fill = file)) +
    geom_bar(data = combined_data[combined_data$Type %in% c("DEL", "INS"), ], stat = "identity", position = "dodge", width = 0.7) +
    scale_fill_manual(values = colors) +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(size = 16, family = "Times New Roman"),
          axis.text.y = element_text(size = 16, family = "Times New Roman")) +
    scale_y_continuous(name = "", limits = c(0, y1_max)) +  # 移除Y轴标题
    
    # 设置图例
    guides(fill = "none") +  # 设置图例标题
    
    # 移除X轴标题
    labs(x = "")  # 移除X轴标题

# 绘制第二个图表（INV和DUP）
p2 <- ggplot(combined_data, aes(x = Type, y = Freq, fill = file)) +
    geom_bar(data = combined_data[combined_data$Type %in% c("INV", "DUP"), ], stat = "identity", position = "dodge", width = 0.7) +
    scale_fill_manual(values = colors) +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(size = 16, family = "Times New Roman"),
          axis.text.y = element_text(size = 16, family = "Times New Roman")) +
    scale_y_continuous(name = "", limits = c(0, y2_max)) +  # 移除Y轴标题
    guides(fill = guide_legend(title = "")) +
    labs(x = "")  # 移除X轴标题

# 将两个图表左右拼接
combined_plot <- p1 + p2 + plot_layout(guides = 'collect')

# 去除标题和其他可能会导致重叠的元素
combined_plot <- combined_plot + theme(plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank())

# 设置字体
combined_plot <- combined_plot + theme(text = element_text(family = "Times New Roman"))

# 输出图形
ggsave("combined_plot.tiff", plot = combined_plot, width = 12, height = 6,  dpi = 600,bg = "white")
