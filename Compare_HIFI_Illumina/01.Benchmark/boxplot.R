
####boxplot, show length of each SV types
H <- read.delim("Hockett.summary", header = FALSE)
R <- read.delim("RGT.summary", header = FALSE)
O <- read.delim("OUN333.summary", header = FALSE)
V <- read.delim("Vlamingh.summary", header = FALSE)
I <- read.delim("Igri.summary", header = FALSE)

# 添加文件标识
H$file <- "Hockett"
R$file <- "RGT"
O$file <- "OUN333"
V$file <- "Vlamingh"
I$file <- "Igri"

# 合并数据
data <- rbind(H, R, O, V, I)

# 取绝对值
data$V3 <- abs(data$V3)

# 确保类型的顺序为DEL, INS, INV, DUP
data$V4 <- factor(data$V4, levels = c("DEL", "INS", "INV", "DUP"))

# 加载绘图所需的库
library(ggplot2)
theme_set(theme_minimal(base_family = "Times New Roman"))

# 创建 ggplot 图
p<-ggplot(data, aes(x = file, y = V3, fill = file)) +
    geom_boxplot() + 
    facet_wrap(~ V4, scales = "free", labeller = labeller(V4 = label_parsed)) +  
    labs(x = "", y = "Length (bp)", title = "") +
    theme(axis.text.y = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 18, face = "bold"),  # 设置子标题的字号和加粗
          text = element_text(size = 16)) +  # 设置其他文本的字号为16
    scale_fill_brewer(palette = "Set3") +
    guides(fill = "none")  # 移除图例
ggsave("boxplot_comparison.tiff", plot = p, width = 10, height = 6, dpi = 600,bg = "white")
