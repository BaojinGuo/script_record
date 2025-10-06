library(ggplot2)
library(dplyr)
library(stringr)

# -----------------------------
# 1. 获取所有 bed 文件
# -----------------------------
files <- list.files(pattern = "*.bed")

all_data <- lapply(files, function(f) {
  df <- read.table(f, header = FALSE, sep = "\t")
  colnames(df) <- c("chr", "start", "end", "depth")
  
  # 样本名取文件名第一个点之前
  sample_name <- str_split(basename(f), "\\.")[[1]][1]
  df$sample <- sample_name
  
  # -----------------------------
  # 2. 去除异常值 (>3*mean)
  # -----------------------------
  m <- mean(df$depth)
  df <- df[df$depth <= 3*m, ]
  
  # -----------------------------
  # 3. 0–1 归一化 (每个文件单独)
  # -----------------------------
  df$depth_norm <- (df$depth - min(df$depth)) / (max(df$depth) - min(df$depth))
  
  # 取区间中点作为 X 轴位置
  df$pos <- (df$start + df$end)/2
  df$pos_Mb <- df$pos / 1e6  # 转换为 Mb
  
  return(df)
})

# -----------------------------
# 4. 合并所有样本
# -----------------------------
plot_data <- do.call(rbind, all_data)

# -----------------------------
# 5. 绘制全基因组 heatmap (不显示数字)
# -----------------------------
p_heatmap_all <- ggplot(plot_data, aes(x=pos_Mb, y=sample, fill=depth_norm)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0.5) +
  theme_bw() +
  labs(x="Genomic position (Mb)", 
       y="Sample", 
       fill="Normalized depth (0-1)",
       title="Depth Heatmap (full genome)")

ggsave("depth_heatmap_full_genome_BWR.png", p_heatmap_all, width=12, height=4, dpi=300)

# -----------------------------
# 6. 绘制 zoom in heatmap (445Mb+) 并显示格子数字
# -----------------------------
zoom_data <- plot_data %>%
  filter(pos >= 445000000)

p_heatmap_zoom <- ggplot(zoom_data, aes(x=pos_Mb, y=sample, fill=depth_norm)) +
  geom_tile() +
  geom_text(aes(label=round(depth_norm, 2)), size=2) +  # 只在 zoom in 显示数字
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0.5) +
  theme_bw() +
  labs(x="Genomic position (Mb)", 
       y="Sample", 
       fill="Normalized depth (0-1)",
       title="Depth Heatmap (zoom in 445Mb+)") +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggsave("depth_heatmap_zoom_445Mb_BWR_numbers.png", p_heatmap_zoom, width=12, height=4, dpi=300)

# -----------------------------
# 7. 显示图形
# -----------------------------
p_heatmap_all
p_heatmap_zoom

write.table(zoom_data, "zoom_445Mb_depth_data.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)
