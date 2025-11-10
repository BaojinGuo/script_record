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

ggsave("depth_heatmap_full_genome_BWR.pdf", p_heatmap_all, width=12, height=24, dpi=300)

# -----------------------------
# 6. 绘制 zoom in heatmap (445Mb+) 并显示格子数字
# -----------------------------
zoom_data <- plot_data %>%
  filter(pos >= 455000000)

p_heatmap_zoom <- ggplot(zoom_data, aes(x=pos_Mb, y=sample, fill=depth_norm)) +
  geom_tile() +
  geom_text(aes(label=round(depth_norm, 2)), size=2) +  # 只在 zoom in 显示数字
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0.5) +
  theme_bw() +
  labs(x="Genomic position (Mb)", 
       y="Sample", 
       fill="Normalized depth (0-1)",
       title="Depth Heatmap (zoom in 455Mb+)") +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggsave("depth_heatmap_zoom_455Mb_BWR_numbers.pdf", p_heatmap_zoom, width=12, height=24, dpi=300)

# -----------------------------
# 7. 显示图形
# -----------------------------
p_heatmap_all
p_heatmap_zoom

write.table(zoom_data, "zoom_455Mb_depth_data.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)

############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
#Add group info
# -----------------------------
# 加载所需包
# -----------------------------
library(ggplot2)
library(dplyr)
library(stringr)

# -----------------------------
# 1. 读取所有 bed 文件
# -----------------------------
files <- list.files(pattern = "*.bed")

all_data <- lapply(files, function(f) {
    df <- read.table(f, header = FALSE, sep = "\t")
    colnames(df) <- c("chr", "start", "end", "depth")
    
    # 样本名取文件名第一个点之前
    sample_name <- str_split(basename(f), "\\.")[[1]][1]
    df$sample <- sample_name
    
    # 去除异常值 (>3*mean)
    m <- mean(df$depth)
    df <- df[df$depth <= 3.5*m, ]
    
    # 0–1归一化
    df$depth_norm <- (df$depth - min(df$depth)) / (max(df$depth) - min(df$depth))
    
    # 区间中点作为位置
    df$pos <- (df$start + df$end)/2
    df$pos_Mb <- df$pos / 1e6
    
    return(df)
})

# -----------------------------
# 2. 合并所有样本
# -----------------------------
plot_data <- do.call(rbind, all_data)

# -----------------------------
# 3. 加入Hull类型信息 (Hull.txt)
# -----------------------------
trait_info <- read.table("Hull.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 只保留存在的样本
trait_info <- trait_info %>% filter(Sample %in% plot_data$sample)

# 将1/0转为Hulless/Covered
trait_info <- trait_info %>%
    mutate(Hull = ifelse(Hull == 1, "Hulless", "Covered"),
           Hull = factor(Hull, levels = c("Hulless", "Covered")))

# 合并
plot_data <- plot_data %>%
    left_join(trait_info, by = c("sample" = "Sample")) %>%
    arrange(desc(Hull), sample)

# 确保y轴顺序
plot_data$sample <- factor(plot_data$sample, levels = unique(plot_data$sample))

# -----------------------------
# 4. 计算分割线位置
# -----------------------------
group_boundaries <- plot_data %>%
    group_by(Hull) %>%
    summarise(y = max(as.numeric(factor(sample)))) %>%
    arrange(desc(Hull))

boundary_lines <- head(group_boundaries$y, -1) + 0.5

# -----------------------------
# 5. 全基因组热图
# -----------------------------
p_all <- ggplot(plot_data, aes(x = pos_Mb, y = sample, fill = depth_norm)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
    geom_hline(yintercept = boundary_lines, color = "black", size = 0.8) +
    annotate(
        "text",
        x = min(plot_data$pos_Mb) - 20,
        y = group_boundaries$y - diff(c(0, group_boundaries$y))/2,
        label = group_boundaries$Hull,
        angle = 90, vjust = 0.5, size = 5
    ) +
    theme_bw() +
    labs(
        x = "Genomic position (Mb)",
        y = "Sample",
        fill = "Normalized depth (0-1)",
        title = "Depth Heatmap (Chr4D)"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(10, 20, 10, 70)
    )

ggsave("depth_heatmap_grouped_full.pdf", p_all, width = 12, height = 24, dpi = 300)

# -----------------------------
# 6. Zoom in (455Mb+)
# -----------------------------
zoom_data <- plot_data %>% filter(pos >= 455000000)

# 再计算zoom分割线
zoom_group_boundaries <- zoom_data %>%
    group_by(Hull) %>%
    summarise(y = max(as.numeric(factor(sample)))) %>%
    arrange(desc(Hull))
zoom_boundary_lines <- head(zoom_group_boundaries$y, -1) + 0.5

p_zoom <- ggplot(zoom_data, aes(x = pos_Mb, y = sample, fill = depth_norm)) +
    geom_tile() +
    geom_text(aes(label = round(depth_norm, 2)), size = 2) +  # 显示数字
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
    geom_hline(yintercept = zoom_boundary_lines, color = "black", size = 0.8) +
    annotate(
        "text",
        x = min(zoom_data$pos_Mb) - 0.1,
        y = zoom_group_boundaries$y - diff(c(0, zoom_group_boundaries$y))/2,
        label = zoom_group_boundaries$Hull,
        angle = 90, vjust = 0.5, size = 5
    ) +
    theme_bw() +
    labs(
        x = "Genomic position (Mb)",
        y = "Sample",
        fill = "Normalized depth (0-1)",
        title = "Depth Heatmap (Chr4D:455Mb+)"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(10, 20, 10, 70)
    )

ggsave("depth_heatmap_grouped_zoom_455Mb.pdf", p_zoom, width = 12, height = 24, dpi = 300)

# -----------------------------
# 7. 输出数据表
# -----------------------------
write.table(zoom_data, "zoom_455Mb_grouped_data.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# 8. 显示图形
# -----------------------------
p_all
p_zoom






###boxplot

# -----------------------------
# 加载包
# -----------------------------
library(ggplot2)
library(dplyr)
library(ggpubr)

# -----------------------------
# 1. 选定区间：455.4–456 Mb
# -----------------------------
focus_data <- plot_data %>%
  filter(pos_Mb >= 455.4 & pos_Mb <= 456) %>%
  filter(!is.na(Hull))   # 去掉没有分组信息的行

# -----------------------------
# 2. 按每100 kb分箱
# -----------------------------
focus_data <- focus_data %>%
  mutate(bin = cut(pos_Mb, breaks = seq(455.4, 456, by = 0.1), include.lowest = TRUE))

# -----------------------------
# 3. 每个bin内比较两组差异并计算p值（Wilcoxon）
# -----------------------------
pvals <- focus_data %>%
  group_by(bin) %>%
  summarise(
    p = ifelse(length(unique(Hull)) == 2,
               wilcox.test(depth_norm ~ Hull)$p.value,
               NA_real_)
  ) %>%
  mutate(sig_label = case_when(
    is.na(p) ~ "ns",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

# -----------------------------
# 4. 绘制箱线图 + jitter + p值标注
# -----------------------------
p_box <- ggplot(focus_data, aes(x = bin, y = depth_norm, fill = Hull)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6, position = position_dodge(width = 0.7)) +
  geom_jitter(
    aes(color = Hull),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
    size = 1.2,
    alpha = 0.7
  ) +
  scale_fill_manual(values = c("Hulless" = "#E64B35", "Covered" = "#4DBBD5")) +
  scale_color_manual(values = c("Hulless" = "#E64B35", "Covered" = "#4DBBD5")) +
  geom_text(
    data = pvals,
    inherit.aes = FALSE,
    aes(x = bin, y = 1.05, label = sig_label),
    color = "black", size = 4, vjust = 0
  ) +
  theme_bw() +
  labs(
    x = "Genomic position (Mb, 100kb bins)",
    y = "Normalized depth (0–1)",
    fill = "Hull type",
    color = "Hull type",
    title = "Depth comparison between Hulless and Covered (455.4–456 Mb)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  ylim(0, 1.1)  # 留出空间放 p 值标注

# -----------------------------
# 5. 保存图形和 p 值表格
# -----------------------------
ggsave("depth_boxplot_455.4_456Mb_by_Hull_jitter.pdf", p_box, width = 10, height = 6, dpi = 300)

write.table(pvals, "pvalues_455.4_456Mb.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

p_box

###################################################################################################################
####################################三组分组

library(ggplot2)
library(dplyr)
library(stringr)

# -----------------------------
# 1. 读取所有 bed 文件
# -----------------------------
files <- list.files(pattern = "*.bed")

all_data <- lapply(files, function(f) {
    df <- read.table(f, header = FALSE, sep = "\t")
    colnames(df) <- c("chr", "start", "end", "depth")
    
    # 样本名取文件名第一个点之前
    sample_name <- str_split(basename(f), "\\.")[[1]][1]
    df$sample <- sample_name
    
    # 去除异常值 (>3.5 * mean)
    m <- mean(df$depth)
    df <- df[df$depth <= 3.5 * m, ]
    
    # 归一化 (0–1)
    df$depth_norm <- (df$depth - min(df$depth)) / (max(df$depth) - min(df$depth))
    
    # 区间中点作为位置
    df$pos <- (df$start + df$end) / 2
    df$pos_Mb <- df$pos / 1e6
    
    return(df)
})

# -----------------------------
# 2. 合并所有样本
# -----------------------------
plot_data <- do.call(rbind, all_data)

# -----------------------------
# 3. 读取并整合 Hull 信息
# -----------------------------
trait_info <- read.table("Hull.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 保留存在于数据中的样本
trait_info <- trait_info %>% filter(Sample %in% plot_data$sample)

# 将数字编码转为分组名称
# 假设 Hull 列为：2=Wild, 1=Hulless, 0=Hulled
trait_info <- trait_info %>%
    mutate(Hull = case_when(
        Hull == 2 ~ "Wild",
        Hull == 1 ~ "Hulless",
        Hull == 0 ~ "Hulled"
    ),
    Hull = factor(Hull, levels = c("Wild", "Hulled", "Hulless")))  # Wild最上，Hulless最下

# 合并 Hull 信息
plot_data <- plot_data %>%
    left_join(trait_info, by = c("sample" = "Sample"))

if(!"Hull" %in% colnames(plot_data)) stop("❌ Hull column missing after join!")

# -----------------------------
# 4. 排序 & Y轴数值化
# -----------------------------
plot_data <- plot_data %>%
    arrange(desc(Hull), sample) %>%
    mutate(y_plot = as.numeric(factor(sample, levels = unique(sample))),
           sample = factor(sample, levels = unique(sample)))

# -----------------------------
# 5. 计算分组边界和标签位置
# -----------------------------
group_boundaries <- plot_data %>%
    group_by(Hull) %>%
    summarise(ymin = min(y_plot, na.rm = TRUE),
              ymax = max(y_plot, na.rm = TRUE)) %>%
    arrange(desc(Hull)) %>%
    mutate(y_label = (ymin + ymax) / 2)

boundary_lines <- group_boundaries$ymax[-nrow(group_boundaries)] + 0.5

# -----------------------------
# 6. 全基因组热图
# -----------------------------
p_all <- ggplot(plot_data, aes(x = pos_Mb, y = sample, fill = depth_norm)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
    geom_hline(yintercept = boundary_lines, color = "black", size = 0.8) +
    annotate(
        "text",
        x = min(plot_data$pos_Mb) - 20,
        y = group_boundaries$y_label,
        label = group_boundaries$Hull,
        angle = 90,
        vjust = 0.5,
        size = 5
    ) +
    theme_bw() +
    labs(
        x = "Genomic position (Mb)",
        y = "Sample",
        fill = "Normalized depth (0–1)",
        title = "Depth Heatmap (Chr4D)"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.margin = margin(10, 20, 10, 70),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.spacing = unit(0, "lines")
    ) +
    coord_cartesian(ylim = c(min(plot_data$y_plot) - 0.5, max(plot_data$y_plot) + 0.5))  # 去除上下留白

ggsave("depth_heatmap_grouped_full.pdf", p_all, width = 12, height = 56, dpi = 300,limitsize = FALSE)

# -----------------------------
# 7. Zoom in (455Mb+)
# -----------------------------
zoom_data <- plot_data %>% filter(pos >= 455000000)

zoom_group_boundaries <- zoom_data %>%
    group_by(Hull) %>%
    summarise(ymin = min(y_plot, na.rm = TRUE),
              ymax = max(y_plot, na.rm = TRUE)) %>%
    arrange(desc(Hull)) %>%
    mutate(y_label = (ymin + ymax) / 2)

zoom_boundary_lines <- zoom_group_boundaries$ymax[-nrow(zoom_group_boundaries)] + 0.5

p_zoom <- ggplot(zoom_data, aes(x = pos_Mb, y = sample, fill = depth_norm)) +
    geom_tile() +
    geom_text(aes(label = round(depth_norm, 2)), size = 2) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
    geom_hline(yintercept = zoom_boundary_lines, color = "black", size = 0.8) +
    annotate(
        "text",
        x = min(zoom_data$pos_Mb) - 0.1,
        y = zoom_group_boundaries$y_label,
        label = zoom_group_boundaries$Hull,
        angle = 90,
        vjust = 0.5,
        size = 5
    ) +
    theme_bw() +
    labs(
        x = "Genomic position (Mb)",
        y = "Sample",
        fill = "Normalized depth (0–1)",
        title = "Depth Heatmap (Chr4D:455Mb+)"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.margin = margin(10, 20, 10, 70),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.spacing = unit(0, "lines")
    ) +
    coord_cartesian(ylim = c(min(zoom_data$y_plot) - 0.5, max(zoom_data$y_plot) + 0.5))  # 紧贴上下边界

ggsave("depth_heatmap_grouped_zoom_455Mb.pdf", p_zoom, width = 12, height =56, dpi = 300,limitsize = FALSE)

# -----------------------------
# 8. 输出 zoom 数据
# -----------------------------
write.table(zoom_data, "zoom_455Mb_grouped_data.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
