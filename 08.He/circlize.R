rm(list = ls())
graphics.off()

# 安装并加载必要包
if (!require("circlize")) install.packages("circlize", dependencies = TRUE)
if (!require("dplyr")) install.packages("dplyr", dependencies = TRUE)
if (!require("readr")) install.packages("readr", dependencies = TRUE)
library(circlize)
library(dplyr)
library(readr)
#chr_size_file <- "chr_size.csv"      # 染色体大小文件
regions_file <- "regions.csv"        # 区域数据文件
# ============================================
# 2. 读取和预处理数据
# ============================================

# 读取染色体大小数据
# 定义输出文件
output_pdf <- "barley_genome_circos_final.pdf"
output_png <- "barley_genome_circos_final.png"
chr_info <- data.frame(
    chr = c("1H", "2H", "3H", "4H", "5H", "6H", "7H"),
    length = c(516.51, 665.59, 621.52, 610.33, 588.22, 561.79, 632.54)
)
chr_order <- paste0(1:7, "H")
chr_info$chr <- factor(chr_info$chr, levels = chr_order)
chr_info <- chr_info[order(chr_info$chr), ]

# 读取区域数据
cat("读取区域数据...\n")
regions <- read_csv(regions_file, col_types = cols(
  chr = col_character(),
  start = col_double(),
  end = col_double(),
  type = col_character(),
  name = col_character()
))
# 确保染色体顺序和数据类型正确
regions$chr <- factor(regions$chr, levels = chr_order)
regions <- regions[order(regions$chr, regions$start), ]

# 检查数据
cat("\n区域数据类型统计:\n")
type_counts <- table(regions$type)
print(type_counts)

# ============================================
# 3. 特别提取1H染色体数据
# ============================================

cat("\n1H染色体数据统计:\n")
regions_1H <- regions %>% filter(chr == "1H")
type_counts_1H <- table(regions_1H$type)
print(type_counts_1H)

# 显示1H染色体详细数据
cat("\n1H染色体区域详情:\n")
for (type in c("Selected region", "Known genes", "QTL")) {
  type_data <- regions_1H %>% filter(type == !!type)
  if (nrow(type_data) > 0) {
    cat(sprintf("\n%s (%d个):\n", type, nrow(type_data)))
    print(type_data[, c("start", "end", "name")])
  }
}

# ============================================
# 4. 处理小区域：扩大显示范围
# ============================================

# 复制原始数据用于处理
regions_enhanced <- regions

# 扩大QTL和Known genes的显示范围
min_display_size <- 2.0  # 最小显示尺寸（Mb）

for (i in 1:nrow(regions_enhanced)) {
  region_width <- regions_enhanced$end[i] - regions_enhanced$start[i]
  
  # 对QTL和Known genes进行扩展
  if (regions_enhanced$type[i] %in% c("QTL", "Known genes")) {
    if (region_width < min_display_size) {
      center <- (regions_enhanced$start[i] + regions_enhanced$end[i]) / 2
      # 扩展为最小显示尺寸
      regions_enhanced$start[i] <- center - min_display_size / 2
      regions_enhanced$end[i] <- center + min_display_size / 2
      
      # 确保不超出染色体边界
      chr_max <- chr_info$length[chr_info$chr == regions_enhanced$chr[i]]
      if (regions_enhanced$start[i] < 0) regions_enhanced$start[i] <- 0
      if (regions_enhanced$end[i] > chr_max) regions_enhanced$end[i] <- chr_max
    }
  }
}
# ============================================
# 5. 设置图形参数和配色
# ============================================

# 配色方案
type_colors <- c(
  "Selected region" = "#E41A1C",    # 红色
  "Known genes" = "#377EB8",        # 蓝色  
  "QTL" = "#4DAF4A"                 # 绿色
)

# 设置图形参数
circos.clear()

circos.par(
  start.degree = 90,                    # 从12点位置开始
  gap.degree = c(rep(2, 6), 8),         # 染色体间空隙
  track.margin = c(0.005, 0.005),       # 轨道边距
  cell.padding = c(0, 0, 0, 0),         # 避免轨道高度错误
  points.overflow.warning = FALSE
)

# ============================================
# 6. 初始化Circos图
# ============================================

cat("\n初始化Circos图...\n")

circos.initialize(
  sectors = chr_info$chr,
  xlim = cbind(rep(0, nrow(chr_info)), chr_info$length)
)

# ============================================
# 7. 绘制染色体轨道
# ============================================

cat("绘制染色体轨道...\n")

# 轨道1：染色体名称
circos.track(
  ylim = c(0, 1),
  bg.border = "gray30",
  bg.col = "white",
  track.height = 0.12,
  panel.fun = function(x, y) {
    # 添加染色体名称
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[2] - 0.3,
      CELL_META$sector.index,
      cex = 1.5,
      col = "black",
      facing = "bending.inside",
      niceFacing = TRUE,
      font = 2
    )
  }
)
# ============================================
# 8. 绘制刻度轨道
# ============================================

cat("绘制刻度轨道...\n")

circos.track(
  ylim = c(0, 1),
  bg.border = NA,
  bg.col = NA,
  track.height = 0.05,
  panel.fun = function(x, y) {
    # 每100Mb一个主刻度
    major_at <- seq(0, CELL_META$xlim[2], by = 100)
    if (max(major_at) > CELL_META$xlim[2]) {
      major_at <- major_at[major_at <= CELL_META$xlim[2]]
    }
    
    # 绘制刻度
    circos.axis(
      h = "top",
      major.at = major_at,
      labels = ifelse(major_at %% 200 == 0, major_at, ""),
      labels.cex = 0.7,
      minor.ticks = 4,
      labels.facing = "clockwise",
      col = "gray30",
      labels.col = "black"
    )
  }
)

# ============================================
# 9. 绘制区域数据轨道（不显示类型名称）
# ============================================

cat("绘制区域数据轨道...\n")

# 定义轨道高度
track_heights <- c(
  "Selected region" = 0.12,
  "QTL" = 0.15,
  "Known genes" = 0.15
)

# 为每种类型创建单独的轨道
for (type_idx in seq_along(names(type_colors))) {
  type <- names(type_colors)[type_idx]
  sub_data <- regions_enhanced %>% filter(type == !!type)
  
  if (nrow(sub_data) > 0) {
    cat(sprintf("  - 绘制 %s: %d 个区域\n", type, nrow(sub_data)))
    
    # 使用预定义的轨道高度
    track_height <- track_heights[type]
    
    # 创建轨道（不添加类型名称标签）
    circos.track(
      ylim = c(0, 1),
      bg.border = "gray30",
      bg.col = "white",
      track.height = track_height,
      panel.fun = function(x, y) {
        current_chr <- CELL_META$sector.index
        chr_data <- sub_data %>% filter(chr == current_chr)
        
        if (nrow(chr_data) > 0) {
          # 绘制每个区域
          for (i in 1:nrow(chr_data)) {
            # 计算区域位置
            x_left <- chr_data$start[i]
            x_right <- chr_data$end[i]
            
            # 确保区域在染色体范围内
            if (x_left >= 0 && x_right <= CELL_META$xlim[2] && x_left < x_right) {
              # 根据类型设置不同的y范围
              if (type == "Selected region") {
                y_bottom <- 0.25
                y_top <- 0.75
              } else {
                y_bottom <- 0.15
                y_top <- 0.85
              }
              
              # 绘制矩形
              circos.rect(
                xleft = x_left,
                xright = x_right,
                ybottom = y_bottom,
                ytop = y_top,
                col = type_colors[type],
                border = "gray30",
                lwd = 0.8
              )
              
              # 只显示Known genes的名称（字号变大）
              if (type == "Known genes") {
                label_name <- chr_data$name[i]
                
                # 简化长名称
                if (nchar(label_name) > 15) {
                  label_name <- paste0(substr(label_name, 1, 12), "...")
                }
                
                # 只在区域足够宽时添加标签
                if ((x_right - x_left) > 0.8) {  # 降低阈值
                  circos.text(
                    x = (x_left + x_right) / 2,
                    y = 0.5,
                    labels = label_name,
                    cex = 0.7,  # 增大字号
                    col = "black",  # 黑色文字
                    facing = "clockwise",
                    niceFacing = TRUE,
                    font = 2
                  )
                }
              }
            }
          }
        }
      }
    )
    
    # 注释掉类型名称标签部分（不显示）
    # 这里完全移除类型名称标签，不占用1H空间
  }
}

# ============================================
# 10. 添加图例（在图形外部）
# ============================================

cat("添加图例...\n")

# 在主图外部添加图例
legend(
  x = -1.2,
  y = -1.1,
  legend = c(
    paste("Selected regions:", type_counts["Selected region"]),
    paste("QTL regions:", type_counts["QTL"]),
    paste("Known genes:", type_counts["Known genes"])
  ),
  fill = type_colors,
  border = "gray30",
  bty = "n",
  cex = 0.9,
  xpd = NA,
  text.col = "black",
  title = "Region Types",
  title.col = "black"
)
# 保存为PDF（白色背景）
pdf(output_pdf, width = 12, height = 12, bg = "white")
# 重新运行绘图代码（从第5步开始）
# 这里需要重新执行绘图代码
dev.off()
cat(sprintf("Circos图PDF保存为: %s\n", output_pdf))

# 保存为PNG（白色背景，高分辨率）
png(output_png, width = 2000, height = 2000, res = 300, bg = "white")
# 重新运行绘图代码（从第5步开始）
dev.off()
#


