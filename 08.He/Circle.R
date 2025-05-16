# 安装 circlize 包（如未安装）
# install.packages("circlize")

library(circlize)

# ======== 输入路径 ========
fai_file <- "NLL.plot.fai"      # FAI 文件
atac_dir <- "./"                # ATAC 文件所在目录

# ======== 读取染色体长度信息 ========
fai <- read.table(fai_file, header = FALSE, stringsAsFactors = FALSE)
colnames(fai) <- c("chr", "length")
chr_lengths <- fai[, c("chr", "length")]

# ======== 读取所有 ATAC 文件 ========
atac_files <- list.files(atac_dir, pattern = "\\.bed$", full.names = TRUE)
atac_list <- lapply(atac_files, function(file) {
  df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("chr", "start", "end", "id", "score")
  df <- df[, c("chr", "start", "end", "score")]  # 保留4列
  return(df)
})
names(atac_list) <- basename(atac_files)

# ======== 统一 y 轴最大值 ========
global_ymax <- max(sapply(atac_list, function(df) max(df$score, na.rm = TRUE)), na.rm = TRUE)
global_ymax <- global_ymax * 1.05  # 增加一点缓冲

# ======== 绘图开始 ========
circos.clear()
circos.par(start.degree = 90, gap.after = c(rep(1, nrow(chr_lengths) - 1), 10))

# 初始化染色体位置
circos.initialize(factors = chr_lengths$chr,
                  xlim = cbind(rep(0, nrow(chr_lengths)), chr_lengths$length))

# 染色体名称 track
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  chr <- CELL_META$sector.index
  xlim <- CELL_META$xlim
  circos.text(mean(xlim), 0.5, chr, cex = 0.6, facing = "inside", niceFacing = TRUE)
}, bg.border = NA, track.height = 0.05)

# ======== 为每个 ATAC 文件画线图轨道 ========
colors <- rainbow(length(atac_list))  # 颜色设置

for (i in seq_along(atac_list)) {
  data <- atac_list[[i]]
  color <- colors[i]

  circos.genomicTrack(data,
                      ylim = c(0, global_ymax),
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region, value,
                                            col = color, lwd = 0.8, ...)
                      },
                      track.height = 0.10, bg.border = NA)
}

# ======== 添加图例 ========
legend("right", legend = names(atac_list), fill = colors, cex = 0.6, bty = "n")

# ======== 结束绘图 ========
circos.clear()
