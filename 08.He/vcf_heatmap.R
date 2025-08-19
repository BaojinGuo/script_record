#######单个运行
library(vcfR)
library(pheatmap)

vcf<-read.vcfR("Oat_Pan.DEL.2D_140744075_144750065.vcf.gz",verbose = FALSE)

gt <- extract.gt(vcf, element = "GT", IDtoRowNames = FALSE)

row.names(gt) <- paste(vcf@fix[,1], vcf@fix[,2], sep = "_")

gt_num <- apply(gt, 2, function(x) {
    sapply(x, function(y) {
        if (is.na(y)) return(NA)
        if (y == "0/0") return(0)
        if (y %in% c("0/1", "1/0")) return(1)
        if (y == "1/1") return(2)
        return(NA)
    })
})


gt_num <- as.matrix(gt_num)

focus_samples <- c("PI182478", "Gehl")
other_samples <- setdiff(colnames(gt_num), focus_samples)
gt_num <- gt_num[, c(focus_samples, other_samples)]

order_idx <- order(-gt_num[,"PI182478"], -gt_num[,"Gehl"], na.last = TRUE)
gt_num_sorted <- gt_num[order_idx, ]

col_palette <- c("white", "lightblue", "blue")
gt_num_filled <- gt_num_sorted
gt_num_filled[is.na(gt_num_filled)] <- 0

pheatmap(gt_num_filled,
         cluster_rows = FALSE,      # 可以选择聚类行
         cluster_cols = FALSE,     # 样本列保持顺序
         color = col_palette,
         show_rownames = FALSE,
         fontsize_col = 10,
         main = "Oat Variants Heatmap")


#####运行文件夹下所有文件
library(vcfR)
library(pheatmap)

# 设置文件夹路径
vcf_dir <- "."  # 当前目录，可改为你的目录
out_dir <- "heatmap_output"
dir.create(out_dir, showWarnings = FALSE)

# 获取所有 vcf.gz 文件
vcf_files <- list.files(vcf_dir, pattern = "\\.vcf\\.gz$", full.names = TRUE)

# 循环处理
for (vcf_file in vcf_files) {
  
  cat("Processing:", vcf_file, "\n")
  
  tryCatch({
    # 读取 VCF
    vcf <- read.vcfR(vcf_file, verbose = FALSE)
    
    # 提取 GT
    gt <- extract.gt(vcf, element = "GT", IDtoRowNames = FALSE)
    
    # 设置行名为 chr_pos
    row.names(gt) <- paste(vcf@fix[,1], vcf@fix[,2], sep = "_")
    
    # 转换为数字矩阵
    gt_num <- apply(gt, 2, function(x) {
      sapply(x, function(y) {
        if (is.na(y)) return(NA)
        if (y == "0/0") return(0)
        if (y %in% c("0/1", "1/0")) return(1)
        if (y == "1/1") return(2)
        return(NA)
      })
    })
    
    gt_num <- as.matrix(gt_num)
    
    # 排列样本顺序
    focus_samples <- c("PI182478", "Gehl")
    other_samples <- setdiff(colnames(gt_num), focus_samples)
    gt_num <- gt_num[, c(focus_samples, other_samples)]
    
    # 按 PI182478 和 Gehl 排序行
    order_idx <- order(-gt_num[,"PI182478"], -gt_num[,"Gehl"], na.last = TRUE)
    gt_num_sorted <- gt_num[order_idx, ]
    
    # NA 填 0
    gt_num_filled <- gt_num_sorted
    gt_num_filled[is.na(gt_num_filled)] <- 0
    
    # 保存 CSV
    out_csv <- file.path(out_dir, paste0(basename(vcf_file), "_gt_matrix.csv"))
    write.csv(gt_num_filled, out_csv, row.names = TRUE)
    
    # 绘制热图
    out_png <- file.path(out_dir, paste0(basename(vcf_file), "_heatmap.png"))
    col_palette <- c("white", "lightblue", "blue")
    png(out_png, width = 1200, height = 800)
    pheatmap(gt_num_filled,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             color = col_palette,
             show_rownames = FALSE,
             fontsize_col = 10,
             main = paste0("Heatmap: ", basename(vcf_file)))
    dev.off()
    
    cat("Finished:", vcf_file, "\n")
    
  }, error = function(e){
    cat("Error in file:", vcf_file, " -> skipped\n")
  })
}








