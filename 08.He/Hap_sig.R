library(tidyverse)
library(ggpubr)

# 创建输出文件夹
if(!dir.exists("violin_plots")) dir.create("violin_plots")

# 读取数据
hap <- read.csv("HAPmatrix-two/HvAGLG1_HvCBF8A.hap.matrix.csv", check.names = FALSE)
trait <- read.table("Trait1-GY.txt", header = TRUE, check.names = FALSE)

# 合并数据并创建组合
data <- merge(hap, trait, by.x = "Accession", by.y = "Taxa")
data$hap_combo <- paste(data$HvAGLG1, data$HvCBF8A, sep = "_")

# 获取性状列
trait_cols <- grep("^T[0-9]{3}$", colnames(data), value = TRUE)


for(trait_name in trait_cols) {
    # 过滤数据
    plot_data <- data[!is.na(data[[trait_name]]), ]
    
    # 计算每组样本数
    group_counts <- table(plot_data$hap_combo)
    valid_groups <- names(group_counts)[group_counts >= 5]
    
    if(length(valid_groups) >= 2) {
        plot_data <- plot_data[plot_data$hap_combo %in% valid_groups, ]
        
        # 设置因子顺序与图例一致
        plot_data$hap_combo <- factor(plot_data$hap_combo, levels = valid_groups)
        
        # 计算Y轴范围
        y_max <- max(plot_data[[trait_name]], na.rm = TRUE)
        y_min <- min(plot_data[[trait_name]], na.rm = TRUE)
        
        # 创建标题
        title_text <- paste("HvAGLG1_HvCBF8A -", trait_name)
        
        # 筛选显著比较
        comparisons <- combn(valid_groups, 2, simplify = FALSE)
        sig_comparisons <- list()
        
        for(comp in comparisons) {
            group1 <- plot_data[plot_data$hap_combo == comp[1], trait_name]
            group2 <- plot_data[plot_data$hap_combo == comp[2], trait_name]
            
            if(t.test(group1, group2)$p.value < 0.05) {
                sig_comparisons <- c(sig_comparisons, list(comp))
            }
        }
        
        # 绘图
        p <- ggviolin(plot_data, 
                      x = "hap_combo", 
                      y = trait_name,
                      fill = "hap_combo",
                      palette = "npg",
                      add = "boxplot",
                      legend = "right") +
            labs(title = title_text, x = "") +
            guides(fill = guide_legend(override.aes = list(color = NA))) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                  axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.title = element_blank()) +
            scale_fill_discrete(labels = paste0(valid_groups, "(", group_counts[valid_groups], ")"))
        
        # 只添加显著比较，并调整位置
        if(length(sig_comparisons) > 0) {
            # 确定每个比较的垂直位置
            n_comparisons <- length(sig_comparisons)
            step_height <- (y_max - y_min) * 0.08
            
            p <- p + stat_compare_means(
                comparisons = sig_comparisons, 
                label = "p.format",
                hide.ns = TRUE,
                tip.length = 0.01,
                step.increase = 0.08,  # 增加比较线之间的间距
                vjust = 1.8  # 垂直调整
            ) +
                # 扩展Y轴上限
                scale_y_continuous(limits = c(y_min, y_max + step_height * n_comparisons * 1.2))
        }
        
        # 保存
        ggsave(paste0("violin_plots/HvAGLG1_HvCBF8A.", trait_name, ".violin.pdf"),
               p, width = 10, height = 7)  # 增加画布尺寸
        
        message(paste("已保存:", trait_name))
    } else {
        message(paste("跳过", trait_name))
    }
}






############循环################
library(tidyverse)
library(ggpubr)

# 主循环：对每个hap文件
for(hap_file in list.files("HAPmatrix-two", pattern = "\\.hap\\.matrix\\.csv$", full.names = TRUE)) {
  # 从文件名提取基因名称
  gene_name <- gsub("\\.hap\\.matrix\\.csv$", "", basename(hap_file))
  
  # 读取hap矩阵文件
  hap <- read.csv(hap_file, check.names = FALSE)
  
  # 内循环：对每个表型文件
  for(trait_file in c("Trait1-GY.txt", "Trait2-GP.txt", "Trait5-Heat2.txt")) {
    # 读取表型文件
    trait <- read.table(trait_file, header = TRUE, check.names = FALSE)
    
    # 从表型文件名提取性状类别
    trait_category <- gsub("\\.txt$", "", basename(trait_file))
    
    # 创建对应的输出文件夹
    output_dir <- file.path("violin_plots", gene_name, trait_category)
    if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # 合并数据
    data <- merge(hap, trait, by.x = "Accession", by.y = "Taxa")
    
    # 获取基因列名
    gene_cols <- colnames(hap)[-1]
    
    # 创建组合
    if(length(gene_cols) >= 2) {
      data$hap_combo <- paste(data[[gene_cols[1]]], data[[gene_cols[2]]], sep = "_")
      combo_name <- paste(gene_cols, collapse = "_")
    } else {
      data$hap_combo <- data[[gene_cols[1]]]
      combo_name <- gene_cols[1]
    }
    
    # 获取性状列
    trait_cols <- grep("^T[0-9]{3}$", colnames(data), value = TRUE)
    
    # 对每个性状绘图
    for(trait_name in trait_cols) {
      # 过滤数据
      plot_data <- data[!is.na(data[[trait_name]]), ]
      
      # 计算每组样本数
      group_counts <- table(plot_data$hap_combo)
      valid_groups <- names(group_counts)[group_counts >= 5]
      
      if(length(valid_groups) >= 2) {
        plot_data <- plot_data[plot_data$hap_combo %in% valid_groups, ]
        plot_data$hap_combo <- factor(plot_data$hap_combo, levels = valid_groups)
        
        y_max <- max(plot_data[[trait_name]], na.rm = TRUE)
        y_min <- min(plot_data[[trait_name]], na.rm = TRUE)
        
        # 筛选显著比较
        comparisons <- combn(valid_groups, 2, simplify = FALSE)
        sig_comparisons <- list()
        
        for(comp in comparisons) {
          group1 <- plot_data[plot_data$hap_combo == comp[1], trait_name]
          group2 <- plot_data[plot_data$hap_combo == comp[2], trait_name]
          
          # 添加错误处理：检查数据是否有变化
          if(length(unique(group1)) > 1 && length(unique(group2)) > 1) {
            test_result <- tryCatch({
              t.test(group1, group2)
            }, error = function(e) NULL)
            
            if(!is.null(test_result) && test_result$p.value < 0.05) {
              sig_comparisons <- c(sig_comparisons, list(comp))
            }
          }
        }
        
        # 绘图
        p <- ggviolin(plot_data, 
                      x = "hap_combo", 
                      y = trait_name,
                      fill = "hap_combo",
                      palette = "npg",
                      add = "boxplot",
                      legend = "right") +
          labs(title = paste(combo_name, "-", trait_name), x = "") +
          guides(fill = guide_legend(override.aes = list(color = NA))) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                axis.text.x = element_text(angle = 45, hjust = 1),
                legend.title = element_blank()) +
          scale_fill_discrete(labels = paste0(valid_groups, "(", group_counts[valid_groups], ")"))
        
        if(length(sig_comparisons) > 0) {
          n_comparisons <- length(sig_comparisons)
          step_height <- (y_max - y_min) * 0.08
          
          p <- p + stat_compare_means(
            comparisons = sig_comparisons, 
            label = "p.format",
            hide.ns = TRUE,
            tip.length = 0.01,
            step.increase = 0.08,
            vjust = 1.8
          ) +
            scale_y_continuous(limits = c(y_min, y_max + step_height * n_comparisons * 1.2))
        }
        
        # 保存到对应的文件夹
        output_path <- file.path(output_dir, paste0(trait_name, ".violin.pdf"))
        ggsave(output_path, p, width = 10, height = 7)
        
        message(paste("已保存:", gene_name, "/", trait_category, "/", trait_name))
      }
    }
  }
}

##################################################################
library(tidyverse)
library(ggpubr)

# 要处理的文件夹列表
hap_folders <- c("HAPmatrix-two", "HAPmatrix-THREE", "HAPmatrix-FOUR")
trait_files <- c("Trait1-GY.txt", "Trait2-GP.txt", "Trait5-Heat2.txt")

# 主循环：对每个hap文件夹
for(hap_folder in hap_folders) {
  # 获取当前文件夹下的所有hap矩阵文件
  hap_files <- list.files(hap_folder, pattern = "\\.hap\\.matrix\\.csv$", full.names = TRUE)
  
  # 对每个hap文件
  for(hap_file in hap_files) {
    # 从文件名提取基因名称
    gene_name <- gsub("\\.hap\\.matrix\\.csv$", "", basename(hap_file))
    
    # 读取hap矩阵文件
    hap <- read.csv(hap_file, check.names = FALSE)
    
    # 内循环：对每个表型文件
    for(trait_file in trait_files) {
      # 读取表型文件
      trait <- read.table(trait_file, header = TRUE, check.names = FALSE)
      
      # 从表型文件名提取性状类别
      trait_category <- gsub("\\.txt$", "", basename(trait_file))
      
      # 创建对应的输出文件夹
      output_dir <- file.path("violin_plots", hap_folder, gene_name, trait_category)
      if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      
      # 合并数据
      data <- merge(hap, trait, by.x = "Accession", by.y = "Taxa")
      
      # 获取基因列名（排除Accession列）
      gene_cols <- setdiff(colnames(hap), "Accession")
      
      # 创建组合：将所有基因列用"_"连接
      data$hap_combo <- apply(data[, gene_cols, drop = FALSE], 1, 
                              function(x) paste(x, collapse = "_"))
      combo_name <- paste(gene_cols, collapse = "_")
      
      # 获取性状列
      trait_cols <- grep("^T[0-9]{3}$", colnames(data), value = TRUE)
      
      # 对每个性状绘图
      for(trait_name in trait_cols) {
        # 过滤数据
        plot_data <- data[!is.na(data[[trait_name]]), ]
        
        # 计算每组样本数
        group_counts <- table(plot_data$hap_combo)
        valid_groups <- names(group_counts)[group_counts >= 5]
        
        if(length(valid_groups) >= 2) {
          plot_data <- plot_data[plot_data$hap_combo %in% valid_groups, ]
          plot_data$hap_combo <- factor(plot_data$hap_combo, levels = valid_groups)
          
          y_max <- max(plot_data[[trait_name]], na.rm = TRUE)
          y_min <- min(plot_data[[trait_name]], na.rm = TRUE)
          
          # 筛选显著比较
          comparisons <- combn(valid_groups, 2, simplify = FALSE)
          sig_comparisons <- list()
          
          for(comp in comparisons) {
            group1 <- plot_data[plot_data$hap_combo == comp[1], trait_name]
            group2 <- plot_data[plot_data$hap_combo == comp[2], trait_name]
            
            # 添加错误处理
            if(length(unique(group1)) > 1 && length(unique(group2)) > 1) {
              test_result <- tryCatch({
                t.test(group1, group2)
              }, error = function(e) NULL)
              
              if(!is.null(test_result) && test_result$p.value < 0.05) {
                sig_comparisons <- c(sig_comparisons, list(comp))
              }
            }
          }
          
          # 绘图
          p <- ggviolin(plot_data, 
                        x = "hap_combo", 
                        y = trait_name,
                        fill = "hap_combo",
                        palette = "npg",
                        add = "boxplot",
                        legend = "right") +
            labs(title = paste(combo_name, "-", trait_name), x = "") +
            guides(fill = guide_legend(override.aes = list(color = NA))) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
                  axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 8)) +
            scale_fill_discrete(labels = paste0(valid_groups, "(", group_counts[valid_groups], ")"))
          
          # 调整图例位置和大小以适应更多组合
          if(length(valid_groups) > 6) {
            p <- p + theme(legend.position = "bottom",
                           legend.direction = "horizontal")
          }
          
          if(length(sig_comparisons) > 0) {
            n_comparisons <- length(sig_comparisons)
            step_height <- (y_max - y_min) * 0.08
            
            p <- p + stat_compare_means(
              comparisons = sig_comparisons, 
              label = "p.format",
              hide.ns = TRUE,
              tip.length = 0.01,
              step.increase = 0.08,
              vjust = 1.8,
              size = 3  # 减小显著性标记字体
            ) +
              scale_y_continuous(limits = c(y_min, y_max + step_height * n_comparisons * 1.2))
          }
          
          # 根据组数调整图形尺寸
          n_groups <- length(valid_groups)
          plot_width <- max(10, n_groups * 1.5)  # 动态调整宽度
          plot_height <- ifelse(length(valid_groups) > 6, 8, 7)  # 更多组时增加高度
          
          # 保存到对应的文件夹
          output_path <- file.path(output_dir, paste0(trait_name, ".violin.pdf"))
          ggsave(output_path, p, width = plot_width, height = plot_height)
          
          message(paste("已保存:", hap_folder, "/", gene_name, "/", 
                       trait_category, "/", trait_name))
        }
      }
    }
  }
}
                                  
