# 加载必要的R包
install.packages("kohonen")  # 如果未安装，请先安装
install.packages("ggplot2")  # 用于可视化
library(viridis)

library(kohonen)
library(ggplot2)

# 读取数据（假设数据格式为基因表达矩阵，行为基因，列为样本）
data <- read.csv("expression_data.csv", row.names = 1)  # 需要替换为你的输入文件

# 数据标准化（SOM 需要标准化数据）
data_scaled <- scale(data)

# 设置SOM网格（5x5, 可调整）
som_grid <- somgrid(xdim = 5, ydim = 5, topo = "hexagonal")  # 5x5网格

# 训练SOM模型
som_model <- som(data_scaled, grid = som_grid, rlen = 100)  # 迭代100次

custom_palette <- function(n) {
  colors <- c(
    colorRampPalette(c("#006400", "#90EE90"))(2),  # 深绿 → 浅绿（Leaf）
    colorRampPalette(c("#00008B", "#87CEFA"))(2),  # 深蓝 → 浅蓝（Phloem）
    colorRampPalette(c("#FF4500", "#FFD700"))(2),  # 深橙 → 浅黄（Pod）
    colorRampPalette(c("#800080", "#DA70D6"))(2)   # 深紫 → 浅紫（Seed）
  )
  return(colors)
}
# 可视化SOM结果
plot(som_model, type = "codes")  # 查看神经元权重向量
plot(som_model, type = "counts") # 查看每个单元中的样本数量
plot(som_model, type = "mapping") # 显示样本的映射情况
plot(som_model, type = "dist.neighbours") # 显示神经元间距离

# 聚类分析（对SOM的神经元进行聚类）
som_cluster <- cutree(hclust(dist(som_model$codes[[1]])), k = 4) # 设定4个聚类
add.cluster.boundaries(som_model, som_cluster)

# 可视化聚类结果
plot(som_model, type = "property", property = som_cluster, main = "SOM Clustering")

# 保存结果
write.csv(som_model$codes[[1]], "som_cluster_centroids.csv")
