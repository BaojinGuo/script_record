library(ggplot2)
library(ggrepel)
library(factoextra)

HTI <- read.csv("Trait5-Hti.csv", row.names = 1, check.names = FALSE)
HSI <- read.csv("Trait5-Hsi.csv", row.names = 1, check.names = FALSE)

HSI_num <- as.data.frame(lapply(HSI, as.numeric))
rownames(HSI_num) <- rownames(HSI)
HSI_inv <- -HSI_num
HTI_z <- scale(HTI)
HSI_z <- scale(HSI_inv)

heat_mat <- cbind(HTI_z, HSI_z)
heat_mat_clean <- heat_mat[complete.cases(heat_mat), ]
pca_res <- prcomp(heat_mat_clean, center = TRUE, scale. = FALSE)
heat_score <- pca_res$x[, 1]
heat_score_df <- data.frame(
    Sample = rownames(heat_mat_clean),
    HeatScore_PC1 = heat_score
)

# 排序
heat_score_df <- heat_score_df[order(heat_score_df$HeatScore_PC1,
                                     decreasing = TRUE), ]

write.csv(heat_score_df,
          file = "Heat_tolerance_score_PC1.csv",
          row.names = FALSE)
pca_plot_df <- data.frame(
    Sample = rownames(heat_mat_clean),
    PC1 = pca_res$x[,1],
    PC2 = pca_res$x[,2]
)

ggplot(pca_plot_df, aes(PC1, PC2)) +
    geom_point(size = 2) +
    theme_bw() +
    labs(title = "PCA of multi-trait heat tolerance",
         x = "PC1 (overall heat tolerance)",
         y = "PC2")
dist_mat <- dist(heat_mat_clean, method = "euclidean")
hc <- hclust(dist_mat, method = "ward.D2")

plot(hc, cex = 0.6, main = "Hierarchical clustering of heat tolerance")
cluster <- cutree(hc, k = 3)

cluster_df <- data.frame(
    Sample = names(cluster),
    Cluster = cluster,
    HeatScore_PC1 = heat_score
)

write.csv(cluster_df,
          file = "Heat_tolerance_clusters.csv",
          row.names = FALSE)
pca_plot_df$Cluster <- factor(cluster[rownames(pca_plot_df)])

ggplot(pca_plot_df, aes(PC1, PC2, color = Cluster)) +
    geom_point(size = 2) +
    theme_bw() +
    labs(title = "PCA and clustering of heat tolerance")
top20 <- rownames(pca_plot_df)[
    order(pca_plot_df$PC1, decreasing = TRUE)[1:20]
]
ggplot(pca_plot_df, aes(PC1, PC2, color = Cluster)) +
    geom_point(size = 2) +
    geom_text_repel(
        data = pca_plot_df[top20, ],
        aes(label = Sample),
        size = 3,
        max.overlaps = Inf
    ) +
    theme_bw() +
    labs(
        title = "PCA and clustering of heat tolerance",
        subtitle = "Top 20 heat-tolerant genotypes highlighted"
    )
fviz_nbclust(heat_mat_clean, FUN = hcut, method = "silhouette")


#################GPToptimize###############################
############################
## 0. Libraries
############################
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(factoextra)
library(dplyr)

############################
## 1. Read data
############################
HTI <- read.csv("Trait5-Hti.csv",
                row.names = 1, check.names = FALSE)

HSI <- read.csv("Trait5-Hsi.csv",
                row.names = 1, check.names = FALSE)

############################
## 2. Numeric conversion
############################
HTI_num <- as.data.frame(lapply(HTI, as.numeric))
HSI_num <- as.data.frame(lapply(HSI, as.numeric))

rownames(HTI_num) <- rownames(HTI)
rownames(HSI_num) <- rownames(HSI)

############################
## 3. Direction unification
##    Define: larger = better heat tolerance
############################
HSI_inv <- -HSI_num   # HSI 越小越耐热 → 取反

############################
## 4. Z-score standardization
############################
HTI_z <- scale(HTI_num)
HSI_z <- scale(HSI_inv)

############################
## 5. Combine traits
############################
heat_mat <- cbind(HTI_z, HSI_z)

# remove missing
heat_mat_clean <- heat_mat[complete.cases(heat_mat), ]

write.csv(heat_mat_clean,
          "Heat_matrix_zscore.csv")

############################
## 6. PCA
############################
pca_res <- prcomp(heat_mat_clean,
                  center = TRUE,
                  scale. = FALSE)

############################
## 7. Force PC1 biological direction
##    PC1 positively correlated with HTI
############################
hti_mean <- rowMeans(
  HTI_z[rownames(heat_mat_clean), ],
  na.rm = TRUE
)

if (cor(pca_res$x[,1], hti_mean) < 0) {
  pca_res$x[,1] <- -pca_res$x[,1]
  pca_res$rotation[,1] <- -pca_res$rotation[,1]
}

############################
## 8. Heat tolerance score
############################
heat_score <- pca_res$x[,1]

heat_score_df <- data.frame(
  Sample = rownames(heat_mat_clean),
  HeatScore_PC1 = heat_score
)

heat_score_df <- heat_score_df %>%
  arrange(desc(HeatScore_PC1))

write.csv(heat_score_df,
          "Heat_tolerance_score_PC1.csv",
          row.names = FALSE)

############################
## 9. Top20 & Bottom20
############################
top20 <- heat_score_df$Sample[1:20]
bottom20 <- heat_score_df$Sample[
  (nrow(heat_score_df)-19):nrow(heat_score_df)
]

write.csv(data.frame(Top20 = top20),
          "Top20_heat_tolerant_genotypes.csv",
          row.names = FALSE)

############################
## 10. PCA scatter plot
############################
pca_plot_df <- data.frame(
  Sample = rownames(heat_mat_clean),
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2]
)

p1 <- ggplot(pca_plot_df, aes(PC1, PC2)) +
  geom_point(size = 2) +
  geom_text_repel(
    data = pca_plot_df[pca_plot_df$Sample %in% top20, ],
    aes(label = Sample),
    size = 3,
    max.overlaps = Inf
  ) +
  theme_bw() +
  labs(
    title = "PCA of multi-trait heat tolerance",
    subtitle = "Top 20 heat-tolerant genotypes highlighted",
    x = paste0("PC1 (", round(
      summary(pca_res)$importance[2,1]*100,1), "%)"),
    y = paste0("PC2 (", round(
      summary(pca_res)$importance[2,2]*100,1), "%)")
  )

ggsave("PCA_heat_tolerance_top20.pdf",
       p1, width = 7, height = 6)

############################
## 11. Hierarchical clustering
############################
dist_mat <- dist(heat_mat_clean, method = "euclidean")
hc <- hclust(dist_mat, method = "ward.D2")

pdf("Heat_tolerance_hclust.pdf", width = 8, height = 6)
plot(hc, cex = 0.6,
     main = "Hierarchical clustering of heat tolerance")
dev.off()

############################
## 12. Silhouette-based K selection
############################
pdf("Silhouette_K_selection.pdf", width = 7, height = 5)
fviz_nbclust(heat_mat_clean,
             FUN = hcut,
             method = "silhouette")
dev.off()

############################
## 13. Final clustering (K = 4)
############################
cluster <- cutree(hc, k = 4)

cluster_df <- data.frame(
  Sample = names(cluster),
  Cluster = factor(cluster),
  HeatScore_PC1 = heat_score
)

write.csv(cluster_df,
          "Heat_tolerance_clusters_K4.csv",
          row.names = FALSE)

############################
## 14. PCA + cluster visualization
############################
pca_plot_df$Cluster <- factor(cluster[rownames(pca_plot_df)])

p2 <- ggplot(pca_plot_df,
             aes(PC1, PC2, color = Cluster)) +
  geom_point(size = 2) +
  geom_text_repel(
    data = pca_plot_df[pca_plot_df$Sample %in% top20, ],
    aes(label = Sample),
    size = 3,
    max.overlaps = Inf
  ) +
  theme_bw() +
  labs(
    title = "PCA and clustering of heat tolerance",
    subtitle = "K = 4 clusters; Top 20 highlighted"
  )

ggsave("PCA_cluster_K4_top20.pdf",
       p2, width = 7, height = 6)

############################
## 15. PC1 loading barplot
############################
loading_df <- data.frame(
  Trait = rownames(pca_res$rotation),
  Loading_PC1 = pca_res$rotation[,1]
)

loading_df <- loading_df %>%
  arrange(Loading_PC1)

p3 <- ggplot(loading_df,
             aes(x = reorder(Trait, Loading_PC1),
                 y = Loading_PC1)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "PC1 loadings of heat tolerance traits",
    x = "Trait",
    y = "PC1 loading"
  )

ggsave("PC1_loading_barplot.pdf",
       p3, width = 7, height = 8)

#########Deepseek modified version#######
library(ggplot2)
library(ggrepel)
library(dplyr)
library(factoextra)

# 1. 读取数据
HTI <- read.csv("Trait5-Hti.csv", row.names = 1, check.names = FALSE)
HSI <- read.csv("Trait5-Hsi.csv", row.names = 1, check.names = FALSE)

# 2. 转换为数值型
HTI_num <- as.data.frame(lapply(HTI, as.numeric))
HSI_num <- as.data.frame(lapply(HSI, as.numeric))

rownames(HTI_num) <- rownames(HTI)
rownames(HSI_num) <- rownames(HSI)

# 3. 处理缺失值（重要！）
cat("HTI缺失值统计:\n")
print(colSums(is.na(HTI_num)))
cat("HSI缺失值统计:\n")
print(colSums(is.na(HSI_num)))

# 建议：如果缺失值过多，考虑插补或删除
# HTI_num <- na.omit(HTI_num)
# HSI_num <- na.omit(HSI_num)

# 4. 明确处理逻辑
# 假设：HTI越高越耐热，HSI越高越敏感
HSI_inv <- -HSI_num  # 转换后：值越高表示越耐热

# 5. 标准化
HTI_z <- scale(HTI_num)
HSI_z <- scale(HSI_inv)

# 6. 合并矩阵（HTI和-HSI都是越高越耐热）
heat_mat <- cbind(HTI_z, HSI_z)
heat_mat <- heat_mat[complete.cases(heat_mat), ]

# 检查矩阵维度
cat("\n最终矩阵维度:", dim(heat_mat), "\n")

# 7. PCA分析
pca_res <- prcomp(heat_mat, center = TRUE, scale. = FALSE)

# 8. 更合理的PCA方向判断方法
# 方法1：基于HTI的平均值（但需要确保HTI是正向指标）
hti_original <- rowMeans(HTI_num[rownames(heat_mat), ], na.rm = TRUE)  # 使用原始HTI，不是标准化后的
if (cor(pca_res$x[,1], hti_original) < 0) {
    pca_res$x[,1] <- -pca_res$x[,1]
    pca_res$rotation[,1] <- -pca_res$rotation[,1]
}

# 方法2：基于生物学意义（推荐）
# 计算综合耐热评分：HTI平均值 + (-HSI平均值)
comprehensive_score <- rowMeans(
    cbind(HTI_num[rownames(heat_mat), ], 
          HSI_inv[rownames(heat_mat), ]), 
    na.rm = TRUE
)

if (cor(pca_res$x[,1], comprehensive_score) < 0) {
    pca_res$x[,1] <- -pca_res$x[,1]
    pca_res$rotation[,1] <- -pca_res$rotation[,1]
    cat("PCA方向已根据综合耐热评分调整\n")
}

# 9. 验证方向是否正确
cat("\n验证PCA方向:\n")
cat("PC1与原始HTI平均值的相关性:", cor(pca_res$x[,1], hti_original), "\n")
cat("PC1与综合耐热评分的相关性:", cor(pca_res$x[,1], comprehensive_score), "\n")

# PC1应该与耐热性正相关
heat_score <- pca_res$x[,1]
heat_score_df <- data.frame(
    Sample = rownames(heat_mat),
    HeatScore_PC1 = heat_score,
    HTI_Mean = hti_original,
    Comprehensive_Score = comprehensive_score
) %>% arrange(desc(HeatScore_PC1))

# 查看前10个样本验证
cat("\n前10个最耐热的样本（按PC1排序）:\n")
print(head(heat_score_df, 10))

# 10. 绘制PCA并验证
# 用颜色表示原始HTI平均值
hti_quantile <- cut(hti_original, 
                    breaks = quantile(hti_original, probs = seq(0, 1, 0.2)),
                    labels = c("Very Low", "Low", "Medium", "High", "Very High"))

pca_df <- data.frame(
    Sample = rownames(heat_mat),
    PC1 = pca_res$x[,1],
    PC2 = pca_res$x[,2],
    HTI_Mean = hti_original,
    HTI_Level = hti_quantile
)

# 绘制验证图
p_validation <- ggplot(pca_df, aes(PC1, PC2, color = HTI_Level)) +
    geom_point(size = 2) +
    scale_color_brewer(palette = "RdYlBu", 
                      name = "HTI Level\n(越高越耐热)") +
    theme_bw() +
    labs(title = "PCA of Heat Tolerance (按HTI平均值着色)",
         subtitle = "验证PC1方向是否正确",
         x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 1), "%)"),
         y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 1), "%)"))

ggsave("PCA_validation.pdf", p_validation, width = 8, height = 6)

# 11. 继续您的后续分析（聚类等）
top20 <- heat_score_df$Sample[1:20]

p1 <- ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(size = 2) +
    geom_text_repel(
        data = pca_df[pca_df$Sample %in% top20, ],
        aes(label = Sample),
        size = 3
    ) +
    theme_bw() +
    labs(title = "PCA of heat tolerance",
         subtitle = "Top 20 heat-tolerant genotypes",
         x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 1), "%)"),
         y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 1), "%)"))

ggsave("PCA_Top20_corrected.pdf", p1, width = 7, height = 6)

# 12. 加载载荷分析
loading_df <- data.frame(
    Trait = rownames(pca_res$rotation),
    Loading_PC1 = pca_res$rotation[,1]
) %>% arrange(Loading_PC1)

# 检查载荷方向
cat("\nPC1载荷最高的10个性状:\n")
print(tail(loading_df, 10))
cat("\nPC1载荷最低的10个性状:\n")
print(head(loading_df, 10))

p4 <- ggplot(loading_df,
             aes(reorder(Trait, Loading_PC1), Loading_PC1,
                 fill = Loading_PC1 > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("#F8766D", "#00BFC4"),
                      labels = c("Negative", "Positive"),
                      name = "Effect on\nHeat Tolerance") +
    theme_bw() +
    labs(title = "PC1 loadings of heat tolerance traits",
         subtitle = "Positive loadings contribute to heat tolerance",
         x = "Trait", y = "Loading on PC1")

ggsave("PC1_loadings_corrected.pdf", p4, width = 8, height = 10)

# 13. 输出结果
write.csv(heat_mat, "Heat_matrix_zscore_corrected.csv")
write.csv(heat_score_df, "HeatScore_PC1_corrected.csv", row.names = FALSE)
write.csv(loading_df, "PC1_loadings_values.csv", row.names = FALSE)
##############最终修改版####################################
###########################
## 0. Libraries
###########################
library(ggplot2)
library(ggrepel)
library(dplyr)
library(factoextra)

###########################
## 1. Read Data
###########################
HTI <- read.csv("Trait5-Hti.csv", row.names = 1, check.names = FALSE)
HSI <- read.csv("Trait5-Hsi.csv", row.names = 1, check.names = FALSE)

HTI_num <- as.data.frame(lapply(HTI, as.numeric))
HSI_num <- as.data.frame(lapply(HSI, as.numeric))
rownames(HTI_num) <- rownames(HTI)
rownames(HSI_num) <- rownames(HSI)

###########################
## 2. Z-score standardization
###########################
HTI_z <- scale(HTI_num)
HSI_z <- scale(-HSI_num)  # HSI 越小越耐热

write.csv(HTI_z, "HTI_zscore.csv")
write.csv(HSI_z, "HSI_zscore.csv")

###########################
## 3. Remove rows with NA / Inf
###########################
HTI_z <- HTI_z[complete.cases(HTI_z) & !apply(HTI_z,1,function(x) any(is.infinite(x))), ]
HSI_z <- HSI_z[complete.cases(HSI_z) & !apply(HSI_z,1,function(x) any(is.infinite(x))), ]

###########################
## 4. HTI PCA & Cluster
###########################

# Silhouette analysis decide K
p_sil_HTI <- fviz_nbclust(HTI_z, FUN=hcut, method="silhouette")
ggsave("HTI_Silhouette.pdf", p_sil_HTI, width=7, height=6)
p_sil_HSI <- fviz_nbclust(HSI_z, FUN=hcut, method="silhouette")                                             
ggsave("HSI_Silhouette.pdf", p_sil_HTI, width=7, height=6)
pca_HTI <- prcomp(HTI_z, center = TRUE, scale. = FALSE)
if(cor(pca_HTI$x[,1], rowMeans(HTI_z)) < 0){
  pca_HTI$x[,1] <- -pca_HTI$x[,1]
  pca_HTI$rotation[,1] <- -pca_HTI$rotation[,1]
}

# 聚类 K=3
dist_HTI <- dist(HTI_z)
hc_HTI <- hclust(dist_HTI, method="ward.D2")
cluster_HTI_k3 <- cutree(hc_HTI, k=3)

# PCA 数据框
pca_df_HTI <- data.frame(
  Sample = rownames(HTI_z),
  PC1 = pca_HTI$x[,1],
  PC2 = pca_HTI$x[,2],
  Cluster = factor(cluster_HTI_k3)
)

# 给 cluster 设置英文耐热等级
# 假设 cluster 1 = sensitive, 2 = moderate, 3 = tolerant （可以根据 PC1 均值调整）
cluster_mean_PC1 <- tapply(pca_df_HTI$PC1, pca_df_HTI$Cluster, mean)
cluster_rank <- order(cluster_mean_PC1)  # 从低到高
level_HTI <- rep(NA, 3)
level_HTI[cluster_rank] <- c("Sensitive", "Moderate", "Tolerant")
names(level_HTI) <- as.character(1:3)

pca_df_HTI$HeatLevel <- level_HTI[as.character(pca_df_HTI$Cluster)]

# 标注 Top20
top20_HTI <- pca_df_HTI %>% arrange(desc(PC1)) %>% slice(1:20)

# 绘图
p_HTI <- ggplot(pca_df_HTI, aes(x=PC1, y=PC2, color=HeatLevel)) +
  geom_point(size=2) +
  geom_text_repel(data=top20_HTI, aes(label=Sample), size=3, max.overlaps = Inf) +
  scale_color_manual(values = c("Sensitive"="red", "Moderate"="orange", "Tolerant"="green")) +
  theme_bw() +
  labs(
    title="HTI: PCA and K=3 Cluster",
    x=paste0("PC1 (", round(summary(pca_HTI)$importance[2,1]*100,1),"%)"),
    y=paste0("PC2 (", round(summary(pca_HTI)$importance[2,2]*100,1),"%)"),
    color="Heat Level"
  )

ggsave("HTI_PCA_K3_Cluster.pdf", p_HTI, width=8, height=6)

###########################
# 2. HSI PCA + Cluster K=2
###########################
pca_HSI <- prcomp(HSI_z, center = TRUE, scale. = FALSE)
if(cor(pca_HSI$x[,1], rowMeans(HSI_z)) < 0){
  pca_HSI$x[,1] <- -pca_HSI$x[,1]
  pca_HSI$rotation[,1] <- -pca_HSI$rotation[,1]
}

# 聚类 K=2
dist_HSI <- dist(HSI_z)
hc_HSI <- hclust(dist_HSI, method="ward.D2")
cluster_HSI_k2 <- cutree(hc_HSI, k=2)

# PCA 数据框
pca_df_HSI <- data.frame(
  Sample = rownames(HSI_z),
  PC1 = pca_HSI$x[,1],
  PC2 = pca_HSI$x[,2],
  Cluster = factor(cluster_HSI_k2)
)

# 设置英文耐热等级（PC1 高 → tolerant, 低 → sensitive）
cluster_mean_PC1_HSI <- tapply(pca_df_HSI$PC1, pca_df_HSI$Cluster, mean)
cluster_rank_HSI <- order(cluster_mean_PC1_HSI)  # 低到高
level_HSI <- rep(NA, 2)
level_HSI[cluster_rank_HSI] <- c("Sensitive", "Tolerant")
names(level_HSI) <- as.character(1:2)

pca_df_HSI$HeatLevel <- level_HSI[as.character(pca_df_HSI$Cluster)]

# 标注 Top20
top20_HSI <- pca_df_HSI %>% arrange(desc(PC1)) %>% slice(1:20)

# 绘图
p_HSI <- ggplot(pca_df_HSI, aes(x=PC1, y=PC2, color=HeatLevel)) +
  geom_point(size=2) +
  geom_text_repel(data=top20_HSI, aes(label=Sample), size=3, max.overlaps = Inf) +
  scale_color_manual(values=c("Sensitive"="red", "Tolerant"="green")) +
  theme_bw() +
  labs(
    title="HSI: PCA and K=2 Cluster",
    x=paste0("PC1 (", round(summary(pca_HSI)$importance[2,1]*100,1),"%)"),
    y=paste0("PC2 (", round(summary(pca_HSI)$importance[2,2]*100,1),"%)"),
    color="Heat Level"
  )

ggsave("HSI_PCA_K2_Cluster.pdf", p_HSI, width=8, height=6)
                                              

