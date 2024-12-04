###first, build index
srun -A pawsey0399 -p highmem -c 8 -n 1 -t 24:00:00 hisat2-build CS2Z559P1.fa CS2Z559P1.fa
samtools faidx CS2Z559P1.fa


####second, align
ls *_1.clean.fq.gz|while read line; do base_name=$(basename "$line" "_1.clean.fq.gz");  echo '#!/bin/bash
#SBATCH --job-name=hisat2
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=12:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate hisat
srun --export=all -n 1 -c 64 hisat2 --dta -x 0.Index/CS2Z559P1.fa -1 '${base_name}'_1.clean.fq.gz -2 '${base_name}'_2.clean.fq.gz -S '${base_name}'.sam'  > ${base_name}.hisat.sh;  done

####third, sort to bam file
ls *.sam|cut -f1 -d"."|while read line; do  echo '#!/bin/bash
#SBATCH --job-name=sort
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --account=pawsey0399
module load samtools/1.15--h3843a85_0
srun --export=all -n 1 -c 16 samtools sort -@ 16 -o '${line}'.sort.bam '${line}'.sam' >$line.sort.sh ; done

###Forth,featurecounts coount reads

srun -A Pawsey0399 -c 64 -n 1 -p work featureCounts -T 64 -t mRNA -g Parent -p -B -C -a 0.Index/CS2.1_Ac.gff3 -o 5113_2305_CS2_Ac.counts *.sort.bam

####anyway,if you mainly focus on special genes, maybe needs unique reads
ls *.sam|cut -f1 -d"."|while read line; do srun -A pawsey0399 -c 32 -n 1 -p work -t 2:00:00 grep "NH:i:1" $line.sam|grep "YT:Z:CP" > unique_reads/$line.unique.sam  & done
ls *.sam|cut -f1 -d"."|while read line; do srun -A pawsey0399 -c 32 -n 1 -p work -t 2:00:00 samtools view -H $line.sort.bam >unique_reads/$line.header; done 
cat $line.header $line.unique.sam |samtools sort -@ 32 -o $line.unique.sort.bam


####DEG
install.packages("pheatmap")
if (!require("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
BiocManager::install("DESeq2",force=T)


library(DESeq2)
library(pheatmap)

col_data<-read.csv("GROUP.csv",row.names = 1)
count_data<-read.csv("5113_2305_CS2_Ac.counts.csv",row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Condition)
rld <- rlog(dds)
sample_dist_matrix <- dist(t(assay(rld)))
sample_cluster <- hclust(sample_dist_matrix)
pheatmap(as.matrix(sample_dist_matrix), clustering_method = "complete", main = "Sample Clustering")

plotPCA(rld, intgroup = "Condition")

results_list <- list()
for (condition in unique(col_data$Condition)) {
    # 提取当前条件的样本
    condition_samples <- rownames(col_data[col_data$Condition == condition, ])
    
    # 提取当前条件的 count_data 和 col_data
    count_data_condition <- count_data[, condition_samples]
    col_data_condition <- col_data[condition_samples, ]
    
    # 确保 Group 是因子
    col_data_condition$Group <- factor(col_data_condition$Group, levels = c("2305", "5113"))
    
    # 创建 DESeqDataSet 对象
    dds_condition <- DESeqDataSetFromMatrix(countData = count_data_condition, colData = col_data_condition, design = ~ Group)
    
    # 运行 DESeq2 分析
    dds_condition <- DESeq(dds_condition)
    
    # 提取 5113 vs 2305 的比较结果
    res_condition <- results(dds_condition, contrast = c("Group", "5113", "2305"))
    
    # 将当前条件的结果存入列表
    results_list[[condition]] <- res_condition
    
    # 将结果保存为 CSV 文件
    output_filename <- paste0("DEG_results_", condition, ".csv")
    write.csv(as.data.frame(res_condition), file = output_filename)
    
    # 输出结果概述
    cat("Summary of Condition:", condition, "\n")
    print(summary(res_condition))
}




#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
####WGCNA

# Step 1: 数据准备与过滤
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# 导入数据
WGCNA.fpkm = read.table("spike.mrna.fpkm.txt", header = TRUE, check.names = FALSE)
dim(WGCNA.fpkm)  # 检查数据维度
names(WGCNA.fpkm)

# 数据转置，基因为列，样本为行
datExpr0 = as.data.frame(t(WGCNA.fpkm[,-1]))
names(datExpr0) = WGCNA.fpkm$sample
rownames(datExpr0) = names(WGCNA.fpkm[,-1])

# 检查样本和基因的质量
gsg = goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# 基因表达量过滤
meanFPKM = 1
n = nrow(datExpr0)
datExpr0[n+1, ] = apply(datExpr0[1:n, ], 2, mean)
datExpr0 = datExpr0[1:n, datExpr0[n+1, ] > meanFPKM]

# 保存过滤后的数据
filtered_fpkm = t(datExpr0)
filtered_fpkm = data.frame(rownames(filtered_fpkm), filtered_fpkm)
names(filtered_fpkm)[1] = "sample"
write.table(filtered_fpkm, file = "spike.mRNA.fpkm1.filter.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Step 2: 样本聚类与离群值检测
sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "1_spike_fpkm1_sampleClustering.pdf", width = 15, height = 8)
par(cex = 0.6)
par(mar = c(0, 6, 6, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 2, cex.axis = 1.5, cex.main = 2)
dev.off()

# 设置剪切高度，移除异常样本（可根据具体数据调整cutHeight）
cutHeight = 15000
clust = cutreeStatic(sampleTree, cutHeight = cutHeight, minClusterSize = 10)
keepSamples = (clust == 1)
datExpr0 = datExpr0[keepSamples, ]

# Step 3: 性状数据关联
traitData = read.csv("Spike_Trait.csv", row.names = 1)
fpkmSamples = rownames(datExpr0)
traitSamples = rownames(traitData)
traitRows = match(fpkmSamples, traitSamples)
datTraits = traitData[traitRows, ,drop=FALSE]
rownames(datTraits) = fpkmSamples

# 绘制样本聚类与性状热图
sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
pdf(file = "3_Sample_spike_fpkm1_dendrogram_and_trait_heatmap.pdf", width = 8, height = 6)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap", cex.colorLabels = 1.5, cex.dendroLabels = 1, cex.rowText = 2)
dev.off()

# Step 4: 软阈值选择
powers = c(1:30)
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# 绘制软阈值图
pdf(file = "4_spike_fpkm1_SF_value_select.pdf", width = 12, height = 8)
par(mfrow = c(1, 2))
cex1 = 0.85
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, cex = cex1, col = "red")
abline(h = 0.85, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

softPower = 6  # 选择最佳软阈值（需根据图表调整）

# Step 5: 网络构建与模块识别
adjacency = adjacency(datExpr0, power = softPower)
TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM

geneTree = hclust(as.dist(dissTOM), method = "average")
pdf(file = "5_spike_fpkm1_Gene_clustering_on_TOM-based_dissimilarity.pdf", width = 8, height = 6)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods)

# Step 6: 模块合并
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1 - cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")

pdf(file = "6_spike_fpkm1_Clustering_of_module_eigengenes.pdf", width = 8, height = 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h = 0.6, col = "red")  # 调整剪切高度
dev.off()

merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = 0.6, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

# 绘制合并模块结果
pdf(file = "7_spike_fpkm1_merged_dynamic.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

# Step 7: 模块-性状相关性分析，找到相关性模块
moduleTraitCor = cor(mergedMEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr0))

pdf(file = "8_spike_fpkm1_Module-trait_relationships.pdf", width = 8, height = 6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits), yLabels = names(mergedMEs), ySymbols = names(mergedMEs), colorLabels = FALSE, colors = greenWhiteRed(50), textMatrix = textMatrix, cex.text = 0.4, zlim = c(-1, 1), main = "Module-trait relationships")
dev.off()


# Step 8: 提取每个模块的基因信息并保存
modules = unique(mergedColors)
geneModules = data.frame(Gene = names(datExpr0), Module = mergedColors)
for (mod in modules) {
  genesInModule = geneModules[geneModules$Module == mod, "Gene"]
  expressionInModule = datExpr0[, genesInModule, drop = FALSE]
  write.table(expressionInModule, file = paste0("Module_", mod, "_expression.txt"),
              row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
}
write.table(geneModules, file = "All_Module_Genes.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Step 8: 提取相关性模块Hub基因
moduleColors = mergedColors

# --- Hub 基因分析 ---
module = "MElightgreen"
moduleGenes = (moduleColors == "lightgreen") ####be careful 两组的名称有“ME”的区别，mergeME中是MElightgreen，mergedColor中是lightgreen

# 计算 MM 和 GS
MM = cor(datExpr0, mergedMEs[, module], use = "p")
GS = cor(datExpr0, datTraits[, 1], use = "p") # 假设 datTraits 第 1 列是目标性状
hubGeneData = data.frame(Gene = colnames(datExpr0), ModuleMembership = MM, GeneSignificance = GS)
hubGenes = hubGeneData[abs(hubGeneData$ModuleMembership) > 0.8 & abs(hubGeneData$GeneSignificance) > 0.2, ]   #这里的0.8是较严格的标准，0.2是很宽松的标准
write.csv(hubGenes, "spike_fpkm1_lightgreen_hub_genes.csv", row.names = FALSE)


lightgreenGenes = colnames(datExpr0)[moduleGenes]
lightgreenTOM = TOM[moduleGenes, moduleGenes]
# 设置行列名为 lightgreenGenes
rownames(lightgreenTOM) = lightgreenGenes
colnames(lightgreenTOM) = lightgreenGenes
cyt = exportNetworkToCytoscape(lightgreenTOM,
                               edgeFile = "spike_fpkm1_lightgreen_edges.txt",
                               nodeFile = "spike_fpkm1_lightgreen_nodes.txt",
                               weighted = TRUE, threshold = 0.1)


# 用 igraph 可视化网络
library(igraph)
hubGeneNames = hubGenes$Gene
hubTOM = lightgreenTOM[hubGeneNames, hubGeneNames]
hubGraph = graph_from_adjacency_matrix(hubTOM, mode = "undirected", weighted = TRUE)
pdf("9_spike_fpkm1_lightgreen_hub_gene_network.pdf", width = 8, height = 6)
plot(hubGraph, vertex.size = 5, vertex.label.cex = 0.7,
     edge.width = E(hubGraph)$weight * 5,
     main = "Hub Gene Network in Lightgreen Module")
dev.off()

##画感兴趣基因的可视化网络
library(igraph)

# 假设 'hubGenes' 数据框中包含基因号和模块信息
# 设置6P相关基因号，示例为 " Ac6P01G496000"
gene_of_interest = "Ac6P01G496000"

# 查找与该基因相关的基因
# 假设你已经知道 `lightgreenTOM` 中哪些基因与 " Ac6P01G496000" 基因相关
# 你可以根据 TOM（相似性矩阵）中与该基因的相关性筛选
relatedGenes = colnames(lightgreenTOM)[apply(lightgreenTOM, 1, function(x) abs(x[gene_of_interest]) > 0.8)]

# 从 `lightgreenTOM` 矩阵中提取与相关基因子集构建子网络
hubTOM_subset = lightgreenTOM[relatedGenes, relatedGenes]
hubGraph_subset = graph_from_adjacency_matrix(hubTOM_subset, mode = "undirected", weighted = TRUE)

# 绘制子网络图
pdf("Ac6P01G496000_network.pdf", width = 8, height = 6)
plot(hubGraph_subset, vertex.size = 5, vertex.label.cex = 0.7,
     edge.width = E(hubGraph_subset)$weight * 5,
     main = paste("Hub Gene Network Related to", gene_of_interest, "in Lightgreen Module"))
dev.off()


######################################igraph画感兴趣基因loop
library(igraph)

# 读取 Cytoscape 导出的 edge 和 node 文件
edges <- read.table("spike_fpkm1_lightgreen_edges.txt", header = TRUE, sep = "\t")
nodes <- read.table("spike_fpkm1_lightgreen_nodes.txt", header = TRUE, sep = "\t")

# 创建 igraph 对象
graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# 筛选包含 "6P" 的基因
genes_of_interest <- nodes$nodeName[grepl("6P", nodes$nodeName)]  # 筛选节点名称中包含 "6P" 的基因

# 针对每个基因绘制网络图
for (gene in genes_of_interest) {
  # 找到与目标基因直接连接的节点
  neighbors <- neighbors(graph, gene, mode = "all")
  subgraph_nodes <- c(gene, names(neighbors))  # 包括目标基因及其邻居节点
  
  # 构建子图
  subgraph <- induced_subgraph(graph, vids = subgraph_nodes)
  
  # 可视化设置
  V(subgraph)$color <- ifelse(V(subgraph)$name == gene, "red", "lightgreen")  # 目标基因用红色，其余用模块色
  V(subgraph)$size <- ifelse(V(subgraph)$name == gene, 15, 10)               # 目标基因节点更大
  V(subgraph)$label <- V(subgraph)$name
  E(subgraph)$width <- E(subgraph)$weight * 5
  
  # 保存网络图到 PDF 文件
  output_file <- paste0("network_", gene, "_6P.pdf")
  pdf(output_file, width = 8, height = 6)
  plot(subgraph,
       vertex.label.color = "black", # 标签颜色
       edge.curved = 0.2,            # 边曲率
       layout = layout_with_fr,      # Fruchterman-Reingold 布局
       main = paste("Network of Gene:", gene)) # 标题
  dev.off()
}




#########不使用merge的module，用分开的module
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr0))

# 绘制原始模块-性状相关性的热图
pdf(file = "8_spike_fpkm1_Original_Module-trait_relationships.pdf", width = 8, height = 36)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(datTraits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = greenWhiteRed(50),
  textMatrix = textMatrix,
  cex.text = 0.4,
  zlim = c(-1, 1),
  main = "Original Module-Trait Relationships"
)
dev.off()


# Step 1: 提取原始模块的基因表达数据
modules = unique(dynamicColors) # 使用原始模块的颜色
geneModules = data.frame(Gene = names(datExpr0), Module = dynamicColors) # 基因和模块对应表
for (mod in modules) {
  genesInModule = geneModules[geneModules$Module == mod, "Gene"]
  expressionInModule = datExpr0[, genesInModule, drop = FALSE]
  write.table(expressionInModule, file = paste0("Original_Module_", mod, "_expression.txt"),
              row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
}
write.table(geneModules, file = "Original_All_Module_Genes.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Step 2: 提取原始模块的Hub基因
moduleColors = dynamicColors  # 使用原始模块颜色

# --- 选定目标模块 ---
module = "lightgreen"  # 注意：这里直接使用原始模块名称
moduleGenes = (moduleColors == module) # 筛选属于目标模块的基因

# 计算 MM 和 GS（与目标模块和性状的相关性）
MM = cor(datExpr0, MEs[, paste0("ME", module)], use = "p") # 使用原始模块特征值
GS = cor(datExpr0, datTraits[, 1], use = "p") # 假设 datTraits 第 1 列是目标性状
hubGeneData = data.frame(Gene = colnames(datExpr0), ModuleMembership = MM, GeneSignificance = GS)
hubGenes = hubGeneData[abs(hubGeneData$ModuleMembership) > 0.8 & abs(hubGeneData$GeneSignificance) > 0.2, ]   # 设置筛选标准
write.csv(hubGenes, paste0("Original_", module, "_hub_genes.csv"), row.names = FALSE)

# Step 3: 网络文件的生成
# 提取目标模块的基因和 TOM 矩阵
moduleGenesList = colnames(datExpr0)[moduleGenes]
moduleTOM = TOM[moduleGenes, moduleGenes]
rownames(moduleTOM) = moduleGenesList
colnames(moduleTOM) = moduleGenesList

# 导出 Cytoscape 文件
cyt = exportNetworkToCytoscape(moduleTOM,
                               edgeFile = paste0("Original_", module, "_edges.txt"),
                               nodeFile = paste0("Original_", module, "_nodes.txt"),
                               weighted = TRUE, threshold = 0.1) # 根据需求调整 threshold




pipeline:
	1.	数据准备与过滤
		数据读取 (WGCNA.fpkm = read.csv(...))： 导入RNA-Seq数据，以FPKM表达量作为分析的基础。
		数据转置 (datExpr0)： 将基因作为列，样本作为行，这是WGCNA分析的标准数据格式。
		样本和基因的质量控制 (goodSamplesGenes)： 检查数据中是否存在异常样本或基因，并剔除异常值。
		表达量过滤 (meanFPKM)： 设置基因的FPKM均值过滤阈值（如0.5），去除低表达或噪音基因，减少计算量。
	2.	样本聚类与离群值检测
		样本聚类 (sampleTree)： 使用层次聚类检测样本是否有异常。异常样本可能是技术误差或实验误差，需要剔除。
		绘制样本聚类图： 通过图示直观展示哪些样本可能需要过滤。
	3.	性状数据关联
		读取性状数据 (Trait.csv)： 导入样本对应的性状数据，为后续的模块-性状相关性分析做准备。
		性状数据匹配 (match)： 确保性状数据与表达数据一一对应，避免数据错位。
	4.	软阈值选择
		设置软阈值 (pickSoftThreshold)： WGCNA构建加权网络需要一个合适的软阈值（power）。通过选择能够最优地满足无尺度分布的值（如R² > 0.85）作为阈值。
		绘制软阈值图 (4_SF_value_select.pdf)： 帮助用户直观选择最佳软阈值。
	5.	网络构建与基因聚类
		网络构建 (adjacency, TOMsimilarity)： 计算基因间的加权邻接矩阵及拓扑相似性矩阵（TOM）。
		基因聚类 (geneTree)： 基于TOM矩阵聚类基因，生成聚类树。
	6.	模块识别与合并
		模块动态切割 (cutreeDynamic)： 使用动态树切割算法识别基因模块。
		模块合并 (mergeCloseModules)： 基于模块特征基因（ME）相似性合并相近模块，减少冗余模块。
	7.	模块-性状相关性分析
		模块特征基因 (moduleEigengenes)： 提取每个模块的代表性基因表达模式（ME）。
		相关性计算 (cor)： 计算模块ME与性状的相关性及显著性，确定性状相关模块。
		绘制模块-性状相关性热图： 展示各模块与性状间的相关性。
	8.	导出模块基因
		模块基因导出 (write.csv)： 根据模块颜色提取基因，并存储为文件，供后续分析使用。


Hub基因是指模块中与其他基因高度连接的基因，其在生物学功能中可能起核心作用。

	1.	筛选Hub基因的标准：
		基因在模块内部的连接度（kME值）显著较高（通常 > 0.8）。
		基因与目标性状的相关性较高（绝对值 > 0.6）。
	2.	提取Hub基因：
		使用moduleEigengenes提取模块特征基因后，计算每个基因的kME值。
		筛选与目标性状显著相关的模块，进一步提取这些模块中的Hub基因。

画Hub基因网络图

	1.	提取目标模块的邻接矩阵：
  	2.	使用igraph构建和绘制网络：
 	3.	可选：筛选高连接基因子图：
		对于大网络，可以设置连接阈值，仅保留权重较高的边。
	4.	保存网络：
可以将网络数据导出为文件，以便在Cytoscape等可视化工具中进一步分析：




