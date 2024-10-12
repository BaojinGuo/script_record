conda install -c bioconda r-wgcna
#library(reshape2)
#library(stringr)
#library(tidyverse)
#library(magrittr) 
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
WGCNA.fpkm = read.csv("SLZ_RNAseq_FPKM.csv",header=T,check.names=F)
dim(WGCNA.fpkm)
names(WGCNA.fpkm)
datExpr0 = as.data.frame(t(WGCNA.fpkm[,-1]))
names(datExpr0) = WGCNA.fpkm$Geneid
rownames(datExpr0) = names(WGCNA.fpkm[,-1])
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

meanFPKM=0.5  ###--过滤标准，可以修改
n=nrow(datExpr0)
datExpr0[n+1,]=apply(datExpr0[c(1:nrow(datExpr0)),],2,mean)
datExpr0=datExpr0[1:n,datExpr0[n+1,] > meanFPKM]

filtered_fpkm=t(datExpr0)
filtered_fpkm=data.frame(rownames(filtered_fpkm),filtered_fpkm)

names(filtered_fpkm)[1]="sample"
write.table(filtered_fpkm, file="SLZ.mRNA.filter.txt",
            row.names=F, col.names=T,quote=FALSE,sep="\t")


sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "1.sampleClustering.pdf", width = 15, height = 8)
par(cex = 0.6)
par(mar = c(0,6,6,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2,
     cex.axis = 1.5, cex.main = 2)
### Plot a line to show the cut
#abline(h = 180, col = "red")##剪切高度不确定，故无红线
dev.off()



pdf("2_sample clutering_delete_outliers.pdf", width = 6, height = )
par(cex = 0.6)
par(mar = c(0,6,6,0))
cutHeight = 15000  ### cut高度
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1.5, cex.main = 2) +
  abline(h = cutHeight, col = "red")    ###'@"h = 1500"参数为你需要过滤的参数高度
dev.off()

traitData = read.csv("Trait.csv",row.names=1)
allTraits = traitData
dim(allTraits)
names(allTraits)
fpkmSamples = rownames(datExpr0)
traitSamples =rownames(allTraits)
traitRows = match(fpkmSamples, traitSamples)
datTraits = allTraits[traitRows,]
datTraits <- as.data.frame(datTraits)
rownames(datTraits) <- fpkmSamples
rownames(datTraits)


sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)                
pdf(file="3_Sample_dendrogram_and_trait_heatmap.pdf",width=8 ,height= 6)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap",cex.colorLabels = 1.5, cex.dendroLabels = 1, cex.rowText = 2)
dev.off()

enableWGCNAThreads()
# 设置soft-thresholding powers的数量
powers = c(1:30)
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)  


#'@绘制soft power plot 
pdf(file="4_SF_value_select.pdf",width=12, height = 8)
par(mfrow = c(1,2))
cex1 = 0.85
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

softPower = 7                    
adjacency = adjacency(datExpr0, power = softPower)
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM                   

geneTree = hclust(as.dist(dissTOM), method = "average")
pdf(file="4_Gene_clustering_on_TOM-based_dissimilarity.pdf",width=6,height=4)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()


minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
pdf(file="5_Dynamic_Tree_Cut.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
#sizeGrWindow(7, 6)
pdf(file="6_Clustering_of_module_eigengenes.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.6######剪切高度可修改
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs
table(mergedColors)

#输出所有modules
color<-unique(mergedColors)
for (i  in 1:length(color)) {
  y=t(assign(paste(color[i],"expr",sep = "."),datExpr0[mergedColors==color[i]]))
  write.csv(y,paste('6',color[i],"csv",sep = "."),quote = F)
}
#save.image(file = "module_splitted.RData")

pdf(file="7_merged_dynamic.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()



moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs


# Define numbers of genes and samples
# 获得基因数和样本数
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);



# Recalculate MEs with color labels
# 用彩色标签重新计算MEs
# 在给定的单个数据集中计算模块的模块本征基因
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
# 对给定的（特征）向量进行重新排序，以使相似的向量（通过相关性度量）彼此相邻
MEs = orderMEs(MEs0)

# 计算module的ME值与表型的相关系数
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);



pdf(file="8_Module-trait_relationships.pdf", width = 8, height = 6)
# sizeGrWindow(10,6)
# 显示相关性及其p值
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot\
# ySymbols 当ylabels使用时所使用的其他标签； colorLabels 应该使用颜色标签吗
# colors 颜色； textMatrix 单元格名字
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(datTraits),yLabels = names(MEs),ySymbols = names(MEs),
               colorLabels = FALSE,colors = greenWhiteRed(50),textMatrix = textMatrix,setStdMargins = FALSE,
               cex.text = 0.4,zlim = c(-1,1),
main = paste("Module-trait relationships"))
dev.off()










