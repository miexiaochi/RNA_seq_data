library(WGCNA)
fpkm = read.table("ExpData_WGCNA.txt",header=T,
                  comment.char = "",
                  check.names=F)
datTraits =read.table("sample.txt",row.names=1,header=T,comment.char = "",check.names=F)

head(fpkm)
RNAseq_voom <- fpkm 
##移除基因名称列
RNAseq_voom <- RNAseq_voom[, sapply(RNAseq_voom, is.numeric)]
## 因为WGCNA针对的是基因进行聚类，而一般我们的聚类是针对样本用hclust即可，所以这个时候需要转置。
WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:5000],])

datExpr0 <- WGCNA_matrix  ## top 5000 mad genes
datExpr <- datExpr0 

## 下面主要是为了防止临床表型与样本名字对不上
sampleNames = rownames(datExpr);
traitRows = match(sampleNames, datTraits$sample)  
rownames(datTraits) = datTraits[traitRows, 1]

powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#设置网络构建参数选择范围，计算无尺度分布拓扑矩阵

# Plot the results:
##sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net = blockwiseModules(
  datExpr,
  power = 10,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = F, 
  verbose = 3
)
table(net$colors) 
print(sft$powerEstimate)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
## assign all of the gene to their corresponding module 
## hclust for the genes.


#明确样本数和基因数
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#首先针对样本做个系统聚类树
datExpr_tree<-hclust(dist(datExpr), method = "average")
par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1, cex.main = 1,cex.lab=1)
## 如果这个时候样本是有性状，或者临床表型的，可以加进去看看是否聚类合理
#针对前面构造的样品矩阵添加对应颜色
sample_colors <- numbers2colors(as.numeric(factor(datTraits$subtype)), 
                                colors = c("blue","red","green"),signed = FALSE)
## 这个给样品添加对应颜色的代码需要自行修改以适应自己的数据分析项目。
#  sample_colors <- numbers2colors( datTraits ,signed = FALSE)
## 如果样品有多种分类情况，而且 datTraits 里面都是分类信息，那么可以直接用上面代码，当然，这样给的颜色不明显，意义不大。
#构造10个样品的系统聚类树及性状热图
par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(datExpr_tree, sample_colors,
                    groupLabels = colnames(sample),
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")


## step 5 (最重要的) 模块和性状的关系
## 这一步主要是针对于连续变量，如果是分类变量，需要转换成连续变量方可使用
table(datTraits$subtype)
if(T){
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  

  datTraits$subtype <- factor(datTraits$subtype, levels = c("8cell", "16cell", "32cell"))
  
  print(levels(datTraits$subtype))
  
  
  
  design=model.matrix(~0+ datTraits$subtype)
  colnames(design)=levels(datTraits$subtype)
  
  moduleColors <- labels2colors(net$colors)
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩 (样本vs模块)
  moduleTraitCor = cor(MEs, design , use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)

  dim(moduleTraitCor)
  dim(design)
  
  png("step5-Module-trait-relationships.png",width = 800,height = 1200,res = 120)
  par(mar = c(6, 8.5, 3, 3));
  print(design)
  print(colnames(design))
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  
  # 除了上面的热图展现形状与基因模块的相关性外
  # 还可以是条形图,但是只能是指定某个形状
  # 或者自己循环一下批量出图。
  Luminal = as.data.frame(design[,3]);
  names(Luminal) = "Luminal"
  y=Luminal
  GS1=as.numeric(cor(y,datExpr, use="p"))
  GeneSignificance=abs(GS1)
  # Next module significance is defined as average gene significance.
  ModuleSignificance=tapply(GeneSignificance,
                            moduleColors, mean, na.rm=T)
  sizeGrWindow(8,7)
  par(mfrow = c(1,1))
  # 如果模块太多，下面的展示就不友好
  # 不过，我们可以自定义出图。
  plotModuleSignificance(GeneSignificance,moduleColors)
  
}

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
## 算出每个模块跟基因的皮尔森相关系数矩阵
## MEs是每个模块在每个样本里面的值
## datExpr是每个基因在每个样本的表达量
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
## 只有连续型性状才能只有计算
## 这里把是否属于 Luminal 表型这个变量用0,1进行数值化。
Luminal = as.data.frame(design[,3]);
names(Luminal) = "Luminal"
geneTraitSignificance = as.data.frame(cor(datExpr, Luminal, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Luminal), sep="");
names(GSPvalue) = paste("p.GS.", names(Luminal), sep="");
module = "red"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Luminal",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# 主要是可视化 TOM矩阵，WGCNA的标准配图
# 然后可视化不同 模块 的相关性 热图
# 不同模块的层次聚类图
# 还有模块诊断，主要是 intramodular connectivity
if(T){
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  geneTree = net$dendrograms[[1]]; 
  dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6); 
  plotTOM = dissTOM^7; 
  diag(plotTOM) = NA; 
  #TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
  nSelect = 400
  # For reproducibility, we set the random seed
  set.seed(10);
  select = sample(nGenes, size = nSelect);
  selectTOM = dissTOM[select, select];
  # There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
  selectTree = hclust(as.dist(selectTOM), method = "average")
  selectColors = moduleColors[select];
  # Open a graphical window
  sizeGrWindow(9,9)
  # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
  # the color palette; setting the diagonal to NA also improves the clarity of the plot
  plotDiss = selectTOM^7;
  diag(plotDiss) = NA;
  
  png("step7-Network-heatmap.png",width = 800,height = 600)
  TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
  dev.off()
  
  # Recalculate module eigengenes
  MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
  ## 只有连续型性状才能只有计算
  ## 这里把是否属 Luminal 表型这个变量0,1进行数值化
  Luminal = as.data.frame(design[,3]);
  names(Luminal) = "Luminal"
  # Add the weight to existing module eigengenes
  MET = orderMEs(cbind(MEs, Luminal))
  # Plot the relationships among the eigengenes and the trait
  sizeGrWindow(5,7.5);
  
  par(cex = 0.9)
  png("step7-Eigengene-dendrogram.png",width = 800,height = 600)
  plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                        = 90)
  dev.off()
  
  # Plot the dendrogram
  sizeGrWindow(6,6);
  par(cex = 1.0)
  ## 模块的进化树
  png("step7-Eigengene-dendrogram-hclust.png",width = 800,height = 600)
  plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                        plotHeatmaps = FALSE)
  dev.off()
  # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
  par(cex = 1.0)
  ## 性状与模块热
  
  png("step7-Eigengene-adjacency-heatmap.png",width = 800,height = 600)
  plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                        plotDendrograms = FALSE, xLabelsAngle = 90)
  dev.off()
  
}

nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
## 只有连续型性状才能只有计算
## 这里把是否属于 Luminal 表型这个变量用0,1进行数值化。
Luminal = as.data.frame(design[,3]);
names(Luminal) = "Luminal"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, Luminal))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
## 模块的聚类图
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
## 性状与模块热图
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
## step 8 
# 主要是关心具体某个模块内部的基因
if(T){
  # Select module
  module = "red";
  # Select module probes
  probes = colnames(datExpr) ## 我们例子里面的probe就是基因
  inModule = (moduleColors==module);
  modProbes = probes[inModule]; 
  head(modProbes)
  
  # 如果使用WGCNA包自带的热图就很丑。
  which.module="red";
  dat=datExpr[,moduleColors==which.module ] 
  plotMat(t(scale(dat)),nrgcols=30,rlabels=T,
          clabels=T,rcols=which.module,
          title=which.module )
  datExpr[1:4,1:4]
  dat=t(datExpr[,moduleColors==which.module ] )
  library(pheatmap)
  pheatmap(dat ,show_colnames =F,show_rownames = F) #对那些提取出来的1000个基因所在的每一行取出，组合起来为一个新的表达矩阵
  n=t(scale(t(log(dat+1)))) # 'scale'可以对log-ratio数值进行归一化
  n[n>2]=2 
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  group_list=datTraits$subtype
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n) 
  pheatmap(n,show_colnames =F,show_rownames = F,
           annotation_col=ac )
  # 可以很清晰的看到，所有的形状相关的模块基因
  # 其实未必就不是差异表达基因。
}
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 6); 
# Select module
module = "red";
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
## 模块对应的基因关系矩阵 
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0)
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule]
);
nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
filter <- modTOM[top, top]
