# 加载必要的R包
library(DESeq2)
library(ggplot2)
library(ggrepel)
getwd()
# 设置工作目录到数据文件所在的目录
setwd("J:/XY/RNASEQ/5.featurecounts,FPKM,差异分析")  # 请将此路径修改为实际数据所在路径

# 读取counts数据
counts <- read.table("counts.8.16.32.txt", header = TRUE, row.names = 1)

# 提取样本名称
samples <- colnames(counts)[1:ncol(counts)]  # 假设前6列是基因注释信息
print(samples)
# 创建样本信息数据框
sample_info <- data.frame(
  sample = samples,
  group = c(rep("8cell", 4), rep("16cell", 4), rep("32cell", 4)),
  condition = c(rep("incl", 2), rep("ivcl", 2), rep("incl", 2), rep("ivcl", 2), rep("incl", 2), rep("ivcl", 2))
)
rownames(sample_info) <- samples

# 提取counts矩阵
counts_matrix <- as.matrix(counts[, samples])

# 创建DESeq2数据集
dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = sample_info, design = ~ condition + group)

# 运行DESeq2差异表达分析
dds <- DESeq(dds)

# 提取差异表达结果
results_8cell_vs_16cell <- results(dds, contrast = c("group", "8cell", "16cell"))
results_8cell_vs_32cell <- results(dds, contrast = c("group", "8cell", "32cell"))
results_16cell_vs_32cell <- results(dds, contrast = c("group", "16cell", "32cell"))

# 保存结果到文件
write.csv(as.data.frame(results_8cell_vs_16cell), file = "results_8cell_vs_16cell.csv")
write.csv(as.data.frame(results_8cell_vs_32cell), file = "results_8cell_vs_32cell.csv")
write.csv(as.data.frame(results_16cell_vs_32cell), file = "results_16cell_vs_32cell.csv")

# 提取显著差异基因
significant_8cell_vs_16cell <- subset(results_8cell_vs_16cell, padj < 0.05 & abs(log2FoldChange) > 1)
significant_8cell_vs_32cell <- subset(results_8cell_vs_32cell, padj < 0.05 & abs(log2FoldChange) > 1)
significant_16cell_vs_32cell <- subset(results_16cell_vs_32cell, padj < 0.05 & abs(log2FoldChange) > 1)

# 保存显著差异基因到文件
write.csv(as.data.frame(significant_8cell_vs_16cell), file = "significant_8cell_vs_16cell.csv")
write.csv(as.data.frame(significant_8cell_vs_32cell), file = "significant_8cell_vs_32cell.csv")
write.csv(as.data.frame(significant_16cell_vs_32cell), file = "significant_16cell_vs_32cell.csv")
# 上调和下调基因
upregulated_8cell_vs_16cell <- subset(significant_8cell_vs_16cell, log2FoldChange > 0)
downregulated_8cell_vs_16cell <- subset(significant_8cell_vs_16cell, log2FoldChange < 0)

upregulated_8cell_vs_32cell <- subset(significant_8cell_vs_32cell, log2FoldChange > 0)
downregulated_8cell_vs_32cell <- subset(significant_8cell_vs_32cell, log2FoldChange < 0)

upregulated_16cell_vs_32cell <- subset(significant_16cell_vs_32cell, log2FoldChange > 0)
downregulated_16cell_vs_32cell <- subset(significant_16cell_vs_32cell, log2FoldChange < 0)

# 保存上调和下调基因到文件
write.csv(as.data.frame(upregulated_8cell_vs_16cell), file = "upregulated_8cell_vs_16cell.csv")
write.csv(as.data.frame(downregulated_8cell_vs_16cell), file = "downregulated_8cell_vs_16cell.csv")

write.csv(as.data.frame(upregulated_8cell_vs_32cell), file = "upregulated_8cell_vs_32cell.csv")
write.csv(as.data.frame(downregulated_8cell_vs_32cell), file = "downregulated_8cell_vs_32cell.csv")

write.csv(as.data.frame(upregulated_16cell_vs_32cell), file = "upregulated_16cell_vs_32cell.csv")
write.csv(as.data.frame(downregulated_16cell_vs_32cell), file = "downregulated_16cell_vs_32cell.csv")

# 定义绘制火山图的函数
plot_volcano <- function(results, title) {
  results$log2FoldChange <- as.numeric(results$log2FoldChange)
  results$padj <- as.numeric(results$padj)
  results$significance <- ifelse(results$padj < 0.05 & abs(results$log2FoldChange) > 1, "Significant", "Not Significant")
  results$significance[results$padj < 0.05 & results$log2FoldChange > 1] <- "Right Significant"
  
  # 添加标签信息
  results$label <- ifelse(results$padj < 0.05 & abs(results$log2FoldChange) > 2, rownames(results), NA)
  
  ggplot(results, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.4, size = 1.75) +
    scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red", "Right Significant" = "blue")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    ggtitle(title) +
    xlab("Log2 Fold Change") +
    ylab("-Log10 Adjusted P-value") +
    theme_minimal(base_family = "Helvetica", base_size = 12) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "top"
    ) +
    geom_text_repel(aes(label = label), size = 3, box.padding = 0.3, max.overlaps = 15)
}
# 绘制并保存火山图，保存为PDF文件
volcano_8cell_vs_16cell <- plot_volcano(as.data.frame(results_8cell_vs_16cell), "8-cell vs 16-cell")
ggsave("volcano_8cell_vs_16cell.pdf", plot = volcano_8cell_vs_16cell, width = 8, height = 6)

volcano_8cell_vs_32cell <- plot_volcano(as.data.frame(results_8cell_vs_32cell), "8-cell vs 32-cell")
ggsave("volcano_8cell_vs_32cell.pdf", plot = volcano_8cell_vs_32cell, width = 8, height = 6)

volcano_16cell_vs_32cell <- plot_volcano(as.data.frame(results_16cell_vs_32cell), "16-cell vs 32-cell")
ggsave("volcano_16cell_vs_32cell.pdf", plot = volcano_16cell_vs_32cell, width = 8, height = 6)
