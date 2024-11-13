rm(list=ls())

#install.packages("clusterProfiler")
library(clusterProfiler)
#install.packages("AnnotationHub")
library(AnnotationDbi)
library(ggplot2)
library(ggpubr)
library(aPEAR)
library(GOSemSim)
library(ggsci)
library(enrichplot)#用于可视化的包
library(AnnotationHub)
library(dplyr)
library(writexl)

###############################################################################################################
f=getwd()
# Read genelist
setwd(f)
###############################件####################################################################
#gene=read.table("gene.txt",header=FALSE) #单列基因名文
sheep <- loadDb(file='org.Ovis_aries.eg.sqlite')
go_organism <- sheep # org.Hs.eg.db / org.Mm.eg.db / org.Gg.eg.db
kegg_organism <- 'oas' # hsa / mmu / gga
##########################################################################
data <- read.table("community_1_nodes.txt",header=FALSE) #单列基因名文件
data$V1 <- as.character(data$V1) #需要character格式，然后进行ID转化

#################################GO分类#############################################
test = bitr(data$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb= go_organism) #将SYMBOL格式转为ENSEMBL和ENTERZID格式 
test
ego_ALL <- enrichGO(gene = test$ENTREZID, 
                    OrgDb = go_organism, #没有organism="human"，改为OrgDb=org.Hs.eg.db
                    #keytype = 'ENSEMBL',
                    ont = "ALL", #也可以是 CC  BP  MF中的一种
                    pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                    pvalueCutoff = 1, #P值会过滤掉很多，可以全部输出
                    qvalueCutoff = 1,
                    readable = TRUE
                    ) #Gene ID 转成gene Symbol ，易读
head(ego_ALL,2)
kk_df <- enrichKEGG(test$ENTREZID,
                    organism = kegg_organism, 
                    keyType = 'kegg', # one of "kegg", 'ncbi-geneid', 'ncbi-proteinid' and 'uniprot'
                    pvalueCutoff = 1,pAdjustMethod = 'BH',
                    qvalueCutoff = 1,
                
)
head(kk_df,2)

#############################保存结果###########################################################################################
res.go = DOSE::setReadable(ego_ALL, OrgDb = go_organism, keyType = 'ENTREZID') # 按需替换ENTREZID为SYMBOL
res.kegg = DOSE::setReadable(kk_df, OrgDb = go_organism, keyType = 'ENTREZID') # 按需替换ENTREZID为SYMBOL
writexl::write_xlsx(list('GO注释结果' = res.go@result, 'KEGG注释结果' = res.kegg@result), 'community_1.xlsx')






# 读取基因表达数据
genes_data <- read.table("community_1.logfc.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(genes_data) <- c("gene", "log2FC")

# 假设你已经进行了富集分析，并得到了 enrichResult 对象
# result <- enrichGO(...) # 进行富集分析

# 获取富集结果
enrich_results <- ego_ALL@result
head(enrich_results)

# 统计每个通路中上调和下调基因的个数
enrichment_summary <- enrich_results %>%
  rowwise() %>%
  mutate(
    gene_list = list(strsplit(geneID, "/")[[1]]),  # 分割基因ID
    up_gene = sum(genes_data$gene %in% gene_list & genes_data$log2FC > 0),  # 上调基因数量
    down_gene = sum(genes_data$gene %in% gene_list & genes_data$log2FC < 0),  # 下调基因数量
    all = as.numeric(sub(".*/", "", GeneRatio))  # 获取背景基因数量
  ) %>%
  select(ID, ontology = ONTOLOGY, all, pvalue, up_gene, down_gene)  # 选择并重命名列

# 查看统计结果
print(enrichment_summary)

# 保存为 XLSX 格式
write_xlsx(enrichment_summary, "16vs32.draw.xlsx")






#####################################################################################################################
# Treeplot
enrichres2 <- pairwise_termsim(ego_ALL) # calculate pairwise similarities of the enriched terms using Jaccard’s similarity index
treeplot(enrichres2,color="pvalue")
enrichres3 <- pairwise_termsim(kk_df) # calculate pairwise similarities of the enriched terms using Jaccard’s similarity index
treeplot(enrichres3,color="pvalue")
#######################################################################################################################################
# library("DOSE")
# dotplot(ego_ALL, 
#         x = "GeneRatio", color = "pvalue", 
#         showCategory =6, #只显示前10
#         split="ONTOLOGY") + #以ONTOLOGY类型分开
#   facet_grid(ONTOLOGY~., scale='free' )
# 
# dotplot(kk_df,x = "GeneRatio", color = "p.adj", 
#         showCategory =20, #只显示前10
# )
# 
# barplot(ego_ALL, x = "GeneRatio", color = "p.adj", #默认参数（x和color可以根据eG里面的内容更改）
#         showCategory =6, #只显示前10
#         split="ONTOLOGY") + #以ONTOLOGY类型分开
#   facet_grid(ONTOLOGY~., scale='free') #以ONTOLOGY类型分开绘图
# barplot(kk_df, x = "GeneRatio", color = "p.adj", #默认参数（x和color可以根据eG里面的内容更改）
#         showCategory =15) #以ONTOLOGY类型分开NTOLOGY类型分开绘图

##################保存结果#################################




#http://127.0.0.1:26009/graphics/6f4ebacd-1c1e-4eee-844c-075b856f2a8c.png
###############基因与通路关系######################################

cnetplot(ego_ALL,colorEdge = TRUE,showCategory = 10, circular= T,categorySize="p.adjust",color.params = list(foldChange =data$V1))
cnetplot(kk_df,colorEdge = TRUE,showCategory = 10, circular= T,categorySize="pvalue",color.params = list(foldChange =data$V1))




# library(aPEAR)
# top10_ego <- head(ego_ALL, 100)
# # 对结果按p值进行排序，并选择前10个通路
# enrichmentNetwork(top10_ego)
# top10_ego <- head(kk_df, 50)
# enrichmentNetwork(top10_ego)
# emapplot(ego_ALL)
# upsetplot(ego_ALL)
# 
# heatplot(ego_ALL) #热图展示富集功能和基因的包含关系


