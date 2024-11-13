## 读入featurecounts矩阵
expr_df<-read.table( "counts.txt",header=T,row.names=1,check.names=F, sep="\t")
dim(expr_df); names(expr_df)
head(expr_df[,1:7])

#提取基因信息,featurecounts前几列
featurecounts_meta <- expr_df[,1:5] ;
head(featurecounts_meta)
#提取counts
expr_df <- expr_df[,6:ncol(expr_df)]
##基因表达量之和>8
expr_df <- expr_df[rowsums(expr_df)>8,]
head(expr_df[,1:6])
## 保存counts矩阵
write.table(expr_df, "merged.counts.txt",quote=F, sep="\t", row.names=T, col.names=T )
prefix <-"merged_samples"#设置输出文件前缀名
#CPM计算
cpm<-t(t(expr_df)/colSums(expr_df)*1000000)
avg_cpm<- data.frame(avg_cpm=rowMeans(cpm))# 保存
write.table(avg_cpm, paste0(prefix," _avg_cpm.xls"),quote=F, sep="\t", row.names=T, col.names=T )
write.table(cpm, paste0(prefix,"_cpm.xls"),quote=F, sep="\t", row.names=T, col.names=T )
#----- TPM计算
#基因长度，目标基因的外显子长度之和除以1800，单位是Kb，不是bp
kb <-featurecounts_meta$Length / 1000
rpk<-expr_df /kb#每千碱基reads(“per million”scaling factor)长度标准化
tpm<-t(t(rpk)/colSums(rpk)*1000000) # 自万缩放因子(“per million”scaling factor )深度标准化
avg_tpm <- data.frame(avg_tpm=rowMeans(tpm))
# 保存
write.table(avg_tpm, paste0(prefix," avg_tpm,xls"),quote=F, sep="\t",row.names=T, col.names=T )
write.table(tpm, paste0(prefix," tpm.xls"),quote=F, sep="\t", row.names=T, col.names=T )
#----- FPKM计算------
fpkm <-t(t(rpk)/colSums(expr_df)*10^6)
head(fpkm[ ,1:3])
# 保存
write.table(fpkm,file= paste0(prefix, "_fpkm.xls"),quote=F, sep="\t", row.names=T, col.names=T )
#-----FPKM转化为TPM-----
fpkm_to_tpm =t(t(fpkm)/colSums(fpkm))*10^6
head(fpkm_to_tpm)
