######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("colorspace")
#install.packages("stringi")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")

setwd("/home/a2431559261/OLDISK/rich/kegg/cluster")
library("clusterProfiler")
rt=read.table("id.txt",sep="\t",header=T,check.names=F)
rt=rt[is.na(rt[,"entrezID"])==F,]

geneFC=rt$logFC
gene=rt$entrezID
names(geneFC)=gene

#kegg富集分析
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)

#柱状图
tiff(file="barplot.tiff",width = 20,height = 20,units ="cm",compression="lzw",bg="white",res=300)
barplot(kk, drop = TRUE, showCategory = 12)
dev.off()

#点图
tiff(file="dotplot.tiff",width = 20,height = 20,units ="cm",compression="lzw",bg="white",res=300)
dotplot(kk)
dev.off()

#通路图
#library("pathview")
#keggxls=read.table("KEGG.txt",sep="\t",header=T)
#for(i in keggxls$ID){
#  pv.out <- pathview(gene.data = geneFC, pathway.id = i, species = "hsa", out.suffix = "pathview")
#}


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
