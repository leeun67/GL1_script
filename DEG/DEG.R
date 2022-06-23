library(DESeq2)
library(ggplot2)

# 1文件读取
mycounts<- read.csv("gene_count_matrix.csv", row.names = 1, stringsAsFactors = F)
colData <- read.csv( "sample.csv",  stringsAsFactors = T)

#2 Deseq运行
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)

#分组
res <- results(dds, contrast=c("condition", "Red","Gr"))

dds <- estimateSizeFactors(dds)
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1), colData=colData(dds))

#3PCA作图
tiff(file="gr_red.tiff",width=40,height=40,units="cm",compression="lzw",bg="white",res=300)
plotPCA(DESeqTransform(se)) + geom_text(aes(label=name),vjust=2)
dev.off()


#4 筛选差异表达基因
DEG <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1.5)
write.csv(DEG,file= "DEG_genes.cvs",row.names = T)
