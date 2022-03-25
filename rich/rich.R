cp ~/OLDISK/script/gwas/mrna.bed .
cp ~/OLDISK/script/gwas/select_gene.py .

cp   /home/a2431559261/OLDISK/dataset/genome/fl/kegg/kegg/kegg/gene2k .
cp  /home/a2431559261/OLDISK/dataset/genome/fl/kegg/kegg/kegg/kegg_pathway.txt .

python select_gene.py fst.sweep


2，目的基因集
#包括  Gene    logFC   K_num
cat gene2k |grep --file gene_select.txt -w|awk 'BEGIN{OFS=FS="\t"}{print $1,-1,$2}' >ge

#cat gene2k |grep --file gene_select.txt -w|awk 'BEGIN{OFS=FS="\t"}{print $1,$2}' >ge

#l gr_red_DEG_genes.cvs |tr ',' '\t' |sed 's/"//g'|cut -f 1,3|tail -n +2 > gr_red_log2
#awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$1]=$2}ARGIND==2{print $1,a[$1],$2}'  gr_red_log2 ge >gr_red.ge2

sed -i '1 i Gene\tlogFC\tK_num' ge
3，R 进行kegg富集
conda activate rep


library(clusterProfiler)

kegg_anno=read.delim('kegg_pathway.txt',colClasses='character')
gene_select=read.delim('ge',stringsAsFactors=FALSE)$Gene

kk=enricher(gene=gene_select,
TERM2GENE=kegg_anno[c('Pathway','Gene')],
TERM2NAME=kegg_anno[c('Pathway','Description')],
pvalueCutoff=0.5,
qvalueCutoff=1,
maxGSSize=500)

write.table(kk,'kegg_rich.significant.txt',sep='\t',row.names=FALSE,quote=FALSE)



tiff(file="barplot_kegg.tiff",width=20,height=20,units="cm",compression="lzw",bg="white",res=300)
dotplot(kk,color="pvalue"，showCategory=10 )
dev.off()




barplot(kk)
dotplot(kk)
cnetplot(kk)   #网络图展示富集功能和基因的包含关系
emapplot(kk)  #网络图展示各富集功能之间共有基因关系
heatplot(kk)  #热图展示富集功能和基因的包含关系
write.table(kk,'kegg_rich.significant.txt',sep='\t',row.names=FALSE,quote=FALSE)

4，pathview 富集结果文件整理后用于后续pathview分析
python kegg_pathview_file_change.py  ###这个脚本用于整理文件

##R pathview

library(pathview)
rt=read.delim('gr_red.ge2',stringsAsFactors=FALSE)
geneFC=rt$logFC
gene=rt$K_num
names(geneFC)=gene
keggxls=read.table("gr_red_kegg_rich.significant.txt",sep="\t",header=T)
(for i in keggxls$ko00940){
     pv.out<-pathview(gene.data=geneFC,pathway.id=i,species="ko",out.suffix="pathview")
     }
     
     
 GO  富集分析
###复制背景数据文件过来 ###！！！！！！！
cp /home/a2431559261/OLDISK/dataset/genome/fl/go/gene/gene2go.txt .
cp /home/a2431559261/OLDISK/dataset/genome/fl/go/gene/go_annot.txt .
cp /home/a2431559261/OLDISK/dataset/genome/fl/go/gene/go_class.txt .




###挑出目的基因
l red_ye2_DEG_genes.cvs |tr ',' '\t' |sed 's/"//g'|cut -f 1|tail -n +2|sed 's/\..*//g'|sort|uniq >red_ye2_clean.gene

cat gene2go.txt |grep --file gene_select.txt -w | sed '1 i Gene\tGOID' >select_gene.txt   #######！！！！！


###R语言go 富集
library(clusterProfiler)

go_anno <- read.delim('gene2go.txt',header=FALSE,stringsAsFactors=FALSE)
names(go_anno) <- c('Gene', 'GOID')

go_class <- read.delim('go_class.txt',header= FALSE,stringsAsFactors=FALSE)
names(go_class) <- c('GOID','ONTOLOGY','TERM')

go_anno <- merge(go_anno,go_class,by='GOID',all.x=TRUE)


gene_select <- read.delim('select_gene.txt', stringsAsFactors =FALSE)$Gene

gg=enricher(gene=gene_select,
TERM2GENE=go_anno[c('GOID', 'Gene')],
TERM2NAME=go_anno[c('GOID', 'TERM')],  
pvalueCutoff=0.20,
qvalueCutoff=1,  
maxGSSize=500)  # maximal size of genes annotated for testing

write.table(gg, 'red_ye2_go_tmp.txt', sep = '\t',row.names = FALSE,quote = FALSE)
#dotplot(gg,showCategory=10, color="pvalue")


tiff(file="red_ye2_go_dotplot.tiff",width=40,height=40,units="cm",compression="lzw",bg="white",res=300)
dotplot(gg,showCategory=10, color="pvalue")
dev.off()
    
 
