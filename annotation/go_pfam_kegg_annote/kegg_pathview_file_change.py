gene2k=open("gene2k","r")
gene2k=gene2k.read()
gene2k=gene2k.split("\n")
del gene2k[-1]

rich=open("kegg_rich.significant.txt","r")
rich=rich.read()
rich=rich.split("\n")
del rich[-1]

g_k={}
for line in gene2k:
     tmp=line.split("\t")
     g_k[tmp[0]]=tmp[1]

rich_head=rich[0]
del rich[0]

ri=""
for line in rich:
     tmp=line.split("\t") 
     #del tmp[1]
     tmp2=tmp[7]
     tmp3=tmp2.split("/")
     for i in range(len(tmp3)):
          tmp3[i]=g_k[tmp3[i]]
     tmp2="~".join(k for k in tmp3)
     tmp[7]=tmp2
     tmp="\t".join(k for k in tmp)
     ri+=tmp+"\n"

ri_fi=rich_head+"\n"+ri

filename="kegg_rich.txt"
m=open(filename,'w')
m.write(ri_fi)

#import os
#sed -i "s/ /\t/1" kegg_rich.txt
#sed -i "s/~/\//g" kegg_rich.txt


##pathview 注释
#geneFC=rt$logFC
#gene=rt$K_num
#names(geneFC)=gene

#keggxls=read.table("kegg_rich.txt",sep="\t",header=T)
#for(i in keggxls$ID){
 #    pv.out <- pathview( gene.data=geneFC,pathway.id = i, species = "ko", out.suffix = "pathview")
  #   }
