
import sys
import re
from operator import itemgetter
import re

# awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"~":$0 }' youzong.pep |grep -v "gene" |sed 's/\./%/1'|sed 's/\.//g' |sed 's/%/\./1'|tr "~" "\n" >youzong2.pep  

#awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"~":$0 }' youzong.cds |grep -v "gene" >max_youzong.cds

# cat wuyou.cds |sed 's/\.[0-9]//1' |sed 's/\.//g' |tr "agct" "ATCT">max_wuyou.cds

#seqtk subseq fl.cds long_mranID.txt |sed 's/\.[0-9]//1' |sed 's/\.//g' >max_fl.cds
#seqtk subseq fl.cds long_mranID.txt |sed 's/\.[0-9]/\t/1' |cut -f1|sed 's/\.//g'| tr "agct" "ATCT" >max_fl.cds



#打开文件
fasta=open("gou.tr","r")
seq={}
#转化为字典的形式
for line in fasta.readlines():
    content=line.strip()
    if content.startswith(">"):
        #name=content[1:]
        #seq_name=re.sub(" .*$","",name) #修改序列名
        nameS=content[1:]
        name=re.sub(" .*$","",nameS)
        seq[name]=''
    else:
        seq[name]+=content

print(seq.keys())

fasta.close()

#转换为转录本名，转录本长度，基因名的三列表
table=""
for i in seq.keys():
    table+=i+"\t"+str(len(seq[i]))+"\t"+i.split(".")[0]+"\n"

#排序
#生成列表
pre_sort=table.split("\n")
del pre_sort[-1]
#生成嵌套式列表
pre_sort2=[]
for line in pre_sort:
    col=line.split("\t")
    col[1]=int(col[1])
    pre_sort2.append(col)

table_sorted=sorted(pre_sort2,key=itemgetter(2,1),reverse=True)

#得到最长转录本
longest_Mrna=""
dict={}
for line  in table_sorted:
    if line[2] not in dict:
        dict[line[2]]=line[0]
    else:
        pass

a=""
longest_Mrna=list(dict.values())
for line in longest_Mrna:
    a+=line +"\n"

filename="id.txt"
f=open(filename,"w")
f.write(a)
f.close






