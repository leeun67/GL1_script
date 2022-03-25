# -*- coding:utf-8 -*-
import sys
##  用于gff文件中总共的序列数量统计，提供的gff文件需要预排序，先排染色体，再排第四列
##  传递两个参数 第一个 gff文件 第二个 输出文件名





modeler_repeatgff=sys.argv[1]
gff_outname=sys.argv[2]


def replace_last(source_string, replace_what, replace_with):
    head, _sep, tail = source_string.rpartition(replace_what)
    return head + replace_with + tail



f=open(modeler_repeatgff,"r")
a=f.read()
f.close()
a=a.split("\n")
del a[-1]

chrr=""
sta_gff=""
for line in a:
	if line[0]=="#":
		continue
	else:
		line=line.split("\t")
		if line[0]==chrr:
			#sta_gff=sta_gff.split("\t")
			if int(line[3])<int(right) or int(line[3])==int(right):
				if  int(line[4])>int(right):
					sta_gff=replace_last(sta_gff,right,line[4])
#					sta_gff='\t'.join(str(x) for x in sta_gff)
					right=line[4]
				else:
					pass
					#sta_gff='\t'.join(str(x) for x in sta_gff)
				
			else:
				#sta_gff='\t'.join(str(x) for x in sta_gff)
				sta_gff+=chrr+"\t"+line[3]+"\t"+line[4]+"\n"
				right=line[4]

		else:
			chrr=""	
			chrr=line[0]
			sta_gff+=chrr+"\t"+line[3]+"\t"+line[4]+"\n"
			right=line[4]


a=sta_gff
repeat_len=0
sta_gff=sta_gff.split("\n")
del sta_gff[-1]
for line in sta_gff:
	line=line.split("\t")
	b=int(line[2])
	c=int(line[1])
	repeat_len=repeat_len+b-c+1
	c=0
	b=0

a=a+str(repeat_len)+"\n"
f=open(gff_outname,"w")
f.write(a)
f.close()
print(repeat_len)






