


##1 加载top fst位点文件
fst=open('maxfst_site','r')
fst=fst.read()
fst=fst.split('\n')

del fst[-1]
fstt=[]
for line in fst:
        line=line.split('\t')
        line.append(int(line[1])+99999)
        line[0]=int(line[0])
        line[1]=int(line[1])
        #line[3]=tmp
        fstt.append(line)


## 2加载mrna.bed文件
genebed=open('mrna.bed','r')
genebed=genebed.read()
genebed=genebed.split('\n')

del genebed[-1]
gbed=[]

for line in genebed:
        line=line.split('\t')
        line[0]=int(line[0])
        line[1]=int(line[1])
	line[2]=int(line[2])
        #line[3]=tmp
        gbed.append(line)


#默认计算25条染色体
gene=''
for chr in range(1,26):
        t_fst=[]
        t_gb=[]
        for line in fstt:
                if line[0]==chr:
                        t_fst.append(line)
	print(t_fst)
        for line in gbed:
                if line[0]==chr:
                        t_gb.append(line)
	print(t_gb)
        for line_gb in t_gb:
                for line in t_fst:
                        #print(line_ta1)
                        if line_gb[1]<line[1] and line_gb[2]>line[1]:
				gene+=line_gb[3]+'\n'
				continue
			elif line_gb[1]<line[3] and line_gb[2]>line[3]:
				gene+=line_gb[3]+'\n'
				continue
			elif line_gb[1]>line[1] and line_gb[2]<line[3]:
				gene+=line_gb[3]+'\n'
				continue
			elif line_gb[1]==line[1] and line_gb[2]>line[1]:
				gene+=line_gb[3]+'\n'
				continue
			elif line_gb[2]==line[3] and line_gb[1]<line[3]:
				gene+=line_gb[3]+'\n'
				continue


f=open('gene','w')
f.write(gene)
f.close()

import os
os.system('cat gene |sort|uniq|sed "/^$/d"|sort -k 1n -k2n >gene_select.txt')
os.system('rm gene')

