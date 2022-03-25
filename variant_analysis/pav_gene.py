import os
import sys


pav_merge_win_file=sys.argv[1]
gff3=sys.argv[2]



os.system("bedtools sort -i {gff3} |cut -f 1,3,4,5,9|grep -P '\tCDS\t'|tr ';' '\t'|cut -f 1,3,4,6|sed 's/Parent=//g' |sed 's/Name=//g'   >2tmp_cds.bed".format(gff3=gff3))


f=open('2tmp_cds.bed','r')
f=f.read()
f=f.split('\n')
del f[-1]

ff=''
for line in f:
	line = line.split('\t')
	line[3] = line[3].split('.')[0]
	line = '\t'.join(k for k in line) 
	ff += line+'\n' 
	
f=open('tmp_cds.bed','w')
f.write(ff)
f.close()

f=open('tmp_cds.bed','r')
f=f.read()
f=f.split('\n')
del f[-1]


gene_list=[]
for line in f:
	line=line.split('\t')
	gene_list.append(line[3])

gene_list=list(set(gene_list))



cds_list=[]
last_chr=''
last_sta=''
last_stop=''
last_name=''


for line in f:
	sin_list=[]
	cds_po=[]
	line=line.split('\t')
	name=line[3]
	chr=line[0]
	sta=int(line[1])
	stop=int(line[2])
	if name==last_name:
		if sta < last_stop or sta==last_stop:
			if sta>=last_sta and stop >last_stop:
				cds_list[-1][-1][1]=stop
				last_stop=stop
			elif sta <last_sta and stop<=last_stop:
				cds_list[-1][-1][0]=sta
				last_sta=sta
			elif sta<last_sta and stop>last_stop:
				cds_list[-1][-1][1]=stop
				cds_list[-1][-1][0]=sta
				last_sta=sta
				last_stop=stop
		else:
			cds_po.append(sta)
			cds_po.append(stop)
			cds_list[-1].append(cds_po)
			last_sta=sta
			last_stop=stop
	else:
		sin_list.append(name)
		sin_list.append(chr)		
		cds_po.append(sta)
		cds_po.append(stop)
		sin_list.append(cds_po)
		last_name=name
		last_sta=sta
		last_stop=stop
		cds_list.append(sin_list)	
	
#print(cds_list)
print(len(cds_list))	



noredun_cds_list=[]
cd_dict={}

for gene in gene_list:
	tmp_list=[]
	for cd in cds_list:
		if cd[0]==gene and len(tmp_list)<1:
			tmp_list.append(cd)
		elif  cd[0]==gene and len(tmp_list)>=1:
			for i in range(2,len(cd)):
				tmp_list[0].append(cd[i])
				tmp_list.append(1)
	num=0
	for i in range(2,len(tmp_list[0])):
		num+=tmp_list[0][i][1]-tmp_list[0][i][0]+1
	noredun_cds_list.append(num)
	noredun_cds_list.append(tmp_list[0])
	cd_dict[gene]=[]
	cd_dict[gene].append(num)
	cd_dict[gene].append(tmp_list[0])

print(len(noredun_cds_list))
#print(cd_dict)


pav_win=open(pav_merge_win_file,'r')
pav_win=pav_win.read()
pav_win=pav_win.split('\n')
del pav_win[-1]

chr_list=[]
for i in range(len(pav_win)):
        pav_win[i]=pav_win[i].split('\t')
        chr_list.append(pav_win[i][0])


chr_list=list(set(chr_list))


pav_dict={}
for chr in chr_list:
        pav_dict[chr]=[]


for line in pav_win:
	chr=line[0]
	tmp_list=[]
	t=[]
	tmp_list.append(int(line[1]))
	tmp_list.append(int(line[2]))
	#t=pav_dict[chr]
	#t.append(tmp_list)
	pav_dict[chr].append(tmp_list)


gene_list=list(cd_dict.values())

pav_gene=''
for gene in gene_list:
	pav_len=0
	chr=gene[1][1]
	for i in range(2,len(gene[1])):
		cd_win=gene[1][i]
		left=cd_win[0]
		right=cd_win[1]
		if chr in chr_list:
			for p_win in pav_dict[chr]:
				if left <= p_win[1] and left >= p_win[0]:
					if right <= p_win[1]:	
						pav_len+=right-left+1
					else:
						pav_len+=p_win[1]-left+1
				elif left < p_win[0]:
					if right >= p_win[0] and right <= p_win[1]:
						pav_len+=right-p_win[0]+1
					elif right >= p_win[1] :
						pav_len+=p_win[1]-p_win[0]+1
		if pav_len/int(gene[0]) >= 0.75:
			pav_gene+=gene[1][0]+'\t'+str(round(pav_len/int(gene[0]),2))+'\n'

f=open('pav.gene','w')
f.write(pav_gene)
f.close()
	
					




























