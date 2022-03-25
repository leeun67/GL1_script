import sys


snp=sys.argv[1]
snp_vcf=sys.argv[2]
indel=sys.argv[3]
indel_vcf=sys.argv[4]



ge_snp=open(snp,'r')
ge_snp=ge_snp.read()
ge_snp=ge_snp.split('\n')
del ge_snp[-1]

genome_snp={}
chr_list=[]
for line in ge_snp:
	line=line.split('\t')
	chr_combi='~'.join(k for k in line)
	if line[0] not in chr_list:
		chr_list.append(line[0])
		genome_snp[line[0]]=[]
		genome_snp[line[0]].append(chr_combi)
	else:
		genome_snp[line[0]].append(chr_combi)






map_snp=open(snp_vcf,'r')
map_snp=map_snp.read()
map_snp=map_snp.split('\n')
del map_snp[-1]

fi_snp=''
for line in map_snp:
	if line[0]=='#':
		fi_snp+=line+'\n'
	else:
		line2=line.split('\t')
		chr=line2[0]
		if chr in chr_list:
			combi=line2[0]+'~'+line2[1]
			if combi in genome_snp[chr] and len(line2[3])<100 and len(line2[4])<100:
				genome_snp[chr].remove(combi)
				fi_snp+=line+'\n'

f=open('snp_fi.vcf','w')
f.write(fi_snp)
f.close()




ge_indel=open(indel,'r')
ge_indel=ge_indel.read()
ge_indel=ge_indel.split('\n')
del ge_indel[-1]

genome_indel={}
chr_list=[]
for line in ge_indel:
	line=line.split('\t')
	chr_combi='~'.join(k for k in line)
	if line[0] not in chr_list:
		chr_list.append(line[0])
		genome_indel[line[0]]=[]
		genome_indel[line[0]].append(chr_combi)
	else:
		genome_indel[line[0]].append(chr_combi)



map_indel=open(indel_vcf,'r')
map_indel=map_indel.read()
map_indel=map_indel.split('\n')
del map_indel[-1]

fi_indel=''

for line in map_indel:
	if line[0]=='#':
		fi_indel+=line+'\n'
	else:
		line2=line.split('\t')
		chr=line2[0]
		if chr in chr_list:
			combi=line2[0]+'~'+line2[1]
			if combi in genome_indel[chr] and len(line2[3])<100 and len(line2[4])<100:	
				genome_indel[chr].remove(combi)
				fi_indel+=line+'\n'

f=open('indel_fi.vcf','w')
f.write(fi_indel)
f.close()
