import sys

file=sys.argv[1]
out_conserve=sys.argv[2]

seq_all=sys.argv[3]

#spec_name


f=open(file,'r')
f=f.read()
f=f.split('\n')
del f[-1]

sam_0_16=''
for line in f:
	line=line.split('\t')
	if int(line[1])==0 or int(line[1])==16:
		sam_0_16+='\t'.join(k for k in line)+'\n'

#print(sam_0_16)
#print(f)		



sam_0_16=sam_0_16.split('\n')
del sam_0_16[-1]



out=''
exclu_seq=[]
for line in sam_0_16:
	line=line.split('\t')
	#print(line)
	num=0
	for i in range(2,len(line)-1):	
		if line[i][-1]=='M':
			line[i]=line[i].rstrip('M')
			n=0
			n=int(line[i])
			num+=n
			#print(num)
	if num > 125 :
		out+=line[0]+'\t'+line[1]+'\t'+str(num)+'\n'
		exclu_seq.append(line[0])
#print(out)	
f=open(out_conserve,'w')
f.write(out)
f.close()


f=open(seq_all,'r')
f=f.read()
f=f.split('\n')
del f[-1]


print(len(exclu_seq))


term_seq=''

#for line in f:
#	if line  in exclu_seq:
#		exclu_seq.remove(line)
#	else:
		#print(line)
#		term_seq+=line+'\n'



#chr_list=open('li','r')
#chr_list=chr_list.read()
#chr_list=chr_list.split('\n')
#del chr_list[-1]


allseq={}
conserve_seq={}

c1_list=[]
c2_list=[]

for line in f:
	chr=line.split('~')[0]
	if chr not in c1_list:
		allseq[chr]=[]
		c1_list.append(chr)
		allseq[chr].append(line)
	else:
		allseq[chr].append(line)
	

for line in exclu_seq:
	chr=line.split('~')[0]
	if chr not in c2_list:
		conserve_seq[chr]=[]
		c2_list.append(chr)
		conserve_seq[chr].append(line)
	else:
		conserve_seq[chr].append(line)






#print(conserve_seq)
for chr in c1_list:
	for line in allseq[chr]:
		if chr not in c2_list:
			term_seq+=line+'\n'
		elif line in conserve_seq[chr]:
			conserve_seq[chr].remove(line)
		else:
			term_seq+=line+'\n'
#			print(line)
















#print(f)
f=open('specific_seq.txt','w')
f.write(term_seq)
f.close()



f=open('specific_seq.txt','r')
f=f.read()
f=f.split('\n')
del f[-1]

q=[]
for line in f:
	line=line.split('~')
	line[1]=int(line[1])
	line[2]=int(line[2])
	q.append(line)


#print(len(q))
o1=str(len(q))


fi=[]
last_chr='begin'
last_right=''
last_left=''
i=-1


for line in q:
	chr=''
	left=''
	right=''
	chr=str(line[0])
	left=line[1]
	right=line[2]
	#print(left)
	#print(right)
	if chr==last_chr:
		if left < last_right  or left==(last_right+1):
			fi[i][2]=right
			last_right=right
			#print(left)	
		else:
			fi.append(line)
			last_left=left
			last_right=right
			i+=1
	else:
		fi.append(line)
		last_left=left
		last_right=right
		last_chr=chr
		i+=1


#print(len(fi))
o2=str(len(fi))


final=''
sum=0

for line in fi:
	single_win_len=int(line[2])-int(line[1])+1
	sum=sum+single_win_len
	line.append(single_win_len)
	for i in range(0,4):
		line[i]=str(line[i])
	line='\t'.join(k for k in line)
	final+=line+'\n'
#	print(line)

o3=str(sum)
final=final
#print(o3)

f=open('merge_win','w')
f.write(final)
f.close()


sta='spe_win\t'+o1+'\n'+'merge_win\t'+o2+'\n'+'len_win'+'\t'+o3
f=open('sta','w')
f.write(sta)
f=f.close()




