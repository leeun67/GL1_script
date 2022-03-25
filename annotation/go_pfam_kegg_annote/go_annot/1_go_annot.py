
import os
import sys

f=open('tmp','r')
f=f.read()
f=f.split('\n')
del f[-1]

f2=''
for line in f:
	line=line.split('\t')
	length=len(line)
	for i in range(1,length):
		f2 += line[0]+'\t'+line[i]+'\n'

f2=f2.rstrip('\n')
#print(f2)	
os.system('rm tmp')

q=open('tmp2','w')
q.write(f2)
q.close()

os.system('sort -k1,2 tmp2 |uniq >tmp3')
os.system('rm tmp2')


#go_2_gene annot
os.system('cat tmp3|sort -k2 >tmp4')


f3=open('tmp4','r')
f3=f3.read()
f3=f3.split('\n')
del f3[-1]


num=0
go2ge=''
go_name=''
gene=''
for line in f3:
	line= line.split('\t')
	if line[1] != go_name:
		go2ge += go_name+'\t'+gene+'\t'+str(num)+'\n'
		
		go_name=line[1]
		gene=''
		gene += line[0]
		num=1

	elif line[1] == go_name:
		gene += ','+line[0]
		num += 1

go2ge += go_name+'\t'+gene+str(num)
q=go2ge.lstrip('\t').lstrip('0').lstrip('\n')
#print(q)

f=open('ge2ge.tmp','w')
f.write(q)
f.close()


#gene2go.txt and others
os.system('cat tmp3 |sed "1 i Gene\tGOID" >gene2go.txt')
os.system('bash /home/a2431559261/OLDISK/script/rich/go_annot/2_go2ge.sh')



	










