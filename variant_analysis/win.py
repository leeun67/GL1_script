import sys
import os

fi=sys.argv[1]


f=open(fi,'r')
f=f.read()



windowSize=500
step=100
i=0
fasta=f.split('\n')

#del fasta[-1]


#fasta='>chr1\nqqqwwweee\n>chr2\nqqqwwweee\n'
#fasta=fasta.split('\n')
#print(fasta)

outfi=''
for line in fasta:
	currentWindow=0
	if line=='':
		break

	elif line[0]=='>':
		pre=line
	else:
		for left_i in range(0,len(line),step):
			if (left_i+windowSize) < len(line):
				currentWindow=pre+'~'+str(left_i+1)+'~'+str(left_i+windowSize)+'\n'+line[left_i:left_i+windowSize]+'\n'
				outfi+=currentWindow
			else:
				currentWindow=pre+'~'+str(left_i+1)+'~'+str(len(line))+'\n'+line[left_i:]+'\n'
				outfi+=currentWindow
				break
#print(currentWindow)
print(outfi)
