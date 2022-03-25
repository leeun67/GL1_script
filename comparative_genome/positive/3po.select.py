f=open('nnn2','r')

f=f.read()
f=f.split('\n')
del f[-1]

import sys
import os

nnn3=''
os.system('rm nnn3')

for line in f:
	line=line.split('\t')
	a=''
	b=''
	a=line[0]
	b=abs(float(line[1]))

	os.system('echo {a} >>nnn3'.format(a=a))
	os.system('chi2 1.5 {b} >>nnn3'.format(b=b))


	
