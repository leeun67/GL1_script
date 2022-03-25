f=open('a','r')

a=f.read()
a=a.split('\n')
del a[-1]

for i in a:
	msa=''
	f=open(i,'r')
	f=f.read()
	f=f.split('\n')

	length=len(f[2].split('\t')[1])
	a=length%3
	length=length-a
	
	f[0]='13'+'   '+str(length)
	
	for num in range(14):
		if num==0:
			continue
		else:
			f[num]=f[num].split('\t')
			#q=''
			#q=f[num][1][:length]+'TGA'
			#f[num][1]=q
			f[num][1]=f[num][1][:length]			
			f[num]='   '.join(k for k in f[num])

	fi='\n'.join(k for k in f)
	b=open('{i}.phy'.format(i=i),'w')
	b.write(fi)
	b.close()


