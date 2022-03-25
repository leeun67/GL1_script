import sys
filename=sys.argv[1]
outfile=sys.argv[2]


f=open(filename,'r')
fi=f.read()
f.close()
fi=fi.split("\n")
del fi[-1]

vcf=''
for line in fi:
	if line[0]=="#":
		vcf+=line+"\n"
	else:
		tmp_line=line.split('\t')
		if len(tmp_line[4])==1:
			vcf+=line+"\n"


aa=open(outfile,'w')
aa.write(vcf)
aa.close()
