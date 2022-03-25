import sys

file_prefix=sys.argv[1]
chr_prefix=sys.argv[2]
chr_num=sys.argv[3]


def mcs_chr(file_prefix,chr_prefix,chr_num):
	chr_num=int(chr_num)
	f=open("{file_prefix}.chr".format(file_prefix=file_prefix),"w")
	for i in range(chr_num):
		i=i+1
		i=str(i)
		f.write(chr_prefix+i+"\n")

	f.close()

mcs_chr(file_prefix,chr_prefix,chr_num)
