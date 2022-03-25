import sys
import os
### trf的.dat文件 整理正gff文件，有待改善 
###提供参数   第一个 。dat文件预整理  一行染色体名跟着一堆 重复序列位置 第二个参数输出文件名


dat_file=sys.argv[1]
outfile_name=sys.argv[2]

os.system("tail {dat_file} -n +8|sed 's/Parameters:.*//g'|sed 's/Sequence: //g'|sed '/^$/d'|tr ' ' '\t'|cut -f 1,2,3,4,8,14 >tmp".format(dat_file=dat_file))


ff=open("tmp","r")
aa=ff.read()
aa=aa.split("\n")

del aa[-1]
trf_gff=""
charac_gff=""
tmp_gff="\n"
num_temgff=""
center=""
for line in aa:
	if line[0:3]=="gro" or line[0:3]=="tig":
		tmp_gff=tmp_gff.split("\n")
		if len(tmp_gff)>10:
			charac_gff+=tmp_gff[-11]+"\n"+tmp_gff[-10]+"\n"+tmp_gff[-9]+"\n"+tmp_gff[-8]+"\n"+tmp_gff[-7]+"\n"+tmp_gff[-6]+"\n"+tmp_gff[-5]+"\n"+tmp_gff[-4]+"\n"+tmp_gff[-3]+"\n"+tmp_gff[-2]+"\n"
		chr_name=""
		chr_name=line
		tmp_gff=""
	else:
		line=line.split("\t")
		a=int(line[2])
		b=float(line[3])
		leng=a*b
		tt=chr_name+"\t"+"trf"+"\t"+"tandemrepeat"+"\t"+line[0]+"\t"+line[1]+"\t"+line[4]+"\t"+"."+"\t"+"."+"\t"+"ID="+chr_name+"_"+line[0]+"_"+line[1]+";"+"NOTE="+line[2]+"_"+line[3]+"_"+str(leng)+";motif="+line[5]+";\n"
		trf_gff+=chr_name+"\t"+"trf"+"\t"+"tandemrepeat"+"\t"+line[0]+"\t"+line[1]+"\t"+line[4]+"\t"+"."+"\t"+"."+"\t"+"ID="+chr_name+"_"+line[0]+"_"+line[1]+";"+"NOTE="+line[2]+"_"+line[3]+"_"+str(leng)+";motif="+line[5]+";\n"

#                tt=chr_name+"\t"+"trf"+"\t"+"tandemrepeat"+"\t"+line[0]+"\t"+line[1]+"\t"+line[4]+"\t"+"."+"\t"+line[2]+"_"+str(leng)+"\t"+"ID="+chr_name+"_"+line[0][1:3]+";\n"#+"NOTE="+line[2]+"_"+line[3]+"_"+line[4]+";\n"#motif="+line[5]+";\n"
#                trf_gff+=chr_name+"\t"+"trf"+"\t"+"tandemrepeat"+"\t"+line[0]+"\t"+line[1]+"\t"+line[4]+"\t"+"."+"\t"+line[2]+"_"+str(leng)+"\t"+"ID="+chr_name+"_"+line[0][1:3]+";\n"#+"NOTE="+line[2]+"_"+line[3]+"_"+line[4]+";\n"#motif="+line[5]+";\n"

	#	print(trf_gff)
		a=0
		b=0
		leng=0
		if float(line[4])>20000:
			charac_gff+=tt
			center+=tt

		tmp_gff+=tt
		num_temgff=tmp_gff.split("\n")
		if len(num_temgff)<100:
			charac_gff+=tt



ff.close()

a=open(outfile_name,"w")
a.write(trf_gff)
a.close()

a=open("charac.gff","w")
a.write(charac_gff)
a.close()

a=open("center_candi.gff","w")
a.write(center)
a.close()


