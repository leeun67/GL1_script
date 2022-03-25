import os
import sys
import subprocess

pep1=sys.argv[1]
pep2=sys.argv[2]
gff1=sys.argv[3]
gff2=sys.argv[4]
threads=sys.argv[5]

os.system('cat  {pep2}.pep >allpep'.format(pep1=pep1,pep2=pep2))

def genomeblast(pep1,pep2,threads):
	os.system('makeblastdb -in allpep -out all -parse_seqids -dbtype prot')
#	os.system('blastp -query {pep1}.pep -db all -outfmt 6 -out {pep1} -num_threads {threads} '.format(pep1=pep1,threads=threads))
	os.system('blastp -query {pep1}.pep -db all -outfmt 6 -out {pep1}_{pep2}.blast -num_threads {threads} -num_alignments 5 '.format(pep1=pep1,pep2=pep2,threads=threads))

genomeblast(pep1,pep2,threads)

os.system("cat {gff1} |awk '$3==\"gene\"'|cut -f 1,4,5,9|sed 's/ID=//g'|tr ';' '\t'|cut -f 1,2,3,4  >tt1.gff ".format(gff1=gff1))
os.system("cat {gff2} |awk '$3==\"gene\"'|cut -f 1,4,5,9|sed 's/ID=//g'|tr ';' '\t'|cut -f 1,2,3,4  >tt2.gff ".format(gff2=gff2))

with open ('tt1.gff','r') as f:
	tt1=f.readlines()

tmp1=''
for line in tt1:
	line=line.rstrip("\n").split("\t")
	tmp1+=line[0]+'\t'+line[3]+'\t'+line[1]+'\t'+line[2]+"\n"	
tmp1_gff=open('tmp1.gff','w')
tmp1_gff.write(tmp1)
tmp1_gff.close()
f.close()


with open ('tt2.gff','r') as f:
        tt2=f.readlines()

tmp2=''
for line in tt2:
        line=line.rstrip("\n").split("\t")
        tmp2+=line[0]+'\t'+line[3]+'\t'+line[1]+'\t'+line[2]+"\n"
tmp2_gff=open('tmp2.gff','w')
tmp2_gff.write(tmp2)
tmp2_gff.close()
f.close()


os.system("cat tmp1.gff tmp2.gff >{pep1}_{pep2}.gff".format(pep1=pep1,pep2=pep2))

os.system("/home/a2431559261/OLDISK/soft/MCScanX/MCScanX {pep1}_{pep2}".format(pep1=pep1,pep2=pep2))


chr1=''
with open("{pep1}.chr".format(pep1=pep1),"r") as f:
        pep1_chr=f.readlines()

for line in pep1_chr:
	line=line.rstrip("\n")
	chr1+=line+","

chr1=chr1.rstrip(",")
f.close()




chr2=''
with open("{pep2}.chr".format(pep2=pep2),"r") as f:
        pep2_chr=f.readlines()
        
for line in pep2_chr:                
	line=line.rstrip("\n")                
	chr2+=line+","
     
chr2=chr2.rstrip(",")
f.close()

mcsplot_PATH="/home/a2431559261/OLDISK/soft/MCScanX/downstream_analyses/"
circle=open("circos_{pep1}.ctl".format(pep1=pep1),"w")
circle.write("800\n")
circle.write(chr1)
circle.close()

circle=open("circos_{pep2}.ctl".format(pep2=pep2),"w")
circle.write("800\n")
circle.write(chr2)
circle.close()


dot=open("dot.ctl","w")
dimension="800\n800\n"
#dian.wrire("800     //dimension (in pixels) of x axis\n800     //dimension (in pixels) of y axis\n")
dot.write(dimension)
dot.write(chr1+"\n")
dot.write(chr2)
dot.close()


dual=open("dual.ctl","w")
dual.write("600\n800\n")
dual.write(chr1+"\n")
dual.write(chr2)
dual.close()


pwd=subprocess.getoutput('pwd') 
os.chdir('/home/a2431559261/OLDISK/soft/MCScanX/downstream_analyses')

##做图
os.system('java circle_plotter -g {pwd}/{pep1}_{pep2}.gff -s {pwd}/{pep1}_{pep2}.collinearity -c {pwd}/circos_{pep1}.ctl -o {pwd}/{pep1}.cir.PNG &> /dev/null'.format(pwd=pwd,pep1=pep1,pep2=pep2))

os.system('java circle_plotter -g {pwd}/{pep1}_{pep2}.gff -s {pwd}/{pep1}_{pep2}.collinearity -c {pwd}/circos_{pep2}.ctl -o {pwd}/{pep2}.cir.PNG &> /dev/null'.format(pwd=pwd,pep1=pep1,pep2=pep2))

os.system('java dot_plotter -g {pwd}/{pep1}_{pep2}.gff -s {pwd}/{pep1}_{pep2}.collinearity -c {pwd}/dot.ctl -o {pwd}/dot.PNG &> /dev/null'.format(pwd=pwd,pep1=pep1,pep2=pep2))
os.system('java dual_synteny_plotter -g {pwd}/{pep1}_{pep2}.gff -s {pwd}/{pep1}_{pep2}.collinearity -c {pwd}/dual.ctl -o {pwd}/dual.PNG &> /dev/null'.format(pwd=pwd,pep1=pep1,pep2=pep2))







