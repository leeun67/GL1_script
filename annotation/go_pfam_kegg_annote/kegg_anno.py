###############################################################
## 输入文件一 ko orthology得到map 文件后删除空行和前面空格
import os
os.system("sed -i '/^$/d' bak")
os.system("sed -i 's/^  //' bak")


############################################################
##输入文件 ，同上
f=open("bak","r")
f1=f.read()
f1=f1.split("\n")
del f1[-1]
bak=""
for line in f1:
    if "ko:" not in line:
        b=line
    elif "ko:" in line:
        a=line.split(" ")
        a=a[0]
        a=a[3:10]
        bak+=a+"\t"+b+"\n"


bak=bak.split("\n")
del bak[-1]

bak2=""
for line in bak:
    a=line.split(" ")
    del a[-1]
    a=" ".join(k for k in a)
    bak2+=a+"\n"



filename="bak2"
m=open(filename,'w')
m.write(bak2)
m.close()
f.close()



##########################################################################
###输入文件2  gene k_num  ,这个kaas 运行完后需要整理一下，去除空格
f=open("ko","r")
f2=f.read()
f2=f2.split("\n")
f.close()

del f2[-1]



bak2=bak2.split("\n")
del bak2[-1]

bak2="\n".join(k for k in bak2)
filename="bak2"
m=open(filename,"w")
m.write(bak2)
m.close()

os.system("cat bak2 |cut -f 1|sort|uniq >k")
os.system("cat ko |grep --file k >bak3")
os.system("cat bak3|sort -k2|tr '\t' '~' >bak4")



f=open("bak4")
a=f.read()
a=a.split("\n")
del a[-1]
f.close()

b=""
c=""
row1=""
row2=""

for line in a:
    d=line.split("~")
    last_row1=row1
    last_row2=row2
    row1=d[0]
    row2=d[1]
    print(row1)
    if row2==last_row2:
        c=c+"~"+row1
    else:
        c=c+"\n"+row2+"\t"+row1

print(c)
filename="unique_ko"
m=open(filename,'w')
m.write(c)



###
bak5=c
b=open("bak2","r")
bak2=b.read()
bak2=bak2.split("\n")
bak5=bak5.split("\n")
b.close()
del bak5[0]


k={}
bak6=""
for i in bak5:
    i=i.split("\t")
    tmp=i[0]
    k[tmp]=i[1]

for ii in bak2:
    ii=ii.split("\t")
    kk=ii[0]
    if ii[1] not in bak6:
        bak6+="\n" + ii[1]+"\t"+k[kk]
    else:
        bak6+="~"+k[kk]


filename="bak6"
m=open(filename,'w')
m.write(bak6)
m.close()

os.system("sed -i '/^$/d' bak6")

###
f=open("bak6")
bak6=f.read()
bak6=bak6.split("\n")
f.close()



b=""
c=""
for line in bak6:
    d=line.split("\t")
    b=d[0]
    na=d[1]
    na=na.split("~")
    for name in na:
        c=c+"\n"+name+"\t"+b 

bak7=c
filename="bak7"
m=open(filename,'w')
m.write(bak7)
m.close()

os.system("sed -i '/^$/d' bak7")
os.system("sed -i 's/ /\t/1' bak7")
os.system("sed -i '1 i Gene\tPathway\tDescription' bak7")


##整理gene2k

os.system("cat bak|grep 'ko:'|sed 's/^ko://g' |sort|uniq|sed 's/^ /\t/1'  >tmp1")

f=open("tmp1","r")
a=f.read()
a=a.split("\n")
f.close()

del a[-1]
k={}
for line in a:
    line=line.replace(" ","\t",1)
    line=line.split("\t")
    k[line[0]]=line[1]

f=open("bak3","r")
bak3=f.read()
f.close()
bak3=bak3.split("\n")
del bak3[-1]

gene2k=""
for line in bak3:
    line=line.split("\t")
    a=line[0]
    b=line[1]
    gene2k+=a+"\t"+b+"\t"+k[b]+"\n"


filename="gene2k"
m=open(filename,'w')
m.write(gene2k)
m.close()



os.system("mv bak7 kegg_pathway.txt")
os.system("rm bak* tmp1 unique_ko k")
