##读取比对文件列表,去掉   .out.fa
#OG1001011
#OG1001012
#OG1001013


f=open("sin_txt","r")
query_list=f.read()
query_list=query_list.split("\n")
del query_list[-1]
f.close()

##每个文件转变为msa
total_msa=""
phylip=""
for list in query_list:
    pre_msa=""
    pre_phylip=""
    tmp1=open("%s.fa.out"%list,"r")
    tmp1=tmp1.read()
    tmp1=tmp1.split("\n")
    del tmp1[-1]
    msa_name=list
    for i in tmp1:
        if i[0]==">" and pre_msa=="":
            pre_msa+=">"+list+"\n"+i[1:]+"\t"
            pre_phylip+=i[1:]+"\t"
        elif  i[0]==">" and pre_msa!="":
            pre_msa+="\n"+i[1:]+"\t"
            pre_phylip+="\n"+i[1:]+"\t"
        else:
            pre_msa+=i
            pre_phylip+=i
    total_msa+=pre_msa+"\n"
    filename="%s.msa" %list
    tmp2=open(filename,"w")
    tmp2.write(pre_msa)
    tmp2.close       
    #filename="%s.phy" %list  

    #print(pre_phylip)
    pre_phylip=pre_phylip.split("\n")
    
    if phylip=="":
        for j in pre_phylip:
            j=j.split("\t")
            seq=j[1]
            name=j[0].split("|")
            name=name[0]
            phylip+=name+"\t"+seq+"\n"
    else:
        k={}
        phy=phylip.split("\n")
        del phy[-1]
      
        for line1 in pre_phylip:
            seq=line1.split("\t")[1]
            sp_name=line1.split("\t")[0].split("|")[0]
            k[sp_name]=seq
        
        phylip=""
        for line2 in phy:
                tmp=line2.split("\t")
                phylip+=line2+k[tmp[0]]+"\n"  
                print(len(line2))             
         
a=phylip.split("\n")
b=a
a=a[1].split("\t")[1]
print(len(a))



file2name="total.msa"
file=open(file2name,"w")
file.write(total_msa)
file.close()

file3name="total.phylip"
file=open(file3name,"w")
file.write(phylip)
file.close()
