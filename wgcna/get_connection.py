import  sys
node = sys.argv[1]
edge= sys.argv[2]

with open(node, 'r') as  f:
    f = f.read().split('\n')
while '' in f:
    f.remove('')

gene_list=[]
f=f[1:]
for line in f:
    line=line.split('\t')
    gene_list.append(line[0])
    
with open(edge, 'r') as  f:
    f = f.read().split('\n')
while '' in f:
    f.remove('')
f=f[1:]
    
connec_dict={}
for i in gene_list:
    connec_dict[i]=0
    
for line in f:
    line=line.split('\t')
    gene1=line[0]
    gene2=line[1]
    weight=float(line[2])
    gene_1_value=connec_dict[gene1]
    gene_2_value=connec_dict[gene2]
    gene_1_value += weight
    gene_2_value += weight
    connec_dict[gene1]=gene_1_value
    connec_dict[gene2]=gene_2_value

for i in connec_dict.keys():
    print('{i}\t{connec}'.format(i=i,connec=str(connec_dict[i])))
