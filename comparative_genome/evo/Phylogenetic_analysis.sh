
#基因家族聚类
srun -c 48 orthofinder -f ./pep/ -S diamond -t 36


#进化
raxmlHPC-PTHREADS -f a -x 12345 -s out.phylip -# 2000  -p 12345 -m PROTGAMMAAUTO -n out -T 30






