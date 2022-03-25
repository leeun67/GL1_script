#opt
#克隆 juicer
#git clone https://github.com/theaidenlab/juicer.git
#ln -s juicer/CPU scripts
#cd scripts/common
#wget http://hicfiles.tc4ga.com.s3.amazonaws.com/public/juicer/juicer_tools.1.7.6_jcuda.0.8.jar
#ln -s juicer_tools.1.7.6_jcuda.0.8.jar juicer_tools.jar
#cd ../..

# 需要bwa,samtools 1.10 bedtools

# 如果没有N卡，把common目录下的post...hiccup 换一下顺序（加上--cpu参数）
#juicerbox 需要 java -jar 来运行


#source activate falcon

# 参考基因组建立索引
#mkdir references
#cd references 
#ln -s /home/a2431559261/FH_L.contigs.fa.gz .
#bwa index curated.fasta
#cd ..


# 添加限制性内切酶位点信息，注意自己的是不是 MboI 酶。
#mkdir restriction_sites
#cd restriction_sites/
#python2 /home/a2431559261/OLDISK/soft/juicer/misc/generate_site_positions.py  MboI  fl /home/a2431559261/OLDISK/genome_assemble/juicer/references/curated.fasta # 生成了 hg38_MboI.txt 文件

#awk 'BEGIN{OFS="\t"}{print $1, $NF}' fl_MboI.txt > fl.chrom.sizes
#cd ..



# 添加 fastq 文件 (官方测试文件)
#mkdir fastq
#cd fastq
#ln -s /home/a2431559261/FL_seq/pineapple/HiC/P-G_1.fq.gz P-G_R1.fastq.gz  
#ln -s /home/a2431559261/FL_seq/pineapple/HiC/P-G_2.fq.gz P-G_R2.fastq.gz
#cd ..

#mkdir fastq && cd fastq
#nohub wget http://juicerawsmirror.s3.amazonaws.com/opt/juicer/work/HIC003/fastq/HIC003_S2_L001_R1_001.fastq.gz &
#nohub wget http://juicerawsmirror.s3.amazonaws.com/opt/juicer/work/HIC003/fastq/HIC003_S2_L001_R2_001.fastq.gz &
#cd ..


# 运行 Juicer


bash /home/a2431559261/OLDISK/soft/juicer/scripts/juicer.sh  -g fl -S postproc -d /home/a2431559261/OLDISK/genome_assemble/juicer   -D /home/a2431559261/OLDISK/soft/juicer  -y restriction_sites/fl_MboI.txt  -z /home/a2431559261/OLDISK/genome_assemble/juicer/references/curated.fasta -p /home/a2431559261/OLDISK/genome_assemble/juicer/restriction_sites/fl.chrom.sizes -s MboI -t 24
# 上面-z -p 设置绝对路径吧，要不然会识别不了。我直接改软件源码了，-d -D 都设置了还需要绝对路径，作者写的不人性化啊


