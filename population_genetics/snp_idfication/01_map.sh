#!/usr/bin/bash
#SBATCH -J gatk                   # 作业名为 test
#SBATCH -o subab               # 屏幕上的输出文件重定向到 test.out
#SBATCH -p big                     # 作业提交的分区
#SBATCH --qos=normal               # 作业使用的 QoS 为 debug ,normal。debug优先级高，可使用资源少
#SBATCH -N 1                      # 作业申请几个节点
#SBATCH -c 10



## 可能需要基因组的染色体顺序按照 1到25 排序排好


# for i in `cat sample.txt` 

export PATH="/home/a2431559261/OLDISK/soft/gatk-4.1.9.0/:$PATH"


alias "gatk=/home/a2431559261/OLDISK/soft/gatk-4.1.9.0/gatk"

genome=fl.fasta
R1=subab_1.fastq
R2=subab_1.fastq
t=10
ID=subab   #输入reads集的ID号 (0,1,2,3,4,5)
SM=subab	#样本名称 (0,1,2,3,4,5)
LB=	#reads 集的文库名
#PL=Illumina     #测序平台 
prefix=subab


### 1 mapping
mkdir -p 1_maping
mkdir -p 2_dedup
mkdir -p 3_snp


cd 1_maping

#01 建立索引
######bwa index ../ref/${genome}

# -a 制定建立索引的算法
# bwtsw  适合基因组大于10m
# is 适合ref < 2G 这是默认选项
# div 不工作 with long geno

samtools faidx ../ref/${genome}
#samtools 建立索引的目的是方便GATK能够快速获取fasta上的序列

##02 生成.dict文件
gatk CreateSequenceDictionary \
-R database/$genome \


#02 mapping sort

bwa mem -M -t $t -R "@RG\tID:${ID}\tSM:${ID}" ../ref/${genome} ../data/$R1 ../data/$R2|samtools view -b -|samtools sort -o ${ID}.sorted.bam  -

#bwa mem -M -t $t -R "@RG\tID:07\tSM:07" ../ref/${genome} ../data/$R1 ../data/$R2 >${ID}.sam
#samtools view -b ${ID}.sam >${ID}.bam

#samtools sort -o ${ID}.sorted.bam  ${ID}.bam
samtools index ${ID}.sorted.bam

# samtools flags 99
cd ..

### 2 dedup  去除重复
mkdir -p 2_dedup
cd 2_dedup

gatk MarkDuplicates \
-I ../1_maping/${ID}.sorted.bam \
-O ${ID}.dedup.bam \
-M ${ID}.mark_dup_metrics.txt

samtools index ${ID}.dedup.bam

cd ..

#NOTE 若原始SAM/BAM文件没有read group和sample信息，可以通过AddOrReplaceReadGroups添加这部分信息
#java -jar picard.jar AddOrReplaceReadGroups \
#      I=input.bam \
#      O=output.bam \
#      RGID=4 \
#      RGLB=lib1 \
#      RGPL=illumina \
#      RGPU=unit1 \
#      RGSM=20






### 3 call snp
mkdir -p 3_snp
cd 3_snp
gatk CreateSequenceDictionary -R ../ref/${genome}
gatk HaplotypeCaller -R ../ref/${genome} -I ../2_dedup/${ID}.dedup.bam -O ${ID}.raw.gvcf.gz -ERC GVCF -ploidy 2  2> log
cd ..


