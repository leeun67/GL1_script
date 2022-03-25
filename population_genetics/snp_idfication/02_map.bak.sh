#!/usr/bin/bash
#SBATCH -J gatk                   # 作业名为 test
#SBATCH -o 8195               # 屏幕上的输出文件重定向到 test.out
#SBATCH -p comp                   # 作业提交的分区
#SBATCH --qos=normal               # 作业使用的 QoS 为 debug ,normal。debug优先级高，可使用资源少
#SBATCH -N 1                      # 作业申请几个节点
#SBATCH -c 16



export PATH="/home/a2431559261/OLDISK/soft/gatk-4.1.9.0/:$PATH"


alias "gatk=/home/a2431559261/OLDISK/soft/gatk-4.1.9.0/gatk"

genome=fl.fasta
#R1=natu_G1.fastq.gz
#R2=natu_G2.fastq.gz
t=12

#ID=ZW   #输入reads集的ID号 (0,1,2,3,4,5)
#SM=     #样本名称 (0,1,2,3,4,5)
#LB=     #reads 集的文库名
#PL=Illumina     #测序平台
prefix=fl

### 1 mapping
for ID in `cat 81_95`
do
#mkdir -p 1_maping
cd 1_maping


#02 mapping sort
bwa mem -M -t $t -R "@RG\tID:${ID}\tSM:${ID}" ../ref/${genome} ../data/SRR59638${ID}_1.fastq.gz ../data/SRR59638${ID}_1.fastq.gz |samtools view -@ 12 -b -|samtools sort -@ 12 -o ${ID}.sorted.bam  -

samtools index ${ID}.sorted.bam


#用picard 中的sortSam 进行排序 他会自动添加一个SO标签说明文件被成功排序
#bwa mem -M -t $t -R "@RG\tID:N\tSM:G" ../ref/${genome} ../data/$R1 ../data/$R2 >${ID}.sam
#gatk SortSam \
#-I ${ID}.sam \
#-O ${ID}.sort.bam \
#-R ../ref/$genome \
#-SO coordinate --CREATE_INDEX
#samtools view -H /path/to/my.bam 检查文件是否成功排序


#bwa mem -M -t $t -R "@RG\tID:N\tSM:G" ../ref/${genome} ../data/$R1 ../data/$R2 >${ID}.sam
#samtools view -b ${ID}.sam >${ID}.bam
#samtools sort -o ${ID}.sorted.bam  ${ID}.bam
#samtools index ${ID}.sorted.bam

# samtools flags 99
cd ..


### 2 dedup  标记重复
#mkdir -p 2_dedup
cd 2_dedup

gatk MarkDuplicates \
-I ../1_maping/${ID}.sorted.bam \
-O ${ID}.dedup.bam \
-M ${ID}.mark_dup_metrics.txt  --CREATE_INDEX

#samtools index ${ID}.dedup.bam

cd ..





### 3 SNP,INDEL 位点识别
#SNP call 策略
#single sample calling：每一个sample的bam file都进行单独的snp calling，然后每个sample单独snp calling结果再合成一个总的snp calling的结果。

#batch calling： 一定数目群集的bamfiles 一起calling snps，然后再merge在一起

#joint calling： 所有samples的BAM files一起call 出一个包含所有samples 变异信息的output

#小编推荐joint calling ,包含两步
#1 单独为每个样本进行从0 到 haplotype 生成 gVCF中间文件的过程
#2 依据gvcf完成群体的joint calling

#mkdir -p 3_snp
cd 3_snp
gatk HaplotypeCaller \
-R ../ref/${genome} \
-I ../2_dedup/${ID}.dedup.bam \
-O ${ID}.raw.gvcf \
-ERC GVCF \
-ploidy 2

#-L 参数指定识别突变位点的区域
cd ..

done


# 在第三步之前，可以有一个质量值矫正的过程，不过需要一个SNP库，如果没有，可以借助samtools 和GATK创建一个
#https://ming-lian.github.io/2019/02/08/call-snp/
#1 建立矫正模型  gatk BaseRecalibrator
#2 质量值矫正 gatk ApplyBQSR
