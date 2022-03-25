#!/usr/bin/bash
#SBATCH -J combine                   # 作业名为 test
#SBATCH -o 222test.out               # 屏幕上的输出文件重定向到 test.out
#SBATCH -p comp                   # 作业提交的分区
#SBATCH --qos=normal               # 作业使用的 QoS 为 debug ,normal。debug优先级高，可使用资源少
#SBATCH -N 1                      # 作业申请几个节点
#SBATCH -c 16






export PATH="/home/a2431559261/OLDISK/soft/gatk-4.1.9.0/:$PATH"


alias "gatk=/home/a2431559261/OLDISK/soft/gatk-4.1.9.0/gatk"

genome=fl.fasta
#R1=natu_G1.fastq.gz
#R2=natu_G2.fastq.gz
t=24

ID=ZW   #输入reads集的ID号 (0,1,2,3,4,5)
SM=     #样本名称 (0,1,2,3,4,5)
LB=     #reads 集的文库名
PL=Illumina     #测序平台
prefix=fl






### 4 combine gvcf
mkdir -p combine5
cd combine5

# 1) 多个样本结合
#find ../3_snp/ -name "*.gvcf" >input.list

gatk CombineGVCFs \
-R ../ref/$genome \
-O combine_variants.raw.gvcf.gz \
--variant input.list \

#--variant ../3_snp/89.raw.gvcf
#--variant ../3_snp/88.raw.gvcf

# 也可以用这个合并，速度快点
#gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
#      -V ../3_snp/90.raw.gvcf \
#      -V ../3_snp/89.raw.gvcf \
#      -V ../3_snp/88.raw.gvcf \
#      --genomicsdb-workspace-path ../data/${genome} \
#      -L group1




gatk GenotypeGVCFs \
-R ../ref/$genome \
-O combine_variants.raw.vcf.gz \
--variant combine_variants.raw.gvcf.gz
cd ..
#或者
#gatk --java-options "-Xmx4g" GenotypeGVCFs \
#   -R Homo_sapiens_assembly38.fasta \
#   -V gendb://my_database \
#   -O output.vcf.gz



