#!/usr/bin/bash
#SBATCH -J sv_svp_map                   # 作业名为 test
#SBATCH -o map_snp.out               # 屏幕上的输出文件重定向到 test.out
#SBATCH -p comp                   # 作业提交的分区
#SBATCH --qos=normal               # 作业使用的 QoS 为 debug ,normal。debug优先级高，可使用资源少
#SBATCH -N 1                      # 作业申请几个节点
#SBATCH -c 22


source /home/a2431559261/OLDISK/soft/miniconda3/etc/profile.d/conda.sh

ref_genome=cb5.fasta
ref_index_dir=~/OLDISK/dataset/genome/2_cb5/index
R1=fengli-L_L4_323323_1.clean.fq.gz
R2=fengli-L_L4_323323_2.clean.fq.gz
t=24
ID=cb5



export PATH="/home/a2431559261/OLDISK/soft/gatk-4.1.9.0/:$PATH"

alias "gatk=/home/a2431559261/OLDISK/soft/gatk-4.1.9.0/gatk"

conda activate gatk

mkdir -p gatk
cd gatk


#samtools faidx    data/${genome}
#bwa index data/${genome}
#gatk CreateSequenceDictionary -R data/$genome 



##1 mapping
#bwa mem -M -t $t   -R "@RG\tID:${ID}\tSM:${ID}"    ${ref_index_dir}/${ref_genome}  ../data/$R1 ../data/$R2 | \
#samtools view -b -|samtools sort -o ${ref_genome}.sorted.bam  -


##2 dedup 去除重复
#gatk MarkDuplicates \
#-I ${ref_genome}.sorted.bam \
#-O ${ref_genome}.dedup.bam \
#-M ${ref_genome}.mark_dup_metrics.txt

#samtools index ${ref_genome}.dedup.bam


##3 call_snp
gatk HaplotypeCaller -R ${ref_index_dir}/${ref_genome} -I ${ref_genome}.dedup.bam -O ${ref_genome}.raw.gvcf -ERC GVCF -ploidy 2  2> log


##4 combine vcf


##5 分型
gatk GenotypeGVCFs \
-R ${ref_index_dir}/${ref_genome} \
-O combine_variants.raw.vcf \
--variant ${ref_genome}.raw.gvcf


##6 提取SNP变异位点并过滤
gatk SelectVariants -R ${ref_index_dir}/${ref_genome}  \
-O combine_SNP.raw.vcf \
--variant combine_variants.raw.vcf \
--select-type-to-include SNP

gatk VariantFiltration \
-R ${ref_index_dir}/${ref_genome}  \
-O combine_SNP.filtered.vcf  \
--variant combine_SNP.raw.vcf \
--filter-name "snp_filter" \
--filter-expression "QD<2.0 || FS>60.0 || SOR >3.0 ||MQ <40.0 || MQRankSum < -12.5 "


gatk SelectVariants  -R ${ref_index_dir}/${ref_genome} \
-V combine_SNP.filtered.vcf \
-O fl.pass.snps.genotype.vcf \
-select "vc.isNotFiltered()"


## 提取INDEL 变异位点并过滤
gatk SelectVariants -R ${ref_index_dir}/${ref_genome}  \
-O combine_INDEL.raw.vcf  \
--variant combine_variants.raw.vcf \
--select-type-to-include INDEL



#过滤
gatk VariantFiltration \
-R ${ref_index_dir}/${ref_genome} \
-O combine_INDEL.filtered.vcf  \
--variant combine_INDEL.raw.vcf \
--filter-name "indel_filter" \
--filter-expression "QD<2.0 || FS>200.0 || SOR >10.0 ||MQ <40.0 || MQRankSum < -12.5 "



gatk SelectVariants  -R ${ref_index_dir}/${ref_genome} \
-V combine_INDEL.filtered.vcf \
-O fl.pass.indel.genotype.vcf \
-select "vc.isNotFiltered()"

vcftools  --vcf combine_SNP.filtered.vcf --recode --recode-INFO-all --stdout --max-missing 0.9 --maf 0.05 --minDP 4 > final.snp.vcf


gatk SelectVariants  -R ${ref_index_dir}/${ref_genome} \
-V final.snp.vcf \
-O 2final.pass.snps.vcf \
-select "vc.isNotFiltered()"




