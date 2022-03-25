#!/usr/bin/bash
#SBATCH -J snp                   # 作业名为 test
#SBATCH -o sv_ge_snp               # 屏幕上的输出文件重定向到 test.out
#SBATCH -p comp                   # 作业提交的分区
#SBATCH --qos=normal               # 作业使用的 QoS 为 debug ,normal。debug优先级高，可使用资源少
#SBATCH -N 1                      # 作业申请几个节点
#SBATCH -c 24




query=da.fna
ref=shamen.fna


source /home/a2431559261/OLDISK/soft/miniconda3/etc/profile.d/conda.sh
conda activate sv


mkdir -p sv_map
cd sv_map

## 1 基于基因组的SNP 鉴定
nucmer -g 1000 -c 90 -l 40 ../data/${ref} ../data/${query}

delta-filter -r -q out.delta >filter.delta


show-snps -ClrT filter.delta >snp.${ref}

show-snps -ClqT filter.delta >snp.${query}




## 2 基于reads 的SNP INDEL 寻找
#conda activate gatk
#sed -i "s/genome=/genome=${ref}/g" sv_gatk.sh  
#sed -i "s/ref_index_dir=/ref_index_dir=${ref_index_dir}/g" sv_gatk.sh
#sed -i "s/R1=/R1=${R1}/g" sv_gatk.sh
#sed -i "s/R1=/R2=${R2}/g" sv_gatk.sh

#bash sv_gatk.sh











