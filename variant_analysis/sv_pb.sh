#!/bin/sh
#SBATCH -J pav_pbmap                   # 作业名为 test
#SBATCH -o pb_map.out               # 屏幕上的输出文件重定向到 test.out
#SBATCH -p comp                   # 作业提交的分区
#SBATCH --qos=normal               # 作业使用的 QoS 为 debug ,normal。debug优先级高，可使用资源少
#SBATCH -N 1                      # 作业申请几个节点
#SBATCH -c 24



source /home/a2431559261/OLDISK/soft/miniconda3/etc/profile.d/conda.sh
conda activate sv

ref_genome=f153.fasta
query_reads=m64033_200214_085627.subreads.fasta
query_pre=gl_f153
work_dir=`pwd`

mkdir -p pav_pb 

#cd /mnt/memorydisk

#ln -s ${work_dir}/data/ .

cd pav_pb

ngmlr -t 24 -r ../data/${ref_genome} -q ../data/${query_reads} -o ${query_pre}.sam


awk -F "\t" '{if ($5 >= 0 || substr ($1, 0, 1)=="@") print $0}' ${query_pre}.sam > ${query_pre}.filter.sam




samtools   view -b ${query_pre}.filter.sam|samtools sort -o ${query_pre}.sorted.bam  - 



sniffles -m ${query_pre}.sorted.bam -v ${query_pre}.vcf

#rm -rf data


#rsync  -auP * ${work_dir}/pav_pb

#rm -rf *
 


